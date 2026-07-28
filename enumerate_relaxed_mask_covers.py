#!/usr/bin/env python3
"""Enumerate the relaxed gain-mask cover space of a finite phase repair.

This exposes backbone gain patterns and their row owners before the more
expensive distinct-row and secondary-loss replay.  A complete enumeration is
an exact statement about the relaxed finite mask system.  A limited or timed
enumeration is exploratory only.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import exact_greedy
from analyze_relaxed_mask_dual import (
    build_gain_mask_support,
    load_point_file,
)


def enumerate_exact_radius_covers(
    masks: list[int],
    point_count: int,
    radius: int,
    solver_name: str,
    model_limit: int,
    time_limit: float,
    save_models: bool = False,
) -> dict:
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    vpool = IDPool()
    variables = {
        mask: vpool.id(("gain_mask", mask)) for mask in masks
    }
    clauses = CardEnc.equals(
        list(variables.values()),
        bound=radius,
        vpool=vpool,
        encoding=EncType.seqcounter,
    ).clauses
    for bit in range(point_count):
        clauses.append(
            [
                variables[mask]
                for mask in masks
                if mask & (1 << bit)
            ]
        )
    solver = Solver(name=solver_name, bootstrap_with=clauses)
    started = time.monotonic()
    models = 0
    frequencies: Counter[int] = Counter()
    common: set[int] | None = None
    saved_models: list[list[int]] = []
    complete = False
    stop_reason = "UNSAT"
    while solver.solve():
        if model_limit and models >= model_limit:
            stop_reason = "MODEL_LIMIT"
            break
        if time_limit > 0 and time.monotonic() - started >= time_limit:
            stop_reason = "TIME_LIMIT"
            break
        positive = {
            literal for literal in solver.get_model() if literal > 0
        }
        selected = {
            mask for mask in masks if variables[mask] in positive
        }
        if len(selected) != radius:
            raise AssertionError("cardinality encoding returned wrong radius")
        models += 1
        frequencies.update(selected)
        common = selected if common is None else common & selected
        if save_models:
            saved_models.append(sorted(selected))
        solver.add_clause([-variables[mask] for mask in selected])
    else:
        complete = True
    solver.delete()
    result = {
        "complete": complete,
        "stop_reason": stop_reason,
        "cover_count": models,
        "observed_common_masks": sorted(common or ()),
        "mask_frequencies": frequencies,
        "elapsed_seconds": time.monotonic() - started,
    }
    if save_models:
        result["models"] = saved_models
    return result


def mask_record(
    mask: int,
    frequency: int,
    support: dict[int, list[tuple[int, int]]],
    rows: list[dict],
    point_count: int,
) -> dict:
    moves = support[mask]
    return {
        "point_indices": [
            bit for bit in range(point_count) if mask & (1 << bit)
        ],
        "point_count": mask.bit_count(),
        "observed_model_frequency": frequency,
        "supporting_move_count": len(moves),
        "supporting_moves": [
            {
                "row_index": row_index,
                "p": int(rows[row_index]["p"]),
                "target": int(target),
            }
            for row_index, target in moves
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path, nargs="+")
    parser.add_argument("--base-phase", type=Path, required=True)
    parser.add_argument("--fixed-primes", default="")
    parser.add_argument("--radius", type=int, required=True)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--model-limit", type=int, default=100_000)
    parser.add_argument("--time-limit", type=float, default=240.0)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--save-models",
        action="store_true",
        help=(
            "store every selected gain-mask tuple; use only when the "
            "expected model family is small enough for the output file"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.radius < 1:
        raise SystemExit("--radius must be positive")
    if args.model_limit < 0 or args.time_limit < 0 or args.top < 0:
        raise SystemExit("limits must be nonnegative")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    candidates = exact_greedy.load_candidates(args.pool, True)
    points = []
    seen = set()
    for path in args.points:
        for point in load_point_file(path):
            if point not in seen:
                seen.add(point)
                points.append(point)
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(
            args.base_phase.read_text()
        ).items()
    }
    fixed_primes = {
        int(value) for value in args.fixed_primes.split(",") if value
    }
    missed_points, support = build_gain_mask_support(
        rows,
        candidates,
        points,
        phases,
        fixed_primes,
        np,
    )
    masks = sorted(support)
    enumeration = enumerate_exact_radius_covers(
        masks,
        len(missed_points),
        args.radius,
        args.solver,
        args.model_limit,
        args.time_limit,
        args.save_models,
    )
    frequencies: Counter[int] = enumeration.pop("mask_frequencies")
    ordered = sorted(
        frequencies,
        key=lambda mask: (
            -frequencies[mask],
            -mask.bit_count(),
            mask,
        ),
    )
    common_masks = enumeration["observed_common_masks"]
    result = {
        "problem": "relaxed exact-radius gain-mask cover space",
        "scope": (
            "the embedded pool, base phase, fixed rows, supplied points, "
            "and exact declared number of selected gain masks"
        ),
        "pool": str(args.pool),
        "base_phase": str(args.base_phase),
        "point_sources": [str(path) for path in args.points],
        "fixed_primes": sorted(fixed_primes),
        "radius": args.radius,
        "input_point_count": len(points),
        "base_miss_count": len(missed_points),
        "gain_mask_count": len(masks),
        **enumeration,
        "observed_common_mask_records": [
            mask_record(
                mask,
                frequencies[mask],
                support,
                rows,
                len(missed_points),
            )
            for mask in common_masks
        ],
        "top_observed_mask_records": [
            mask_record(
                mask,
                frequencies[mask],
                support,
                rows,
                len(missed_points),
            )
            for mask in ordered[: args.top]
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"points={len(points)} misses={len(missed_points)} "
        f"masks={len(masks)} radius={args.radius} "
        f"covers={result['cover_count']} complete={result['complete']} "
        f"common={len(common_masks)} stop={result['stop_reason']} "
        f"seconds={result['elapsed_seconds']:.3f}",
        flush=True,
    )
    return 0 if result["complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
