#!/usr/bin/env python3
"""Replay every row-owner assignment of a saved relaxed mask family."""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from pathlib import Path

import exact_greedy
from analyze_relaxed_mask_dual import (
    build_gain_mask_support,
    load_point_file,
)
from local_phase_cegis import build_targets


def enumerate_owner_assignments(
    models: list[list[int]],
    support: dict[int, list[tuple[int, int]]],
    assignment: list[int],
    targets,
    rows: list[dict],
    np,
) -> dict:
    assignment_array = np.asarray(assignment, dtype=targets.dtype)
    base_cover = np.count_nonzero(
        targets == assignment_array,
        axis=1,
    )
    raw_owner_products = 0
    distinct_row_trials = 0
    exact_repairs = []
    replay_failures = 0
    seen_repairs = set()
    for model_index, model in enumerate(models, 1):
        missing = [mask for mask in model if mask not in support]
        if missing:
            raise RuntimeError(
                f"model {model_index} contains unsupported masks"
            )
        owner_lists = [support[mask] for mask in model]
        for owners in itertools.product(*owner_lists):
            raw_owner_products += 1
            if len({row_index for row_index, _target in owners}) < len(
                owners
            ):
                continue
            distinct_row_trials += 1
            cover = base_cover.astype(np.int32, copy=True)
            for row_index, target in owners:
                cover -= (
                    targets[:, row_index]
                    == int(assignment[row_index])
                )
                cover += targets[:, row_index] == int(target)
            misses = int(np.count_nonzero(cover == 0))
            if misses:
                replay_failures += 1
                continue
            moves = tuple(
                sorted(
                    (
                        int(rows[row_index]["p"]),
                        int(target),
                    )
                    for row_index, target in owners
                )
            )
            if moves in seen_repairs:
                continue
            seen_repairs.add(moves)
            exact_repairs.append(
                {
                    "source_model": model_index,
                    "moves": [
                        {"p": prime, "target": target}
                        for prime, target in moves
                    ],
                }
            )
    exact_repairs.sort(
        key=lambda repair: tuple(
            (move["p"], move["target"]) for move in repair["moves"]
        )
    )
    union_moves = sorted(
        {
            (move["p"], move["target"])
            for repair in exact_repairs
            for move in repair["moves"]
        }
    )
    return {
        "raw_owner_products": raw_owner_products,
        "distinct_row_owner_trials": distinct_row_trials,
        "replay_failures": replay_failures,
        "exact_repair_count": len(exact_repairs),
        "exact_repairs": exact_repairs,
        "union_moves": [
            {"p": prime, "target": target}
            for prime, target in union_moves
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path, nargs="+")
    parser.add_argument("--base-phase", type=Path, required=True)
    parser.add_argument("--fixed-primes", default="")
    parser.add_argument("--models", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

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
    assignment = [
        phases[int(row["p"])] % int(row["h"]) for row in rows
    ]
    targets, build_seconds = build_targets(points, candidates, np)
    missed_points, support = build_gain_mask_support(
        rows,
        candidates,
        points,
        phases,
        fixed_primes,
        np,
    )
    model_payload = json.loads(args.models.read_text())
    if not model_payload.get("complete"):
        raise RuntimeError("saved relaxed model enumeration is incomplete")
    models = [
        [int(mask) for mask in model]
        for model in model_payload.get("models", [])
    ]
    if len(models) != int(model_payload["cover_count"]):
        raise RuntimeError("saved model count does not match cover count")
    enumeration = enumerate_owner_assignments(
        models,
        support,
        assignment,
        targets,
        rows,
        np,
    )
    result = {
        "problem": "exact row-owner replay of saved relaxed mask covers",
        "scope": (
            "the embedded pool, declared base phase and fixed rows, supplied "
            "finite point corpus, and saved complete relaxed model family"
        ),
        "pool": str(args.pool),
        "points": [str(path) for path in args.points],
        "base_phase": str(args.base_phase),
        "fixed_primes": sorted(fixed_primes),
        "models": str(args.models),
        "input_point_count": len(points),
        "base_miss_count": len(missed_points),
        "relaxed_model_count": len(models),
        "gain_mask_count": len(support),
        "build_seconds": build_seconds,
        **enumeration,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"models={len(models)} owner_trials="
        f"{result['distinct_row_owner_trials']} "
        f"exact_repairs={result['exact_repair_count']} "
        f"replay_failures={result['replay_failures']} "
        f"union_moves={len(result['union_moves'])}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
