#!/usr/bin/env python3
"""Minimize a finite radius obstruction already UNSAT in its relaxed mask layer."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import exact_greedy
from certify_anchor_phase_quotient import load_cells
from finite_sample_mask_repair import search_mask_repair
from local_phase_cegis import build_targets


def greedy_unsat_assumption_core(solver, assumptions: list[int]) -> list[int]:
    if solver.solve(assumptions=assumptions):
        raise RuntimeError("assumption model is SAT")
    core_set = set(solver.get_core() or assumptions)
    core = [literal for literal in assumptions if literal in core_set]
    for literal in list(core):
        trial = [value for value in core if value != literal]
        if not solver.solve(assumptions=trial):
            core = trial
    return core


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path, nargs="+")
    parser.add_argument("--base-phase", type=Path, required=True)
    parser.add_argument("--fixed-primes", required=True)
    parser.add_argument("--max-changes", type=int, required=True)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    if args.max_changes < 0:
        raise SystemExit("--max-changes must be nonnegative")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    candidates = exact_greedy.load_candidates(args.pool, True)
    points = []
    seen = set()
    for path in args.points:
        for point in load_cells(path):
            if point not in seen:
                seen.add(point)
                points.append(point)
    initial = {
        int(prime): int(target)
        for prime, target in json.loads(
            args.base_phase.read_text()
        ).items()
    }
    fixed_primes = {
        int(value) for value in args.fixed_primes.split(",") if value
    }
    assignment = []
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        target = int(initial.get(prime, residue)) % h
        if h % modulus or target % modulus != residue:
            raise RuntimeError(f"invalid initial phase for p={prime}")
        assignment.append(target)

    targets, build_seconds = build_targets(points, candidates, np)
    assignment_array = np.asarray(assignment, dtype=targets.dtype)
    base_cover = np.count_nonzero(
        targets == assignment_array,
        axis=1,
    )
    miss_indices = [
        int(index) for index in np.flatnonzero(base_cover == 0)
    ]
    masks = set()
    for row_index, row in enumerate(rows):
        if int(row["p"]) in fixed_primes:
            continue
        current = assignment[row_index]
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        by_target: dict[int, int] = {}
        for bit, point_index in enumerate(miss_indices):
            target = int(targets[point_index, row_index])
            if target == current or target % modulus != residue:
                continue
            by_target[target] = (
                by_target.get(target, 0) | (1 << bit)
            )
        masks.update(mask for mask in by_target.values() if mask)
    ordered_masks = sorted(masks)

    vpool = IDPool()
    mask_vars = [vpool.id(("mask", mask)) for mask in ordered_masks]
    clauses = CardEnc.atmost(
        mask_vars,
        bound=args.max_changes,
        vpool=vpool,
        encoding=EncType.seqcounter,
    ).clauses
    activation_vars = []
    for bit in range(len(miss_indices)):
        activation = vpool.id(("point", bit))
        activation_vars.append(activation)
        clauses.append(
            [
                -activation,
                *[
                    variable
                    for mask, variable in zip(
                        ordered_masks,
                        mask_vars,
                    )
                    if mask & (1 << bit)
                ],
            ]
        )
    solver = Solver(name=args.solver, bootstrap_with=clauses)
    started = time.monotonic()
    core = greedy_unsat_assumption_core(solver, activation_vars)
    solver.delete()
    core_bits = {
        activation_vars.index(literal) for literal in core
    }
    core_points = [
        points[miss_indices[bit]]
        for bit in sorted(core_bits)
    ]
    replay = search_mask_repair(
        rows,
        candidates,
        core_points,
        initial,
        fixed_primes,
        args.max_changes,
        args.solver,
        0,
        0.0,
    )
    if replay["status"] != "UNSAT" or not replay["complete_negative"]:
        raise RuntimeError("minimized relaxed core failed complete replay")
    args.output.write_text(
        json.dumps([[k, l] for k, l in core_points]) + "\n"
    )
    report = {
        "pool": str(args.pool),
        "point_sources": [str(path) for path in args.points],
        "base_phase": str(args.base_phase),
        "fixed_primes": sorted(fixed_primes),
        "max_changes": args.max_changes,
        "input_point_count": len(points),
        "initial_miss_count": len(miss_indices),
        "relaxed_gain_mask_count": len(ordered_masks),
        "minimized_assumption_core_count": len(
            set(activation_vars) & set(core)
        ),
        "minimized_point_count": len(core_points),
        "complete_replay": {
            key: value
            for key, value in replay.items()
            if key != "repaired"
        },
        "build_seconds": build_seconds,
        "elapsed_seconds": time.monotonic() - started,
    }
    args.report_output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"input={len(points)} misses={len(miss_indices)} "
        f"masks={len(ordered_masks)} core={len(core_points)} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
