#!/usr/bin/env python3
"""Exact minimum-phase-change repair of a finite affine-cover critical core."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

from local_phase_cegis import build_targets


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path)
    parser.add_argument("--initial-phases", type=Path, required=True)
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    parser.add_argument("--fixed-primes", default="")
    parser.add_argument(
        "--core-max-coverage",
        type=int,
        default=1,
        help="repair exactly the points covered at most this many times",
    )
    parser.add_argument(
        "--option-max-coverage",
        type=int,
        help=(
            "offer changed phases observed on points up to this initial "
            "coverage; defaults to --core-max-coverage"
        ),
    )
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument(
        "--sat-only",
        action="store_true",
        help="seek any hard-core repair without minimizing phase changes",
    )
    parser.add_argument(
        "--max-changes",
        type=int,
        help=(
            "hard upper bound on changed mutable phases; implies a SAT "
            "feasibility search instead of MaxSAT optimization"
        ),
    )
    args = parser.parse_args()
    if args.core_max_coverage < 0:
        raise SystemExit("--core-max-coverage must be nonnegative")
    if args.option_max_coverage is not None and args.option_max_coverage < 0:
        raise SystemExit("--option-max-coverage must be nonnegative")
    if args.max_changes is not None and args.max_changes < 0:
        raise SystemExit("--max-changes must be nonnegative")
    if args.max_changes is not None:
        args.sat_only = True

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool, WCNF  # type: ignore
    from pysat.examples.rc2 import RC2  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    points = [
        (int(k), int(l)) for k, l in json.loads(args.points.read_text())
    ]
    initial = {
        int(prime): int(target)
        for prime, target in json.loads(args.initial_phases.read_text()).items()
    }
    fixed_primes = {
        int(value) for value in args.fixed_primes.split(",") if value
    }
    candidates = [
        (
            int(row["h"]),
            int(row["p"]),
            int(row["a"]),
            int(row["b"]),
            int(row.get("ord2", row["h"])),
            int(row.get("ord3", row["h"])),
        )
        for row in rows
    ]
    assignment = np.empty(len(rows), dtype=np.uint32)
    for row_index, row in enumerate(rows):
        prime = int(row["p"])
        h = int(row["h"])
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        target = initial.get(prime, residue) % h
        if target % modulus != residue:
            raise RuntimeError(f"invalid initial phase for p={prime}")
        assignment[row_index] = target
    missing_fixed = fixed_primes - {int(row["p"]) for row in rows}
    if missing_fixed:
        raise RuntimeError(f"fixed primes absent from pool: {missing_fixed}")

    started = time.monotonic()
    targets, build_seconds = build_targets(points, candidates, np)
    full_cover = np.count_nonzero(targets == assignment, axis=1)
    core_indices = np.flatnonzero(
        full_cover <= args.core_max_coverage
    )
    if not len(core_indices):
        args.phase_output.write_text(args.initial_phases.read_text())
        result = {
            "pool": str(args.pool),
            "points": str(args.points),
            "initial_phases": str(args.initial_phases),
            "point_count": len(points),
            "core_point_count": 0,
            "initial_misses": 0,
            "result": "ALREADY_COVERED",
            "full_misses": 0,
            "build_seconds": build_seconds,
            "elapsed_seconds": time.monotonic() - started,
        }
        args.result_output.write_text(json.dumps(result, indent=2) + "\n")
        print("result=ALREADY_COVERED full_misses=0", flush=True)
        return 0
    core_targets = targets[core_indices]
    option_max_coverage = (
        args.core_max_coverage
        if args.option_max_coverage is None
        else args.option_max_coverage
    )
    option_indices = np.flatnonzero(full_cover <= option_max_coverage)
    option_targets = targets[option_indices]
    initial_misses = int(np.count_nonzero(full_cover == 0))
    print(
        f"rows={len(rows)} points={len(points)} core={len(core_indices)} "
        f"option_points={len(option_indices)} initial_misses={initial_misses} "
        f"matrix_s={build_seconds:.3f}",
        flush=True,
    )

    vpool = IDPool()
    wcnf = WCNF()
    target_variables: list[dict[int, int]] = []
    hard_clause_count = 0
    option_count = 0
    for row_index, row in enumerate(rows):
        prime = int(row["p"])
        h = int(row["h"])
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        current = int(assignment[row_index])
        observed = np.unique(option_targets[:, row_index])
        valid_targets = {
            int(target)
            for target in observed
            if int(target) % modulus == residue
        }
        valid_targets.add(current)
        if prime in fixed_primes:
            valid_targets = {current}
        ordered_targets = sorted(valid_targets)
        variable_by_target = {
            target: vpool.id(("target", row_index, target))
            for target in ordered_targets
        }
        target_variables.append(variable_by_target)
        variables = list(variable_by_target.values())
        option_count += len(variables)
        if prime in fixed_primes:
            wcnf.append([variable_by_target[current]])
            hard_clause_count += 1
        else:
            encoding = CardEnc.equals(
                variables,
                bound=1,
                vpool=vpool,
                encoding=EncType.seqcounter,
            )
            for clause in encoding.clauses:
                wcnf.append(clause)
            hard_clause_count += len(encoding.clauses)
            # One unit of cost for changing this row's phase.
            wcnf.append([variable_by_target[current]], weight=1)

    for core_row in core_targets:
        clause = [
            target_variables[row_index][int(target)]
            for row_index, target in enumerate(core_row)
            if int(target) in target_variables[row_index]
        ]
        if not clause:
            raise RuntimeError("a critical-core point has no valid phase")
        wcnf.append(clause)
        hard_clause_count += 1
    if args.max_changes is not None:
        changed_literals = [
            -target_variables[row_index][int(assignment[row_index])]
            for row_index, row in enumerate(rows)
            if int(row["p"]) not in fixed_primes
        ]
        change_encoding = CardEnc.atmost(
            changed_literals,
            bound=args.max_changes,
            vpool=vpool,
            encoding=EncType.seqcounter,
        )
        for clause in change_encoding.clauses:
            wcnf.append(clause)
        hard_clause_count += len(change_encoding.clauses)
    encode_seconds = time.monotonic() - started - build_seconds
    print(
        f"options={option_count} variables={vpool.top} "
        f"hard_clauses={hard_clause_count} soft={len(rows)-len(fixed_primes)} "
        f"encode_s={encode_seconds:.3f}",
        flush=True,
    )

    solve_started = time.monotonic()
    if args.sat_only:
        from pysat.solvers import Solver  # type: ignore

        oracle = Solver(name=args.solver, bootstrap_with=wcnf.hard)
        current_hints = [
            target_variables[row_index][int(assignment[row_index])]
            for row_index in range(len(rows))
        ]
        try:
            oracle.set_phases(current_hints)
        except NotImplementedError:
            pass
        model = oracle.get_model() if oracle.solve() else None
        oracle.delete()
        if model is not None:
            preliminary_positive = {
                literal for literal in model if literal > 0
            }
            optimum = sum(
                target_variables[row_index][int(assignment[row_index])]
                not in preliminary_positive
                for row_index in range(len(rows))
            )
        else:
            optimum = None
    else:
        with RC2(wcnf, solver=args.solver, adapt=True) as optimizer:
            model = optimizer.compute()
            optimum = optimizer.cost if model is not None else None
    solve_seconds = time.monotonic() - solve_started
    if model is None:
        result = {
            "pool": str(args.pool),
            "points": str(args.points),
            "initial_phases": str(args.initial_phases),
            "point_count": len(points),
            "core_point_count": len(core_indices),
            "initial_misses": initial_misses,
            "option_count": option_count,
            "variable_count": vpool.top,
            "hard_clause_count": hard_clause_count,
            "solver": args.solver,
            "sat_only": args.sat_only,
            "max_changes": args.max_changes,
            "result": "UNSAT",
            "build_seconds": build_seconds,
            "encode_seconds": encode_seconds,
            "solve_seconds": solve_seconds,
        }
        args.result_output.write_text(json.dumps(result, indent=2) + "\n")
        print(f"result=UNSAT solve_s={solve_seconds:.3f}", flush=True)
        return 2

    positive = {literal for literal in model if literal > 0}
    repaired = assignment.copy()
    selected_count = 0
    for row_index, variable_by_target in enumerate(target_variables):
        selected = [
            target
            for target, variable in variable_by_target.items()
            if variable in positive
        ]
        if len(selected) != 1:
            raise AssertionError("model does not select exactly one phase")
        repaired[row_index] = selected[0]
        selected_count += 1
    if selected_count != len(rows):
        raise AssertionError("incomplete repaired assignment")
    repaired_core_cover = np.count_nonzero(
        core_targets == repaired, axis=1
    )
    if np.any(repaired_core_cover == 0):
        raise AssertionError("MaxSAT model misses a critical-core point")
    repaired_full_cover = np.count_nonzero(targets == repaired, axis=1)
    miss_indices = np.flatnonzero(repaired_full_cover == 0)
    phase_map = {
        str(int(row["p"])): int(repaired[row_index])
        for row_index, row in enumerate(rows)
    }
    args.phase_output.write_text(json.dumps(phase_map) + "\n")
    result = {
        "pool": str(args.pool),
        "points": str(args.points),
        "initial_phases": str(args.initial_phases),
        "point_count": len(points),
        "core_max_coverage": args.core_max_coverage,
        "option_max_coverage": option_max_coverage,
        "option_point_count": len(option_indices),
        "core_point_count": len(core_indices),
        "initial_misses": initial_misses,
        "option_count": option_count,
        "variable_count": vpool.top,
        "hard_clause_count": hard_clause_count,
        "solver": args.solver,
        "sat_only": args.sat_only,
        "max_changes": args.max_changes,
        "result": "SAT",
        "changed_phases": int(optimum),
        "full_misses": int(len(miss_indices)),
        "misses": [
            [points[int(index)][0], points[int(index)][1]]
            for index in miss_indices
        ],
        "build_seconds": build_seconds,
        "encode_seconds": encode_seconds,
        "solve_seconds": solve_seconds,
        "elapsed_seconds": time.monotonic() - started,
    }
    args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"result=SAT changed={optimum} full_misses={len(miss_indices)} "
        f"solve_s={solve_seconds:.3f}",
        flush=True,
    )
    return 0 if not len(miss_indices) else 1


if __name__ == "__main__":
    raise SystemExit(main())
