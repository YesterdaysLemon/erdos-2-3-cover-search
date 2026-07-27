#!/usr/bin/env python3
"""Sparse MILP search for a low-change finite affine-cover repair."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

from certify_anchor_phase_quotient import load_cells
from local_phase_cegis import build_targets


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path, nargs="+")
    parser.add_argument("--initial-phases", type=Path, required=True)
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    parser.add_argument("--fixed-primes", default="")
    parser.add_argument("--core-max-coverage", type=int, default=1)
    parser.add_argument("--option-max-coverage", type=int, default=0)
    parser.add_argument("--max-changes", type=int)
    parser.add_argument(
        "--feasibility-only",
        action="store_true",
        help="use a zero objective and return any integer repair",
    )
    parser.add_argument("--time-limit", type=float, default=600.0)
    args = parser.parse_args()
    if args.core_max_coverage < 0 or args.option_max_coverage < 0:
        raise SystemExit("coverage thresholds must be nonnegative")
    if args.max_changes is not None and args.max_changes < 0:
        raise SystemExit("--max-changes must be nonnegative")
    if args.time_limit <= 0:
        raise SystemExit("--time-limit must be positive")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from scipy.optimize import Bounds, LinearConstraint, milp  # type: ignore
    from scipy.sparse import coo_matrix  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    points = []
    seen_points = set()
    for path in args.points:
        for point in load_cells(path):
            if point not in seen_points:
                seen_points.add(point)
                points.append(point)
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
    max_h = max((int(row["h"]) for row in rows), default=1)
    if max_h <= 1 << 32:
        assignment_dtype = np.uint32
    elif max_h <= 1 << 64:
        assignment_dtype = np.uint64
    else:
        assignment_dtype = object
    assignment = np.empty(len(rows), dtype=assignment_dtype)
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
    option_indices = np.flatnonzero(
        full_cover <= args.option_max_coverage
    )
    core_targets = targets[core_indices]
    option_targets = targets[option_indices]
    initial_misses = int(np.count_nonzero(full_cover == 0))
    print(
        f"rows={len(rows)} points={len(points)} core={len(core_indices)} "
        f"option_points={len(option_indices)} initial_misses={initial_misses} "
        f"matrix_s={build_seconds:.3f}",
        flush=True,
    )

    variable_targets: list[list[int]] = []
    target_variables: list[dict[int, int]] = []
    current_variables = []
    variable_count = 0
    for row_index, row in enumerate(rows):
        prime = int(row["p"])
        modulus = int(row["target_modulus"])
        residue = int(row["target_residue"])
        current = int(assignment[row_index])
        observed = {
            int(target)
            for target in np.unique(option_targets[:, row_index])
            if int(target) % modulus == residue
        }
        observed.add(current)
        if prime in fixed_primes:
            observed = {current}
        ordered = sorted(observed)
        mapping = {}
        for target in ordered:
            mapping[target] = variable_count
            variable_count += 1
        variable_targets.append(ordered)
        target_variables.append(mapping)
        current_variables.append(mapping[current])

    matrix_rows = []
    matrix_columns = []
    matrix_values = []
    lower = []
    upper = []

    # Exactly one phase per source row.
    constraint_index = 0
    for mapping in target_variables:
        for variable in mapping.values():
            matrix_rows.append(constraint_index)
            matrix_columns.append(variable)
            matrix_values.append(1.0)
        lower.append(1.0)
        upper.append(1.0)
        constraint_index += 1

    # Every critical-core point must be covered.
    for core_row in core_targets:
        clause_variables = []
        for row_index, target in enumerate(core_row):
            variable = target_variables[row_index].get(int(target))
            if variable is not None:
                clause_variables.append(variable)
        if not clause_variables:
            raise RuntimeError("a critical-core point has no offered phase")
        for variable in set(clause_variables):
            matrix_rows.append(constraint_index)
            matrix_columns.append(variable)
            matrix_values.append(1.0)
        lower.append(1.0)
        upper.append(np.inf)
        constraint_index += 1

    if args.max_changes is not None:
        for variable in current_variables:
            matrix_rows.append(constraint_index)
            matrix_columns.append(variable)
            matrix_values.append(1.0)
        lower.append(float(len(rows) - args.max_changes))
        upper.append(np.inf)
        constraint_index += 1

    matrix = coo_matrix(
        (matrix_values, (matrix_rows, matrix_columns)),
        shape=(constraint_index, variable_count),
        dtype=np.float64,
    ).tocsr()
    objective = np.zeros(variable_count, dtype=np.float64)
    if not args.feasibility_only:
        # Prefer the smallest number of changes even when a looser hard bound
        # is supplied.  Positive candidates are rechecked exactly below.
        objective[np.asarray(current_variables, dtype=np.int64)] = -1.0
    encode_seconds = time.monotonic() - started - build_seconds
    print(
        f"variables={variable_count} constraints={constraint_index} "
        f"nonzeros={matrix.nnz} encode_s={encode_seconds:.3f}",
        flush=True,
    )

    solve_started = time.monotonic()
    result = milp(
        c=objective,
        integrality=np.ones(variable_count, dtype=np.uint8),
        bounds=Bounds(0.0, 1.0),
        constraints=LinearConstraint(
            matrix,
            np.asarray(lower, dtype=np.float64),
            np.asarray(upper, dtype=np.float64),
        ),
        options={
            "time_limit": args.time_limit,
            "mip_rel_gap": 0.0,
            "presolve": True,
        },
    )
    solve_seconds = time.monotonic() - solve_started
    if result.x is None:
        output = {
            "pool": str(args.pool),
            "points": [str(path) for path in args.points],
            "initial_phases": str(args.initial_phases),
            "result": "NO_INTEGER_MODEL",
            "milp_status": int(result.status),
            "milp_message": str(result.message),
            "initial_misses": initial_misses,
            "core_point_count": len(core_indices),
            "option_point_count": len(option_indices),
            "variable_count": variable_count,
            "constraint_count": constraint_index,
            "nonzero_count": matrix.nnz,
            "build_seconds": build_seconds,
            "encode_seconds": encode_seconds,
            "solve_seconds": solve_seconds,
        }
        args.result_output.write_text(json.dumps(output, indent=2) + "\n")
        print(
            f"result=NO_INTEGER_MODEL status={result.status} "
            f"solve_s={solve_seconds:.3f}",
            flush=True,
        )
        return 2

    selected = np.rint(result.x).astype(np.int8)
    if np.max(np.abs(result.x - selected)) > 1e-6:
        raise RuntimeError("MILP returned a nonintegral incumbent")
    repaired = assignment.copy()
    for row_index, mapping in enumerate(target_variables):
        chosen = [
            target
            for target, variable in mapping.items()
            if selected[variable]
        ]
        if len(chosen) != 1:
            raise RuntimeError("MILP incumbent violates exact-one")
        repaired[row_index] = chosen[0]
    repaired_core_cover = np.count_nonzero(
        core_targets == repaired, axis=1
    )
    if np.any(repaired_core_cover == 0):
        raise RuntimeError("MILP incumbent misses a critical-core point")
    changed = int(np.count_nonzero(repaired != assignment))
    if args.max_changes is not None and changed > args.max_changes:
        raise RuntimeError("MILP incumbent violates the change bound")
    repaired_full_cover = np.count_nonzero(targets == repaired, axis=1)
    miss_indices = np.flatnonzero(repaired_full_cover == 0)
    phase_map = {
        str(int(row["p"])): int(repaired[row_index])
        for row_index, row in enumerate(rows)
    }
    args.phase_output.write_text(json.dumps(phase_map) + "\n")
    output = {
        "pool": str(args.pool),
        "points": [str(path) for path in args.points],
        "initial_phases": str(args.initial_phases),
        "result": "INTEGER_MODEL",
        "milp_status": int(result.status),
        "milp_message": str(result.message),
        "initial_misses": initial_misses,
        "core_point_count": len(core_indices),
        "option_point_count": len(option_indices),
        "max_changes": args.max_changes,
        "feasibility_only": args.feasibility_only,
        "changed_phases": changed,
        "full_misses": int(len(miss_indices)),
        "misses": [
            [points[int(index)][0], points[int(index)][1]]
            for index in miss_indices
        ],
        "variable_count": variable_count,
        "constraint_count": constraint_index,
        "nonzero_count": matrix.nnz,
        "build_seconds": build_seconds,
        "encode_seconds": encode_seconds,
        "solve_seconds": solve_seconds,
        "elapsed_seconds": time.monotonic() - started,
    }
    args.result_output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"result=INTEGER_MODEL changed={changed} "
        f"full_misses={len(miss_indices)} status={result.status} "
        f"solve_s={solve_seconds:.3f}",
        flush=True,
    )
    return 0 if not len(miss_indices) else 1


if __name__ == "__main__":
    raise SystemExit(main())
