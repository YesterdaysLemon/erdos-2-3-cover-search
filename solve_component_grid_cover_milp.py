#!/usr/bin/env python3
"""Sparse HiGHS MILP for an exact affine-line cover of F_q^2."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path


def parse_targets(text: str) -> dict[int, int]:
    result = {}
    for item in text.split(","):
        if not item:
            continue
        row_prime, target = item.split(":", 1)
        result[int(row_prime)] = int(target)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--fixed-targets", default="")
    parser.add_argument("--initial-phases", type=Path, required=True)
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_matrix

    rows = [
        row
        for row in json.loads(args.pool.read_text())["choices"]
        if int(row["h"]) == args.prime
    ]
    fixed = parse_targets(args.fixed_targets)
    row_primes = {int(row["p"]) for row in rows}
    if fixed.keys() - row_primes:
        raise RuntimeError("a fixed target is not a component row")
    row_by_prime = {int(row["p"]): row for row in rows}
    for prime, target in fixed.items():
        row = row_by_prime[prime]
        if target % int(row["target_modulus"]) != int(
            row["target_residue"]
        ):
            raise RuntimeError(
                f"fixed target violates row restriction p={prime}"
            )
        fixed[prime] = target % args.prime
    initial = {
        int(row_prime): int(target)
        for row_prime, target in json.loads(
            args.initial_phases.read_text()
        ).items()
    }
    mutable = [
        row for row in rows if int(row["p"]) not in fixed
    ]
    fixed_rows = [
        row for row in rows if int(row["p"]) in fixed
    ]
    variable_count = len(mutable) * args.prime

    def variable(row_index: int, target: int) -> int:
        return row_index * args.prime + target

    objective = np.ones(variable_count, dtype=float)
    variable_upper = np.zeros(variable_count, dtype=float)
    for row_index, row in enumerate(mutable):
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        for target in range(args.prime):
            if target % modulus == residue:
                variable_upper[variable(row_index, target)] = 1.0
        initial_target = initial[int(row["p"])] % args.prime
        if initial_target % modulus == residue:
            objective[variable(row_index, initial_target)] = 0.0

    matrix_rows = []
    matrix_columns = []
    matrix_values = []
    lower = []
    upper = []
    constraint_index = 0
    for row_index, _row in enumerate(mutable):
        for target in range(args.prime):
            matrix_rows.append(constraint_index)
            matrix_columns.append(variable(row_index, target))
            matrix_values.append(1.0)
        lower.append(1.0)
        upper.append(1.0)
        constraint_index += 1

    uncovered_by_fixed = []
    for k in range(args.prime):
        for l in range(args.prime):
            if any(
                (
                    int(row["a"]) * k
                    + int(row["b"]) * l
                    - fixed[int(row["p"])]
                )
                % args.prime
                == 0
                for row in fixed_rows
            ):
                continue
            for row_index, row in enumerate(mutable):
                target = (
                    int(row["a"]) * k
                    + int(row["b"]) * l
                ) % args.prime
                matrix_rows.append(constraint_index)
                matrix_columns.append(variable(row_index, target))
                matrix_values.append(1.0)
            lower.append(1.0)
            upper.append(math.inf)
            constraint_index += 1
            uncovered_by_fixed.append((k, l))

    matrix = coo_matrix(
        (
            np.asarray(matrix_values, dtype=float),
            (
                np.asarray(matrix_rows, dtype=np.int32),
                np.asarray(matrix_columns, dtype=np.int32),
            ),
        ),
        shape=(constraint_index, variable_count),
    ).tocsc()
    print(
        f"prime={args.prime} rows={len(rows)} fixed={len(fixed_rows)} "
        f"mutable={len(mutable)} cells={len(uncovered_by_fixed)} "
        f"variables={variable_count} constraints={constraint_index} "
        f"nonzeros={matrix.nnz}",
        flush=True,
    )
    started = time.monotonic()
    result = milp(
        c=objective,
        integrality=np.ones(variable_count, dtype=np.uint8),
        bounds=Bounds(
            np.zeros(variable_count),
            variable_upper,
        ),
        constraints=LinearConstraint(
            matrix,
            np.asarray(lower),
            np.asarray(upper),
        ),
        options={
            "time_limit": args.time_limit,
            "mip_rel_gap": 0.0,
            "presolve": True,
        },
    )
    elapsed = time.monotonic() - started
    phases = dict(fixed)
    if result.x is not None:
        for row_index, row in enumerate(mutable):
            selected = [
                target
                for target in range(args.prime)
                if result.x[variable(row_index, target)] > 0.5
            ]
            if len(selected) != 1:
                raise AssertionError("MILP row phase is not integral")
            phases[int(row["p"])] = selected[0]

    misses = []
    if result.x is not None:
        for k in range(args.prime):
            for l in range(args.prime):
                if not any(
                    (
                        int(row["a"]) * k
                        + int(row["b"]) * l
                        - phases[int(row["p"])]
                    )
                    % args.prime
                    == 0
                    for row in rows
                ):
                    misses.append((k, l))
        if misses:
            raise AssertionError("decoded MILP assignment leaves a cell")
        args.phase_output.write_text(
            json.dumps(
                {
                    str(row_prime): target
                    for row_prime, target in sorted(phases.items())
                }
            )
            + "\n"
        )

    output = {
        "pool": str(args.pool),
        "prime": args.prime,
        "status": int(result.status),
        "success": bool(result.success),
        "message": str(result.message),
        "has_incumbent": result.x is not None,
        "objective": (
            float(result.fun) if result.fun is not None else None
        ),
        "fixed_rows": len(fixed_rows),
        "mutable_rows": len(mutable),
        "target_restrictions_enforced": True,
        "cells_not_covered_by_fixed": len(uncovered_by_fixed),
        "variables": variable_count,
        "constraints": constraint_index,
        "nonzeros": int(matrix.nnz),
        "verified_misses": len(misses) if result.x is not None else None,
        "elapsed_seconds": elapsed,
        "phase_file": (
            str(args.phase_output) if result.x is not None else None
        ),
    }
    args.result_output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"status={result.status} incumbent={result.x is not None} "
        f"objective={output['objective']} elapsed_s={elapsed:.3f}",
        flush=True,
    )
    if result.status == 2:
        return 2
    return 0 if result.x is not None else 4


if __name__ == "__main__":
    raise SystemExit(main())
