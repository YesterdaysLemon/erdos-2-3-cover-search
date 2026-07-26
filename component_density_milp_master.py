#!/usr/bin/env python3
"""Sparse HiGHS master for accumulated component-density CEGIS cells."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import exact_uncovered
from component_density_cegis import q_part


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("cells", type=Path)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--initial-phases", type=Path, required=True)
    parser.add_argument("--fixed-targets", default="")
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import coo_matrix

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    if any(int(row["target_modulus"]) != 1 for row in rows):
        raise RuntimeError("MILP master requires target_modulus=1")
    phases = {
        int(row_prime): int(target)
        for row_prime, target in json.loads(
            args.initial_phases.read_text()
        ).items()
    }
    fixed_targets = {}
    for item in args.fixed_targets.split(","):
        if not item:
            continue
        row_prime, target = item.split(":", 1)
        fixed_targets[int(row_prime)] = int(target)
    phases.update(fixed_targets)
    cells = [
        (int(k), int(l))
        for k, l in json.loads(args.cells.read_text())
    ]
    if len(cells) != len(set(cells)):
        raise RuntimeError("cell file contains duplicates")

    row_data = []
    residual_modulus = 1
    for row in rows:
        h = int(row["h"])
        residual = q_part(h, args.prime)
        residual_modulus = max(residual_modulus, residual)
        row_data.append((h // residual, residual))

    options: dict[int, set[int]] = defaultdict(set)
    required_by_cell = []
    for k, l in cells:
        required = []
        for row_index, (row, (other, _residual)) in enumerate(
            zip(rows, row_data)
        ):
            target = (
                int(row["a"]) * k + int(row["b"]) * l
            ) % other
            required.append(target)
            if other > 1 and int(row["p"]) not in fixed_targets:
                options[row_index].add(target)
        required_by_cell.append(required)

    option_list = [
        (row_index, target)
        for row_index, targets in sorted(options.items())
        for target in sorted(targets)
    ]
    option_index = {
        option: index for index, option in enumerate(option_list)
    }
    variable_count = len(option_list)
    objective = np.ones(variable_count, dtype=float)
    for index, (row_index, target) in enumerate(option_list):
        row_prime = int(rows[row_index]["p"])
        other, _residual = row_data[row_index]
        if phases[row_prime] % other == target:
            objective[index] = 0.0

    matrix_rows = []
    matrix_columns = []
    matrix_values = []
    lower = []
    upper = []
    constraint_index = 0

    for row_index, targets in sorted(options.items()):
        for target in targets:
            matrix_rows.append(constraint_index)
            matrix_columns.append(option_index[(row_index, target)])
            matrix_values.append(1.0)
        lower.append(-math.inf)
        upper.append(1.0)
        constraint_index += 1

    cut_count = 0
    for required in required_by_cell:
        constant = 0
        entries = []
        for row_index, (row, (other, residual)) in enumerate(
            zip(rows, row_data)
        ):
            weight = residual_modulus // residual
            row_prime = int(row["p"])
            if other == 1:
                constant += weight
            elif row_prime in fixed_targets:
                if fixed_targets[row_prime] % other == required[row_index]:
                    constant += weight
            else:
                entries.append(
                    (
                        option_index[
                            (row_index, required[row_index])
                        ],
                        float(weight),
                    )
                )
        bound = residual_modulus - constant
        if bound <= 0:
            continue
        if sum(value for _index, value in entries) < bound:
            result = {
                "sat": False,
                "reason": "one cut lacks sufficient compatible weight",
                "cells": len(cells),
            }
            args.result_output.write_text(
                json.dumps(result, indent=2) + "\n"
            )
            print("UNSAT deficient cut", flush=True)
            return 2
        for column, value in entries:
            matrix_rows.append(constraint_index)
            matrix_columns.append(column)
            matrix_values.append(value)
        lower.append(float(bound))
        upper.append(math.inf)
        constraint_index += 1
        cut_count += 1

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
        f"rows={len(rows)} cells={len(cells)} cuts={cut_count} "
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
            np.ones(variable_count),
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
    decoded_phases = None
    if result.x is not None:
        decoded_phases = dict(phases)
        selected_by_row = {}
        for value, (row_index, target) in zip(
            result.x, option_list
        ):
            if value <= 0.5:
                continue
            if row_index in selected_by_row:
                raise AssertionError("MILP selected two row targets")
            selected_by_row[row_index] = target
        for row_index, target in selected_by_row.items():
            row = rows[row_index]
            row_prime = int(row["p"])
            h = int(row["h"])
            other, residual = row_data[row_index]
            if residual == 1:
                decoded_phases[row_prime] = target
            else:
                decoded_phases[row_prime] = exact_uncovered.crt(
                    [
                        (target, other),
                        (
                            decoded_phases[row_prime] % residual,
                            residual,
                        ),
                    ]
                ) % h
        decoded_phases.update(fixed_targets)

        for k, l in cells:
            scaled_density = 0
            for row, (other, residual) in zip(rows, row_data):
                if (
                    int(row["a"]) * k
                    + int(row["b"]) * l
                    - decoded_phases[int(row["p"])]
                ) % other == 0:
                    scaled_density += residual_modulus // residual
            if scaled_density < residual_modulus:
                raise AssertionError("decoded MILP phase violates a cut")
        args.phase_output.write_text(
            json.dumps(
                {
                    str(row_prime): target
                    for row_prime, target in decoded_phases.items()
                }
            )
            + "\n"
        )

    output = {
        "pool": str(args.pool),
        "cells": str(args.cells),
        "prime": args.prime,
        "fixed_targets": {
            str(row_prime): target
            for row_prime, target in fixed_targets.items()
        },
        "variables": variable_count,
        "constraints": constraint_index,
        "nonzeros": int(matrix.nnz),
        "status": int(result.status),
        "success": bool(result.success),
        "message": str(result.message),
        "objective": (
            float(result.fun) if result.fun is not None else None
        ),
        "has_incumbent": decoded_phases is not None,
        "elapsed_seconds": elapsed,
        "phase_file": (
            str(args.phase_output) if decoded_phases is not None else None
        ),
    }
    args.result_output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"status={result.status} incumbent={decoded_phases is not None} "
        f"objective={output['objective']} elapsed_s={elapsed:.3f}",
        flush=True,
    )
    if result.status == 2:
        return 2
    return 0 if decoded_phases is not None else 4


if __name__ == "__main__":
    raise SystemExit(main())
