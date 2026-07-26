#!/usr/bin/env python3
"""Exact rational two-anchor capacity obstruction for a derived line pool."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from fractions import Fraction
from pathlib import Path

import anchor_capacity


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--anchor-indices", default="0,1")
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from scipy.optimize import linprog  # type: ignore

    candidates = json.loads(args.pool.read_text())["choices"]
    anchor_indices = tuple(int(value) for value in args.anchor_indices.split(","))
    if len(anchor_indices) != 2:
        raise SystemExit("exactly two anchor indices are required")
    anchor_rows = [candidates[index] for index in anchor_indices]
    anchors = [
        (int(row["h"]), int(row["a"]), int(row["b"])) for row in anchor_rows
    ]
    period = math.lcm(*(row[0] for row in anchors))
    cells = tuple(
        sorted(
            {
                tuple((a * k + b * l) % h for h, a, b in anchors)
                for k in range(period)
                for l in range(period)
            }
        )
    )
    expected = math.prod(row[0] for row in anchors)
    if len(cells) != expected:
        raise RuntimeError("anchor maps are not jointly surjective")
    required = tuple(cell for cell in cells if all(value != 0 for value in cell))

    rows = []
    excluded = set(anchor_indices)
    for index, row in enumerate(candidates):
        if index in excluded:
            continue
        options = anchor_capacity.option_cells(
            int(row["h"]), int(row["a"]), int(row["b"]), anchors, cells
        )
        rows.append((int(row["p"]), options))
    columns = [
        (row_no, option_no)
        for row_no, (_, options) in enumerate(rows)
        for option_no in range(len(options))
    ]
    a_eq = np.zeros((len(rows), len(columns) + 1))
    for column, (row_no, option_no) in enumerate(columns):
        a_eq[row_no, column] = 1.0
    a_ub = np.zeros((len(required), len(columns) + 1))
    for cell_no, cell in enumerate(required):
        for column, (row_no, option_no) in enumerate(columns):
            compatible, contribution = rows[row_no][1][option_no]
            if cell in compatible:
                a_ub[cell_no, column] = -float(contribution)
        a_ub[cell_no, -1] = 1.0
    objective = np.zeros(len(columns) + 1)
    objective[-1] = -1.0
    result = linprog(
        objective,
        A_ub=a_ub,
        b_ub=np.zeros(len(required)),
        A_eq=a_eq,
        b_eq=np.ones(len(rows)),
        bounds=[(0.0, 1.0)] * len(columns) + [(None, None)],
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)

    weights = [
        Fraction(-float(value)).limit_denominator(1_000_000)
        for value in result.ineqlin.marginals
    ]
    weight_sum = sum(weights, Fraction())
    weights = [weight / weight_sum for weight in weights]
    exact_upper = Fraction()
    for p, options in rows:
        exact_upper += max(
            sum(
                (
                    weights[cell_no] * contribution
                    for cell_no, cell in enumerate(required)
                    if cell in compatible
                ),
                Fraction(),
            )
            for compatible, contribution in options
        )
    threshold = Fraction(1, len(cells))
    impossible = exact_upper < threshold
    print(
        f"pool={args.pool} candidates={len(candidates)} "
        f"anchors={[(row['p'], row['h']) for row in anchor_rows]}"
    )
    print(
        f"anchor_cells={len(cells)} required_cells={len(required)} "
        f"variables={len(columns)}"
    )
    print(f"fractional_optimum={result.x[-1]:.15f}")
    print(
        f"exact_upper_float={float(exact_upper):.15f} "
        f"exact_upper_bits=({exact_upper.numerator.bit_length()},"
        f"{exact_upper.denominator.bit_length()})"
    )
    print(f"threshold={threshold} ({float(threshold):.15f})")
    print(f"result={'IMPOSSIBLE' if impossible else 'INCONCLUSIVE'}")
    return 2 if impossible else 0


if __name__ == "__main__":
    raise SystemExit(main())
