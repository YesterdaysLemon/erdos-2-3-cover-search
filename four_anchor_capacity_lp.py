#!/usr/bin/env python3
"""Four-anchor exact-dual obstruction with a coprime fourth modulus."""

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
    parser.add_argument("--anchor-indices", default="0,1,2,6")
    parser.add_argument("--weight-denominator", type=int, default=100000)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from scipy.optimize import linprog  # type: ignore
    from scipy.sparse import coo_matrix  # type: ignore

    candidates = json.loads(args.pool.read_text())["choices"]
    indices = tuple(int(value) for value in args.anchor_indices.split(","))
    if len(indices) != 4:
        raise SystemExit("exactly four anchor indices are required")
    anchor_rows = [candidates[index] for index in indices]
    first = anchor_rows[:3]
    fourth = anchor_rows[3]
    q = int(fourth["h"])
    if any(math.gcd(q, int(row["h"])) != 1 for row in first):
        raise RuntimeError("fourth anchor modulus must be coprime to the first three")
    if len(__import__("exact_uncovered").factor(q)) != 1:
        raise RuntimeError("fourth anchor modulus must be prime")
    base_anchors = [
        (int(row["h"]), int(row["a"]), int(row["b"])) for row in first
    ]
    base_period = math.lcm(*(item[0] for item in base_anchors))
    base_cells = tuple(
        sorted(
            {
                tuple((a * k + b * l) % h for h, a, b in base_anchors)
                for k in range(base_period)
                for l in range(base_period)
            }
        )
    )
    if len(base_cells) != math.prod(item[0] for item in base_anchors):
        raise RuntimeError("first three anchors are not jointly surjective")
    required_base = [cell for cell in base_cells if all(value != 0 for value in cell)]
    required = tuple(
        (base_cell, fourth_value)
        for base_cell in required_base
        for fourth_value in range(1, q)
    )
    required_index = {item: index for index, item in enumerate(required)}

    fourth_direction = (int(fourth["a"]) % q, int(fourth["b"]) % q)
    excluded = set(indices)
    rows = []
    for index, row in enumerate(candidates):
        if index in excluded:
            continue
        h = int(row["h"])
        a = int(row["a"])
        b = int(row["b"])
        base_options = anchor_capacity.option_cells(
            h, a, b, base_anchors, base_cells
        )
        parallel = (
            h % q == 0
            and (a * fourth_direction[1] - b * fourth_direction[0]) % q == 0
        )
        options = []
        for base_compatible, _old_contribution in base_options:
            if parallel:
                fourth_options = ((value,) for value in range(q))
            else:
                fourth_options = (tuple(range(q)),)
            for fourth_compatible in fourth_options:
                compatible = tuple(
                    (base_cell, value)
                    for base_cell in base_compatible
                    for value in fourth_compatible
                )
                options.append(
                    (compatible, Fraction(1, h * len(compatible)))
                )
        rows.append((int(row["p"]), options))

    columns = [
        (row_index, option_index)
        for row_index, (_p, options) in enumerate(rows)
        for option_index in range(len(options))
    ]
    eq = coo_matrix(
        (
            np.ones(len(columns)),
            (
                [row_index for row_index, _option_index in columns],
                range(len(columns)),
            ),
        ),
        shape=(len(rows), len(columns) + 1),
    ).tocsr()
    ub_rows = []
    ub_columns = []
    ub_values = []
    for column, (row_index, option_index) in enumerate(columns):
        compatible, contribution = rows[row_index][1][option_index]
        for item in compatible:
            cell_index = required_index.get(item)
            if cell_index is not None:
                ub_rows.append(cell_index)
                ub_columns.append(column)
                ub_values.append(-float(contribution))
    for cell_index in range(len(required)):
        ub_rows.append(cell_index)
        ub_columns.append(len(columns))
        ub_values.append(1.0)
    ub = coo_matrix(
        (ub_values, (ub_rows, ub_columns)),
        shape=(len(required), len(columns) + 1),
    ).tocsr()
    objective = np.zeros(len(columns) + 1)
    objective[-1] = -1.0
    result = linprog(
        objective,
        A_ub=ub,
        b_ub=np.zeros(len(required)),
        A_eq=eq,
        b_eq=np.ones(len(rows)),
        bounds=[(0.0, 1.0)] * len(columns) + [(None, None)],
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    threshold = Fraction(1, len(base_cells) * q)
    weights = [
        Fraction(-float(value)).limit_denominator(args.weight_denominator)
        for value in result.ineqlin.marginals
    ]
    weight_sum = sum(weights, Fraction())
    weights = [weight / weight_sum for weight in weights]
    weight_by_item = dict(zip(required, weights))
    exact_upper = Fraction()
    for _p, options in rows:
        exact_upper += max(
            sum(
                (
                    weight_by_item[item] * contribution
                    for item in compatible
                    if item in weight_by_item
                ),
                Fraction(),
            )
            for compatible, contribution in options
        )
    impossible = exact_upper < threshold
    print(
        f"pool={args.pool} candidates={len(candidates)} "
        f"anchors={[(row['p'], row['h']) for row in anchor_rows]} "
        f"cells={len(base_cells) * q} required={len(required)} "
        f"columns={len(columns)}"
    )
    print(
        f"lp={result.x[-1]:.15f} exact_upper={float(exact_upper):.15f} "
        f"threshold={float(threshold):.15f} "
        f"bits=({exact_upper.numerator.bit_length()},"
        f"{exact_upper.denominator.bit_length()})"
    )
    print(f"result={'IMPOSSIBLE' if impossible else 'INCONCLUSIVE'}")
    return 2 if impossible else 0


if __name__ == "__main__":
    raise SystemExit(main())
