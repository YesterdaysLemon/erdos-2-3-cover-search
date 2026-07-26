#!/usr/bin/env python3
"""Memory-efficient dual LP for a coprime fourth anchor."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from fractions import Fraction
from pathlib import Path

import anchor_capacity
import exact_uncovered


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--anchor-indices", required=True)
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
    if len(exact_uncovered.factor(q)) != 1:
        raise RuntimeError("fourth anchor modulus must be prime")
    if any(math.gcd(q, int(row["h"])) != 1 for row in first):
        raise RuntimeError("fourth anchor must be coprime to the first three")
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
    required_base = tuple(cell for cell in base_cells if all(value != 0 for value in cell))
    base_index = {cell: index for index, cell in enumerate(required_base)}
    fourth_direction = (int(fourth["a"]) % q, int(fourth["b"]) % q)

    excluded = set(indices)
    rows = []
    for index, row in enumerate(candidates):
        if index in excluded:
            continue
        h = int(row["h"])
        a = int(row["a"])
        b = int(row["b"])
        options = anchor_capacity.option_cells(
            h, a, b, base_anchors, base_cells
        )
        compact_options = [
            (
                tuple(
                    base_index[cell]
                    for cell in compatible
                    if cell in base_index
                ),
                contribution,
            )
            for compatible, contribution in options
        ]
        parallel = (
            h % q == 0
            and (a * fourth_direction[1] - b * fourth_direction[0]) % q == 0
        )
        rows.append((int(row["p"]), parallel, compact_options))

    base_count = len(required_base)
    scale = len(base_cells) * q
    y_count = base_count * (q - 1)
    w_start = y_count
    z_start = w_start + base_count
    variable_count = z_start + len(rows)

    # Equality rows define w_base=sum_d y_(base,d), followed by sum w=1.
    eq_rows = []
    eq_columns = []
    eq_values = []
    for base in range(base_count):
        eq_rows.append(base)
        eq_columns.append(w_start + base)
        eq_values.append(1.0)
        for d_index in range(q - 1):
            eq_rows.append(base)
            eq_columns.append(base * (q - 1) + d_index)
            eq_values.append(-1.0)
    for base in range(base_count):
        eq_rows.append(base_count)
        eq_columns.append(w_start + base)
        eq_values.append(1.0)
    eq = coo_matrix(
        (eq_values, (eq_rows, eq_columns)),
        shape=(base_count + 1, variable_count),
    ).tocsr()
    b_eq = np.zeros(base_count + 1)
    b_eq[-1] = 1.0

    # For every row option, weighted contribution <= z_row.  Nonparallel
    # rows use aggregate w variables; parallel rows select one fourth value.
    ub_rows = []
    ub_columns = []
    ub_values = []
    constraint = 0
    for row_index, (_p, parallel, options) in enumerate(rows):
        if parallel:
            for bases, contribution in options:
                # fourth value zero gives score zero and is automatically
                # dominated; required fourth values are 1,...,q-1.
                for d_index in range(q - 1):
                    for base in bases:
                        ub_rows.append(constraint)
                        ub_columns.append(base * (q - 1) + d_index)
                        ub_values.append(float(contribution * scale))
                    ub_rows.append(constraint)
                    ub_columns.append(z_start + row_index)
                    ub_values.append(-1.0)
                    constraint += 1
        else:
            for bases, contribution in options:
                for base in bases:
                    ub_rows.append(constraint)
                    ub_columns.append(w_start + base)
                    ub_values.append(float(contribution * scale / q))
                ub_rows.append(constraint)
                ub_columns.append(z_start + row_index)
                ub_values.append(-1.0)
                constraint += 1
    ub = coo_matrix(
        (ub_values, (ub_rows, ub_columns)),
        shape=(constraint, variable_count),
    ).tocsr()
    objective = np.zeros(variable_count)
    objective[z_start:] = 1.0
    result = linprog(
        objective,
        A_ub=ub,
        b_ub=np.zeros(constraint),
        A_eq=eq,
        b_eq=b_eq,
        bounds=[(0.0, None)] * variable_count,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)

    rational_y = [
        Fraction(float(value)).limit_denominator(args.weight_denominator)
        for value in result.x[:y_count]
    ]
    weight_sum = sum(rational_y, Fraction())
    rational_y = [value / weight_sum for value in rational_y]
    rational_w = [
        sum(
            rational_y[base * (q - 1) + d_index]
            for d_index in range(q - 1)
        )
        for base in range(base_count)
    ]
    exact_upper = Fraction()
    for _p, parallel, options in rows:
        if parallel:
            scores = [Fraction()]
            for bases, contribution in options:
                for d_index in range(q - 1):
                    scores.append(
                        contribution
                        * scale
                        * sum(
                            rational_y[base * (q - 1) + d_index]
                            for base in bases
                        )
                    )
        else:
            scores = [
                contribution
                * scale
                * sum(rational_w[base] for base in bases)
                / q
                for bases, contribution in options
            ]
        exact_upper += max(scores)
    impossible = exact_upper < 1
    print(
        f"pool={args.pool} candidates={len(candidates)} "
        f"anchors={[(row['p'], row['h']) for row in anchor_rows]} "
        f"cells={len(base_cells) * q} required={base_count * (q - 1)} "
        f"dual_variables={variable_count} constraints={constraint}"
    )
    print(
        f"lp_ratio={result.fun:.15f} exact_ratio_upper={float(exact_upper):.15f} "
        f"bits=({exact_upper.numerator.bit_length()},"
        f"{exact_upper.denominator.bit_length()})"
    )
    print(f"result={'IMPOSSIBLE' if impossible else 'INCONCLUSIVE'}")
    return 2 if impossible else 0


if __name__ == "__main__":
    raise SystemExit(main())
