#!/usr/bin/env python3
"""Sparse dual LP for three base anchors and three prime anchors.

The first and third extra coordinates may be correlated.  The middle extra
coordinate is conditionally independent given the base cell.  This is a
genuine joint dual measure, and is sufficient whenever no candidate row is
simultaneously restricted in the middle coordinate and another extra
coordinate.  The final result is checked with exact rational arithmetic.
"""

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
    parser.add_argument("--anchor-indices", default="0,1,2,6,11,18")
    parser.add_argument("--weight-denominator", type=int, default=100000)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from scipy.optimize import linprog  # type: ignore
    from scipy.sparse import coo_matrix  # type: ignore

    candidates = json.loads(args.pool.read_text())["choices"]
    indices = tuple(int(value) for value in args.anchor_indices.split(","))
    if len(indices) != 6:
        raise SystemExit("exactly six anchor indices are required")
    anchor_rows = [candidates[index] for index in indices]
    first = anchor_rows[:3]
    extra = anchor_rows[3:]
    q1, q2, q3 = (int(row["h"]) for row in extra)
    if any(len(exact_uncovered.factor(q)) != 1 for q in (q1, q2, q3)):
        raise RuntimeError("extra anchor moduli must be prime")
    if any(
        math.gcd(left, right) != 1
        for pos, left in enumerate((q1, q2, q3))
        for right in (q1, q2, q3)[pos + 1 :]
    ) or any(
        math.gcd(q, int(row["h"])) != 1
        for q in (q1, q2, q3)
        for row in first
    ):
        raise RuntimeError("extra anchors must be pairwise coprime to the base")

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
        raise RuntimeError("base anchors are not jointly surjective")
    required_base = tuple(cell for cell in base_cells if all(value != 0 for value in cell))
    base_index = {cell: index for index, cell in enumerate(required_base)}
    directions = [
        (int(row["a"]) % q, int(row["b"]) % q)
        for row, q in zip(extra, (q1, q2, q3))
    ]

    excluded = set(indices)
    rows = []
    restriction_counts: dict[tuple[bool, bool, bool], int] = {}
    for index, row in enumerate(candidates):
        if index in excluded:
            continue
        h = int(row["h"])
        a = int(row["a"])
        b = int(row["b"])
        options = anchor_capacity.option_cells(h, a, b, base_anchors, base_cells)
        compact = [
            (
                tuple(base_index[cell] for cell in compatible if cell in base_index),
                contribution,
            )
            for compatible, contribution in options
        ]
        restricted = tuple(
            h % q == 0 and (a * direction[1] - b * direction[0]) % q == 0
            for q, direction in zip((q1, q2, q3), directions)
        )
        if restricted[1] and (restricted[0] or restricted[2]):
            raise RuntimeError(
                "sparse factorization is invalid: a row restricts the middle "
                "extra coordinate and another extra coordinate"
            )
        restriction_counts[restricted] = restriction_counts.get(restricted, 0) + 1
        rows.append((int(row["p"]), restricted, compact))

    base_count = len(required_base)
    d1_count, d2_count, d3_count = q1 - 1, q2 - 1, q3 - 1
    x_count = base_count * d1_count * d3_count
    u1_start = x_count
    u1_count = base_count * d1_count
    u2_start = u1_start + u1_count
    u2_count = base_count * d2_count
    u3_start = u2_start + u2_count
    u3_count = base_count * d3_count
    w_start = u3_start + u3_count
    z_start = w_start + base_count
    variable_count = z_start + len(rows)
    scale = len(base_cells) * q1 * q2 * q3

    def x_index(base: int, d1: int, d3: int) -> int:
        return (base * d1_count + d1) * d3_count + d3

    def u1_index(base: int, d1: int) -> int:
        return u1_start + base * d1_count + d1

    def u2_index(base: int, d2: int) -> int:
        return u2_start + base * d2_count + d2

    def u3_index(base: int, d3: int) -> int:
        return u3_start + base * d3_count + d3

    eq_rows: list[int] = []
    eq_columns: list[int] = []
    eq_values: list[float] = []
    equality = 0
    for base in range(base_count):
        for d1 in range(d1_count):
            eq_rows.append(equality)
            eq_columns.append(u1_index(base, d1))
            eq_values.append(1.0)
            for d3 in range(d3_count):
                eq_rows.append(equality)
                eq_columns.append(x_index(base, d1, d3))
                eq_values.append(-1.0)
            equality += 1
    for base in range(base_count):
        for d3 in range(d3_count):
            eq_rows.append(equality)
            eq_columns.append(u3_index(base, d3))
            eq_values.append(1.0)
            for d1 in range(d1_count):
                eq_rows.append(equality)
                eq_columns.append(x_index(base, d1, d3))
                eq_values.append(-1.0)
            equality += 1
    for base in range(base_count):
        eq_rows.append(equality)
        eq_columns.append(w_start + base)
        eq_values.append(1.0)
        for d1 in range(d1_count):
            eq_rows.append(equality)
            eq_columns.append(u1_index(base, d1))
            eq_values.append(-1.0)
        equality += 1
        eq_rows.append(equality)
        eq_columns.append(w_start + base)
        eq_values.append(1.0)
        for d2 in range(d2_count):
            eq_rows.append(equality)
            eq_columns.append(u2_index(base, d2))
            eq_values.append(-1.0)
        equality += 1
    for base in range(base_count):
        eq_rows.append(equality)
        eq_columns.append(w_start + base)
        eq_values.append(1.0)
    equality += 1
    eq = coo_matrix(
        (eq_values, (eq_rows, eq_columns)), shape=(equality, variable_count)
    ).tocsr()
    b_eq = np.zeros(equality)
    b_eq[-1] = 1.0

    ub_rows: list[int] = []
    ub_columns: list[int] = []
    ub_values: list[float] = []
    constraint = 0

    def append_constraint(columns: list[int], coefficient: Fraction, z: int) -> None:
        nonlocal constraint
        for column in columns:
            ub_rows.append(constraint)
            ub_columns.append(column)
            ub_values.append(float(coefficient))
        ub_rows.append(constraint)
        ub_columns.append(z)
        ub_values.append(-1.0)
        constraint += 1

    for row_index, (_p, restricted, options) in enumerate(rows):
        r1, r2, r3 = restricted
        z = z_start + row_index
        for bases, contribution in options:
            if r1 and r3:
                for d1 in range(d1_count):
                    for d3 in range(d3_count):
                        append_constraint(
                            [x_index(base, d1, d3) for base in bases],
                            contribution * scale / q2,
                            z,
                        )
            elif r1:
                for d1 in range(d1_count):
                    append_constraint(
                        [u1_index(base, d1) for base in bases],
                        contribution * scale / (q2 * q3),
                        z,
                    )
            elif r2:
                for d2 in range(d2_count):
                    append_constraint(
                        [u2_index(base, d2) for base in bases],
                        contribution * scale / (q1 * q3),
                        z,
                    )
            elif r3:
                for d3 in range(d3_count):
                    append_constraint(
                        [u3_index(base, d3) for base in bases],
                        contribution * scale / (q1 * q2),
                        z,
                    )
            else:
                append_constraint(
                    [w_start + base for base in bases],
                    contribution * scale / (q1 * q2 * q3),
                    z,
                )
    ub = coo_matrix(
        (ub_values, (ub_rows, ub_columns)), shape=(constraint, variable_count)
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

    rational_x = [
        Fraction(float(value)).limit_denominator(args.weight_denominator)
        if value > 1e-14
        else Fraction()
        for value in result.x[:x_count]
    ]
    weight_sum = sum(rational_x, Fraction())
    rational_x = [value / weight_sum for value in rational_x]
    rational_u1 = [
        sum(rational_x[x_index(base, d1, d3)] for d3 in range(d3_count))
        for base in range(base_count)
        for d1 in range(d1_count)
    ]
    rational_u3 = [
        sum(rational_x[x_index(base, d1, d3)] for d1 in range(d1_count))
        for base in range(base_count)
        for d3 in range(d3_count)
    ]
    rational_w = [
        sum(rational_u1[base * d1_count + d1] for d1 in range(d1_count))
        for base in range(base_count)
    ]
    rational_u2: list[Fraction] = []
    for base in range(base_count):
        raw = [max(0.0, float(result.x[u2_index(base, d2)])) for d2 in range(d2_count)]
        raw_sum = sum(raw)
        if not rational_w[base] or raw_sum <= 1e-14:
            rational_u2.extend([Fraction()] * d2_count)
            continue
        conditional = [
            Fraction(value / raw_sum).limit_denominator(args.weight_denominator)
            if value > 1e-14
            else Fraction()
            for value in raw
        ]
        conditional_sum = sum(conditional, Fraction())
        rational_u2.extend(
            rational_w[base] * value / conditional_sum for value in conditional
        )

    exact_upper = Fraction()
    for _p, restricted, options in rows:
        r1, r2, r3 = restricted
        scores = [Fraction()]
        for bases, contribution in options:
            if r1 and r3:
                for d1 in range(d1_count):
                    for d3 in range(d3_count):
                        scores.append(
                            contribution
                            * scale
                            / q2
                            * sum(rational_x[x_index(base, d1, d3)] for base in bases)
                        )
            elif r1:
                for d1 in range(d1_count):
                    scores.append(
                        contribution
                        * scale
                        / (q2 * q3)
                        * sum(rational_u1[base * d1_count + d1] for base in bases)
                    )
            elif r2:
                for d2 in range(d2_count):
                    scores.append(
                        contribution
                        * scale
                        / (q1 * q3)
                        * sum(rational_u2[base * d2_count + d2] for base in bases)
                    )
            elif r3:
                for d3 in range(d3_count):
                    scores.append(
                        contribution
                        * scale
                        / (q1 * q2)
                        * sum(rational_u3[base * d3_count + d3] for base in bases)
                    )
            else:
                scores.append(
                    contribution
                    * scale
                    / (q1 * q2 * q3)
                    * sum(rational_w[base] for base in bases)
                )
        exact_upper += max(scores)

    impossible = exact_upper < 1
    print(
        f"pool={args.pool} candidates={len(candidates)} "
        f"anchors={[(row['p'], row['h']) for row in anchor_rows]} "
        f"cells={len(base_cells) * q1 * q2 * q3} "
        f"required={base_count * d1_count * d2_count * d3_count} "
        f"variables={variable_count} constraints={constraint}"
    )
    print(f"restriction_counts={restriction_counts}")
    print(
        f"lp_ratio={result.fun:.15f} exact_ratio_upper={float(exact_upper):.15f} "
        f"nonzero_pair_weights={sum(bool(value) for value in rational_x)} "
        f"bits=({exact_upper.numerator.bit_length()},"
        f"{exact_upper.denominator.bit_length()})"
    )
    print(f"result={'IMPOSSIBLE' if impossible else 'INCONCLUSIVE'}")
    return 2 if impossible else 0


if __name__ == "__main__":
    raise SystemExit(main())
