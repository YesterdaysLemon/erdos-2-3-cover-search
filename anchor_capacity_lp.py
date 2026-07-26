#!/usr/bin/env python3
"""Fast LP relaxation of the p=5,p=7 anchor-capacity obstruction."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from fractions import Fraction
from pathlib import Path

import anchor_capacity
import cegis_cover
import local_cover
import search_cover


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("period", type=int)
    parser.add_argument("--candidate-pool", type=Path)
    parser.add_argument("--count", type=int)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from scipy.optimize import linprog  # type: ignore

    if args.candidate_pool:
        source = json.loads(args.candidate_pool.read_text())
        candidates = []
        for row in source["choices"]:
            p = int(row["p"])
            ord2 = search_cover.multiplicative_order(2, p)
            ord3 = search_cover.multiplicative_order(3, p)
            h, a, b, ord2, ord3 = cegis_cover.signature(p, ord2, ord3)
            candidates.append((h, p, a, b, ord2, ord3))
        candidates.sort()
        if args.count is not None:
            candidates = candidates[: args.count]
        common = None
        factorization = {}
    else:
        candidates, common, factorization = local_cover.get_complete_period_candidates(
            args.period
        )
    by_prime = {row[1]: row for row in candidates}
    anchors = [
        (by_prime[p][0], by_prime[p][2], by_prime[p][3])
        for p in anchor_capacity.ANCHOR_PRIMES
    ]
    anchor_period = math.lcm(*(row[0] for row in anchors))
    cells = tuple(
        sorted(
            {
                tuple((a * k + b * l) % h for h, a, b in anchors)
                for k in range(anchor_period)
                for l in range(anchor_period)
            }
        )
    )
    required = tuple(cell for cell in cells if all(value != 0 for value in cell))
    rows = []
    for h, p, a, b, ord2, ord3 in candidates:
        if p in anchor_capacity.ANCHOR_PRIMES:
            continue
        rows.append((p, anchor_capacity.option_cells(h, a, b, anchors, cells)))

    columns = [(row_no, option_no) for row_no, (_, opts) in enumerate(rows) for option_no in range(len(opts))]
    a_eq = np.zeros((len(rows), len(columns)))
    for column, (row_no, option_no) in enumerate(columns):
        a_eq[row_no, column] = 1.0
    b_eq = np.ones(len(rows))
    # Add a final variable z and maximize it subject to every required cell's
    # capacity being at least z.  This gives the strongest fractional version
    # of the capacity test and exposes dual cell weights.
    z_column = len(columns)
    a_ub = np.zeros((len(required), len(columns) + 1))
    for cell_no, cell in enumerate(required):
        for column, (row_no, option_no) in enumerate(columns):
            compatible, contribution = rows[row_no][1][option_no]
            if cell in compatible:
                a_ub[cell_no, column] = -float(contribution)
        a_ub[cell_no, z_column] = 1.0
    b_ub = np.zeros(len(required))
    a_eq_with_z = np.pad(a_eq, ((0, 0), (0, 1)))
    objective = np.zeros(len(columns) + 1)
    objective[z_column] = -1.0
    result = linprog(
        objective,
        A_ub=a_ub,
        b_ub=b_ub,
        A_eq=a_eq_with_z,
        b_eq=b_eq,
        bounds=[(0.0, 1.0)] * len(columns) + [(None, None)],
        method="highs",
    )
    print(
        f"period={args.period} eligible_primes={len(candidates)} "
        f"gcd_bits={common.bit_length() if common is not None else 'n/a'} "
        f"factors={len(factorization)}"
    )
    print(f"variables={len(columns)} required_cells={len(required)}")
    print(f"lp_status={result.status} success={result.success} message={result.message}")
    if not result.success:
        return 1

    optimum = float(result.x[z_column])
    threshold = Fraction(1, len(cells))
    # HiGHS' inequality marginals for z-capacity <= 0 are non-positive.
    # Rationalize and renormalize them, then verify the resulting separating
    # inequality using exact Fraction arithmetic rather than floating point.
    approximate_weights = [-float(value) for value in result.ineqlin.marginals]
    rational_weights = [Fraction(value).limit_denominator(1_000_000) for value in approximate_weights]
    weight_sum = sum(rational_weights, Fraction())
    rational_weights = [weight / weight_sum for weight in rational_weights]
    upper = Fraction()
    for p, options in rows:
        upper += max(
            sum(
                (
                    rational_weights[cell_no] * contribution
                    for cell_no, cell in enumerate(required)
                    if cell in compatible
                ),
                Fraction(),
            )
            for compatible, contribution in options
        )
    exact_result = "IMPOSSIBLE" if upper < threshold else "INCONCLUSIVE"
    print(f"fractional_optimum={optimum:.15f}")
    print(f"required={threshold} ({float(threshold):.15f})")
    print(f"dual_upper_exact={upper} ({float(upper):.15f})")
    print(f"dual_weights={[str(value) for value in rational_weights]}")
    print(f"result={exact_result}")
    return 2 if upper < threshold else 0


if __name__ == "__main__":
    raise SystemExit(main())
