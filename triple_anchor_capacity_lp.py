#!/usr/bin/env python3
"""LP capacity obstruction using p=5, p=7, and p=11 as anchors."""

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
import search_cover


ANCHORS = (5, 7, 11)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("period", type=int)
    parser.add_argument("candidate_pool", type=Path)
    parser.add_argument("--count", type=int)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from scipy.optimize import linprog  # type: ignore

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
    by_prime = {row[1]: row for row in candidates}
    if not all(p in by_prime for p in ANCHORS):
        raise RuntimeError("candidate pool is missing an anchor")
    anchors = [(by_prime[p][0], by_prime[p][2], by_prime[p][3]) for p in ANCHORS]
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
    rows = []
    for h, p, a, b, ord2, ord3 in candidates:
        if p in ANCHORS:
            continue
        rows.append((p, anchor_capacity.option_cells(h, a, b, anchors, cells)))
    columns = [
        (row_no, option_no)
        for row_no, (_, options) in enumerate(rows)
        for option_no in range(len(options))
    ]
    a_eq = np.zeros((len(rows), len(columns) + 1))
    for column, (row_no, option_no) in enumerate(columns):
        a_eq[row_no, column] = 1.0
    b_eq = np.ones(len(rows))
    objective = np.zeros(len(columns) + 1)
    objective[-1] = -1.0

    all_impossible = True
    for third_target in range(by_prime[11][0]):
        required = tuple(
            cell
            for cell in cells
            if cell[0] != 0 and cell[1] != 0 and cell[2] != third_target
        )
        a_ub = np.zeros((len(required), len(columns) + 1))
        for cell_no, cell in enumerate(required):
            for column, (row_no, option_no) in enumerate(columns):
                compatible, contribution = rows[row_no][1][option_no]
                if cell in compatible:
                    a_ub[cell_no, column] = -float(contribution)
            a_ub[cell_no, -1] = 1.0
        result = linprog(
            objective,
            A_ub=a_ub,
            b_ub=np.zeros(len(required)),
            A_eq=a_eq,
            b_eq=b_eq,
            bounds=[(0.0, 1.0)] * len(columns) + [(None, None)],
            method="highs",
        )
        if not result.success:
            raise RuntimeError(result.message)
        threshold = Fraction(1, len(cells))
        weights = [
            Fraction(-float(value)).limit_denominator(1_000_000)
            for value in result.ineqlin.marginals
        ]
        weight_sum = sum(weights, Fraction())
        weights = [weight / weight_sum for weight in weights]
        upper = Fraction()
        for p, options in rows:
            upper += max(
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
        impossible = upper < threshold
        all_impossible &= impossible
        print(
            f"target11={third_target} required_cells={len(required)} "
            f"optimum={result.x[-1]:.15f} exact_upper={upper} "
            f"result={'IMPOSSIBLE' if impossible else 'INCONCLUSIVE'}"
        )
    print(
        f"period={args.period} candidates={len(candidates)} anchor_cells={len(cells)} "
        f"result={'IMPOSSIBLE' if all_impossible else 'INCONCLUSIVE'}"
    )
    return 2 if all_impossible else 0


if __name__ == "__main__":
    raise SystemExit(main())
