#!/usr/bin/env python3
"""Coarse exact-dual power obstruction on normalized p=5,p=7 cells."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
from fractions import Fraction
from pathlib import Path

import exact_uncovered
from power_anchor_capacity_lp import power_target_congruence


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--weight-denominator", type=int, default=100000)
    args = parser.parse_args()
    if args.power < 1:
        raise SystemExit("--power must be positive")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from scipy.optimize import linprog  # type: ignore
    from scipy.sparse import coo_matrix  # type: ignore

    candidates = json.loads(args.pool.read_text())["choices"]
    by_prime = {int(row["p"]): (index, row) for index, row in enumerate(candidates)}
    if 5 not in by_prime or 7 not in by_prime:
        raise RuntimeError("pool must contain p=5 and p=7 anchors")
    anchor_indices = {by_prime[5][0], by_prime[7][0]}
    anchors = [by_prime[5][1], by_prime[7][1]]
    algebraic_primes = tuple(
        prime for prime in exact_uncovered.factor(args.power) if prime % 2
    )
    sophie_germain = args.power % 4 == 0
    period = math.lcm(4, 6, *algebraic_primes, *(4,) if sophie_germain else ())

    axis = np.arange(period, dtype=np.int64)
    k = np.repeat(axis, period)
    l = np.tile(axis, period)
    anchor5 = (k + 3 * l) % 4
    anchor7 = (2 * k + l) % 6
    cell_id = anchor5 * 6 + anchor7
    valid_anchor_targets = []
    for row in anchors:
        residue, modulus = power_target_congruence(
            int(row["h"]), int(row["p"]), args.power
        )
        valid_anchor_targets.append(
            tuple(
                target
                for target in range(int(row["h"]))
                if target % modulus == residue
            )
        )
    valid_phase_space = set(itertools.product(*valid_anchor_targets))
    translation_image = {
        tuple(
            (
                int(row["a"]) * args.power * shift_k
                + int(row["b"]) * args.power * shift_l
            )
            % int(row["h"])
            for row in anchors
        )
        for shift_k in range(12)
        for shift_l in range(12)
    }
    normalized_targets = min(valid_phase_space)
    normalized_orbit = {
        tuple(
            (target + shift) % int(row["h"])
            for target, shift, row in zip(normalized_targets, image, anchors)
        )
        for image in translation_image
    }
    if normalized_orbit != valid_phase_space:
        raise RuntimeError("p=5,p=7 targets are not jointly normalizable")
    algebraic_uncovered = np.ones(period * period, dtype=bool)
    for prime in algebraic_primes:
        algebraic_uncovered &= ~((k % prime == 0) & (l % prime == 0))
    if sophie_germain:
        algebraic_uncovered &= ~((k % 4 == 2) & (l % 4 == 0))
    required_cells = [
        value5 * 6 + value7
        for value5 in range(4)
        if value5 != normalized_targets[0]
        for value7 in range(6)
        if value7 != normalized_targets[1]
    ]
    demand_counts = np.bincount(
        cell_id[algebraic_uncovered], minlength=24
    ).astype(np.int64)
    if any(demand_counts[cell] == 0 for cell in required_cells):
        raise RuntimeError("an anchor residual cell has zero demand")

    pattern_cache = {}
    rows = []
    for index, row in enumerate(candidates):
        if index in anchor_indices:
            continue
        h = int(row["h"])
        p = int(row["p"])
        a = int(row["a"])
        b = int(row["b"])
        shared = math.gcd(h, period)
        key = (shared, a % shared, b % shared)
        counts = pattern_cache.get(key)
        if counts is None:
            targets = (key[1] * k + key[2] * l) % shared
            combined = targets[algebraic_uncovered] * 24 + cell_id[algebraic_uncovered]
            counts = np.bincount(
                combined, minlength=shared * 24
            ).reshape(shared, 24)
            pattern_cache[key] = counts
        try:
            residue, modulus = power_target_congruence(h, p, args.power)
        except RuntimeError:
            continue
        target_period = math.lcm(shared, modulus)
        targets = sorted(
            {
                target % shared
                for target in range(target_period)
                if target % modulus == residue
            }
        )
        options = []
        seen = set()
        for target in targets:
            vector = tuple(int(counts[target, cell]) for cell in required_cells)
            if vector in seen:
                continue
            seen.add(vector)
            options.append((shared, h, vector))
        if not options:
            raise RuntimeError(f"no valid target for prime {p}")
        rows.append((p, options))

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
        shared, h, vector = rows[row_index][1][option_index]
        for cell_index, count in enumerate(vector):
            if count:
                ub_rows.append(cell_index)
                ub_columns.append(column)
                ub_values.append(
                    -count * shared / (h * int(demand_counts[required_cells[cell_index]]))
                )
    for cell_index in range(len(required_cells)):
        ub_rows.append(cell_index)
        ub_columns.append(len(columns))
        ub_values.append(1.0)
    ub = coo_matrix(
        (ub_values, (ub_rows, ub_columns)),
        shape=(len(required_cells), len(columns) + 1),
    ).tocsr()
    objective = np.zeros(len(columns) + 1)
    objective[-1] = -1.0
    result = linprog(
        objective,
        A_ub=ub,
        b_ub=np.zeros(len(required_cells)),
        A_eq=eq,
        b_eq=np.ones(len(rows)),
        bounds=[(0.0, 1.0)] * len(columns) + [(None, None)],
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    weights = [
        Fraction(-float(value)).limit_denominator(args.weight_denominator)
        for value in result.ineqlin.marginals
    ]
    weight_sum = sum(weights, Fraction())
    weights = [weight / weight_sum for weight in weights]
    exact_upper = Fraction()
    for _p, options in rows:
        exact_upper += max(
            sum(
                (
                    weights[cell_index]
                    * Fraction(
                        count * shared,
                        h * int(demand_counts[required_cells[cell_index]]),
                    )
                    for cell_index, count in enumerate(vector)
                ),
                Fraction(),
            )
            for shared, h, vector in options
        )
    impossible = exact_upper < 1
    print(
        f"pool={args.pool} power={args.power} candidates={len(candidates)} "
        f"normalized_targets={normalized_targets} "
        f"algebraic_primes={algebraic_primes} "
        f"sophie_germain={sophie_germain} period={period} "
        f"patterns={len(pattern_cache)} columns={len(columns)}"
    )
    print(
        f"lp_ratio={result.x[-1]:.15f} "
        f"exact_ratio_upper={float(exact_upper):.15f} "
        f"bits=({exact_upper.numerator.bit_length()},"
        f"{exact_upper.denominator.bit_length()})"
    )
    print(f"result={'IMPOSSIBLE' if impossible else 'INCONCLUSIVE'}")
    return 2 if impossible else 0


if __name__ == "__main__":
    raise SystemExit(main())
