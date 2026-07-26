#!/usr/bin/env python3
"""Exact-dual anchor-capacity obstruction including perfect-power factors."""

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


def power_target_congruence(h: int, p: int, power: int) -> tuple[int, int]:
    d = math.gcd(power, p - 1)
    step = (p - 1) // h
    divisor = math.gcd(step, d)
    half = (p - 1) // 2
    if half % divisor:
        raise RuntimeError(f"prime {p} has no power-compatible target")
    modulus = d // divisor
    residue = (
        0
        if modulus == 1
        else (half // divisor) * pow(step // divisor, -1, modulus) % modulus
    )
    return residue, modulus


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--anchor-indices", default="0,1,2")
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
    indices = tuple(int(value) for value in args.anchor_indices.split(","))
    if len(indices) not in (2, 3):
        raise SystemExit("two or three anchor indices are required")
    anchors = [candidates[index] for index in indices]
    algebraic_primes = tuple(
        prime for prime in exact_uncovered.factor(args.power) if prime % 2
    )
    sophie_germain = args.power % 4 == 0
    period = math.lcm(
        *(int(row["h"]) for row in anchors),
        *algebraic_primes,
        *(4,) if sophie_germain else (),
    )
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
    normalization_period = math.lcm(*(int(row["h"]) for row in anchors))
    translation_image = {
        tuple(
            (
                int(row["a"]) * args.power * k
                + int(row["b"]) * args.power * l
            )
            % int(row["h"])
            for row in anchors
        )
        for k in range(normalization_period)
        for l in range(normalization_period)
    }
    valid_phase_space = set(itertools.product(*valid_anchor_targets))
    normalized_targets = min(valid_phase_space)
    normalized_orbit = {
        tuple(
            (target + shift) % int(row["h"])
            for target, shift, row in zip(normalized_targets, image, anchors)
        )
        for image in translation_image
    }
    if normalized_orbit != valid_phase_space:
        raise RuntimeError("anchor targets are not jointly normalizable in power mode")

    def cell(k: int, l: int) -> tuple[int, ...]:
        values = tuple(
            (int(row["a"]) * k + int(row["b"]) * l) % int(row["h"])
            for row in anchors
        )
        power_values = tuple(
            coordinate % prime
            for prime in algebraic_primes
            for coordinate in (k, l)
        )
        sophie_values = (k % 4, l % 4) if sophie_germain else ()
        return values + power_values + sophie_values

    cells = tuple(sorted({cell(k, l) for k in range(period) for l in range(period)}))
    required = tuple(
        item
        for item in cells
        if all(
            value != normalized_targets[index]
            for index, value in enumerate(item[: len(anchors)])
        )
        and not any(
            item[len(anchors) + 2 * index] == 0
            and item[len(anchors) + 2 * index + 1] == 0
            for index in range(len(algebraic_primes))
        )
        and not (
            sophie_germain
            and item[len(anchors) + 2 * len(algebraic_primes)] == 2
            and item[len(anchors) + 2 * len(algebraic_primes) + 1] == 0
        )
    )

    excluded = set(indices)
    rows = []
    for index, row in enumerate(candidates):
        if index in excluded:
            continue
        h = int(row["h"])
        p = int(row["p"])
        a = int(row["a"])
        b = int(row["b"])
        shared = math.gcd(h, period)
        compatible_by_target = [set() for _ in range(shared)]
        for k in range(period):
            for l in range(period):
                compatible_by_target[(a * k + b * l) % shared].add(cell(k, l))
        try:
            residue, modulus = power_target_congruence(h, p, args.power)
        except RuntimeError:
            continue
        target_period = math.lcm(shared, modulus)
        target_residues = {
            target % shared
            for target in range(target_period)
            if target % modulus == residue
        }
        options = []
        seen = set()
        for target in sorted(target_residues):
            compatible = tuple(
                item for item in cells if item in compatible_by_target[target]
            )
            if compatible in seen:
                continue
            seen.add(compatible)
            options.append((compatible, Fraction(1, h * len(compatible))))
        if not options:
            raise RuntimeError(f"no valid options for prime {p}")
        rows.append((p, options))

    columns = [
        (row_index, option_index)
        for row_index, (_p, options) in enumerate(rows)
        for option_index in range(len(options))
    ]
    a_eq = coo_matrix(
        (
            np.ones(len(columns)),
            (
                [row_index for row_index, _option_index in columns],
                range(len(columns)),
            ),
        ),
        shape=(len(rows), len(columns) + 1),
    ).tocsr()
    required_index = {item: index for index, item in enumerate(required)}
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
    a_ub = coo_matrix(
        (ub_values, (ub_rows, ub_columns)),
        shape=(len(required), len(columns) + 1),
    ).tocsr()
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

    threshold = Fraction(1, len(cells))
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
                    weights[cell_index] * contribution
                    for cell_index, required_cell in enumerate(required)
                    if required_cell in compatible
                ),
                Fraction(),
            )
            for compatible, contribution in options
        )
    impossible = exact_upper < threshold
    print(
        f"pool={args.pool} power={args.power} candidates={len(candidates)} "
        f"anchors={[(row['p'], row['h']) for row in anchors]} "
        f"normalized_targets={normalized_targets} "
        f"algebraic_primes={algebraic_primes} "
        f"sophie_germain={sophie_germain} cells={len(cells)} "
        f"required={len(required)}"
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
