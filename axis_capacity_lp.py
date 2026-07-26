#!/usr/bin/env python3
"""Exact-dual capacity screen for a perfect-power construction on one axis.

The axis is partitioned into cells modulo a small period L.  For every prime
fibre and every power-compatible target, its restriction to the axis is either
empty or one ordinary residue class.  A sparse LP chooses a fractional target
mixture for each fibre and maximizes the minimum normalized conditional
capacity across the residual cells.  A rationally checked ratio below one
proves that no phase assignment from the declared finite pool covers that
axis (with the requested anchors fixed).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from fractions import Fraction
from pathlib import Path

import exact_uncovered
from power_anchor_capacity_lp import power_target_congruence
from round_fractional_phases import combine_congruences


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--axis", choices=("k", "l"), required=True)
    parser.add_argument("--cell-period", type=int, default=60)
    parser.add_argument("--fixed-primes", required=True)
    parser.add_argument("--max-component", type=int, default=0)
    parser.add_argument("--weight-denominator", type=int, default=100000)
    args = parser.parse_args()
    if args.power < 1 or args.cell_period < 1:
        raise SystemExit("--power and --cell-period must be positive")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from scipy.optimize import linprog  # type: ignore
    from scipy.sparse import coo_matrix  # type: ignore

    source = json.loads(args.pool.read_text())
    phases = json.loads(args.phase_file.read_text())
    fixed_primes = {
        int(value) for value in args.fixed_primes.split(",") if value
    }
    rows = []
    for raw in source["choices"]:
        h = int(raw["h"])
        if args.max_component and max(
            (
                prime**exponent
                for prime, exponent in exact_uncovered.factor(h).items()
            ),
            default=1,
        ) > args.max_component:
            continue
        rows.append(raw)
    by_prime = {int(row["p"]): row for row in rows}
    missing = sorted(fixed_primes - set(by_prime))
    if missing:
        raise RuntimeError(f"fixed primes are absent: {missing}")

    cell_period = args.cell_period
    required = [Fraction(1) for _ in range(cell_period)]
    algebraic_primes = tuple(
        prime
        for prime in exact_uncovered.factor(args.power)
        if prime % 2
    )
    for prime in algebraic_primes:
        if cell_period % prime == 0:
            for cell in range(cell_period):
                if cell % prime == 0:
                    required[cell] = Fraction()
        else:
            for cell in range(cell_period):
                required[cell] *= Fraction(prime - 1, prime)
    if args.axis == "k" and args.power % 4 == 0:
        common = math.gcd(4, cell_period)
        for cell in range(cell_period):
            if (cell - 2) % common == 0:
                required[cell] = Fraction()

    # Fixed fibres are treated by a union bound.  In the intended anchor use,
    # their reduced axis moduli divide L and hence they remove whole cells.
    for prime in fixed_primes:
        row = by_prime[prime]
        h = int(row["h"])
        coefficient = int(row["a" if args.axis == "k" else "b"]) % h
        c = int(phases[str(prime)]) % h
        common = math.gcd(coefficient, h)
        if c % common:
            continue
        modulus = h // common
        residue = (
            0
            if modulus == 1
            else (
                (c // common)
                * pow(coefficient // common, -1, modulus)
            )
            % modulus
        )
        shared = math.gcd(modulus, cell_period)
        contribution = Fraction(shared, modulus)
        for cell in range(cell_period):
            if (cell - residue) % shared == 0:
                required[cell] = max(
                    Fraction(), required[cell] - contribution
                )

    active_cells = [
        cell for cell in range(cell_period) if required[cell] > 0
    ]
    active_index = {cell: index for index, cell in enumerate(active_cells)}
    candidate_rows = []
    for raw in rows:
        p = int(raw["p"])
        if p in fixed_primes:
            continue
        h = int(raw["h"])
        coefficient = int(raw["a" if args.axis == "k" else "b"]) % h
        common = math.gcd(coefficient, h)
        modulus = h // common
        shared = math.gcd(modulus, cell_period)
        power_residue, power_modulus = power_target_congruence(
            h, p, args.power
        )
        reduced_coefficient = coefficient // common
        options = []
        seen = set()
        for target_residue in range(shared):
            # The reduced axis class is x=target_residue (mod shared)
            # exactly when c/common = reduced_coefficient*target_residue
            # modulo shared.
            local_c = (
                common * reduced_coefficient * target_residue
            ) % (common * shared)
            try:
                combine_congruences(
                    power_residue,
                    power_modulus,
                    local_c,
                    common * shared,
                )
            except RuntimeError:
                continue
            cells = tuple(
                active_index[cell]
                for cell in active_cells
                if cell % shared == target_residue
            )
            if not cells or cells in seen:
                continue
            seen.add(cells)
            options.append((cells, Fraction(shared, modulus)))
        if options:
            candidate_rows.append((p, options))

    columns = [
        (row_index, option_index)
        for row_index, (_p, options) in enumerate(candidate_rows)
        for option_index in range(len(options))
    ]
    rho_column = len(columns)
    eq_rows = []
    eq_columns = []
    for column, (row_index, _option_index) in enumerate(columns):
        eq_rows.append(row_index)
        eq_columns.append(column)
    a_eq = coo_matrix(
        (
            np.ones(len(columns)),
            (eq_rows, eq_columns),
        ),
        shape=(len(candidate_rows), len(columns) + 1),
    ).tocsr()
    ub_rows = []
    ub_columns = []
    ub_values = []
    for column, (row_index, option_index) in enumerate(columns):
        cells, contribution = candidate_rows[row_index][1][option_index]
        for cell_index in cells:
            ub_rows.append(cell_index)
            ub_columns.append(column)
            ub_values.append(-float(contribution))
    for cell_index, cell in enumerate(active_cells):
        ub_rows.append(cell_index)
        ub_columns.append(rho_column)
        ub_values.append(float(required[cell]))
    a_ub = coo_matrix(
        (ub_values, (ub_rows, ub_columns)),
        shape=(len(active_cells), len(columns) + 1),
    ).tocsr()
    objective = np.zeros(len(columns) + 1)
    objective[rho_column] = -1.0
    result = linprog(
        objective,
        A_ub=a_ub,
        b_ub=np.zeros(len(active_cells)),
        A_eq=a_eq,
        b_eq=np.ones(len(candidate_rows)),
        bounds=[(0.0, 1.0)] * len(columns) + [(None, None)],
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)

    raw_weights = [
        Fraction(max(0.0, -float(value))).limit_denominator(
            args.weight_denominator
        )
        for value in result.ineqlin.marginals
    ]
    denominator = sum(
        (
            raw_weights[index] * required[cell]
            for index, cell in enumerate(active_cells)
        ),
        Fraction(),
    )
    exact_upper = Fraction()
    for _p, options in candidate_rows:
        exact_upper += max(
            sum(
                (
                    raw_weights[cell_index] * contribution
                    for cell_index in cells
                ),
                Fraction(),
            )
            for cells, contribution in options
        )
    exact_ratio = exact_upper / denominator
    print(
        f"axis={args.axis} rows={len(candidate_rows)} "
        f"active_cells={len(active_cells)} options={len(columns)} "
        f"float_ratio={result.x[rho_column]:.15f} "
        f"exact_ratio={float(exact_ratio):.15f} "
        f"exact_numerator={exact_upper} "
        f"exact_denominator={denominator} "
        f"impossible={exact_ratio < 1}",
        flush=True,
    )
    return 0 if exact_ratio >= 1 else 2


if __name__ == "__main__":
    raise SystemExit(main())
