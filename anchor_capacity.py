#!/usr/bin/env python3
"""Joint-anchor capacity bound for periodic Erdos-203 affine covers.

The p=5 and p=7 line maps are jointly surjective, so translation normalizes
both chosen fibres to zero.  Their joint residues partition the exponent
lattice into 4*6=24 equal cells, nine already covered by an anchor.  For every
other eligible prime, this program computes the exact anchor cells a fibre can
meet and optimizes its target choice.  If even the maximum minimum raw capacity
of an uncovered cell is below 1/24, the complete period is impossible.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from fractions import Fraction
from pathlib import Path

import local_cover


ANCHOR_PRIMES = (5, 7)

_OPTION_COMPAT_CACHE: dict[
    tuple[
        tuple[tuple[int, int, int], ...],
        tuple[tuple[int, ...], ...],
        int,
        int,
        int,
    ],
    tuple[tuple[tuple[int, ...], ...], ...],
] = {}


def option_cells(
    h: int,
    a: int,
    b: int,
    anchors: list[tuple[int, int, int]],
    cells: tuple[tuple[int, ...], ...],
) -> list[tuple[tuple[tuple[int, ...], ...], Fraction]]:
    anchor_period = math.lcm(*(modulus for modulus, _, _ in anchors))
    shared = math.gcd(h, anchor_period)
    cache_key = (tuple(anchors), cells, shared, a % shared, b % shared)
    cached = _OPTION_COMPAT_CACHE.get(cache_key)
    if cached is not None:
        return [
            (compatible, Fraction(1, h * len(compatible)))
            for compatible in cached
        ]
    seen: set[tuple[tuple[int, ...], ...]] = set()
    options = []
    # Bucket the torus in one pass.  The previous implementation rescanned
    # all anchor_period**2 points independently for each target residue,
    # adding an unnecessary factor of ``shared``.
    compatible_by_target = [set() for _ in range(shared)]
    for k in range(anchor_period):
        for l in range(anchor_period):
            target = (a * k + b * l) % shared
            cell = tuple(
                (anchor_a * k + anchor_b * l) % modulus
                for modulus, anchor_a, anchor_b in anchors
            )
            compatible_by_target[target].add(cell)
    for compatible_set in compatible_by_target:
        compatible = tuple(cell for cell in cells if cell in compatible_set)
        if compatible in seen:
            continue
        seen.add(compatible)
        options.append((compatible, Fraction(1, h * len(compatible))))
    _OPTION_COMPAT_CACHE[cache_key] = tuple(
        compatible for compatible, contribution in options
    )
    return options


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("period", type=int)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import z3  # type: ignore

    candidates, common, factorization = local_cover.get_complete_period_candidates(
        args.period
    )
    by_prime = {p: row for row in candidates for p in (row[1],)}
    if not all(p in by_prime for p in ANCHOR_PRIMES):
        raise RuntimeError("period does not contain both anchor primes")
    anchors = [
        (by_prime[p][0], by_prime[p][2], by_prime[p][3]) for p in ANCHOR_PRIMES
    ]
    anchor_period = math.lcm(*(modulus for modulus, _, _ in anchors))
    cell_set = {
        tuple(
            (anchor_a * k + anchor_b * l) % modulus
            for modulus, anchor_a, anchor_b in anchors
        )
        for k in range(anchor_period)
        for l in range(anchor_period)
    }
    cells = tuple(sorted(cell_set))
    expected_cells = math.prod(modulus for modulus, _, _ in anchors)
    if len(cells) != expected_cells:
        raise RuntimeError("anchor maps are not jointly surjective")
    required_cells = tuple(cell for cell in cells if all(value != 0 for value in cell))

    rows = []
    for h, p, a, b, ord2, ord3 in candidates:
        if p in ANCHOR_PRIMES:
            continue
        rows.append((p, h, option_cells(h, a, b, anchors, cells)))
    scale = math.lcm(
        *(contribution.denominator for _, _, opts in rows for _, contribution in opts),
        len(cells),
    )

    solver = z3.Solver()
    variables = []
    for p, h, options in rows:
        choices = [z3.Bool(f"choice_{p}_{index}") for index in range(len(options))]
        solver.add(z3.PbEq([(choice, 1) for choice in choices], 1))
        variables.append((choices, options))
    threshold = Fraction(1, len(cells))
    threshold_scaled = int(threshold * scale)
    for cell in required_cells:
        literals = []
        weights = []
        for choices, options in variables:
            for choice, (compatible, contribution) in zip(choices, options):
                if cell in compatible:
                    literals.append(choice)
                    weights.append(int(contribution * scale))
        solver.add(z3.PbGe(list(zip(literals, weights)), threshold_scaled))
    feasible = solver.check() == z3.sat
    print(f"period={args.period}")
    print(
        f"eligible_primes={len(candidates)} gcd_bits={common.bit_length()} "
        f"factor_count={len(factorization)}"
    )
    print(f"anchor_cells={len(cells)} required_cells={len(required_cells)}")
    print(f"required_capacity={threshold} ({float(threshold):.12f})")
    print(f"capacity_feasibility={'SAT' if feasible else 'UNSAT'}")
    print(f"result={'INCONCLUSIVE' if feasible else 'IMPOSSIBLE'}")
    return 0 if feasible else 2


if __name__ == "__main__":
    raise SystemExit(main())
