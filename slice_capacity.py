#!/usr/bin/env python3
"""Exact 2-adic slice-capacity bound for periodic Erdos-203 covers.

Normalize the p=5 fibre to k+3l == 0 (mod 4).  For every other prime,
compute which of the four k+3l slices its chosen fibre can meet and the exact
density contributed to each compatible slice.  An integer optimization then
maximizes the minimum raw capacity of the three slices not covered by p=5.

Capacity below 1/4 proves that the complete prime pool for the period cannot
cover Z^2.  Capacity at least 1/4 is only a necessary condition.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from fractions import Fraction
from pathlib import Path

import local_cover


def v2(n: int) -> int:
    exponent = 0
    while n % 2 == 0:
        exponent += 1
        n //= 2
    return exponent


def slice_options(h: int, a: int, b: int) -> list[tuple[tuple[int, ...], Fraction]]:
    two_part = 1 << v2(h)
    modulus = max(4, two_part)
    seen: set[tuple[int, ...]] = set()
    options = []
    for c2 in range(two_part):
        compatible = tuple(
            residue
            for residue in range(4)
            if any(
                (a * k + b * l - c2) % two_part == 0
                and (k + 3 * l - residue) % 4 == 0
                for k in range(modulus)
                for l in range(modulus)
            )
        )
        if compatible in seen:
            continue
        seen.add(compatible)
        options.append((compatible, Fraction(1, h * len(compatible))))
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
    rows = []
    for h, p, a, b, ord2, ord3 in candidates:
        if p == 5:
            continue
        rows.append((p, h, slice_options(h, a, b)))

    scale = math.lcm(
        *(contribution.denominator for _, _, opts in rows for _, contribution in opts),
        4,
    )
    optimizer = z3.Optimize()
    variables = []
    for p, h, options in rows:
        choice = z3.Int(f"choice_{p}")
        optimizer.add(choice >= 0, choice < len(options))
        variables.append((choice, options))

    totals = []
    for residue in (1, 2, 3):
        pieces = []
        for choice, options in variables:
            pieces.append(
                z3.Sum(
                    *[
                        z3.If(
                            choice == index,
                            int(contribution * scale) if residue in compatible else 0,
                            0,
                        )
                        for index, (compatible, contribution) in enumerate(options)
                    ]
                )
            )
        totals.append(z3.Sum(*pieces))
    minimum = z3.Int("minimum_slice_capacity")
    for total in totals:
        optimizer.add(minimum <= total)
    optimizer.maximize(minimum)
    if optimizer.check() != z3.sat:
        raise RuntimeError("slice-capacity optimization unexpectedly failed")
    model = optimizer.model()
    optimum = Fraction(model.eval(minimum).as_long(), scale)
    threshold = Fraction(1, 4)
    print(f"period={args.period}")
    print(f"eligible_primes={len(candidates)} gcd_bits={common.bit_length()}")
    print(f"factor_count={len(factorization)}")
    print(f"minimum_slice_capacity={optimum} ({float(optimum):.12f})")
    print(f"required={threshold} result={'IMPOSSIBLE' if optimum < threshold else 'INCONCLUSIVE'}")
    return 2 if optimum < threshold else 0


if __name__ == "__main__":
    raise SystemExit(main())
