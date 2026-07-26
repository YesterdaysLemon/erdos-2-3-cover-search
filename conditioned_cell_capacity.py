#!/usr/bin/env python3
"""Exact per-cell union-bound capacity for a derived affine-fibre pool.

Fix (k,l) = (r,s) modulo A.  A row

    a*k + b*l = c (mod h),  c = target_residue (mod target_modulus)

either misses the whole A-cell for every allowed c, or can cover a fraction
g/h of it, where g = gcd(h, A*a, A*b).  Summing those exact fractions gives
an upper bound even though every row must use one global phase.  Therefore an
admissible cell with total below one is a rigorous obstruction to a cover by
the declared pool.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def congruences_intersect(
    residue_a: int,
    modulus_a: int,
    residue_b: int,
    modulus_b: int,
) -> bool:
    return (residue_a - residue_b) % math.gcd(modulus_a, modulus_b) == 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--cell-modulus", type=int)
    parser.add_argument("--cell-k", type=int, required=True)
    parser.add_argument("--cell-l", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    algebraic_primes = tuple(
        int(value) for value in payload.get("algebraic_primes", ())
    )
    cell_modulus = (
        args.cell_modulus
        if args.cell_modulus is not None
        else math.prod(algebraic_primes)
    )
    if cell_modulus < 1:
        raise SystemExit("--cell-modulus must be positive")
    r = args.cell_k % cell_modulus
    s = args.cell_l % cell_modulus
    excluded_by = [
        prime
        for prime in algebraic_primes
        if r % prime == 0 and s % prime == 0
    ]

    compatible = []
    capacity = Fraction(0)
    for row in rows:
        h = int(row["h"])
        a = int(row["a"])
        b = int(row["b"])
        target_residue = int(row["target_residue"])
        target_modulus = int(row["target_modulus"])
        if h % target_modulus:
            raise RuntimeError(
                f"target modulus {target_modulus} does not divide h={h}"
            )
        if math.gcd(h, a, b) != 1:
            raise RuntimeError(
                f"row p={row['p']} is not surjective modulo h={h}"
            )

        image_modulus = math.gcd(
            h,
            cell_modulus * a,
            cell_modulus * b,
        )
        image_residue = (a * r + b * s) % image_modulus
        if not congruences_intersect(
            image_residue,
            image_modulus,
            target_residue,
            target_modulus,
        ):
            continue
        contribution = Fraction(image_modulus, h)
        capacity += contribution
        compatible.append(
            {
                "p": int(row["p"]),
                "h": h,
                "a": a,
                "b": b,
                "target_residue": target_residue,
                "target_modulus": target_modulus,
                "image_residue": image_residue,
                "image_modulus": image_modulus,
                "contribution_numerator": contribution.numerator,
                "contribution_denominator": contribution.denominator,
            }
        )

    result = {
        "pool": str(args.pool),
        "row_count": len(rows),
        "cell_modulus": cell_modulus,
        "cell": [r, s],
        "algebraic_primes": list(algebraic_primes),
        "excluded_by_algebraic_identity": excluded_by,
        "required_cell": not excluded_by,
        "compatible_row_count": len(compatible),
        "capacity_numerator": capacity.numerator,
        "capacity_denominator": capacity.denominator,
        "capacity_decimal": float(capacity),
        "capacity_below_one": capacity < 1,
        "compatible_rows": compatible,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"cell=({r},{s}) mod {cell_modulus} required={not excluded_by} "
        f"compatible={len(compatible)}/{len(rows)} "
        f"capacity={capacity.numerator}/{capacity.denominator} "
        f"decimal={float(capacity):.15f} below_one={capacity < 1}",
        flush=True,
    )
    if excluded_by:
        return 3
    return 0 if capacity < 1 else 2


if __name__ == "__main__":
    raise SystemExit(main())
