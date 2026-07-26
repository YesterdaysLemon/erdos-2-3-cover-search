#!/usr/bin/env python3
"""Condition a derived affine-fibre pool on a rectangular lattice coset.

Write

    k = k0 + K*x,    l = l0 + L*y.

For each derived row, intersect its saved target congruence with the
congruence required to meet this coset, divide by the exact common factor,
and emit the resulting row in ``(x, y)`` coordinates.  This generalizes
``condition_derived_cell.py`` to different coordinate periods.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import exact_uncovered
from local_phase_cegis import merge_congruences


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--k-period", type=int, required=True)
    parser.add_argument("--l-period", type=int, required=True)
    parser.add_argument("--k-residue", type=int, required=True)
    parser.add_argument("--l-residue", type=int, required=True)
    parser.add_argument("--exclude-primes", default="")
    parser.add_argument("--canonicalize-algebraic", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.k_period < 1 or args.l_period < 1:
        raise SystemExit("coordinate periods must be positive")

    source = json.loads(args.pool.read_text())
    excluded = {
        int(value) for value in args.exclude_primes.split(",") if value
    }
    k0 = args.k_residue % args.k_period
    l0 = args.l_residue % args.l_period
    algebraic_primes = tuple(
        int(prime) for prime in source.get("algebraic_primes", ())
    )

    canonical_primes = []
    shift_x_congruences = []
    shift_y_congruences = []
    for prime in algebraic_primes:
        if (
            math.gcd(args.k_period, prime) != 1
            or math.gcd(args.l_period, prime) != 1
        ):
            if args.canonicalize_algebraic:
                raise RuntimeError(
                    "algebraic canonicalization requires both coordinate "
                    f"periods to be invertible modulo {prime}"
                )
            continue
        canonical_primes.append(prime)
        shift_x_congruences.append(
            (
                -k0 * pow(args.k_period, -1, prime) % prime,
                prime,
            )
        )
        shift_y_congruences.append(
            (
                -l0 * pow(args.l_period, -1, prime) % prime,
                prime,
            )
        )
    if args.canonicalize_algebraic and canonical_primes:
        shift_x = exact_uncovered.crt(shift_x_congruences)
        shift_y = exact_uncovered.crt(shift_y_congruences)
    else:
        shift_x = 0
        shift_y = 0
        canonical_primes = []

    rows = []
    incompatible = 0
    full_cell_primes = []
    for raw in source["choices"]:
        prime = int(raw["p"])
        if prime in excluded:
            continue
        h = int(raw["h"])
        a = int(raw["a"]) % h
        b = int(raw["b"]) % h
        residue = int(raw["target_residue"]) % h
        modulus = int(raw["target_modulus"])
        if h % modulus:
            raise RuntimeError(
                f"target modulus does not divide h for prime {prime}"
            )

        base = (a * k0 + b * l0) % h
        common = math.gcd(
            args.k_period * a,
            args.l_period * b,
            h,
        )
        merged = merge_congruences(
            residue,
            modulus,
            base,
            common,
        )
        if merged is None:
            incompatible += 1
            continue
        combined, combined_modulus = merged

        new_h = h // common
        if new_h == 1:
            full_cell_primes.append(prime)
            continue
        new_a = (args.k_period * a // common) % new_h
        new_b = (args.l_period * b // common) % new_h
        new_residue = ((combined - base) // common) % new_h
        new_modulus = combined_modulus // common
        target_shift = (new_a * shift_x + new_b * shift_y) % new_h
        new_residue = (new_residue - target_shift) % new_modulus
        if new_h % new_modulus:
            raise AssertionError("child target modulus does not divide h")
        if math.gcd(new_a, new_b, new_h) != 1:
            raise AssertionError("child affine map is not surjective")

        rows.append(
            {
                "h": new_h,
                "p": prime,
                "a": new_a,
                "b": new_b,
                "ord2": new_h // math.gcd(new_a, new_h),
                "ord3": new_h // math.gcd(new_b, new_h),
                "c": new_residue,
                "target_residue": new_residue,
                "target_modulus": new_modulus,
                "parent_h": h,
                "parent_a": a,
                "parent_b": b,
                "parent_base": base,
                "parent_common": common,
                "coordinate_shift": target_shift,
            }
        )

    rows.sort(key=lambda row: (int(row["h"]), int(row["p"])))
    result = {
        "source": str(args.pool),
        "power": int(source.get("power", 1)),
        "rectangle_condition": {
            "k_period": args.k_period,
            "l_period": args.l_period,
            "k_residue": k0,
            "l_residue": l0,
            "coordinate_shift": {"x": shift_x, "y": shift_y},
        },
        "algebraic_primes": canonical_primes,
        "sophie_germain": bool(source.get("sophie_germain", False)),
        "excluded_primes": sorted(excluded),
        "incompatible_rows": incompatible,
        "full_cell_primes": sorted(full_cell_primes),
        "choices": rows,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"rows={len(rows)} incompatible={incompatible} "
        f"full_cell={len(full_cell_primes)} "
        f"density={sum(1 / int(row['h']) for row in rows):.12f} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
