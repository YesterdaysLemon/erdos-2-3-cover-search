#!/usr/bin/env python3
"""Condition a derived affine-fibre pool on an affine lattice coset.

The parent coordinates are parameterized by

    k = k0 + A*x + B*y
    l = l0 + C*x + D*y.

Each saved target restriction is intersected with the constant congruence
forced by the lattice, then divided by the exact common factor of the two new
coefficients and the row modulus.
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
    parser.add_argument("--matrix", required=True, help="A,B,C,D")
    parser.add_argument("--offset", required=True, help="k0,l0")
    parser.add_argument("--exclude-primes", default="")
    parser.add_argument("--canonicalize-algebraic", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    matrix_values = tuple(int(value) for value in args.matrix.split(","))
    offset_values = tuple(int(value) for value in args.offset.split(","))
    if len(matrix_values) != 4 or len(offset_values) != 2:
        raise SystemExit("--matrix needs A,B,C,D and --offset needs k0,l0")
    matrix_a, matrix_b, matrix_c, matrix_d = matrix_values
    k0, l0 = offset_values
    determinant = matrix_a * matrix_d - matrix_b * matrix_c
    if determinant == 0:
        raise SystemExit("lattice matrix must have nonzero determinant")

    source = json.loads(args.pool.read_text())
    excluded = {
        int(value) for value in args.exclude_primes.split(",") if value
    }
    algebraic_primes = tuple(
        int(prime) for prime in source.get("algebraic_primes", ())
    )
    canonical_primes = []
    shift_x_congruences = []
    shift_y_congruences = []
    for prime in algebraic_primes:
        if math.gcd(determinant, prime) != 1:
            if args.canonicalize_algebraic:
                raise RuntimeError(
                    f"lattice determinant is not invertible modulo {prime}"
                )
            continue
        inverse_det = pow(determinant, -1, prime)
        shift_x_congruences.append(
            (
                inverse_det
                * (-matrix_d * k0 + matrix_b * l0)
                % prime,
                prime,
            )
        )
        shift_y_congruences.append(
            (
                inverse_det
                * (matrix_c * k0 - matrix_a * l0)
                % prime,
                prime,
            )
        )
        canonical_primes.append(prime)
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
        coeff_x = a * matrix_a + b * matrix_c
        coeff_y = a * matrix_b + b * matrix_d
        common = math.gcd(coeff_x, coeff_y, h)
        merged = merge_congruences(
            residue, modulus, base, common
        )
        if merged is None:
            incompatible += 1
            continue
        combined, combined_modulus = merged
        new_h = h // common
        if new_h == 1:
            full_cell_primes.append(prime)
            continue
        new_a = (coeff_x // common) % new_h
        new_b = (coeff_y // common) % new_h
        new_modulus = combined_modulus // common
        target_shift = (new_a * shift_x + new_b * shift_y) % new_h
        new_residue = (
            (combined - base) // common - target_shift
        ) % new_modulus
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
                "parent_coeff_x": coeff_x,
                "parent_coeff_y": coeff_y,
                "parent_common": common,
                "coordinate_shift": target_shift,
            }
        )

    rows.sort(key=lambda row: (int(row["h"]), int(row["p"])))
    result = {
        "source": str(args.pool),
        "power": int(source.get("power", 1)),
        "lattice_condition": {
            "matrix": {
                "a": matrix_a,
                "b": matrix_b,
                "c": matrix_c,
                "d": matrix_d,
            },
            "offset": {"k": k0, "l": l0},
            "determinant": determinant,
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
        f"full_cell={len(full_cell_primes)} determinant={determinant} "
        f"density={sum(1 / int(row['h']) for row in rows):.12f} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
