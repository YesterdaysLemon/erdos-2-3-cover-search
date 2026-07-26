#!/usr/bin/env python3
"""Condition a power-compatible affine-fibre pool on one lattice cell.

Write k = k0 + L*x and l = l0 + L*y.  For every original fibre

    a*k + b*l = c (mod h),

intersect the allowed power targets for c with the congruence needed for the
fibre to meet the selected cell, divide out the common gcd, and emit the
resulting affine fibre in (x,y).  The extra target_residue/target_modulus
fields record the exact remaining phase restriction; component_core.py may
ignore them when applying target-independent necessary conditions.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import exact_uncovered
from power_anchor_capacity_lp import power_target_congruence
from round_fractional_phases import combine_congruences


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--cell-period", type=int, required=True)
    parser.add_argument("--k-residue", type=int, required=True)
    parser.add_argument("--l-residue", type=int, required=True)
    parser.add_argument("--exclude-primes", default="")
    parser.add_argument(
        "--max-component",
        type=int,
        default=0,
        help="if nonzero, retain rows whose original prime-power components fit",
    )
    parser.add_argument(
        "--canonicalize-algebraic",
        action="store_true",
        help=(
            "translate the cell coordinates so every surviving odd "
            "power-factor identity excludes the origin modulo that prime"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.power < 1 or args.cell_period < 1:
        raise SystemExit("--power and --cell-period must be positive")

    excluded = {
        int(value) for value in args.exclude_primes.split(",") if value
    }
    k0 = args.k_residue % args.cell_period
    l0 = args.l_residue % args.cell_period
    source = json.loads(args.pool.read_text())
    algebraic_primes = tuple(
        prime
        for prime in exact_uncovered.factor(args.power)
        if prime % 2
    )
    canonical_primes = []
    shift_x_congruences = []
    shift_y_congruences = []
    for prime in algebraic_primes:
        common = math.gcd(args.cell_period, prime)
        if common != 1:
            if k0 % prime == 0 and l0 % prime == 0:
                raise RuntimeError(
                    f"the whole selected cell is algebraically covered mod {prime}"
                )
            continue
        canonical_primes.append(prime)
        inverse = pow(args.cell_period, -1, prime)
        shift_x_congruences.append(((-k0 * inverse) % prime, prime))
        shift_y_congruences.append(((-l0 * inverse) % prime, prime))
    if args.canonicalize_algebraic and canonical_primes:
        shift_x = exact_uncovered.crt(shift_x_congruences)
        shift_y = exact_uncovered.crt(shift_y_congruences)
    else:
        shift_x = 0
        shift_y = 0
        canonical_primes = []

    rows = []
    incompatible = 0
    full_cell_rows = 0
    for raw in source["choices"]:
        p = int(raw["p"])
        if p in excluded:
            continue
        h = int(raw["h"])
        if args.max_component and max(
            (
                prime**exponent
                for prime, exponent in exact_uncovered.factor(h).items()
            ),
            default=1,
        ) > args.max_component:
            continue
        a = int(raw["a"]) % h
        b = int(raw["b"]) % h
        power_residue, power_modulus = power_target_congruence(
            h, p, args.power
        )
        base = (a * k0 + b * l0) % h
        common = math.gcd(args.cell_period * a, args.cell_period * b, h)
        try:
            combined, combined_modulus = combine_congruences(
                power_residue,
                power_modulus,
                base,
                common,
            )
        except RuntimeError:
            incompatible += 1
            continue

        new_h = h // common
        new_a = (args.cell_period * a // common) % new_h
        new_b = (args.cell_period * b // common) % new_h
        target_residue = ((combined - base) // common) % new_h
        target_modulus = combined_modulus // common
        target_shift = (new_a * shift_x + new_b * shift_y) % new_h
        target_residue = (target_residue - target_shift) % target_modulus
        if new_h % target_modulus:
            raise AssertionError("derived target modulus does not divide h")
        if math.gcd(new_a, new_b, new_h) != 1:
            raise AssertionError("derived affine map is not surjective")
        if new_h == 1:
            full_cell_rows += 1
        rows.append(
            {
                "h": new_h,
                "p": p,
                "a": new_a,
                "b": new_b,
                "ord2": new_h // math.gcd(new_a, new_h),
                "ord3": new_h // math.gcd(new_b, new_h),
                "c": target_residue,
                "target_residue": target_residue,
                "target_modulus": target_modulus,
                "original_h": h,
                "original_a": a,
                "original_b": b,
                "coordinate_shift": target_shift,
            }
        )
    rows.sort(key=lambda row: (int(row["h"]), int(row["p"])))
    payload = {
        "source": str(args.pool),
        "source_scan_status": {
            "complete_order_interval": bool(
                source.get("complete_order_interval", False)
            ),
            "unresolved": int(source.get("unresolved", 0)),
        },
        "power": args.power,
        "cell_condition": {
            "period": args.cell_period,
            "k_residue": k0,
            "l_residue": l0,
            "excluded_primes": sorted(excluded),
            "coordinate_shift": {
                "x": shift_x,
                "y": shift_y,
            },
        },
        "algebraic_primes": canonical_primes,
        "incompatible_rows": incompatible,
        "full_cell_rows": full_cell_rows,
        "choices": rows,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"rows={len(rows)} incompatible={incompatible} "
        f"full_cell_rows={full_cell_rows} "
        f"density={sum(1 / int(row['h']) for row in rows):.12f} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
