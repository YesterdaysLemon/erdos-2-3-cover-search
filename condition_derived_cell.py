#!/usr/bin/env python3
"""Condition an already-derived affine-fibre pool on a further lattice cell.

Unlike condition_power_cell.py, this composes the row-specific
target_residue/target_modulus restrictions already present in the source.
An optional parent phase assignment is projected to the child coordinates.
Fixed parent phases that do not meet the selected cell are rigorously inactive
and are omitted; mutable phases may be retargeted to any compatible child
fibre because only coverage inside this cell is being synthesized.
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
    parser.add_argument("--cell-period", type=int, required=True)
    parser.add_argument("--k-residue", type=int, required=True)
    parser.add_argument("--l-residue", type=int, required=True)
    parser.add_argument("--phase-file", type=Path)
    parser.add_argument("--phase-output", type=Path)
    parser.add_argument("--fixed-primes", default="")
    parser.add_argument("--canonicalize-algebraic", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.cell_period < 1:
        raise SystemExit("--cell-period must be positive")
    if bool(args.phase_output) and not args.phase_file:
        raise SystemExit("--phase-output requires --phase-file")

    source = json.loads(args.pool.read_text())
    parent_phases = (
        {
            int(prime): int(target)
            for prime, target in json.loads(args.phase_file.read_text()).items()
        }
        if args.phase_file
        else {}
    )
    fixed_primes = {
        int(value) for value in args.fixed_primes.split(",") if value
    }
    k0 = args.k_residue % args.cell_period
    l0 = args.l_residue % args.cell_period
    algebraic_primes = tuple(
        int(prime) for prime in source.get("algebraic_primes", ())
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
    projected_phases = {}
    incompatible = 0
    fixed_inactive = []
    full_cell_primes = []
    for raw in source["choices"]:
        h = int(raw["h"])
        p = int(raw["p"])
        a = int(raw["a"]) % h
        b = int(raw["b"]) % h
        target_residue = int(raw["target_residue"])
        target_modulus = int(raw["target_modulus"])
        if h % target_modulus:
            raise RuntimeError(f"invalid parent target modulus for prime {p}")
        base = (a * k0 + b * l0) % h
        common = math.gcd(args.cell_period * a, args.cell_period * b, h)
        merged = merge_congruences(
            target_residue,
            target_modulus,
            base,
            common,
        )
        if merged is None:
            incompatible += 1
            continue
        combined, combined_modulus = merged

        parent_phase = parent_phases.get(p)
        parent_phase_meets = (
            parent_phase is not None
            and (parent_phase - base) % common == 0
        )
        if p in fixed_primes and not parent_phase_meets:
            fixed_inactive.append(p)
            continue

        new_h = h // common
        if new_h == 1:
            full_cell_primes.append(p)
            continue
        new_a = (args.cell_period * a // common) % new_h
        new_b = (args.cell_period * b // common) % new_h
        new_residue = ((combined - base) // common) % new_h
        new_modulus = combined_modulus // common
        target_shift = (new_a * shift_x + new_b * shift_y) % new_h
        new_residue = (new_residue - target_shift) % new_modulus
        if new_h % new_modulus:
            raise AssertionError("child target modulus does not divide h")
        if math.gcd(new_a, new_b, new_h) != 1:
            raise AssertionError("child affine map is not surjective")

        if parent_phase_meets:
            child_phase = (
                (parent_phase - base) // common - target_shift
            ) % new_h
            if child_phase % new_modulus != new_residue:
                raise AssertionError("projected phase violates child target")
        else:
            child_phase = new_residue
        projected_phases[p] = child_phase
        rows.append(
            {
                "h": new_h,
                "p": p,
                "a": new_a,
                "b": new_b,
                "ord2": new_h // math.gcd(new_a, new_h),
                "ord3": new_h // math.gcd(new_b, new_h),
                "c": child_phase,
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
    payload = {
        "source": str(args.pool),
        "power": int(source.get("power", 1)),
        "cell_condition": {
            "period": args.cell_period,
            "k_residue": k0,
            "l_residue": l0,
            "coordinate_shift": {"x": shift_x, "y": shift_y},
        },
        "algebraic_primes": canonical_primes,
        "sophie_germain": bool(source.get("sophie_germain", False)),
        "incompatible_rows": incompatible,
        "fixed_inactive_primes": sorted(fixed_inactive),
        "full_cell_primes": sorted(full_cell_primes),
        "choices": rows,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    if args.phase_output:
        args.phase_output.write_text(
            json.dumps(
                {
                    str(row["p"]): int(projected_phases[row["p"]])
                    for row in rows
                }
            )
            + "\n"
        )
    print(
        f"rows={len(rows)} incompatible={incompatible} "
        f"fixed_inactive={len(fixed_inactive)} "
        f"full_cell={len(full_cell_primes)} "
        f"density={sum(1 / int(row['h']) for row in rows):.12f} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
