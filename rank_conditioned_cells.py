#!/usr/bin/env python3
"""Rank coarse open cells by the strength of their conditioned candidate pool.

This performs the modulus bookkeeping of condition_derived_cell.py without
materializing thousands of full JSON pools.  It is useful for choosing the
next cell in a monotone certificate search: higher conditioned density and
fewer active protected rows generally give the phase synthesizer more room.
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
    parser.add_argument("phases", type=Path)
    parser.add_argument("locked_primes", type=Path)
    parser.add_argument("cells", type=Path)
    parser.add_argument("--cell-period", type=int, default=432)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    phases = {
        int(prime): int(value)
        for prime, value in json.loads(args.phases.read_text()).items()
    }
    lock_payload = json.loads(args.locked_primes.read_text())
    locked = {
        int(prime)
        for prime in lock_payload.get("primes", lock_payload)
    }
    cell_payload = json.loads(args.cells.read_text())
    raw_cells = cell_payload.get("cells", cell_payload)

    ranked = []
    for raw_cell in raw_cells:
        if len(raw_cell) != 4:
            raise RuntimeError(f"expected four-coordinate cell, got {raw_cell}")
        k16, l16, k27, l27 = map(int, raw_cell)
        k0 = exact_uncovered.crt([(k16 % 16, 16), (k27 % 27, 27)])
        l0 = exact_uncovered.crt([(l16 % 16, 16), (l27 % 27, 27)])
        density = 0.0
        rows = 0
        incompatible = 0
        fixed_active = []
        fixed_inactive = []
        mutable_full_cell = []
        fixed_full_cell = []
        max_component = 1
        for raw in source["choices"]:
            h = int(raw["h"])
            prime = int(raw["p"])
            a = int(raw["a"]) % h
            b = int(raw["b"]) % h
            target_residue = int(raw["target_residue"])
            target_modulus = int(raw["target_modulus"])
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
            parent_meets = (phases[prime] - base) % common == 0
            if prime in locked and not parent_meets:
                fixed_inactive.append(prime)
                continue
            new_h = h // common
            if new_h == 1:
                if prime in locked:
                    fixed_full_cell.append(prime)
                else:
                    mutable_full_cell.append(prime)
                continue
            if prime in locked:
                fixed_active.append(prime)
            rows += 1
            density += 1 / new_h
            for component in exact_uncovered.factor(new_h):
                exponent = 0
                value = new_h
                while value % component == 0:
                    exponent += 1
                    value //= component
                max_component = max(max_component, component**exponent)

        ranked.append(
            {
                "cell": [k16, l16, k27, l27],
                "k_residue": k0,
                "l_residue": l0,
                "rows": rows,
                "density": density,
                "max_component": max_component,
                "incompatible": incompatible,
                "fixed_active_count": len(fixed_active),
                "fixed_active_primes": sorted(fixed_active),
                "fixed_inactive_count": len(fixed_inactive),
                "mutable_full_cell_primes": sorted(mutable_full_cell),
                "fixed_full_cell_primes": sorted(fixed_full_cell),
            }
        )

    ranked.sort(
        key=lambda item: (
            -len(item["mutable_full_cell_primes"]),
            -float(item["density"]),
            int(item["fixed_active_count"]),
            item["cell"],
        )
    )
    payload = {
        "pool": str(args.pool),
        "phases": str(args.phases),
        "locked_primes": str(args.locked_primes),
        "cells": str(args.cells),
        "cell_period": args.cell_period,
        "ranked": ranked,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"cells={len(ranked)} "
        f"best_density={ranked[0]['density']:.12f} "
        f"best_cell={ranked[0]['cell']} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
