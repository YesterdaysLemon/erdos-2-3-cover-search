#!/usr/bin/env python3
"""Project a saved phase map onto the valid targets for a new power."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from power_anchor_capacity_lp import power_target_congruence


def cyclic_distance(left: int, right: int, modulus: int) -> int:
    delta = (left - right) % modulus
    return min(delta, modulus - delta)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--normalize-primes", default="")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.power < 1:
        raise SystemExit("--power must be positive")
    normalizers = {
        int(value) for value in args.normalize_primes.split(",") if value
    }
    source = {
        int(prime): int(target)
        for prime, target in json.loads(args.source.read_text()).items()
    }
    projected: dict[str, int] = {}
    retained = 0
    changed = 0
    absent = 0
    incompatible = 0
    density = 0.0
    present_primes = set()
    for raw in json.loads(args.pool.read_text())["choices"]:
        h = int(raw["h"])
        p = int(raw["p"])
        present_primes.add(p)
        try:
            residue, target_modulus = power_target_congruence(
                h, p, args.power
            )
        except RuntimeError:
            incompatible += 1
            continue
        if p in normalizers:
            if residue:
                raise RuntimeError(
                    f"normalization target zero is invalid for prime {p}"
                )
            target = 0
        elif p not in source:
            target = residue
            absent += 1
        else:
            old = source[p] % h
            lower_delta = (old - residue) % target_modulus
            lower = (old - lower_delta) % h
            upper = (old + (-lower_delta) % target_modulus) % h
            target = min(
                (lower, upper),
                key=lambda value: (
                    cyclic_distance(value, old, h),
                    value,
                ),
            )
            retained += int(target == old)
            changed += int(target != old)
        if target % target_modulus != residue:
            raise AssertionError("projected target violates power congruence")
        projected[str(p)] = target
        density += 1 / h

    missing_normalizers = normalizers - present_primes
    if missing_normalizers:
        raise RuntimeError(
            f"normalization primes missing from pool: {missing_normalizers}"
        )
    args.output.write_text(json.dumps(projected) + "\n")
    print(
        f"rows={len(projected)} density={density:.12f} "
        f"retained={retained} changed={changed} absent={absent} "
        f"incompatible={incompatible} power={args.power} "
        f"output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
