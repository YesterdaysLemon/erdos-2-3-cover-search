#!/usr/bin/env python3
"""Exact complement check for the concurrent zero-phase power construction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import exact_uncovered
from power_anchor_capacity_lp import power_target_congruence


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument("--max-component", type=int, default=256)
    parser.add_argument("--limit", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.power < 1:
        raise SystemExit("--power must be positive")
    if args.period < 1:
        raise SystemExit("--period must be positive")

    rows = []
    incompatible = 0
    outside_period = 0
    for source in json.loads(args.pool.read_text())["choices"]:
        h = int(source["h"])
        p = int(source["p"])
        if args.period % h:
            outside_period += 1
            continue
        residue, modulus = power_target_congruence(h, p, args.power)
        if residue % modulus:
            incompatible += 1
            continue
        rows.append(
            {
                "h": h,
                "p": p,
                "a": int(source["a"]),
                "b": int(source["b"]),
                "c": 0,
                "ord2": int(source["ord2"]),
                "ord3": int(source["ord3"]),
            }
        )
    if not rows:
        raise RuntimeError("zero-phase pool is empty")

    algebraic_primes = tuple(
        prime for prime in exact_uncovered.factor(args.power) if prime % 2
    )
    sophie_germain = args.power % 4 == 0
    missed, meta = exact_uncovered.find_uncovered(
        rows,
        max_component=args.max_component,
        limit=args.limit,
        algebraic_primes=algebraic_primes,
        sophie_germain=sophie_germain,
    )
    payload = {
        "pool": str(args.pool),
        "power": args.power,
        "period": args.period,
        "rows": len(rows),
        "density": sum(1 / int(row["h"]) for row in rows),
        "outside_period": outside_period,
        "zero_incompatible": incompatible,
        "checker": meta,
        "uncovered": missed,
    }
    if args.output:
        args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 1 if missed else 0


if __name__ == "__main__":
    raise SystemExit(main())
