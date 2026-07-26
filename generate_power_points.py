#!/usr/bin/env python3
"""Generate reproducible non-algebraic exponent pairs for CEGIS retention."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import exact_uncovered


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--bits", type=int, default=64)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.count < 1:
        raise SystemExit("--count must be positive")
    if args.bits < 16:
        raise SystemExit("--bits must be at least 16")

    algebraic_primes = tuple(
        prime for prime in exact_uncovered.factor(args.power) if prime % 2
    )
    sophie_germain = args.power % 4 == 0
    rng = random.Random(args.seed)
    points = []
    seen = set()
    while len(points) < args.count:
        point = (rng.getrandbits(args.bits), rng.getrandbits(args.bits))
        if point in seen:
            continue
        if any(
            point[0] % prime == 0 and point[1] % prime == 0
            for prime in algebraic_primes
        ):
            continue
        if sophie_germain and point[0] % 4 == 2 and point[1] % 4 == 0:
            continue
        seen.add(point)
        points.append(point)
    args.output.write_text(json.dumps(points) + "\n")
    print(
        f"points={len(points)} bits={args.bits} power={args.power} "
        f"algebraic_primes={algebraic_primes} "
        f"sophie_germain={sophie_germain} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
