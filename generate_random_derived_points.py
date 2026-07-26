#!/usr/bin/env python3
"""Generate reproducible admissible random points for a derived pool."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path


def integer(text: str) -> int:
    return int(text, 0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--bits", type=int, default=64)
    parser.add_argument("--seed", type=integer, default=0x203)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.count < 1 or args.bits < 1:
        raise SystemExit("--count and --bits must be positive")

    payload = json.loads(args.pool.read_text())
    algebraic_primes = tuple(
        int(prime) for prime in payload.get("algebraic_primes", ())
    )
    sophie_germain = bool(payload.get("sophie_germain", False))
    rng = random.Random(args.seed)
    points = []
    seen = set()
    draws = 0
    while len(points) < args.count:
        k = rng.getrandbits(args.bits)
        l = rng.getrandbits(args.bits)
        draws += 1
        if any(
            k % prime == 0 and l % prime == 0
            for prime in algebraic_primes
        ):
            continue
        if sophie_germain and k % 4 == 2 and l % 4 == 0:
            continue
        point = (k, l)
        if point in seen:
            continue
        seen.add(point)
        points.append([k, l])
    args.output.write_text(json.dumps(points) + "\n")
    print(
        f"points={len(points)} draws={draws} bits={args.bits} "
        f"seed={args.seed} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
