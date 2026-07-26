#!/usr/bin/env python3
"""Discover primes of small <2,3> subgroup order by factoring common powers."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import cegis_cover
import search_cover


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-h", type=int, required=True)
    parser.add_argument("--factor-limit", type=int, default=100_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import gmpy2  # type: ignore
    import sympy  # type: ignore

    primes: set[int] = set()
    unresolved = 0
    two = gmpy2.mpz(1)
    three = gmpy2.mpz(1)
    for exponent in range(1, args.max_h + 1):
        two *= 2
        three *= 3
        common = int(gmpy2.gcd(two - 1, three - 1))
        for p in tuple(primes):
            while common % p == 0:
                common //= p
        if common > 1:
            factors = sympy.factorint(common, limit=args.factor_limit)
            for factor in factors:
                value = int(factor)
                if sympy.isprime(value):
                    primes.add(value)
                else:
                    unresolved += 1
        if exponent % 100 == 0:
            print(
                f"exponent={exponent} primes={len(primes)} unresolved={unresolved}",
                flush=True,
            )

    candidates = []
    for p in primes:
        if p in (2, 3):
            continue
        ord2 = search_cover.multiplicative_order(2, p)
        ord3 = search_cover.multiplicative_order(3, p)
        h, a, b, ord2, ord3 = cegis_cover.signature(p, ord2, ord3)
        if h <= args.max_h:
            candidates.append((h, p, a, b, ord2, ord3))
    candidates.sort()
    payload = {
        "max_h": args.max_h,
        "factor_limit": args.factor_limit,
        "unresolved": unresolved,
        "choices": [
            {
                "h": h,
                "p": p,
                "a": a,
                "b": b,
                "ord2": ord2,
                "ord3": ord3,
                "c": 0,
            }
            for h, p, a, b, ord2, ord3 in candidates
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"wrote={args.output} candidates={len(candidates)} "
        f"density={sum(1 / row[0] for row in candidates):.12f} "
        f"max_prime={max(row[1] for row in candidates)} unresolved={unresolved}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
