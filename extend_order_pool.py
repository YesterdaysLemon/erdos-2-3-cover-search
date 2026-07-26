#!/usr/bin/env python3
"""Extend a completely factored order-first Erdos-203 candidate pool."""

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
    parser.add_argument("source", type=Path)
    parser.add_argument("--max-h", type=int, required=True)
    parser.add_argument("--factor-limit", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import gmpy2  # type: ignore
    import sympy  # type: ignore

    source = json.loads(args.source.read_text())
    start = int(source["max_h"])
    if int(source.get("unresolved", 0)):
        raise RuntimeError("source pool has unresolved cofactors")
    if args.max_h <= start:
        raise SystemExit("--max-h must exceed the source max_h")

    primes = {int(row["p"]) for row in source["choices"]}
    unresolved = 0
    two = gmpy2.mpz(2) ** start
    three = gmpy2.mpz(3) ** start
    for exponent in range(start + 1, args.max_h + 1):
        two *= 2
        three *= 3
        common = int(gmpy2.gcd(two - 1, three - 1))
        for prime in tuple(primes):
            while common % prime == 0:
                common //= prime
        if common > 1:
            factorization = sympy.factorint(common, limit=args.factor_limit)
            for factor in factorization:
                value = int(factor)
                if sympy.isprime(value):
                    primes.add(value)
                else:
                    unresolved += 1
        if exponent % 100 == 0:
            print(
                f"exponent={exponent} primes={len(primes)} "
                f"unresolved={unresolved}",
                flush=True,
            )

    candidates = []
    for prime in primes:
        if prime in (2, 3):
            continue
        ord2 = search_cover.multiplicative_order(2, prime)
        ord3 = search_cover.multiplicative_order(3, prime)
        h, a, b, ord2, ord3 = cegis_cover.signature(prime, ord2, ord3)
        if h <= args.max_h:
            candidates.append((h, prime, a, b, ord2, ord3))
    candidates.sort()
    payload = {
        "max_h": args.max_h,
        "factor_limit": args.factor_limit,
        "unresolved": unresolved,
        "extended_from": str(args.source),
        "choices": [
            {
                "h": h,
                "p": prime,
                "a": a,
                "b": b,
                "ord2": ord2,
                "ord3": ord3,
                "c": 0,
            }
            for h, prime, a, b, ord2, ord3 in candidates
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"wrote={args.output} candidates={len(candidates)} "
        f"density={sum(1 / row[0] for row in candidates):.12f} "
        f"max_prime={max(row[1] for row in candidates)} "
        f"unresolved={unresolved}",
        flush=True,
    )
    return 0 if unresolved == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
