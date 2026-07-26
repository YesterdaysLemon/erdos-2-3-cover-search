#!/usr/bin/env python3
"""Scan one independent subgroup-order interval for Erdos-203 fibres."""

from __future__ import annotations

import argparse
import collections
import json
import os
import sys
from pathlib import Path

import cegis_cover
import search_cover


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start-h", type=int, required=True)
    parser.add_argument("--end-h", type=int, required=True)
    parser.add_argument("--known-pool", type=Path)
    parser.add_argument("--factor-limit", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.start_h <= args.end_h:
        raise SystemExit("require 1 <= start-h <= end-h")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import gmpy2  # type: ignore
    import sympy  # type: ignore

    known_primes: set[int] = set()
    strip_schedule: dict[int, list[int]] = collections.defaultdict(list)
    if args.known_pool:
        source = json.loads(args.known_pool.read_text())
        if int(source.get("unresolved", 0)):
            raise RuntimeError("known pool has unresolved cofactors")
        for row in source["choices"]:
            prime = int(row["p"])
            h = int(row["h"])
            known_primes.add(prime)
            first = ((args.start_h + h - 1) // h) * h
            for exponent in range(first, args.end_h + 1, h):
                strip_schedule[exponent].append(prime)
    discovered: set[int] = set()
    discovered_rows: dict[int, tuple[int, int, int, int, int]] = {}
    unresolved_values: set[int] = set()
    two = gmpy2.mpz(2) ** (args.start_h - 1)
    three = gmpy2.mpz(3) ** (args.start_h - 1)
    for exponent in range(args.start_h, args.end_h + 1):
        two *= 2
        three *= 3
        common = int(gmpy2.gcd(two - 1, three - 1))
        # A prime p divides both 2^e-1 and 3^e-1 exactly when its stored
        # subgroup order h divides e.  Scheduling only those multiples is
        # equivalent to testing every known/discovered prime at every
        # exponent, but reduces the expected number of trial divisions in an
        # interval of width W from W*|P| to about
        # W*sum(1/h for p in P).
        for prime in strip_schedule.get(exponent, ()):
            while common % prime == 0:
                common //= prime
        if common > 1:
            factorization = sympy.factorint(common, limit=args.factor_limit)
            for factor in factorization:
                value = int(factor)
                if sympy.isprime(value):
                    if value in known_primes or value in discovered:
                        continue
                    discovered.add(value)
                    ord2 = search_cover.multiplicative_order(2, value)
                    ord3 = search_cover.multiplicative_order(3, value)
                    h, a, b, ord2, ord3 = cegis_cover.signature(
                        value, ord2, ord3
                    )
                    discovered_rows[value] = (h, a, b, ord2, ord3)
                    next_exponent = ((exponent // h) + 1) * h
                    for future in range(
                        next_exponent, args.end_h + 1, h
                    ):
                        strip_schedule[future].append(value)
                else:
                    unresolved_values.add(value)
        if exponent % 100 == 0:
            print(
                f"exponent={exponent} discovered={len(discovered)} "
                f"unresolved={len(unresolved_values)}",
                flush=True,
            )

    candidates = []
    for prime, (h, a, b, ord2, ord3) in discovered_rows.items():
        if prime in (2, 3):
            continue
        if args.start_h <= h <= args.end_h:
            candidates.append((h, prime, a, b, ord2, ord3))
    candidates.sort()
    payload = {
        "start_h": args.start_h,
        "end_h": args.end_h,
        "factor_limit": args.factor_limit,
        "known_pool": str(args.known_pool) if args.known_pool else None,
        "unresolved": len(unresolved_values),
        "unresolved_values": sorted(unresolved_values),
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
        f"unresolved={len(unresolved_values)}",
        flush=True,
    )
    return 0 if not unresolved_values else 2


if __name__ == "__main__":
    raise SystemExit(main())
