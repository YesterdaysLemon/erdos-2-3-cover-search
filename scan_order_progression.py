#!/usr/bin/env python3
"""Discover Erdos-203 fibres on a selected progression of exponents.

This is a construction-oriented scan, not a complete subgroup-order interval
scan.  For example, scanning only exponents divisible by 60 preferentially
finds rows whose modulus shrinks substantially after conditioning on a
period-60 cell.  Every emitted row is independently checkable, but absence
from the output is not an exclusion theorem.
"""

from __future__ import annotations

import argparse
import collections
import json
import math
import os
import sys
from pathlib import Path

import cegis_cover
import search_cover


def first_multiple_at_least(value: int, modulus: int) -> int:
    return ((value + modulus - 1) // modulus) * modulus


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start-exponent", type=int, required=True)
    parser.add_argument("--end-exponent", type=int, required=True)
    parser.add_argument("--step", type=int, required=True)
    parser.add_argument("--min-order", type=int, required=True)
    parser.add_argument("--max-order", type=int, required=True)
    parser.add_argument("--known-pool", type=Path, required=True)
    parser.add_argument("--factor-limit", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not (
        1 <= args.start_exponent <= args.end_exponent
        and args.step >= 1
        and 1 <= args.min_order <= args.max_order
    ):
        raise SystemExit("invalid exponent range, step, or order range")

    first_exponent = first_multiple_at_least(
        args.start_exponent, args.step
    )
    if first_exponent > args.end_exponent:
        raise SystemExit("the exponent interval contains no selected multiple")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import gmpy2  # type: ignore
    import sympy  # type: ignore

    source = json.loads(args.known_pool.read_text())
    if int(source.get("unresolved", 0)):
        raise RuntimeError("known pool has unresolved cofactors")

    known_primes: set[int] = set()
    strip_schedule: dict[int, list[int]] = collections.defaultdict(list)

    def schedule_prime(prime: int, order: int, after: int) -> None:
        period = math.lcm(order, args.step)
        first = first_multiple_at_least(after, period)
        for exponent in range(first, args.end_exponent + 1, period):
            strip_schedule[exponent].append(prime)

    for row in source["choices"]:
        prime = int(row["p"])
        order = int(row["h"])
        known_primes.add(prime)
        schedule_prime(prime, order, first_exponent)

    discovered: set[int] = set()
    discovered_rows: dict[int, tuple[int, int, int, int, int]] = {}
    unresolved_values: set[int] = set()
    safe_factorization_fallbacks = 0
    two_multiplier = gmpy2.mpz(2) ** args.step
    three_multiplier = gmpy2.mpz(3) ** args.step
    two = gmpy2.mpz(2) ** (first_exponent - args.step)
    three = gmpy2.mpz(3) ** (first_exponent - args.step)
    selected_count = 0

    for exponent in range(
        first_exponent, args.end_exponent + 1, args.step
    ):
        selected_count += 1
        two *= two_multiplier
        three *= three_multiplier
        common = int(gmpy2.gcd(two - 1, three - 1))
        for prime in strip_schedule.get(exponent, ()):
            while common % prime == 0:
                common //= prime
        if common > 1:
            try:
                factorization = sympy.factorint(
                    common, limit=args.factor_limit
                )
            except ValueError:
                # SymPy 1.14 can pass an incompletely factored composite
                # returned by a recursive Pollard call into its prime-only
                # factor cache.  Retry with deterministic trial division;
                # any leftover composite is recorded as unresolved below.
                sympy.ntheory.factor_.factor_cache.cache_clear()
                factorization = sympy.factorint(
                    common,
                    limit=args.factor_limit,
                    use_pm1=False,
                    use_rho=False,
                    use_ecm=False,
                )
                safe_factorization_fallbacks += 1
            for factor in factorization:
                value = int(factor)
                if not sympy.isprime(value):
                    unresolved_values.add(value)
                    continue
                if value in known_primes or value in discovered:
                    continue
                discovered.add(value)
                ord2 = search_cover.multiplicative_order(2, value)
                ord3 = search_cover.multiplicative_order(3, value)
                h, a, b, ord2, ord3 = cegis_cover.signature(
                    value, ord2, ord3
                )
                discovered_rows[value] = (h, a, b, ord2, ord3)
                schedule_prime(
                    value,
                    h,
                    exponent + args.step,
                )
        if selected_count % 25 == 0:
            print(
                f"exponent={exponent} selected={selected_count} "
                f"discovered={len(discovered)} "
                f"unresolved={len(unresolved_values)}",
                flush=True,
            )

    candidates = []
    for prime, (h, a, b, ord2, ord3) in discovered_rows.items():
        if (
            prime not in (2, 3)
            and args.min_order <= h <= args.max_order
        ):
            candidates.append((h, prime, a, b, ord2, ord3))
    candidates.sort()
    payload = {
        "scan_kind": "selected_exponent_progression",
        "complete_order_interval": False,
        "start_h": args.min_order,
        "end_h": args.max_order,
        "start_exponent": args.start_exponent,
        "end_exponent": args.end_exponent,
        "first_exponent": first_exponent,
        "step": args.step,
        "selected_exponents": selected_count,
        "factor_limit": args.factor_limit,
        "safe_factorization_fallbacks": safe_factorization_fallbacks,
        "known_pool": str(args.known_pool),
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
        f"wrote={args.output} selected={selected_count} "
        f"candidates={len(candidates)} "
        f"density={sum(1 / row[0] for row in candidates):.12f} "
        f"unresolved={len(unresolved_values)} "
        f"safe_factorization_fallbacks={safe_factorization_fallbacks}",
        flush=True,
    )
    return 0 if not unresolved_values else 2


if __name__ == "__main__":
    raise SystemExit(main())
