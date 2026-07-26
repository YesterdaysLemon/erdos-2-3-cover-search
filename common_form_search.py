#!/usr/bin/env python3
"""Search for a common-form reduction of Erdos problem 203.

If p divides 2^B-3 then 2^k 3^l == 2^(k+B*l) (mod p).  A set of
congruence classes in the single variable n=k+B*l that covers every integer
therefore gives a complete two-dimensional CRT construction.
"""

from __future__ import annotations

import argparse
import math
from array import array

import search_cover


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-limit", type=int, default=2_000_000)
    parser.add_argument("--max-order", type=int, default=20_000)
    parser.add_argument("--max-b", type=int, default=2_000_000)
    parser.add_argument("--top", type=int, default=20)
    args = parser.parse_args()

    relations = []
    for p in search_cover.sieve(args.prime_limit):
        order = search_cover.multiplicative_order(2, p)
        if order > args.max_order:
            continue
        value = 1
        log3 = None
        for exponent in range(order):
            if value == 3:
                log3 = exponent
                break
            value = value * 2 % p
        if log3 is not None:
            relations.append((order, log3, p))

    score = array("d", [0.0]) * (args.max_b + 1)
    for order, log3, p in relations:
        start = log3 if log3 > 0 else order
        for B in range(start, args.max_b + 1, order):
            score[B] += 1.0 / order
    best = sorted(
        ((score[B], B) for B in range(1, args.max_b + 1)), reverse=True
    )[: args.top]
    print(f"relations={len(relations)}")
    for total, B in best:
        compatible = [item for item in relations if B % item[0] == item[1]]
        print(
            f"B={B} density={total:.12f} primes={len(compatible)} "
            f"rows={compatible}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
