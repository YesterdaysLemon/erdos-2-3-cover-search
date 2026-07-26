#!/usr/bin/env python3
"""Search for a finite covering-congruence certificate for Erdos problem 203.

For an odd prime p, fixing m modulo p covers precisely one fibre of

    (k, l) |-> 2**k * 3**l (mod p).

If ord_p(2) divides A and ord_p(3) divides B, this fibre is periodic on the
A-by-B exponent torus.  The SAT instance below chooses at most one fibre per
prime and requires every torus point to be covered.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path


def sieve(limit: int) -> list[int]:
    is_prime = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        is_prime[0] = 0
    if limit >= 1:
        is_prime[1] = 0
    for n in range(2, math.isqrt(limit) + 1):
        if is_prime[n]:
            is_prime[n * n : limit + 1 : n] = b"\x00" * (
                (limit - n * n) // n + 1
            )
    return [n for n in range(5, limit + 1) if is_prime[n]]


def divisors(n: int) -> list[int]:
    small: list[int] = []
    large: list[int] = []
    for d in range(1, math.isqrt(n) + 1):
        if n % d == 0:
            small.append(d)
            if d * d != n:
                large.append(n // d)
    return small + large[::-1]


def multiplicative_order(a: int, p: int) -> int:
    # Factor the fixed group order before shrinking ``order``.  Using the
    # changing order as the trial-division bound can skip a remaining prime
    # factor (for example ord_34511(3) is 17, not 17*29).
    order = p - 1
    remaining = p - 1
    factors = []
    d = 2
    while d * d <= remaining:
        if remaining % d == 0:
            factors.append(d)
            while remaining % d == 0:
                remaining //= d
        d += 1 if d == 2 else 2
    if remaining > 1:
        factors.append(remaining)
    for prime in factors:
        while order % prime == 0 and pow(a, order // prime, p) == 1:
            order //= prime
    return order


@dataclass(frozen=True)
class PrimeData:
    p: int
    ord2: int
    ord3: int
    values: tuple[int, ...]


def candidate_primes(A: int, B: int, limit: int) -> list[PrimeData]:
    out: list[PrimeData] = []
    for p in sieve(limit):
        if pow(2, A, p) != 1 or pow(3, B, p) != 1:
            continue
        ord2 = multiplicative_order(2, p)
        ord3 = multiplicative_order(3, p)
        # F_p^* is cyclic, so <2, 3> has order lcm(ord_p(2), ord_p(3)).
        # Generate it without the potentially quadratic Cartesian product.
        vals = {1}
        frontier = [1]
        while frontier:
            value = frontier.pop()
            for generator in (2, 3):
                nxt = value * generator % p
                if nxt not in vals:
                    vals.add(nxt)
                    frontier.append(nxt)
        out.append(PrimeData(p, ord2, ord3, tuple(sorted(vals))))
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("A", type=int, help="period in the exponent k")
    parser.add_argument("B", type=int, help="period in the exponent l")
    parser.add_argument("--prime-limit", type=int, default=2_000_000)
    parser.add_argument("--timeout-ms", type=int, default=60_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from z3 import Bool, If, Or, PbEq, Solver, Sum, sat  # type: ignore

    A, B = args.A, args.B
    data = candidate_primes(A, B, args.prime_limit)
    print(
        f"period={A}x{B} prime_limit={args.prime_limit} "
        f"candidate_primes={len(data)}"
    )
    print(
        "primes="
        + ",".join(
            f"{d.p}[{d.ord2},{d.ord3};{len(d.values)}]" for d in data
        )
    )

    solver = Solver()
    solver.set(timeout=args.timeout_ms)
    choose = {
        (d.p, value): Bool(f"p{d.p}_v{value}")
        for d in data
        for value in d.values
    }
    for d in data:
        solver.add(PbEq([(choose[d.p, value], 1) for value in d.values], 1))

    # "Exactly one" above is harmless: an otherwise unused prime can choose an
    # arbitrary fibre.  It also reduces symmetry compared with an optional prime.
    for k in range(A):
        two = pow(2, k)
        for l in range(B):
            covering = []
            for d in data:
                value = (pow(2, k, d.p) * pow(3, l, d.p)) % d.p
                covering.append(choose[d.p, value])
            if not covering:
                print("UNSAT: no candidate prime covers even one torus point")
                return 1
            solver.add(Or(covering))

    result = solver.check()
    print(f"result={result}")
    if result != sat:
        return 2

    model = solver.model()
    selected = []
    for d in data:
        for value in d.values:
            if bool(model.eval(choose[d.p, value])):
                residue = (-pow(value, -1, d.p)) % d.p
                selected.append(
                    {
                        "prime": d.p,
                        "target": value,
                        "m_residue": residue,
                        "ord2": d.ord2,
                        "ord3": d.ord3,
                    }
                )
                break

    certificate = {
        "A": A,
        "B": B,
        "prime_limit": args.prime_limit,
        "congruences": selected,
    }
    print(json.dumps(certificate, indent=2))
    if args.output:
        args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
