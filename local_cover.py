#!/usr/bin/env python3
"""Local-search optimizer for periodic Erdos-203 covering systems.

The objective is the number of uncovered sampled exponent pairs.  A zero on a
sample is only a candidate; exact verification is a separate required step.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

import cegis_cover
import search_cover


def get_candidates(period: int | None, prime_limit: int, count: int | None):
    raw = []
    for p in search_cover.sieve(prime_limit):
        ord2 = search_cover.multiplicative_order(2, p)
        ord3 = search_cover.multiplicative_order(3, p)
        h = math.lcm(ord2, ord3)
        if period is None or period % h == 0:
            raw.append((h, p, ord2, ord3))
    raw.sort()
    if count is not None:
        raw = raw[:count]
    out = []
    for h, p, ord2, ord3 in raw:
        h, a, b, ord2, ord3 = cegis_cover.signature(p, ord2, ord3)
        out.append((h, p, a, b, ord2, ord3))
    out.sort()
    return out


def get_complete_period_candidates(period: int):
    """Return every prime whose 2- and 3-orders both divide ``period``.

    Such primes are exactly the prime divisors of
    gcd(2**period - 1, 3**period - 1).  This is preferable to a prime cutoff
    whenever that gcd can be factored completely: a resulting UNSAT answer
    then rules out the whole declared period, not merely the scanned prefix.
    """
    import sympy  # Loaded from the task-local dependency directory by callers.

    try:
        import gmpy2  # type: ignore
    except ImportError:
        common = math.gcd((1 << period) - 1, pow(3, period) - 1)
    else:
        common = int(
            gmpy2.gcd(gmpy2.mpz(2) ** period - 1, gmpy2.mpz(3) ** period - 1)
        )

    factorization = sympy.factorint(common)
    reconstructed = math.prod(int(p) ** int(e) for p, e in factorization.items())
    if reconstructed != common or not all(sympy.isprime(p) for p in factorization):
        raise RuntimeError("period gcd was not factored completely")

    out = []
    for p_value in factorization:
        p = int(p_value)
        if p in (2, 3):
            continue
        ord2 = search_cover.multiplicative_order(2, p)
        ord3 = search_cover.multiplicative_order(3, p)
        h, a, b, ord2, ord3 = cegis_cover.signature(p, ord2, ord3)
        out.append((h, p, a, b, ord2, ord3))
    out.sort()
    return out, common, {int(p): int(e) for p, e in factorization.items()}


def optimize(
    candidates,
    points,
    rng: random.Random,
    restarts: int,
    sweeps: int,
    power: int,
):
    bucket_maps = []
    for h, p, a, b, ord2, ord3 in candidates:
        buckets = {}
        if power != 1:
            root = cegis_cover.primitive_root(p)
            generator = pow(root, (p - 1) // h, p)
            power_gcd = math.gcd(power, p - 1)
        for index, (k, l) in enumerate(points):
            c = (a * k + b * l) % h
            # The target g^c means m == -g^(-c) (mod p).  Since the
            # eventual m is M^power, retain only fibres for which this residue
            # is a power-th power modulo p.
            if power != 1:
                target = pow(generator, c, p)
                m_residue = -pow(target, -1, p) % p
                if pow(m_residue, (p - 1) // power_gcd, p) != 1:
                    continue
            buckets.setdefault(c, []).append(index)
        bucket_maps.append(buckets)

    best = None
    best_uncovered = len(points) + 1
    order = list(range(len(candidates)))
    for restart in range(restarts):
        # Random initial fibres, with the p=5 translation normalized to zero.
        choices = []
        cover_count = [0] * len(points)
        for item, buckets in zip(candidates, bucket_maps):
            h, p, a, b, ord2, ord3 = item
            if not buckets:
                choices.append(None)
                continue
            c = rng.choice(list(buckets))
            choices.append(c)
            for index in buckets.get(c, ()):
                cover_count[index] += 1

        for sweep in range(sweeps):
            before = sum(value == 0 for value in cover_count)
            rng.shuffle(order)
            for i in order:
                h, p, a, b, ord2, ord3 = candidates[i]
                if not bucket_maps[i]:
                    continue
                buckets = bucket_maps[i]
                old = choices[i]
                for index in buckets.get(old, ()):
                    cover_count[index] -= 1
                # After removing this prime, choose the fibre that covers the
                # largest number of currently uncovered sample points.
                best_score = -1
                best_targets = []
                for c, indices in buckets.items():
                    score = sum(cover_count[index] == 0 for index in indices)
                    if score > best_score:
                        best_score = score
                        best_targets = [c]
                    elif score == best_score:
                        best_targets.append(c)
                new = rng.choice(best_targets)
                choices[i] = new
                for index in buckets[new]:
                    cover_count[index] += 1
            after = sum(value == 0 for value in cover_count)
            if after < best_uncovered:
                best_uncovered = after
                best = choices[:]
                print(
                    f"restart={restart} sweep={sweep} "
                    f"uncovered={after}/{len(points)}"
                )
            if after == 0 or after >= before:
                break
        if best_uncovered == 0:
            break
    return best_uncovered, best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--period", type=int)
    parser.add_argument("--count", type=int)
    parser.add_argument("--sample-range", type=int, default=10**12)
    parser.add_argument("--prime-limit", type=int, default=400_000)
    parser.add_argument("--points", type=int, default=20_000)
    parser.add_argument("--restarts", type=int, default=20)
    parser.add_argument("--sweeps", type=int, default=50)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--quiet-choices", action="store_true")
    parser.add_argument(
        "--power",
        type=int,
        default=1,
        help="search for m=M^power; power should be odd and divide the period",
    )
    args = parser.parse_args()
    rng = random.Random(203)
    if args.period is None and args.count is None:
        args.period = 55_440
    candidates = get_candidates(args.period, args.prime_limit, args.count)
    if args.power % 2 == 0 or (
        args.period is not None and args.period % args.power != 0
    ):
        raise SystemExit("--power must be odd and, when given, divide --period")
    sample_range = args.period or args.sample_range
    points = []
    while len(points) < args.points:
        point = (rng.randrange(sample_range), rng.randrange(sample_range))
        # If gcd(k,l,power)>1, an odd divisor d gives X^d+1 directly.
        if math.gcd(point[0], point[1], args.power) > 1:
            continue
        points.append(point)
    print(
        f"period={args.period} primes={len(candidates)} "
        f"density={sum(1/item[0] for item in candidates):.9f} power={args.power}"
    )
    uncovered, choices = optimize(
        candidates, points, rng, args.restarts, args.sweeps, args.power
    )
    print(f"best_uncovered={uncovered}")
    if choices is not None:
        if not args.quiet_choices:
            print("choices")
            for item, c in zip(candidates, choices):
                print((*item, c))
        if args.output:
            payload = {
                "period": args.period,
                "power": args.power,
                "sample_range": sample_range,
                "sample_points": len(points),
                "uncovered": uncovered,
                "choices": [
                    {
                        "h": item[0],
                        "p": item[1],
                        "a": item[2],
                        "b": item[3],
                        "ord2": item[4],
                        "ord3": item[5],
                        "c": c,
                    }
                    for item, c in zip(candidates, choices)
                ],
            }
            args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0 if uncovered == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
