#!/usr/bin/env python3
"""Memory-light exact greedy covering search for Erdos problem 203.

The exact checker returns genuine exponent classes missed by all lines chosen
so far.  Each round greedily spends unused prime fibres to hit the whole batch,
then rebuilds the checker on the enlarged CRT system.  Checker UNSAT is a real
finite cover; every intermediate result is only exploratory.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import cegis_cover
import exact_uncovered
import search_cover

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def load_candidates(
    path: Path, derived_pool: bool = False
) -> list[tuple[int, int, int, int, int, int]]:
    source = json.loads(path.read_text())
    candidates = []
    for raw in source["choices"]:
        p = int(raw["p"])
        if derived_pool:
            h = int(raw["h"])
            a = int(raw["a"]) % h
            b = int(raw["b"]) % h
            ord2 = int(raw.get("ord2", h))
            ord3 = int(raw.get("ord3", h))
            if math.gcd(a, b, h) != 1:
                raise RuntimeError(f"non-surjective derived row for prime {p}")
        else:
            ord2 = search_cover.multiplicative_order(2, p)
            ord3 = search_cover.multiplicative_order(3, p)
            h, a, b, ord2, ord3 = cegis_cover.signature(p, ord2, ord3)
        candidates.append((h, p, a, b, ord2, ord3))
    candidates.sort()
    return candidates


def as_row(item: tuple[int, int, int, int, int, int], c: int) -> dict:
    h, p, a, b, ord2, ord3 = item
    return {
        "h": h,
        "p": p,
        "a": a,
        "b": b,
        "ord2": ord2,
        "ord3": ord3,
        "c": c,
    }


def power_compatible(
    item: tuple[int, int, int, int, int, int], c: int, power: int
) -> bool:
    h, p, a, b, ord2, ord3 = item
    divisor = math.gcd(power, p - 1)
    exponent = (p - 1) // 2 - ((p - 1) // h) * c
    return exponent % divisor == 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_pool", type=Path)
    parser.add_argument("--batch", type=int, default=100)
    parser.add_argument("--rounds", type=int, default=10000)
    parser.add_argument("--max-component", type=int, default=100000)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--power", type=int, default=1)
    parser.add_argument("--derived-pool", action="store_true")
    parser.add_argument("--normalize-primes")
    parser.add_argument("--checker-solver", default="cadical195")
    args = parser.parse_args()
    if args.power < 1 or args.power % 2 == 0:
        raise SystemExit("--power must be a positive odd integer")
    algebraic_primes = tuple(exact_uncovered.factor(args.power))

    candidates = load_candidates(args.candidate_pool, args.derived_pool)
    by_prime = {item[1]: item for item in candidates}
    if args.normalize_primes:
        normalize_primes = tuple(
            int(value) for value in args.normalize_primes.split(",") if value
        )
    elif args.derived_pool:
        normalize_primes = ()
    else:
        normalize_primes = (5, 7)
    if normalize_primes:
        anchors = [(by_prime[p][0], by_prime[p][2], by_prime[p][3]) for p in normalize_primes]
        period = math.lcm(*(row[0] for row in anchors))
        image_size = len(
            {
                tuple((a * k + b * l) % h for h, a, b in anchors)
                for k in range(period)
                for l in range(period)
            }
        )
        if image_size != math.prod(row[0] for row in anchors):
            raise RuntimeError("declared normalization maps are not jointly surjective")
    selected: list[dict] = []
    used: set[int] = set()
    start_round = 0
    if args.checkpoint and args.checkpoint.exists():
        saved = json.loads(args.checkpoint.read_text())
        start_round = int(saved.get("round", 0))
        for row in saved["choices"]:
            p = int(row["p"])
            selected.append(as_row(by_prime[p], int(row["c"])))
            used.add(p)
    else:
        for p in normalize_primes:
            if p in by_prime:
                selected.append(as_row(by_prime[p], 0))
                used.add(p)

    for round_no in range(start_round + 1, args.rounds + 1):
        missed, meta = exact_uncovered.find_uncovered(
            selected,
            max_component=args.max_component,
            limit=args.batch,
            algebraic_primes=algebraic_primes,
            solver_name=args.checker_solver,
        )
        print(
            f"round={round_no} selected={len(selected)} "
            f"checker={'SAT' if missed else 'UNSAT'} misses={len(missed)}",
            flush=True,
        )
        if not missed:
            payload = {
                "candidate_pool": str(args.candidate_pool),
                "derived_pool": args.derived_pool,
                "choices": selected,
                "power": args.power,
                "checker": meta,
            }
            if args.output:
                args.output.write_text(json.dumps(payload, indent=2) + "\n")
            else:
                print(json.dumps(payload, indent=2))
            return 0

        remaining = set(range(len(missed)))
        while remaining:
            best = None
            for item in candidates:
                h, p, a, b, ord2, ord3 = item
                if p in used:
                    continue
                buckets: dict[int, list[int]] = {}
                for index in remaining:
                    k, l = missed[index]
                    c = (a * k + b * l) % h
                    if not power_compatible(item, c, args.power):
                        continue
                    buckets.setdefault(c, []).append(index)
                if not buckets:
                    continue
                c, hit = max(buckets.items(), key=lambda pair: len(pair[1]))
                score = (len(hit), -h, -p)
                if best is None or score > best[0]:
                    best = (score, item, c, hit)
            if best is None:
                print("EXHAUSTED candidate fibres without a cover", flush=True)
                return 2
            _, item, c, hit = best
            selected.append(as_row(item, c))
            used.add(item[1])
            remaining.difference_update(hit)

        if args.checkpoint:
            args.checkpoint.write_text(
                json.dumps({"round": round_no, "choices": selected}, indent=2)
                + "\n"
            )

    print(f"ROUND_LIMIT selected={len(selected)}", flush=True)
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
