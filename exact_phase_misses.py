#!/usr/bin/env python3
"""Extract exact CRT counterexamples for a saved power-compatible phase map."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import exact_uncovered
from power_anchor_capacity_lp import power_target_congruence


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument(
        "--row-max-component",
        type=int,
        default=0,
        help="if nonzero, retain rows whose prime-power components are at most this",
    )
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument("--max-component", type=int, default=256)
    parser.add_argument("--diversity-primes", default="")
    parser.add_argument("--diversity-quota", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    phases = json.loads(args.phase_file.read_text())
    rows = []
    for raw in source["choices"]:
        h = int(raw["h"])
        p = int(raw["p"])
        if args.period % h:
            continue
        if args.row_max_component and max(
            (
                prime**exponent
                for prime, exponent in exact_uncovered.factor(h).items()
            ),
            default=1,
        ) > args.row_max_component:
            continue
        try:
            residue, modulus = power_target_congruence(h, p, args.power)
        except RuntimeError:
            continue
        c = int(phases.get(str(p), residue)) % h
        if c % modulus != residue:
            raise RuntimeError(f"saved target for prime {p} is incompatible")
        row = {
            key: int(raw[key])
            for key in ("h", "p", "a", "b", "ord2", "ord3")
        }
        row["c"] = c
        rows.append(row)

    algebraic_primes = tuple(
        prime
        for prime in exact_uncovered.factor(args.power)
        if prime % 2
    )
    diversity_primes = tuple(
        int(value)
        for value in args.diversity_primes.split(",")
        if value
    )
    misses, meta = exact_uncovered.find_uncovered(
        rows,
        max_component=args.max_component,
        limit=args.limit,
        diversity_primes=diversity_primes,
        diversity_quota=args.diversity_quota,
        algebraic_primes=algebraic_primes,
        sophie_germain=args.power % 4 == 0,
    )
    args.output.write_text(json.dumps(misses) + "\n")
    print(
        f"rows={len(rows)} misses={len(misses)} "
        f"sat={meta['sat']} components={len(meta['components'])} "
        f"period={meta['period']} output={args.output}"
    )
    return 1 if misses else 0


if __name__ == "__main__":
    raise SystemExit(main())
