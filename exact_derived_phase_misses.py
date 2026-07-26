#!/usr/bin/env python3
"""Run the exact modular checker on a saved derived-pool phase map."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import exact_greedy
import exact_uncovered
import exact_uncovered_z3_bv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument(
        "--period",
        type=int,
        default=0,
        help="if nonzero, retain only rows whose modulus divides this period",
    )
    parser.add_argument("--max-component", type=int, default=256)
    parser.add_argument("--limit", type=int, default=1)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument(
        "--independent-z3",
        action="store_true",
        help="also replay the assignment with the independent Z3-BV checker",
    )
    parser.add_argument("--diversity-primes", default="")
    parser.add_argument("--diversity-coordinate-moduli", default="")
    parser.add_argument("--diversity-quota", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    candidates = exact_greedy.load_candidates(args.pool, True)
    if args.period:
        candidates = [
            candidate
            for candidate in candidates
            if args.period % int(candidate[0]) == 0
        ]
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(args.phase_file.read_text()).items()
    }
    row_by_prime = {
        int(row["p"]): row for row in payload["choices"]
    }
    rows = []
    for candidate in candidates:
        h, prime, *_rest = candidate
        source = row_by_prime[prime]
        residue = int(source["target_residue"])
        modulus = int(source["target_modulus"])
        target = phases.get(prime, residue) % h
        if target % modulus != residue:
            raise RuntimeError(f"phase for p={prime} violates target restriction")
        rows.append(exact_greedy.as_row(candidate, target))

    algebraic_primes = tuple(
        int(value) for value in payload.get("algebraic_primes", ())
    )
    diversity_primes = tuple(
        int(value)
        for value in args.diversity_primes.split(",")
        if value
    )
    diversity_coordinate_moduli = tuple(
        int(value)
        for value in args.diversity_coordinate_moduli.split(",")
        if value
    )
    misses, meta = exact_uncovered.find_uncovered(
        rows,
        max_component=args.max_component,
        limit=args.limit,
        algebraic_primes=algebraic_primes,
        solver_name=args.solver,
        diversity_primes=diversity_primes,
        diversity_coordinate_moduli=diversity_coordinate_moduli,
        diversity_quota=args.diversity_quota,
        sophie_germain=bool(payload.get("sophie_germain", False)),
    )
    independent_misses = []
    independent_meta = None
    if args.independent_z3:
        independent_misses, independent_meta = (
            exact_uncovered_z3_bv.find_uncovered(
                rows,
                max_component=args.max_component,
                limit=args.limit,
                algebraic_primes=algebraic_primes,
                sophie_germain=bool(
                    payload.get("sophie_germain", False)
                ),
            )
        )
    result = {
        "pool": str(args.pool),
        "phase_file": str(args.phase_file),
        "row_count": len(rows),
        "checker": meta,
        "misses": [[k, l] for k, l in misses],
        "independent_checker": independent_meta,
        "independent_misses": [
            [k, l] for k, l in independent_misses
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"rows={len(rows)} misses={len(misses)} sat={meta['sat']} "
        f"independent_misses={len(independent_misses)} "
        f"output={args.output}",
        flush=True,
    )
    return 1 if misses or independent_misses else 0


if __name__ == "__main__":
    raise SystemExit(main())
