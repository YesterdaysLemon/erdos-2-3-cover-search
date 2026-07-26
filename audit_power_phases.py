#!/usr/bin/env python3
"""Fast reproducible random audit of a perfect-power phase assignment."""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from pathlib import Path

import exact_greedy
import exact_uncovered
from power_anchor_capacity_lp import power_target_congruence


def integer(text: str) -> int:
    return int(text, 0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("phases", type=Path)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--period", type=int, default=0)
    parser.add_argument("--max-component", type=int, default=0)
    parser.add_argument("--draws", type=int, default=100000)
    parser.add_argument("--seed", type=integer, default=0x203)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    candidates = exact_greedy.load_candidates(args.pool, True)
    if args.period:
        candidates = [row for row in candidates if args.period % row[0] == 0]
    if args.max_component:
        candidates = [
            row
            for row in candidates
            if max(
                (
                    prime**exponent
                    for prime, exponent in exact_uncovered.factor(row[0]).items()
                ),
                default=1,
            )
            <= args.max_component
        ]
    compatible = []
    power_congruences = {}
    for row in candidates:
        h, p = row[:2]
        try:
            power_congruences[p] = power_target_congruence(h, p, args.power)
        except RuntimeError:
            continue
        compatible.append(row)
    candidates = compatible
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(args.phases.read_text()).items()
    }
    missing = [p for _h, p, *_rest in candidates if p not in phases]
    if missing:
        raise RuntimeError(f"phase assignment is missing {len(missing)} primes")
    invalid = [
        p
        for h, p, *_rest in candidates
        if phases[p] % h % power_congruences[p][1]
        != power_congruences[p][0]
    ]
    if invalid:
        raise RuntimeError(
            f"phase assignment has {len(invalid)} power-incompatible targets"
        )

    rng = random.Random(args.seed)
    ks = np.fromiter(
        (rng.getrandbits(64) for _ in range(args.draws)),
        dtype=np.uint64,
        count=args.draws,
    )
    ls = np.fromiter(
        (rng.getrandbits(64) for _ in range(args.draws)),
        dtype=np.uint64,
        count=args.draws,
    )
    algebraic_primes = tuple(
        prime for prime in exact_uncovered.factor(args.power) if prime % 2
    )
    eligible = np.ones(args.draws, dtype=np.bool_)
    for prime in algebraic_primes:
        eligible &= (ks % prime != 0) | (ls % prime != 0)
    if args.power % 4 == 0:
        eligible &= (ks % 4 != 2) | (ls % 4 != 0)

    uncovered = eligible.copy()
    started = time.monotonic()
    for h, p, a, b, _ord2, _ord3 in candidates:
        indices = np.flatnonzero(uncovered)
        if not len(indices):
            break
        values = (
            a * (ks[indices] % h) + b * (ls[indices] % h)
        ) % h
        uncovered[indices[values == phases[p]]] = False
    tested = int(np.count_nonzero(eligible))
    missed = int(np.count_nonzero(uncovered))
    print(
        f"rows={len(candidates)} draws={args.draws} tested={tested} "
        f"missed={missed} rate={missed / tested:.15f} "
        f"seed={args.seed} seconds={time.monotonic() - started:.3f}"
    )
    return 1 if missed else 0


if __name__ == "__main__":
    raise SystemExit(main())
