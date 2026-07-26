#!/usr/bin/env python3
"""Reproducible random coverage histogram for a derived affine pool."""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from collections import Counter
from pathlib import Path


def integer(text: str) -> int:
    return int(text, 0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("phases", type=Path)
    parser.add_argument("--draws", type=int, default=100_000)
    parser.add_argument("--batch", type=int, default=100_000)
    parser.add_argument("--seed", type=integer, default=0x203)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.draws < 1 or args.batch < 1:
        raise SystemExit("--draws and --batch must be positive")

    dependency_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dependency_path))
    import numpy as np  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    phases = {
        int(prime): int(value)
        for prime, value in json.loads(args.phases.read_text()).items()
    }
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        if prime not in phases:
            raise RuntimeError(f"phase file is missing prime {prime}")
        phase = phases[prime] % h
        if (
            phase % int(row.get("target_modulus", 1))
            != int(row.get("target_residue", 0))
        ):
            raise RuntimeError(f"illegal phase for prime {prime}")

    algebraic_primes = tuple(
        int(prime) for prime in payload.get("algebraic_primes", ())
    )
    sophie_germain = bool(payload.get("sophie_germain", False))
    rng = random.Random(args.seed)
    histogram = Counter()
    generated = 0
    eligible_total = 0
    started = time.monotonic()
    while generated < args.draws:
        count = min(args.batch, args.draws - generated)
        ks = np.fromiter(
            (rng.getrandbits(64) for _ in range(count)),
            dtype=np.uint64,
            count=count,
        )
        ls = np.fromiter(
            (rng.getrandbits(64) for _ in range(count)),
            dtype=np.uint64,
            count=count,
        )
        eligible = np.ones(count, dtype=np.bool_)
        for prime in algebraic_primes:
            eligible &= (ks % prime != 0) | (ls % prime != 0)
        if sophie_germain:
            eligible &= (ks % 4 != 2) | (ls % 4 != 0)
        indices = np.flatnonzero(eligible)
        eligible_total += len(indices)
        coverage = np.zeros(len(indices), dtype=np.uint16)
        batch_k = ks[indices]
        batch_l = ls[indices]
        for row in rows:
            h = int(row["h"])
            values = (
                int(row["a"]) * (batch_k % h)
                + int(row["b"]) * (batch_l % h)
            ) % h
            coverage += (
                values == phases[int(row["p"])] % h
            )
        unique, counts = np.unique(coverage, return_counts=True)
        for value, frequency in zip(unique, counts):
            histogram[int(value)] += int(frequency)
        generated += count

    uncovered = histogram[0]
    result = {
        "pool": str(args.pool),
        "phase_file": str(args.phases),
        "draws": args.draws,
        "eligible": eligible_total,
        "seed": args.seed,
        "row_count": len(rows),
        "coverage_histogram": {
            str(value): histogram[value] for value in sorted(histogram)
        },
        "uncovered": uncovered,
        "uncovered_rate": (
            uncovered / eligible_total if eligible_total else None
        ),
        "mean_coverage": (
            sum(value * count for value, count in histogram.items())
            / eligible_total
            if eligible_total
            else None
        ),
        "seconds": time.monotonic() - started,
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"rows={len(rows)} draws={args.draws} eligible={eligible_total} "
        f"uncovered={uncovered} "
        f"rate={result['uncovered_rate']:.12g} "
        f"mean={result['mean_coverage']:.12f} "
        f"seconds={result['seconds']:.3f}",
        flush=True,
    )
    return 1 if uncovered else 0


if __name__ == "__main__":
    raise SystemExit(main())
