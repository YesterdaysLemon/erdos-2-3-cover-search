#!/usr/bin/env python3
"""Generate reproducible exponent points conditioned on weak anchor cells."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("profile", type=Path)
    parser.add_argument("--coverage-below", type=float, default=1.0)
    parser.add_argument("--per-cell", type=int, default=10000)
    parser.add_argument("--bits", type=int, default=64)
    parser.add_argument("--seed", type=int, default=203)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.per_cell < 1:
        raise SystemExit("--per-cell must be positive")
    if args.bits < 16 or args.bits > 64:
        raise SystemExit("--bits must lie in [16,64]")

    profile = json.loads(args.profile.read_text())
    pool = json.loads(Path(profile["pool"]).read_text())
    by_prime = {int(row["p"]): row for row in pool["choices"]}
    anchors = [by_prime[int(prime)] for prime in profile["anchor_primes"]]
    period = int(profile["anchor_period"])
    algebraic_primes = tuple(int(p) for p in profile["algebraic_primes"])
    selected = [
        record
        for record in profile["cells"]
        if float(record["coverage"]) < args.coverage_below
    ]
    if not selected:
        raise RuntimeError("profile has no cells below the requested threshold")

    preimages = {tuple(int(v) for v in record["cell"]): [] for record in selected}
    for k in range(period):
        for l in range(period):
            cell = tuple(
                (
                    int(row["a"]) * k + int(row["b"]) * l
                )
                % int(row["h"])
                for row in anchors
            )
            if cell in preimages:
                preimages[cell].append((k, l))
    missing = [cell for cell, values in preimages.items() if not values]
    if missing:
        raise RuntimeError(f"anchor cells have no period preimage: {missing}")

    rng = random.Random(args.seed)
    limit = 1 << args.bits
    points = []
    seen = set()
    counts = {}
    for record in selected:
        cell = tuple(int(v) for v in record["cell"])
        cell_points = 0
        while cell_points < args.per_cell:
            base_k, base_l = preimages[cell][
                rng.randrange(len(preimages[cell]))
            ]
            k = base_k + period * rng.randrange(
                (limit - 1 - base_k) // period + 1
            )
            l = base_l + period * rng.randrange(
                (limit - 1 - base_l) // period + 1
            )
            point = (k, l)
            if point in seen:
                continue
            if any(
                k % prime == 0 and l % prime == 0
                for prime in algebraic_primes
            ):
                continue
            seen.add(point)
            points.append(point)
            cell_points += 1
        counts[str(cell)] = cell_points
    args.output.write_text(json.dumps(points) + "\n")
    print(
        f"cells={len(selected)} points={len(points)} per_cell={args.per_cell} "
        f"threshold={args.coverage_below} counts={counts} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
