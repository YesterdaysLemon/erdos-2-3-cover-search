#!/usr/bin/env python3
"""Generate reproducible random power-compatible phases for a period pool."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from power_anchor_capacity_lp import power_target_congruence


class Lcg32:
    """Small cross-language reproducible unsigned 32-bit generator."""

    def __init__(self, seed: int):
        self.state = seed & 0xFFFFFFFF

    def randrange(self, stop: int) -> int:
        if stop <= 0:
            raise ValueError("empty range")
        self.state = (
            1664525 * self.state + 1013904223
        ) & 0xFFFFFFFF
        return self.state % stop


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument(
        "--engine",
        choices=("python", "lcg32"),
        default="python",
        help="random-number engine (lcg32 is reproducible across languages)",
    )
    parser.add_argument("--normalize-primes", default="")
    parser.add_argument(
        "--skip-normalizer-draws",
        action="store_true",
        help=(
            "do not advance the random engine for primes whose target is "
            "forced to zero"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    normalizers = {
        int(value) for value in args.normalize_primes.split(",") if value
    }
    rng = (
        random.Random(args.seed)
        if args.engine == "python"
        else Lcg32(args.seed)
    )
    phases = {}
    density = 0.0
    rows = 0
    for row in json.loads(args.pool.read_text())["choices"]:
        h = int(row["h"])
        p = int(row["p"])
        if args.period % h:
            continue
        residue, modulus = power_target_congruence(h, p, args.power)
        count = (h - 1 - residue) // modulus + 1
        if p in normalizers and args.skip_normalizer_draws:
            if residue:
                raise RuntimeError(
                    f"normalization target zero is invalid for prime {p}"
                )
            target = 0
        else:
            target = residue + modulus * rng.randrange(count)
            if p in normalizers:
                if residue:
                    raise RuntimeError(
                        f"normalization target zero is invalid for prime {p}"
                    )
                target = 0
        phases[str(p)] = target
        density += 1 / h
        rows += 1
    missing = normalizers - {int(prime) for prime in phases}
    if missing:
        raise RuntimeError(f"normalization primes missing from pool: {missing}")
    args.output.write_text(json.dumps(phases) + "\n")
    print(
        f"rows={rows} density={density:.12f} seed={args.seed} "
        f"engine={args.engine} "
        f"skip_normalizer_draws={args.skip_normalizer_draws} "
        f"normalizers={sorted(normalizers)} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
