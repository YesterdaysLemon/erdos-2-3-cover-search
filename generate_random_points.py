#!/usr/bin/env python3
"""Generate deterministic random nonnegative lattice points for SAT screens."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--bits", type=int, default=128)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.count < 1:
        raise SystemExit("--count must be positive")
    if args.bits < 1:
        raise SystemExit("--bits must be positive")

    rng = random.Random(args.seed)
    points = [
        [rng.getrandbits(args.bits), rng.getrandbits(args.bits)]
        for _ in range(args.count)
    ]
    args.output.write_text(json.dumps(points) + "\n")
    print(
        f"points={len(points)} bits={args.bits} seed={args.seed} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
