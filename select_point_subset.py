#!/usr/bin/env python3
"""Write a deterministic prefix or seeded sample of a JSON point list."""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    points = json.loads(args.source.read_text())
    if not 0 <= args.count <= len(points):
        raise SystemExit("--count is outside the point list")
    if args.seed is None:
        selected = points[: args.count]
    else:
        selected = random.Random(args.seed).sample(points, args.count)
    args.output.write_text(json.dumps(selected) + "\n")
    print(
        f"source={len(points)} selected={len(selected)} "
        f"seed={args.seed} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
