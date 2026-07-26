#!/usr/bin/env python3
"""Concatenate JSON exponent-pair lists with stable deduplication."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    points = []
    seen = set()
    source_counts = []
    for path in args.inputs:
        added = 0
        for raw_k, raw_l in json.loads(path.read_text()):
            point = (int(raw_k), int(raw_l))
            if point in seen:
                continue
            seen.add(point)
            points.append(point)
            added += 1
        source_counts.append((str(path), added))
    args.output.write_text(json.dumps(points) + "\n")
    print(
        f"points={len(points)} sources={source_counts} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
