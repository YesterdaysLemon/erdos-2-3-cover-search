#!/usr/bin/env python3
"""Extract a top-level JSON point-list key as a standalone point list."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("--key", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.source.read_text())
    points = [[int(k), int(l)] for k, l in payload[args.key]]
    args.output.write_text(json.dumps(points) + "\n")
    print(
        f"source={args.source} key={args.key} points={len(points)} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
