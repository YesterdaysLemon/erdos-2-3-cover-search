#!/usr/bin/env python3
"""Translate a JSON affine-coordinate point list exactly."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("points", type=Path)
    parser.add_argument("--shift-k", type=int)
    parser.add_argument("--shift-l", type=int)
    parser.add_argument(
        "--normalization-audit",
        type=Path,
        help="read shift_k and shift_l from a normalization audit",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.normalization_audit:
        if args.shift_k is not None or args.shift_l is not None:
            raise SystemExit(
                "explicit shifts and --normalization-audit are exclusive"
            )
        audit = json.loads(args.normalization_audit.read_text())
        shift_k = int(audit["shift_k"])
        shift_l = int(audit["shift_l"])
    else:
        if args.shift_k is None or args.shift_l is None:
            raise SystemExit("both --shift-k and --shift-l are required")
        shift_k = args.shift_k
        shift_l = args.shift_l
    translated = []
    seen = set()
    for raw_k, raw_l in json.loads(args.points.read_text()):
        point = (int(raw_k) + shift_k, int(raw_l) + shift_l)
        if point not in seen:
            seen.add(point)
            translated.append([point[0], point[1]])
    args.output.write_text(json.dumps(translated) + "\n")
    print(
        f"input={args.points} points={len(translated)} "
        f"shift=({shift_k},{shift_l}) output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
