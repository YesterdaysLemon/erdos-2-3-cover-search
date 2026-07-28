#!/usr/bin/env python3
"""Merge exact integer point lists from JSON files without numeric coercion."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from certify_anchor_phase_quotient import load_cells


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "inputs",
        nargs="+",
        help="PATH or PATH:KEY for a point list stored under a top-level key",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    merged: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    source_counts = []
    for spec in args.inputs:
        path = Path(spec)
        separator = ""
        key = ""
        if not path.exists():
            path_text, separator, key = spec.rpartition(":")
            path = Path(path_text)
            if not separator or not path.exists():
                raise RuntimeError(f"point source does not exist: {spec}")
        if separator:
            payload = json.loads(path.read_text())
            if not isinstance(payload, dict) or key not in payload:
                raise RuntimeError(f"{path} has no top-level key {key!r}")
            payload = payload[key]
        else:
            payload = load_cells(path)
        count_before = len(merged)
        for raw_k, raw_l in payload:
            point = (int(raw_k), int(raw_l))
            if point in seen:
                continue
            seen.add(point)
            merged.append(point)
        source_counts.append(
            {
                "source": spec,
                "new_points": len(merged) - count_before,
            }
        )

    args.output.write_text(
        json.dumps([[k, ell] for k, ell in merged]) + "\n"
    )
    print(
        f"sources={len(args.inputs)} points={len(merged)} "
        f"output={args.output}"
    )
    for record in source_counts:
        print(
            f"source={record['source']} new_points={record['new_points']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
