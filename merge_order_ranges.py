#!/usr/bin/env python3
"""Merge a complete base pool with adjacent exact range scans."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("base", type=Path)
    parser.add_argument("ranges", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    base = json.loads(args.base.read_text())
    if int(base.get("unresolved", 0)):
        raise RuntimeError("base pool has unresolved cofactors")
    rows = {int(row["p"]): row for row in base["choices"]}
    expected_start = int(base["max_h"]) + 1
    intervals = []
    for path in args.ranges:
        payload = json.loads(path.read_text())
        if int(payload.get("unresolved", 0)):
            raise RuntimeError(f"range has unresolved cofactors: {path}")
        intervals.append((int(payload["start_h"]), int(payload["end_h"]), payload))
    intervals.sort()
    for start, end, payload in intervals:
        if start != expected_start:
            raise RuntimeError(
                f"range gap or overlap: expected {expected_start}, got {start}"
            )
        for row in payload["choices"]:
            prime = int(row["p"])
            if prime in rows and rows[prime] != row:
                raise RuntimeError(f"conflicting signature for prime {prime}")
            rows[prime] = row
        expected_start = end + 1
    choices = sorted(rows.values(), key=lambda row: (int(row["h"]), int(row["p"])))
    result = {
        "max_h": expected_start - 1,
        "factor_limit": min(
            [int(base.get("factor_limit", 0))]
            + [int(payload.get("factor_limit", 0)) for _, _, payload in intervals]
        ),
        "unresolved": 0,
        "merged_from": [str(args.base)] + [str(path) for path in args.ranges],
        "choices": choices,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"wrote={args.output} candidates={len(choices)} "
        f"density={sum(1 / int(row['h']) for row in choices):.12f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
