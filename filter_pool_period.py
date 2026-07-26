#!/usr/bin/env python3
"""Retain derived-pool rows whose moduli divide a declared period."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.period < 1:
        raise SystemExit("--period must be positive")
    payload = json.loads(args.pool.read_text())
    input_rows = payload["choices"]
    rows = [
        row for row in input_rows
        if args.period % int(row["h"]) == 0
    ]
    result = dict(payload)
    result["source"] = str(args.pool)
    result["period_filter"] = args.period
    result["choices"] = rows
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"input={len(input_rows)} output={len(rows)} "
        f"density={sum(1 / int(row['h']) for row in rows):.12f} "
        f"period={args.period}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
