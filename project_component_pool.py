#!/usr/bin/env python3
"""Project compatible CRT fibre rows onto one component modulus."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--component", type=int, required=True)
    parser.add_argument(
        "--period",
        type=int,
        help="retain source rows whose modulus divides this period",
    )
    parser.add_argument(
        "--all-incident",
        action="store_true",
        help="project every source row whose modulus contains the component",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    q = args.component
    if q < 2:
        raise SystemExit("component must be at least two")
    if args.all_incident == (args.period is not None):
        raise SystemExit("choose exactly one of --period and --all-incident")
    if args.period is not None and args.period % q:
        raise SystemExit("component must divide the requested period")

    payload = json.loads(args.pool.read_text())
    projected = []
    for row in payload["choices"]:
        h = int(row["h"])
        if h <= 1 or h % q:
            continue
        if args.period is not None and args.period % h:
            continue
        modulus = int(row["target_modulus"])
        if modulus != 1:
            raise RuntimeError(
                f"row {row['p']} has constrained target modulus {modulus}"
            )
        item = dict(row)
        item["source_h"] = h
        item["h"] = q
        item["a"] = int(row["a"]) % q
        item["b"] = int(row["b"]) % q
        item["c"] = int(row.get("c", 0)) % q
        item["target_residue"] = 0
        item["target_modulus"] = 1
        projected.append(item)

    result = {
        "source_pool": str(args.pool),
        "projection_component": q,
        "source_period": args.period,
        "all_incident": args.all_incident,
        "choices": projected,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"component={q} period={args.period} "
        f"all_incident={args.all_incident} rows={len(projected)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
