#!/usr/bin/env python3
"""Filter derived affine-cover rows by component divisibility."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--component", type=int, required=True)
    parser.add_argument(
        "--mode",
        choices=("independent", "incident"),
        required=True,
        help=(
            "retain rows whose modulus is respectively not divisible or "
            "divisible by the component"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.component < 2:
        raise SystemExit("component must be at least two")
    payload = json.loads(args.pool.read_text())
    source_rows = payload["choices"]
    if args.mode == "independent":
        retained = [
            row
            for row in source_rows
            if int(row["h"]) % args.component != 0
        ]
    else:
        retained = [
            row
            for row in source_rows
            if int(row["h"]) % args.component == 0
        ]

    result = dict(payload)
    result["source_pool"] = str(args.pool)
    result["component_filter"] = {
        "component": args.component,
        "mode": args.mode,
    }
    result["source_choice_count"] = len(source_rows)
    result["retained_choice_count"] = len(retained)
    result["choices"] = retained
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"component={args.component} mode={args.mode} "
        f"source={len(source_rows)} retained={len(retained)} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
