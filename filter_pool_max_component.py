#!/usr/bin/env python3
"""Retain affine-fibre rows below an exact prime-power component bound."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import exact_uncovered


def max_prime_power(value: int) -> int:
    return max(
        (
            prime**exponent
            for prime, exponent in exact_uncovered.factor(value).items()
        ),
        default=1,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--max-component", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.max_component < 1:
        raise SystemExit("--max-component must be positive")

    payload = json.loads(args.pool.read_text())
    source_rows = payload["choices"]
    retained = [
        row
        for row in source_rows
        if max_prime_power(int(row["h"])) <= args.max_component
    ]
    result = dict(payload)
    result["choices"] = retained
    result["component_filter"] = {
        "source": str(args.pool),
        "max_component": args.max_component,
        "input_rows": len(source_rows),
        "output_rows": len(retained),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"input={len(source_rows)} retained={len(retained)} "
        f"max_component={args.max_component} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
