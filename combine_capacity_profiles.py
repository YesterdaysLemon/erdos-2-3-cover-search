#!/usr/bin/env python3
"""Union weak cells from compatible exact conditional-capacity profiles."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("profiles", type=Path, nargs="+")
    parser.add_argument("--below", type=float, default=1.1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.below <= 0:
        raise SystemExit("--below must be positive")

    payloads = [json.loads(path.read_text()) for path in args.profiles]
    compatibility_keys = (
        "power",
        "max_component",
        "anchor_primes",
        "anchor_targets",
        "anchor_period",
        "algebraic_primes",
    )
    reference = payloads[0]
    for path, payload in zip(args.profiles[1:], payloads[1:]):
        for key in compatibility_keys:
            if payload.get(key) != reference.get(key):
                raise RuntimeError(f"{path} has incompatible {key}")

    selected = {}
    source_counts = []
    for path, payload in zip(args.profiles, payloads):
        added = 0
        for record in payload["cells"]:
            if float(record["coverage"]) >= args.below:
                continue
            cell = tuple(int(value) for value in record["cell"])
            demand = int(record["demand"])
            if cell in selected:
                if selected[cell] != demand:
                    raise RuntimeError(f"inconsistent demand for cell {cell}")
                continue
            selected[cell] = demand
            added += 1
        source_counts.append((str(path), added))
    if not selected:
        raise RuntimeError("no weak cells selected")

    output = {
        key: reference[key]
        for key in (
            "pool",
            "power",
            "max_component",
            "anchor_primes",
            "anchor_targets",
            "anchor_period",
            "algebraic_primes",
        )
    }
    output["profile_union"] = {
        "below": args.below,
        "sources": [str(path) for path in args.profiles],
    }
    output["cells"] = [
        {"cell": list(cell), "coverage": 0.0, "demand": demand}
        for cell, demand in sorted(selected.items())
    ]
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"cells={len(selected)} below={args.below} "
        f"sources={source_counts} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
