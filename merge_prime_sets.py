#!/usr/bin/env python3
"""Merge integer-prime lists from JSON files and command-line additions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load_values(spec: str) -> list[int]:
    path_text, separator, key = spec.partition(":")
    payload = json.loads(Path(path_text).read_text())
    if separator:
        for component in key.split("."):
            payload = payload[component]
    if isinstance(payload, dict):
        payload = payload.get("primes", ())
    return [int(value) for value in payload]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="*", help="PATH or PATH:KEY")
    parser.add_argument("--add", default="", help="comma-separated integers")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    values = {
        value
        for spec in args.inputs
        for value in load_values(spec)
    }
    values.update(
        int(value) for value in args.add.split(",") if value
    )
    ordered = sorted(values)
    args.output.write_text(json.dumps(ordered, indent=2) + "\n")
    print(f"primes={len(ordered)} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
