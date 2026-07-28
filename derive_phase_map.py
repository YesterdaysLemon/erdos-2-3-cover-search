#!/usr/bin/env python3
"""Derive one phase map by applying explicit prime-target updates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_update(raw: str) -> tuple[int, int]:
    try:
        prime_text, target_text = raw.split("=", 1)
        prime = int(prime_text)
        target = int(target_text)
    except (ValueError, TypeError) as error:
        raise argparse.ArgumentTypeError(
            "updates must have the form PRIME=TARGET"
        ) from error
    if prime < 2 or target < 0:
        raise argparse.ArgumentTypeError(
            "PRIME must be at least 2 and TARGET must be nonnegative"
        )
    return prime, target


def derive_phase(
    base: dict[str, int],
    updates: list[tuple[int, int]],
) -> dict[str, int]:
    result = {str(int(prime)): int(target) for prime, target in base.items()}
    seen = set()
    for prime, target in updates:
        key = str(prime)
        if key not in result:
            raise KeyError(f"phase map omits prime {prime}")
        if prime in seen:
            raise ValueError(f"prime {prime} is updated more than once")
        seen.add(prime)
        result[key] = target
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("base", type=Path)
    parser.add_argument("--set", dest="updates", action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    updates = [parse_update(raw) for raw in args.updates]
    base = json.loads(args.base.read_text())
    derived = derive_phase(base, updates)
    args.output.write_text(json.dumps(derived) + "\n")
    print(
        f"base={args.base} updates={len(updates)} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
