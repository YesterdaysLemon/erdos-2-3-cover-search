#!/usr/bin/env python3
"""Copy a phase map while replacing explicitly named prime targets."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("--set", required=True, dest="assignments")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    phases = json.loads(args.source.read_text())
    replacements = {}
    for item in args.assignments.split(","):
        prime, target = item.split("=", 1)
        replacements[str(int(prime))] = int(target)
    missing = [prime for prime in replacements if prime not in phases]
    if missing:
        raise RuntimeError(f"source phase map is missing primes {missing}")
    phases.update(replacements)
    args.output.write_text(json.dumps(phases) + "\n")
    print(
        f"source={args.source} replacements={replacements} "
        f"output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
