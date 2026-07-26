#!/usr/bin/env python3
"""Filter exponent-pair samples already covered by power identities."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import exact_uncovered


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--split-index", type=int)
    parser.add_argument("--first-output", type=Path)
    parser.add_argument("--second-output", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.power < 1:
        raise SystemExit("--power must be positive")
    source = [
        (int(raw_k), int(raw_l))
        for raw_k, raw_l in json.loads(args.input.read_text())
    ]
    if args.split_index is not None and not 0 <= args.split_index <= len(source):
        raise SystemExit("--split-index is outside the input point list")
    if bool(args.first_output) != bool(args.second_output):
        raise SystemExit("--first-output and --second-output must be used together")
    if args.first_output and args.split_index is None:
        raise SystemExit("partition outputs require --split-index")
    algebraic_primes = tuple(
        prime for prime in exact_uncovered.factor(args.power) if prime % 2
    )
    sophie_germain = args.power % 4 == 0

    kept = []
    kept_before_split = 0
    for index, point in enumerate(source):
        k, l = point
        if any(k % prime == 0 and l % prime == 0 for prime in algebraic_primes):
            continue
        if sophie_germain and k % 4 == 2 and l % 4 == 0:
            continue
        kept.append(point)
        if args.split_index is not None and index < args.split_index:
            kept_before_split += 1
    args.output.write_text(json.dumps(kept) + "\n")
    if args.first_output and args.second_output:
        args.first_output.write_text(json.dumps(kept[:kept_before_split]) + "\n")
        args.second_output.write_text(json.dumps(kept[kept_before_split:]) + "\n")
    print(
        f"input={len(source)} kept={len(kept)} "
        f"removed={len(source) - len(kept)} "
        f"algebraic_primes={algebraic_primes} "
        f"sophie_germain={sophie_germain} "
        f"split_index={kept_before_split if args.split_index is not None else None} "
        f"output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
