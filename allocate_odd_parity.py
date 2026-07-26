#!/usr/bin/env python3
"""Allocate unused odd-order prime fibres disjointly among parity slices."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import cegis_cover
import search_cover


def induced(parent: str, h: int, a: int, b: int) -> tuple[int, int]:
    if parent == "k_even":
        return 2 * a % h, b % h
    if parent == "l_even":
        return a % h, 2 * b % h
    if parent == "sum_even":
        return (a + b) % h, 2 * b % h
    raise ValueError(parent)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("--k-base", type=Path, required=True)
    parser.add_argument("--l-base", type=Path, required=True)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument(
        "--pattern",
        default="kl",
        help="cyclic allocation pattern using k, l, and optionally s",
    )
    args = parser.parse_args()
    if not args.pattern or any(value not in "kls" for value in args.pattern):
        raise SystemExit("--pattern must contain only k, l, s")

    pools = {
        "k": json.loads(args.k_base.read_text()),
        "l": json.loads(args.l_base.read_text()),
        "s": {"choices": []},
    }
    odd = []
    for old in json.loads(args.source.read_text())["choices"]:
        p = int(old["p"])
        ord2 = search_cover.multiplicative_order(2, p)
        ord3 = search_cover.multiplicative_order(3, p)
        h, a, b, ord2, ord3 = cegis_cover.signature(p, ord2, ord3)
        if h % 2:
            odd.append((h, p, a, b, ord2, ord3))
    odd.sort()
    names = {"k": "k_even", "l": "l_even", "s": "sum_even"}
    for index, (h, p, a, b, ord2, ord3) in enumerate(odd):
        key = args.pattern[index % len(args.pattern)]
        parent = names[key]
        aa, bb = induced(parent, h, a, b)
        if math.gcd(aa, bb, h) != 1:
            raise RuntimeError((parent, p, h, aa, bb))
        pools[key]["choices"].append(
            {
                "h": h,
                "p": p,
                "a": aa,
                "b": bb,
                "ord2": ord2,
                "ord3": ord3,
                "c": 0,
                "parent": parent,
                "original_h": h,
                "original_a": a,
                "original_b": b,
                "target_scale": 1,
            }
        )
    for key in sorted(set(args.pattern)):
        payload = pools[key]
        payload["choices"].sort(key=lambda row: (row["h"], row["p"]))
        payload["source"] = str(args.source)
        payload["odd_allocation_pattern"] = args.pattern
        output = args.prefix.with_name(f"{args.prefix.name}_{names[key]}.json")
        output.write_text(json.dumps(payload, indent=2) + "\n")
        print(
            f"{names[key]} output={output} candidates={len(payload['choices'])} "
            f"density={sum(1 / row['h'] for row in payload['choices']):.12f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
