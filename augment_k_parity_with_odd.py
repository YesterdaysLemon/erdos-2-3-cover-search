#!/usr/bin/env python3
"""Allocate every unused odd-order prime fibre to the k-even subproblem."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import cegis_cover
import search_cover


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("base", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    payload = json.loads(args.base.read_text())
    choices = list(payload["choices"])
    used = {int(row["p"]) for row in choices}
    for old in json.loads(args.source.read_text())["choices"]:
        p = int(old["p"])
        ord2 = search_cover.multiplicative_order(2, p)
        ord3 = search_cover.multiplicative_order(3, p)
        h, a, b, ord2, ord3 = cegis_cover.signature(p, ord2, ord3)
        if h % 2 == 0 or p in used:
            continue
        # On k=2x, l=y, the original line becomes
        # (2a)x + by = c (mod h).  Since h is odd, 2 is a unit and the
        # induced map remains surjective with the same index h.
        aa = 2 * a % h
        bb = b % h
        if math.gcd(aa, bb, h) != 1:
            raise RuntimeError((p, h, aa, bb))
        choices.append(
            {
                "h": h,
                "p": p,
                "a": aa,
                "b": bb,
                "ord2": ord2,
                "ord3": ord3,
                "c": 0,
                "parent": "k_even",
                "original_h": h,
                "original_a": a,
                "original_b": b,
                "target_scale": 1,
            }
        )
    choices.sort(key=lambda row: (row["h"], row["p"]))
    payload.update(
        {
            "source": str(args.source),
            "parent": "k_even_plus_odd",
            "choices": choices,
        }
    )
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"output={args.output} candidates={len(choices)} "
        f"density={sum(1 / row['h'] for row in choices):.12f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
