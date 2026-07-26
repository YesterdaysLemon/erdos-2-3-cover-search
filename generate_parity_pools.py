#!/usr/bin/env python3
"""Split corrected prime fibres into three induced parity-line problems."""

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
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument(
        "--share-odd",
        action="store_true",
        help=(
            "include every odd-order fibre in all three parent pools; its "
            "target remains shared in any eventual combined construction"
        ),
    )
    args = parser.parse_args()

    raw = json.loads(args.source.read_text())
    groups: dict[str, list[dict]] = {"k_even": [], "l_even": [], "sum_even": []}
    for old in raw["choices"]:
        p = int(old["p"])
        ord2 = search_cover.multiplicative_order(2, p)
        ord3 = search_cover.multiplicative_order(3, p)
        h, a, b, ord2, ord3 = cegis_cover.signature(p, ord2, ord3)
        if h % 2:
            if args.share_odd:
                transforms = {
                    "k_even": (2 * a, b),
                    "l_even": (a, 2 * b),
                    "sum_even": (a + b, 2 * b),
                }
                for name, (aa, bb) in transforms.items():
                    aa %= h
                    bb %= h
                    if math.gcd(aa, bb, h) != 1:
                        raise RuntimeError((p, h, a, b, name, aa, bb))
                    groups[name].append(
                        {
                            "h": h,
                            "p": p,
                            "a": aa,
                            "b": bb,
                            "ord2": ord2,
                            "ord3": ord3,
                            "c": 0,
                            "parent": name,
                            "original_h": h,
                            "original_a": a,
                            "original_b": b,
                            "target_scale": 1,
                            "odd_shared": True,
                        }
                    )
            continue
        if (a % 2, b % 2) == (1, 0):
            name, aa, bb = "k_even", a, b // 2
        elif (a % 2, b % 2) == (0, 1):
            name, aa, bb = "l_even", a // 2, b
        elif (a % 2, b % 2) == (1, 1):
            name, aa, bb = "sum_even", (a + b) // 2, b
        else:
            continue
        modulus = h // 2
        aa %= modulus
        bb %= modulus
        if math.gcd(aa, bb, modulus) != 1:
            raise RuntimeError((p, h, a, b, modulus, aa, bb))
        groups[name].append(
            {
                "h": modulus,
                "p": p,
                "a": aa,
                "b": bb,
                "ord2": ord2,
                "ord3": ord3,
                "c": 0,
                "parent": name,
                "original_h": h,
                "original_a": a,
                "original_b": b,
                "target_scale": 2,
            }
        )

    for name, choices in groups.items():
        choices.sort(key=lambda row: (row["h"], row["p"]))
        output = args.prefix.with_name(f"{args.prefix.name}_{name}.json")
        output.write_text(
            json.dumps(
                {
                    "source": str(args.source),
                    "parent": name,
                    "choices": choices,
                },
                indent=2,
            )
            + "\n"
        )
        print(
            f"{name}: output={output} candidates={len(choices)} "
            f"density={sum(1 / row['h'] for row in choices):.12f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
