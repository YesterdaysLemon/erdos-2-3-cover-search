#!/usr/bin/env python3
"""Certify a no-cover result from one unavoidable pair intersection."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def joint_index(left: dict, right: dict) -> int:
    h1, a1, b1 = (int(left[key]) for key in ("h", "a", "b"))
    h2, a2, b2 = (int(right[key]) for key in ("h", "a", "b"))
    return math.gcd(
        a1 * b2 - a2 * b1,
        a2 * h1,
        b2 * h1,
        a1 * h2,
        b1 * h2,
        h1 * h2,
    )


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--first-prime", type=int, required=True)
    parser.add_argument("--second-prime", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.first_prime == args.second_prime:
        raise SystemExit("the two fibre primes must be distinct")

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    if len(by_prime) != len(rows):
        raise RuntimeError("pool contains a repeated fibre prime")
    missing = {
        args.first_prime,
        args.second_prime,
    } - by_prime.keys()
    if missing:
        raise RuntimeError(f"anchor primes are absent: {sorted(missing)}")
    for row in rows:
        h = int(row["h"])
        if math.gcd(int(row["a"]), int(row["b"]), h) != 1:
            raise RuntimeError(f"row p={row['p']} is not surjective")
        if int(row.get("target_modulus", 1)) != 1:
            raise RuntimeError(
                "certificate currently requires unrestricted row targets"
            )

    left = by_prime[args.first_prime]
    right = by_prime[args.second_prime]
    index = joint_index(left, right)
    joint_surjective = index == 1
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    forced_overlap = (
        Fraction(1, int(left["h"]) * int(right["h"]))
        if joint_surjective
        else Fraction(0)
    )
    union_upper_bound = total_density - forced_overlap
    proved = joint_surjective and union_upper_bound < 1

    result = {
        "pool": str(args.pool),
        "row_count": len(rows),
        "anchor_primes": [args.first_prime, args.second_prime],
        "anchor_rows": [
            {
                key: int(row[key])
                for key in ("p", "h", "a", "b")
            }
            for row in (left, right)
        ],
        "joint_image_index": index,
        "joint_target_map_surjective": joint_surjective,
        "total_reciprocal_density": fraction_payload(total_density),
        "forced_pair_overlap_density": fraction_payload(forced_overlap),
        "union_density_upper_bound": fraction_payload(union_upper_bound),
        "proved_no_cover": proved,
        "argument": (
            "The two anchor target maps are jointly surjective, so every "
            "pair of anchor fibres intersects with the stated density. "
            "The union of all rows is at most the sum of their individual "
            "densities minus this unavoidable intersection."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"rows={len(rows)} total={total_density} "
        f"forced_overlap={forced_overlap} upper={union_upper_bound} "
        f"proved_no_cover={proved}",
        flush=True,
    )
    return 0 if proved else 1


if __name__ == "__main__":
    raise SystemExit(main())
