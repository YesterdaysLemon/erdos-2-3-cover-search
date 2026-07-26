#!/usr/bin/env python3
"""Certify no cover using disjoint unavoidable pair intersections."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from certify_forced_pair_overlap import joint_index


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    if len({int(row["p"]) for row in rows}) != len(rows):
        raise RuntimeError("pool contains a repeated fibre prime")
    if any(int(row.get("target_modulus", 1)) != 1 for row in rows):
        raise RuntimeError("certificate requires unrestricted row targets")
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    edges = sorted(
        (
            (
                Fraction(1, int(left["h"]) * int(right["h"])),
                min(int(left["p"]), int(right["p"])),
                max(int(left["p"]), int(right["p"])),
                left,
                right,
            )
            for index, left in enumerate(rows)
            for right in rows[:index]
            if joint_index(left, right) == 1
        ),
        key=lambda edge: (-edge[0], edge[1], edge[2]),
    )
    used = set()
    selected_pairs = []
    forced_overlap = Fraction(0)
    for weight, _first, _second, left, right in edges:
        left_prime = int(left["p"])
        right_prime = int(right["p"])
        if left_prime in used or right_prime in used:
            continue
        used.update((left_prime, right_prime))
        forced_overlap += weight
        selected_pairs.append(
            {
                "primes": [left_prime, right_prime],
                "rows": [
                    {
                        key: int(row[key])
                        for key in ("p", "h", "a", "b")
                    }
                    for row in (left, right)
                ],
                "forced_intersection_density": fraction_payload(weight),
            }
        )
        if total_density - forced_overlap < 1:
            break

    union_upper_bound = total_density - forced_overlap
    proved = union_upper_bound < 1
    result = {
        "pool": str(args.pool),
        "row_count": len(rows),
        "total_reciprocal_density": fraction_payload(total_density),
        "selected_disjoint_pairs": selected_pairs,
        "forced_overlap_sum": fraction_payload(forced_overlap),
        "union_density_upper_bound": fraction_payload(union_upper_bound),
        "proved_no_cover": proved,
        "argument": (
            "Within each selected pair, the two target maps are jointly "
            "surjective and hence the chosen fibres always intersect with "
            "the stated density. The pairs use disjoint rows, so unioning "
            "each pair first subtracts all listed intersections from the "
            "sum of individual row densities."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"rows={len(rows)} pairs={len(selected_pairs)} "
        f"total={total_density} overlap={forced_overlap} "
        f"upper={union_upper_bound} proved_no_cover={proved}",
        flush=True,
    )
    return 0 if proved else 1


if __name__ == "__main__":
    raise SystemExit(main())
