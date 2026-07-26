#!/usr/bin/env python3
"""Resolve prior scan exceptions with exact single-anchor star forests."""

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
    parser.add_argument("prior_scan", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    prior = json.loads(args.prior_scan.read_text())
    rows = source["choices"]
    moduli = [int(row["h"]) for row in rows]
    compatibility = [bytearray(len(rows)) for _ in rows]
    for first in range(len(rows)):
        for second in range(first):
            if joint_index(rows[first], rows[second]) == 1:
                compatibility[first][second] = 1
                compatibility[second][first] = 1
    results = []
    proved = 0
    for family in prior["unresolved"]:
        period = int(family["period"])
        selected_indices = [
            index
            for index, modulus in enumerate(moduli)
            if period % modulus == 0
        ]
        selected = [rows[index] for index in selected_indices]
        total_density = Fraction(
            sum(period // moduli[index] for index in selected_indices),
            period,
        )
        choices = []
        for anchor_index in selected_indices:
            overlap = Fraction(
                sum(
                    period // moduli[outside_index]
                    for outside_index in selected_indices
                    if outside_index != anchor_index
                    and compatibility[anchor_index][outside_index]
                ),
                period * moduli[anchor_index],
            )
            choices.append(
                (
                    total_density - overlap,
                    int(rows[anchor_index]["p"]),
                    overlap,
                )
            )
        if choices:
            upper, anchor_prime, overlap = min(choices)
        else:
            upper = total_density
            anchor_prime = None
            overlap = Fraction(0)
        no_cover = upper < 1
        proved += int(no_cover)
        results.append(
            {
                "period": period,
                "components": family["components"],
                "rows": len(selected),
                "total_density": fraction_payload(total_density),
                "anchor_prime": anchor_prime,
                "star_overlap_sum": fraction_payload(overlap),
                "union_upper_bound": fraction_payload(upper),
                "proved_no_cover": no_cover,
            }
        )

    unresolved = [
        result for result in results if not result["proved_no_cover"]
    ]
    unresolved.sort(
        key=lambda item: (
            item["union_upper_bound"]["decimal"],
            item["period"],
        )
    )
    output = {
        "pool": str(args.pool),
        "prior_scan": str(args.prior_scan),
        "families_checked": len(results),
        "proved_no_cover": proved,
        "unresolved_count": len(unresolved),
        "unresolved": unresolved,
        "results": results,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"families={len(results)} proved={proved} "
        f"unresolved={len(unresolved)}",
        flush=True,
    )
    for result in unresolved[:20]:
        print(
            f"period={result['period']} rows={result['rows']} "
            f"upper={result['union_upper_bound']['decimal']:.12f}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
