#!/usr/bin/env python3
"""Unverified preview of the four-way correction on current survivors.

This is a discovery scan only.  Any apparent elimination requires a separate
replay before it is promoted.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from certify_block_plus_overlap_forest import (
    read_fraction,
    recorded_anchor_rows,
)
from scan_all_ranked_pairanchor_star import (
    best_anchor_lower,
    fraction_payload,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("aggregate_scan", type=Path)
    parser.add_argument("survivor_scan", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    aggregate = json.loads(args.aggregate_scan.read_text())
    survivors = json.loads(args.survivor_scan.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    aggregate_by_period = {
        int(record["period"]): record for record in aggregate["results"]
    }
    survivor_periods = {
        int(record["period"]) for record in survivors["unresolved"]
    }
    block_cache = {}
    row_cache = {}
    results = []
    improved = 0
    newly_proved = 0
    for period in sorted(survivor_periods):
        prior = aggregate_by_period[period]
        block_path = prior["block_certificate"]
        if block_path is None:
            results.append(
                {
                    "period": period,
                    "block_certificate": None,
                    "prior_upper": prior["union_upper_bound"],
                    "strengthened_upper": prior["union_upper_bound"],
                    "star_improvement": fraction_payload(Fraction(0)),
                    "newly_proved_no_cover": False,
                }
            )
            continue
        if block_path not in block_cache:
            block = json.loads(Path(block_path).read_text())
            recorded = recorded_anchor_rows(block)
            anchors = [by_prime[int(row["p"])] for row in recorded]
            block_cache[block_path] = (
                anchors,
                {int(row["p"]) for row in anchors},
                read_fraction(block["forced_overlap_loss"]),
            )
        anchors, anchor_primes, block_loss = block_cache[block_path]
        selected = [
            row for row in rows if period % int(row["h"]) == 0
        ]
        total_density = sum(
            (Fraction(1, int(row["h"])) for row in selected),
            Fraction(0),
        )
        strengthened_star = Fraction(0)
        for row in selected:
            prime = int(row["p"])
            if prime in anchor_primes:
                continue
            key = (block_path, prime)
            if key not in row_cache:
                row_cache[key] = best_anchor_lower(
                    row,
                    anchors,
                    fourway_triples=True,
                )
            strengthened_star += row_cache[key]
        prior_star = read_fraction(prior["pair_anchor_star_loss"])
        improvement = strengthened_star - prior_star
        upper = total_density - block_loss - strengthened_star
        improved += int(improvement > 0)
        newly_proved += int(upper < 1)
        results.append(
            {
                "period": period,
                "block_certificate": block_path,
                "prior_upper": prior["union_upper_bound"],
                "strengthened_upper": fraction_payload(upper),
                "star_improvement": fraction_payload(improvement),
                "newly_proved_no_cover": upper < 1,
            }
        )
        print(
            f"period={period} improvement={float(improvement):.12g} "
            f"upper={float(upper):.12f} proved={upper < 1}",
            flush=True,
        )

    output = {
        "status": "UNVERIFIED_PREVIEW",
        "pool": str(args.pool),
        "aggregate_scan": str(args.aggregate_scan),
        "survivor_scan": str(args.survivor_scan),
        "families_checked": len(results),
        "families_with_strict_improvement": improved,
        "newly_proved_no_cover_unverified": newly_proved,
        "results": results,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"families={len(results)} improved={improved} "
        f"newly_proved_unverified={newly_proved}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
