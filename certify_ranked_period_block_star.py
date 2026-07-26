#!/usr/bin/env python3
"""Certify a block-star upper bound for one divisor-period family."""

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
    parser.add_argument("block_certificate", type=Path)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument(
        "--fourway-triples",
        action="store_true",
        help=(
            "restore the positive four-fold term for jointly-surjective "
            "outside-plus-three-anchor target maps"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    block = json.loads(args.block_certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    if len(by_prime) != len(rows):
        raise RuntimeError("pool has repeated fibre primes")
    recorded = recorded_anchor_rows(block)
    anchors = []
    for row in recorded:
        prime = int(row["p"])
        if prime not in by_prime or any(
            int(by_prime[prime][key]) != int(row[key])
            for key in ("h", "a", "b")
        ):
            raise RuntimeError("block anchor differs from the pool")
        anchors.append(by_prime[prime])
    anchor_primes = {int(row["p"]) for row in anchors}
    selected = [
        row for row in rows if args.period % int(row["h"]) == 0
    ]
    selected_primes = {int(row["p"]) for row in selected}
    if not anchor_primes <= selected_primes:
        raise RuntimeError("block does not lie in the requested family")

    total_density = sum(
        (Fraction(1, int(row["h"])) for row in selected),
        Fraction(0),
    )
    row_lowers = []
    star_loss = Fraction(0)
    for row in selected:
        prime = int(row["p"])
        if prime in anchor_primes:
            continue
        lower = best_anchor_lower(
            row,
            anchors,
            fourway_triples=args.fourway_triples,
        )
        star_loss += lower
        row_lowers.append(
            {
                "outside_prime": prime,
                "block_intersection_lower_bound": fraction_payload(lower),
            }
        )
    block_loss = read_fraction(block["forced_overlap_loss"])
    upper = total_density - block_loss - star_loss
    result = {
        "schema": (
            "ranked_period_block_star_v2"
            if args.fourway_triples
            else "ranked_period_block_star_v1"
        ),
        "pool": str(args.pool),
        "period": args.period,
        "fourway_triples": args.fourway_triples,
        "row_count": len(selected),
        "block_certificate": str(args.block_certificate),
        "block_anchor_primes": sorted(anchor_primes),
        "block_overlap_loss": fraction_payload(block_loss),
        "outside_row_lowers": row_lowers,
        "star_overlap_sum": fraction_payload(star_loss),
        "total_density": fraction_payload(total_density),
        "union_density_upper_bound": fraction_payload(upper),
        "proved_no_cover": upper < 1,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"period={args.period} rows={len(selected)} "
        f"anchors={len(anchors)} block_loss={block_loss} "
        f"star_loss={star_loss} upper={upper} "
        f"proved_no_cover={upper < 1}",
        flush=True,
    )
    return 0 if upper < 1 else 1


if __name__ == "__main__":
    raise SystemExit(main())
