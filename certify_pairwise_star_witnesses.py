#!/usr/bin/env python3
"""Promote selected pairwise star subsets to an exact witness certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from certify_block_plus_overlap_forest import recorded_anchor_rows
from certify_forced_pair_overlap import joint_index
from search_mitm_pairwise_conditional_bounds import (
    exact_pairwise_subset_lower,
    fraction_payload,
)


def recorded_row(row: dict) -> dict:
    return {key: int(row[key]) for key in ("p", "h", "a", "b")}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("block_certificate", type=Path)
    parser.add_argument("discovery", type=Path)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    block = json.loads(args.block_certificate.read_text())
    discovery = json.loads(args.discovery.read_text())
    if discovery.get("schema") != "mitm_pairwise_conditional_search_v1":
        raise RuntimeError("unexpected discovery schema")
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    anchors = [
        by_prime[int(row["p"])]
        for row in recorded_anchor_rows(block)
    ]
    anchor_primes = {int(row["p"]) for row in anchors}
    family_primes = {
        int(row["p"])
        for row in rows
        if args.period % int(row["h"]) == 0
    }
    outside_primes = family_primes - anchor_primes
    discovery_by_prime = {
        int(record["outside_prime"]): record
        for record in discovery["results"]
    }
    if (
        len(discovery_by_prime) != len(discovery["results"])
        or set(discovery_by_prime) != outside_primes
    ):
        raise RuntimeError("discovery rows do not match the period outside set")

    records = []
    for outside_prime in sorted(outside_primes):
        outside = by_prime[outside_prime]
        selected_primes = tuple(
            int(prime)
            for prime in discovery_by_prime[outside_prime][
                "selected_anchor_primes"
            ]
        )
        if (
            len(set(selected_primes)) != len(selected_primes)
            or not set(selected_primes) <= anchor_primes
        ):
            raise RuntimeError(
                f"invalid selected anchor subset for p={outside_prime}"
            )
        selected = tuple(by_prime[prime] for prime in selected_primes)
        if any(joint_index(outside, anchor) != 1 for anchor in selected):
            raise RuntimeError(
                f"incompatible selected anchor for p={outside_prime}"
            )
        records.append(
            {
                "outside_prime": outside_prime,
                "outside_row": recorded_row(outside),
                "selected_anchor_primes": list(selected_primes),
                "pairwise_intersection_lower_bound": fraction_payload(
                    exact_pairwise_subset_lower(outside, selected)
                ),
            }
        )

    result = {
        "schema": "pairwise_star_witnesses_v1",
        "pool": str(args.pool),
        "block_certificate": str(args.block_certificate),
        "period": args.period,
        "row_count": len(
            [row for row in rows if args.period % int(row["h"]) == 0]
        ),
        "anchor_primes": [int(row["p"]) for row in anchors],
        "outside_count": len(records),
        "records": records,
        "argument": (
            "Every record supplies one explicit compatible anchor subset. "
            "The recorded exact pairwise Bonferroni value is a valid lower "
            "bound for that outside fibre. Discovery optimality is neither "
            "claimed nor used."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"period={args.period} anchors={len(anchors)} "
        f"outside_witnesses={len(records)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
