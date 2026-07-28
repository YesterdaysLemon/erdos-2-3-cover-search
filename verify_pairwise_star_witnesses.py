#!/usr/bin/env python3
"""Independent replay of explicit pairwise star-subset witnesses."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from certify_block_plus_overlap_forest import recorded_anchor_rows
from verify_all_ranked_pairanchor_star import pair_index, read_fraction
from verify_ranked_period_conditional_star import witnessed_pairwise_lower


def recorded_row(row: dict) -> dict:
    return {key: int(row[key]) for key in ("p", "h", "a", "b")}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    cert = json.loads(args.certificate.read_text())
    block = json.loads(Path(cert["block_certificate"]).read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    anchors = [
        by_prime[int(row["p"])]
        for row in recorded_anchor_rows(block)
    ]
    anchor_primes = {int(row["p"]) for row in anchors}
    period = int(cert["period"])
    family = [
        row for row in rows if period % int(row["h"]) == 0
    ]
    outside_primes = {
        int(row["p"]) for row in family
    } - anchor_primes
    record_by_prime = {
        int(record["outside_prime"]): record
        for record in cert["records"]
    }

    records_valid = (
        len(record_by_prime) == len(cert["records"])
        and set(record_by_prime) == outside_primes
    )
    checked = 0
    total_lower = Fraction(0)
    if records_valid:
        for outside_prime in sorted(outside_primes):
            record = record_by_prime[outside_prime]
            outside = by_prime[outside_prime]
            selected_primes = tuple(
                int(prime)
                for prime in record["selected_anchor_primes"]
            )
            valid = (
                record["outside_row"] == recorded_row(outside)
                and len(set(selected_primes)) == len(selected_primes)
                and set(selected_primes) <= anchor_primes
                and all(
                    pair_index(outside, by_prime[prime]) == 1
                    for prime in selected_primes
                )
            )
            if not valid:
                records_valid = False
                break
            value = witnessed_pairwise_lower(
                outside,
                tuple(by_prime[prime] for prime in selected_primes),
            )
            if read_fraction(
                record["pairwise_intersection_lower_bound"]
            ) != value:
                records_valid = False
                break
            total_lower += value
            checked += 1

    verified = (
        cert.get("schema") == "pairwise_star_witnesses_v1"
        and cert.get("pool") == str(args.pool)
        and len(by_prime) == len(rows)
        and [int(row["p"]) for row in anchors] == cert["anchor_primes"]
        and int(cert["row_count"]) == len(family)
        and int(cert["outside_count"]) == len(outside_primes)
        and records_valid
    )
    result = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "period": period,
        "row_count": len(family),
        "anchor_count": len(anchors),
        "outside_witnesses_checked": checked,
        "pairwise_lower_sum": {
            "numerator": total_lower.numerator,
            "denominator": total_lower.denominator,
        },
        "verified": verified,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"period={period} outside_witnesses={checked} "
        f"verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
