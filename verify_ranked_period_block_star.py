#!/usr/bin/env python3
"""Independent replay of a one-period block-star certificate."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from verify_all_ranked_pairanchor_star import (
    best_lower,
    read_fraction,
    recorded_anchor_rows,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("block_verification", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    cert = json.loads(args.certificate.read_text())
    block_path = Path(cert["block_certificate"])
    block = json.loads(block_path.read_text())
    block_verification = json.loads(args.block_verification.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    recorded = recorded_anchor_rows(block)
    anchors_match = all(
        int(row["p"]) in by_prime
        and all(
            int(by_prime[int(row["p"])][key]) == int(row[key])
            for key in ("h", "a", "b")
        )
        for row in recorded
    )
    anchors = [by_prime[int(row["p"])] for row in recorded]
    anchor_primes = {int(row["p"]) for row in anchors}
    period = int(cert["period"])
    fourway_triples = bool(cert.get("fourway_triples", False))
    selected = [row for row in rows if period % int(row["h"]) == 0]
    selected_primes = {int(row["p"]) for row in selected}
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in selected),
        Fraction(0),
    )
    block_loss = read_fraction(block["forced_overlap_loss"])
    expected_lowers = {}
    star_loss = Fraction(0)
    for row in selected:
        prime = int(row["p"])
        if prime in anchor_primes:
            continue
        lower = best_lower(
            row,
            anchors,
            fourway_triples=fourway_triples,
        )
        expected_lowers[prime] = lower
        star_loss += lower
    recorded_lowers = {
        int(record["outside_prime"]): read_fraction(
            record["block_intersection_lower_bound"]
        )
        for record in cert["outside_row_lowers"]
    }
    upper = total_density - block_loss - star_loss
    block_verification_valid = (
        bool(block_verification.get("verified"))
        and block_verification.get("certificate") == str(block_path)
        and read_fraction(block_verification["forced_overlap_loss"])
        == block_loss
    )
    verified = (
        cert.get("schema")
        == (
            "ranked_period_block_star_v2"
            if fourway_triples
            else "ranked_period_block_star_v1"
        )
        and str(args.pool) == cert["pool"]
        and bool(cert.get("fourway_triples", False)) == fourway_triples
        and len(by_prime) == len(rows)
        and anchors_match
        and anchor_primes <= selected_primes
        and sorted(anchor_primes) == cert["block_anchor_primes"]
        and block_verification_valid
        and int(cert["row_count"]) == len(selected)
        and recorded_lowers == expected_lowers
        and read_fraction(cert["block_overlap_loss"]) == block_loss
        and read_fraction(cert["star_overlap_sum"]) == star_loss
        and read_fraction(cert["total_density"]) == total_density
        and read_fraction(cert["union_density_upper_bound"]) == upper
        and bool(cert["proved_no_cover"]) == (upper < 1)
    )
    result = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "block_verification": str(args.block_verification),
        "period": period,
        "fourway_triples": fourway_triples,
        "row_count": len(selected),
        "anchor_count": len(anchors),
        "outside_rows_checked": len(expected_lowers),
        "block_overlap_loss": {
            "numerator": block_loss.numerator,
            "denominator": block_loss.denominator,
        },
        "star_overlap_sum": {
            "numerator": star_loss.numerator,
            "denominator": star_loss.denominator,
        },
        "union_density_upper_bound": {
            "numerator": upper.numerator,
            "denominator": upper.denominator,
        },
        "verified": verified,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"period={period} rows={len(selected)} anchors={len(anchors)} "
        f"outside={len(expected_lowers)} upper={upper} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
