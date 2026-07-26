#!/usr/bin/env python3
"""Certify a period block-star bound with stronger conditional fibre edges."""

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
        "--conditional-certificate",
        action="append",
        type=Path,
        default=[],
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    block = json.loads(args.block_certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
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

    conditionals = {}
    for path in args.conditional_certificate:
        cert = json.loads(path.read_text())
        outside_prime = int(cert["outside_prime"])
        if outside_prime in conditionals:
            raise RuntimeError("repeated conditional outside prime")
        if (
            outside_prime in anchor_primes
            or outside_prime not in selected_primes
        ):
            raise RuntimeError(
                "conditional outside prime is not an outside family row"
            )
        conditional_anchor_primes = {
            int(prime) for prime in cert["anchor_primes"]
        }
        if not conditional_anchor_primes <= anchor_primes:
            raise RuntimeError(
                "conditional certificate uses a non-block anchor"
            )
        for row in cert["anchor_rows"]:
            prime = int(row["p"])
            if prime not in by_prime or any(
                int(by_prime[prime][key]) != int(row[key])
                for key in ("h", "a", "b")
            ):
                raise RuntimeError(
                    "conditional anchor differs from the pool"
                )
        outside_record = cert["outside_row"]
        if any(
            int(by_prime[outside_prime][key])
            != int(outside_record[key])
            for key in ("p", "h", "a", "b")
        ):
            raise RuntimeError(
                "conditional outside row differs from the pool"
            )
        conditionals[outside_prime] = (
            path,
            read_fraction(cert["forced_intersection_density"]),
        )

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
        baseline = best_anchor_lower(row, anchors)
        path = None
        enhanced = baseline
        if prime in conditionals:
            path, conditional_lower = conditionals[prime]
            enhanced = max(baseline, conditional_lower)
        star_loss += enhanced
        row_lowers.append(
            {
                "outside_prime": prime,
                "baseline_intersection_lower_bound":
                fraction_payload(baseline),
                "conditional_certificate": (
                    str(path) if path is not None else None
                ),
                "used_intersection_lower_bound":
                fraction_payload(enhanced),
            }
        )
    block_loss = read_fraction(block["forced_overlap_loss"])
    upper = total_density - block_loss - star_loss
    result = {
        "schema": "ranked_period_conditional_star_v1",
        "pool": str(args.pool),
        "period": args.period,
        "row_count": len(selected),
        "block_certificate": str(args.block_certificate),
        "block_anchor_primes": sorted(anchor_primes),
        "block_overlap_loss": fraction_payload(block_loss),
        "conditional_certificates": [
            str(path) for path in args.conditional_certificate
        ],
        "outside_row_lowers": row_lowers,
        "star_overlap_sum": fraction_payload(star_loss),
        "total_density": fraction_payload(total_density),
        "union_density_upper_bound": fraction_payload(upper),
        "proved_no_cover": upper < 1,
        "argument": (
            "The exact anchor block is one node. Each outside fibre is a "
            "leaf of a star and contributes a certified lower bound for its "
            "intersection with the block union. Conditional certificates may "
            "replace the generic Bonferroni edge by an exact stronger edge. "
            "Subtracting every star edge is pointwise valid."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"period={args.period} rows={len(selected)} anchors={len(anchors)} "
        f"conditional_edges={len(conditionals)} block_loss={block_loss} "
        f"star_loss={star_loss} upper={upper} "
        f"proved_no_cover={upper < 1}",
        flush=True,
    )
    return 0 if upper < 1 else 1


if __name__ == "__main__":
    raise SystemExit(main())
