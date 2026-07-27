#!/usr/bin/env python3
"""Independent replay of a conditional period block-star certificate."""

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
    parser.add_argument(
        "--conditional-verification",
        action="append",
        type=Path,
        default=[],
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    cert = json.loads(args.certificate.read_text())
    fourway_triples = bool(cert.get("fourway_triples", False))
    fourterm_quads = bool(cert.get("fourterm_quads", False))
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
    selected = [row for row in rows if period % int(row["h"]) == 0]
    selected_primes = {int(row["p"]) for row in selected}
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in selected),
        Fraction(0),
    )
    block_loss = read_fraction(block["forced_overlap_loss"])

    verification_by_certificate = {}
    for path in args.conditional_verification:
        verification = json.loads(path.read_text())
        certificate_path = str(verification["certificate"])
        if certificate_path in verification_by_certificate:
            raise RuntimeError("repeated conditional verification")
        verification_by_certificate[certificate_path] = verification
    conditional_values = {}
    conditional_valid = True
    for certificate_path in cert["conditional_certificates"]:
        path = Path(certificate_path)
        conditional = json.loads(path.read_text())
        verification = verification_by_certificate.get(certificate_path)
        outside_prime = int(conditional["outside_prime"])
        anchor_subset = {
            int(prime) for prime in conditional["anchor_primes"]
        }
        rows_match = (
            outside_prime in by_prime
            and all(
                int(by_prime[outside_prime][key])
                == int(conditional["outside_row"][key])
                for key in ("p", "h", "a", "b")
            )
            and all(
                int(row["p"]) in by_prime
                and all(
                    int(by_prime[int(row["p"])][key]) == int(row[key])
                    for key in ("h", "a", "b")
                )
                for row in conditional["anchor_rows"]
            )
        )
        value = read_fraction(
            conditional["forced_intersection_density"]
        )
        valid = (
            verification is not None
            and bool(verification.get("verified"))
            and verification.get("certificate") == certificate_path
            and int(verification["outside_prime"]) == outside_prime
            and read_fraction(
                verification["forced_intersection_density"]
            ) == value
            and outside_prime in selected_primes - anchor_primes
            and anchor_subset <= anchor_primes
            and rows_match
        )
        conditional_valid = conditional_valid and valid
        conditional_values[outside_prime] = (
            certificate_path,
            value,
        )

    expected_lowers = {}
    star_loss = Fraction(0)
    for row in selected:
        prime = int(row["p"])
        if prime in anchor_primes:
            continue
        baseline = best_lower(
            row,
            anchors,
            fourway_triples=fourway_triples,
            fourterm_quads=fourterm_quads,
        )
        path = None
        used = baseline
        if prime in conditional_values:
            path, value = conditional_values[prime]
            used = max(baseline, value)
        expected_lowers[prime] = (baseline, path, used)
        star_loss += used
    recorded_lowers = {
        int(record["outside_prime"]): (
            read_fraction(
                record["baseline_intersection_lower_bound"]
            ),
            record["conditional_certificate"],
            read_fraction(record["used_intersection_lower_bound"]),
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
        cert.get("schema") == "ranked_period_conditional_star_v1"
        and str(args.pool) == cert["pool"]
        and len(by_prime) == len(rows)
        and anchors_match
        and anchor_primes <= selected_primes
        and sorted(anchor_primes) == cert["block_anchor_primes"]
        and block_verification_valid
        and conditional_valid
        and set(verification_by_certificate)
        == set(cert["conditional_certificates"])
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
        "conditional_verifications": [
            str(path) for path in args.conditional_verification
        ],
        "period": period,
        "fourway_triples": fourway_triples,
        "fourterm_quads": fourterm_quads,
        "row_count": len(selected),
        "anchor_count": len(anchors),
        "outside_rows_checked": len(expected_lowers),
        "conditional_edges_checked": len(conditional_values),
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
        f"conditional={len(conditional_values)} upper={upper} "
        f"verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
