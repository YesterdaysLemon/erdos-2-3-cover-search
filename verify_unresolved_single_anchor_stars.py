#!/usr/bin/env python3
"""Independently replay single-anchor star certificates."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def pair_index_independent(first: dict, second: dict) -> int:
    h1, a1, b1 = (int(first[key]) for key in ("h", "a", "b"))
    h2, a2, b2 = (int(second[key]) for key in ("h", "a", "b"))
    minors = (
        h1 * h2,
        h1 * a2,
        h1 * b2,
        h2 * a1,
        h2 * b1,
        a1 * b2 - a2 * b1,
    )
    return math.gcd(*(abs(value) for value in minors))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    certificate = json.loads(args.certificate.read_text())
    rows = source["choices"]
    checks = []
    all_valid = True
    proved = 0
    for record in certificate["results"]:
        period = int(record["period"])
        selected = [
            row for row in rows if period % int(row["h"]) == 0
        ]
        by_prime = {int(row["p"]): row for row in selected}
        anchor_prime = record["anchor_prime"]
        anchor = (
            by_prime.get(int(anchor_prime))
            if anchor_prime is not None
            else None
        )
        total_density = sum(
            (Fraction(1, int(row["h"])) for row in selected),
            Fraction(0),
        )
        if anchor is None:
            overlap = Fraction(0)
            anchor_valid = not selected
            edge_count = 0
        else:
            compatible = [
                outside
                for outside in selected
                if int(outside["p"]) != int(anchor["p"])
                and pair_index_independent(anchor, outside) == 1
            ]
            overlap = sum(
                (
                    Fraction(
                        1,
                        int(anchor["h"]) * int(outside["h"]),
                    )
                    for outside in compatible
                ),
                Fraction(0),
            )
            anchor_valid = True
            edge_count = len(compatible)
        upper = total_density - overlap
        no_cover = upper < 1
        valid = (
            int(record["rows"]) == len(selected)
            and read_fraction(record["total_density"]) == total_density
            and anchor_valid
            and read_fraction(record["star_overlap_sum"]) == overlap
            and read_fraction(record["union_upper_bound"]) == upper
            and bool(record["proved_no_cover"]) == no_cover
        )
        all_valid &= valid
        proved += int(no_cover)
        checks.append(
            {
                "period": period,
                "anchor_prime": anchor_prime,
                "edge_count": edge_count,
                "overlap_sum": {
                    "numerator": overlap.numerator,
                    "denominator": overlap.denominator,
                },
                "union_upper_bound": {
                    "numerator": upper.numerator,
                    "denominator": upper.denominator,
                },
                "proved_no_cover": no_cover,
                "valid": valid,
            }
        )

    verified = (
        str(args.pool) == certificate["pool"]
        and int(certificate["families_checked"]) == len(checks)
        and int(certificate["proved_no_cover"]) == proved
        and int(certificate["unresolved_count"]) == len(checks) - proved
        and all_valid
    )
    output = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "families_checked": len(checks),
        "proved_no_cover": proved,
        "unresolved_count": len(checks) - proved,
        "checks": checks,
        "verified": verified,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"families={len(checks)} proved={proved} "
        f"unresolved={len(checks)-proved} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
