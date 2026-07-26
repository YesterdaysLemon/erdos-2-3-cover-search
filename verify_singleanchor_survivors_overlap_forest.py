#!/usr/bin/env python3
"""Independently verify every claimed survivor-overlap forest proof."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

from verify_block_plus_overlap_forest import CycleChecker


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def pair_index(first: dict, second: dict) -> int:
    h1, a1, b1 = (int(first[key]) for key in ("h", "a", "b"))
    h2, a2, b2 = (int(second[key]) for key in ("h", "a", "b"))
    return math.gcd(
        abs(h1 * h2),
        abs(h1 * a2),
        abs(h1 * b2),
        abs(h2 * a1),
        abs(h2 * b1),
        abs(a1 * b2 - a2 * b1),
    )


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
    invalid_periods = []
    proved = 0
    for record in certificate["results"]:
        period = int(record["period"])
        selected = [
            row for row in rows if period % int(row["h"]) == 0
        ]
        by_prime = {int(row["p"]): row for row in selected}
        total_density = sum(
            (Fraction(1, int(row["h"])) for row in selected),
            Fraction(0),
        )
        if record["proved_no_cover"]:
            cycle_checker = CycleChecker(set(by_prime))
            overlap = Fraction(0)
            edges_valid = True
            seen = set()
            for first, second in record["selected_forest_edges"]:
                first = int(first)
                second = int(second)
                canonical = tuple(sorted((first, second)))
                endpoints_valid = (
                    first != second
                    and first in by_prime
                    and second in by_prime
                    and canonical not in seen
                )
                seen.add(canonical)
                acyclic = (
                    cycle_checker.add_edge(first, second)
                    if endpoints_valid
                    else False
                )
                surjective = (
                    pair_index(by_prime[first], by_prime[second]) == 1
                    if endpoints_valid
                    else False
                )
                edges_valid &= endpoints_valid and acyclic and surjective
                if endpoints_valid and surjective:
                    overlap += Fraction(
                        1,
                        int(by_prime[first]["h"])
                        * int(by_prime[second]["h"]),
                    )
            upper = total_density - overlap
            valid = (
                int(record["rows"]) == len(selected)
                and read_fraction(record["total_density"]) == total_density
                and edges_valid
                and read_fraction(record["forest_overlap_sum"]) == overlap
                and read_fraction(record["union_upper_bound"]) == upper
                and upper < 1
            )
            proved += int(valid)
        else:
            overlap = Fraction(0)
            upper = total_density
            valid = (
                int(record["rows"]) == len(selected)
                and read_fraction(record["total_density"]) == total_density
                and not record["selected_forest_edges"]
            )
        if not valid:
            invalid_periods.append(period)
        checks.append(
            {
                "period": period,
                "claimed_no_cover": bool(record["proved_no_cover"]),
                "verified_no_cover": bool(
                    record["proved_no_cover"] and valid
                ),
                "valid": valid,
            }
        )

    claimed_proved = sum(
        bool(record["proved_no_cover"])
        for record in certificate["results"]
    )
    verified = (
        str(args.pool) == certificate["pool"]
        and not invalid_periods
        and int(certificate["families_checked"]) == len(checks)
        and int(certificate["proved_no_cover"]) == claimed_proved
        and proved == claimed_proved
        and int(certificate["unresolved_count"])
        == len(checks) - claimed_proved
    )
    output = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "families_checked": len(checks),
        "claimed_no_cover": claimed_proved,
        "verified_no_cover": proved,
        "unresolved_count": len(checks) - claimed_proved,
        "invalid_periods": invalid_periods,
        "checks": checks,
        "verified": verified,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"families={len(checks)} proved={proved} "
        f"unresolved={len(checks)-claimed_proved} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
