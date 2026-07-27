#!/usr/bin/env python3
"""Independently replay a conditional block-extension certificate."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from certify_block_plus_overlap_forest import (
    read_fraction,
    recorded_anchor_rows,
)


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
    }


def recorded_row(row: dict) -> dict:
    return {key: int(row[key]) for key in ("p", "h", "a", "b")}


def same_path(left: str | Path, right: str | Path) -> bool:
    return Path(left).resolve() == Path(right).resolve()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("base_verification", type=Path)
    parser.add_argument("conditional_verification", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    certificate = json.loads(args.certificate.read_text())
    base = json.loads(Path(certificate["base_certificate"]).read_text())
    conditional = json.loads(
        Path(certificate["conditional_certificate"]).read_text()
    )
    base_verification = json.loads(args.base_verification.read_text())
    conditional_verification = json.loads(
        args.conditional_verification.read_text()
    )
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}

    base_anchor_rows = recorded_anchor_rows(base)
    base_primes = tuple(int(row["p"]) for row in base_anchor_rows)
    extra_prime = int(conditional["outside_prime"])
    extra = by_prime.get(extra_prime)
    base_loss = read_fraction(base["forced_overlap_loss"])
    base_sum = sum(
        (Fraction(1, int(row["h"])) for row in base_anchor_rows),
        Fraction(0),
    )
    conditional_lower = read_fraction(
        conditional["forced_intersection_density"]
    )
    extended_sum = (
        base_sum + Fraction(1, int(extra["h"]))
        if extra is not None
        else Fraction(0)
    )
    extended_loss = base_loss + conditional_lower
    extended_union_upper = extended_sum - extended_loss
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    pool_upper = total_density - extended_loss

    base_rows_match = (
        base_primes == tuple(int(value) for value in base["anchor_primes"])
        and all(
            int(row["p"]) in by_prime
            and recorded_row(by_prime[int(row["p"])]) == recorded_row(row)
            for row in base_anchor_rows
        )
    )
    conditional_matches = (
        same_path(
            conditional["block_certificate"],
            certificate["base_certificate"],
        )
        and tuple(
            int(value) for value in conditional["anchor_primes"]
        )
        == base_primes
        and [
            recorded_row(row) for row in conditional["anchor_rows"]
        ]
        == [recorded_row(row) for row in base_anchor_rows]
        and extra is not None
        and extra_prime not in base_primes
        and recorded_row(extra) == recorded_row(conditional["outside_row"])
    )
    base_verified = (
        bool(base_verification.get("verified"))
        and same_path(
            base_verification["certificate"],
            certificate["base_certificate"],
        )
        and read_fraction(base_verification["forced_overlap_loss"])
        == base_loss
    )
    conditional_verified = (
        bool(conditional_verification.get("verified"))
        and same_path(
            conditional_verification["certificate"],
            certificate["conditional_certificate"],
        )
        and int(conditional_verification["outside_prime"]) == extra_prime
        and read_fraction(
            conditional_verification["forced_intersection_density"]
        )
        == conditional_lower
    )
    expected_anchor_rows = [
        *[recorded_row(row) for row in base_anchor_rows],
        *((
            recorded_row(extra),
        ) if extra is not None else ()),
    ]
    verified = (
        certificate.get("schema") == "conditional_block_extension_v1"
        and same_path(certificate["pool"], args.pool)
        and len(by_prime) == len(rows)
        and int(certificate["row_count"]) == len(rows)
        and base_rows_match
        and conditional_matches
        and base_verified
        and conditional_verified
        and tuple(int(value) for value in certificate["base_anchor_primes"])
        == base_primes
        and [
            recorded_row(row) for row in certificate["base_anchor_rows"]
        ]
        == [recorded_row(row) for row in base_anchor_rows]
        and recorded_row(certificate["extra_row"]) == recorded_row(extra)
        and tuple(int(value) for value in certificate["anchor_primes"])
        == (*base_primes, extra_prime)
        and [
            recorded_row(row) for row in certificate["anchor_rows"]
        ]
        == expected_anchor_rows
        and read_fraction(certificate["base_individual_density_sum"])
        == base_sum
        and read_fraction(certificate["base_forced_overlap_loss"])
        == base_loss
        and read_fraction(
            certificate["conditional_intersection_lower_bound"]
        )
        == conditional_lower
        and read_fraction(certificate["block_individual_density_sum"])
        == extended_sum
        and read_fraction(certificate["forced_overlap_loss"])
        == extended_loss
        and read_fraction(
            certificate["block_union_density_upper_bound"]
        )
        == extended_union_upper
        and read_fraction(certificate["total_pool_density"])
        == total_density
        and read_fraction(certificate["pool_union_density_upper_bound"])
        == pool_upper
        and bool(certificate["proved_no_cover"]) == (pool_upper < 1)
    )
    result = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "base_verification": str(args.base_verification),
        "conditional_verification": str(args.conditional_verification),
        "base_verified": base_verified,
        "conditional_verified": conditional_verified,
        "base_rows_match": base_rows_match,
        "conditional_matches": conditional_matches,
        "anchor_primes": [*base_primes, extra_prime],
        "block_individual_density_sum": fraction_payload(extended_sum),
        "forced_overlap_loss": fraction_payload(extended_loss),
        "block_union_density_upper_bound": fraction_payload(
            extended_union_upper
        ),
        "pool_union_density_upper_bound": fraction_payload(pool_upper),
        "verified": verified,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base_verified={base_verified} "
        f"conditional_verified={conditional_verified} "
        f"anchors={len(base_primes) + 1} "
        f"loss={extended_loss} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
