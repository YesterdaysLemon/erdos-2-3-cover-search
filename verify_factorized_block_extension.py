#!/usr/bin/env python3
"""Independently verify a coprime factorized block-extension certificate."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("base_verification", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    certificate = json.loads(args.certificate.read_text())
    base_verification = json.loads(args.base_verification.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    base_primes = tuple(
        int(p) for p in certificate["base_anchor_primes"]
    )
    extra_record = certificate["extra_row"]
    extra = by_prime[int(extra_record["p"])]
    base_period = int(certificate["base_period"])
    extra_h = int(extra["h"])
    base_union = read_fraction(
        certificate["base_maximum_union_density"]
    )
    base_sum = read_fraction(
        certificate["base_individual_density_sum"]
    )
    extended_union = base_union + (1 - base_union) / extra_h
    extended_sum = base_sum + Fraction(1, extra_h)
    overlap_loss = extended_sum - extended_union
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    upper_bound = total_density - overlap_loss
    verification_anchor_primes = tuple(
        int(p)
        for p in base_verification.get(
            "anchor_primes",
            [
                *base_verification.get("base_primes", ()),
                base_verification.get("extra_prime"),
            ],
        )
        if p is not None
    )
    verification_period = int(
        base_verification.get(
            "period",
            base_verification.get(
                "enumerated_period",
                int(base_verification.get("base_period", 1))
                * int(
                    base_verification.get(
                        "residual_modulus",
                        base_verification.get("extra_modulus", 1),
                    )
                ),
            ),
        )
    )
    verification_union_payload = base_verification.get(
        "maximum_union_density",
        base_verification.get(
            "maximum_block_union_density",
            base_verification.get("extended_maximum_union_density"),
        ),
    )
    if verification_union_payload is None:
        raise RuntimeError("base verification omits its maximum union")
    base_verified = (
        bool(base_verification.get("verified"))
        and verification_anchor_primes == base_primes
        and verification_period == base_period
        and Fraction(
            int(verification_union_payload["numerator"]),
            int(verification_union_payload["denominator"]),
        )
        == base_union
        and Fraction(
            int(base_verification["forced_overlap_loss"]["numerator"]),
            int(base_verification["forced_overlap_loss"]["denominator"]),
        )
        == base_sum - base_union
    )
    extra_matches = all(
        int(extra[key]) == int(extra_record[key])
        for key in ("p", "h", "a", "b")
    )
    verified = (
        str(args.pool) == certificate["pool"]
        and int(certificate["row_count"]) == len(rows)
        and base_verified
        and set(base_primes) <= by_prime.keys()
        and extra_matches
        and math.gcd(extra_h, base_period) == 1
        and math.gcd(int(extra["a"]), int(extra["b"]), extra_h) == 1
        and int(extra.get("target_modulus", 1)) == 1
        and read_fraction(
            certificate["extended_maximum_union_density"]
        )
        == extended_union
        and read_fraction(
            certificate["extended_individual_density_sum"]
        )
        == extended_sum
        and read_fraction(certificate["forced_overlap_loss"])
        == overlap_loss
        and read_fraction(certificate["total_pool_density"])
        == total_density
        and read_fraction(
            certificate["pool_union_density_upper_bound"]
        )
        == upper_bound
        and bool(certificate["proved_no_cover"]) == (upper_bound < 1)
    )
    result = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "base_verification": str(args.base_verification),
        "base_verified": base_verified,
        "base_period": base_period,
        "extra_prime": int(extra["p"]),
        "extra_modulus": extra_h,
        "coprime": math.gcd(extra_h, base_period) == 1,
        "anchor_primes": [*base_primes, int(extra["p"])],
        "enumerated_period": math.lcm(base_period, extra_h),
        "extended_maximum_union_density": {
            "numerator": extended_union.numerator,
            "denominator": extended_union.denominator,
        },
        "maximum_block_union_density": {
            "numerator": extended_union.numerator,
            "denominator": extended_union.denominator,
        },
        "block_individual_density_sum": {
            "numerator": extended_sum.numerator,
            "denominator": extended_sum.denominator,
        },
        "forced_overlap_loss": {
            "numerator": overlap_loss.numerator,
            "denominator": overlap_loss.denominator,
        },
        "pool_union_upper_bound": {
            "numerator": upper_bound.numerator,
            "denominator": upper_bound.denominator,
        },
        "proved_no_cover": upper_bound < 1,
        "verified": verified,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base_verified={base_verified} coprime="
        f"{math.gcd(extra_h, base_period) == 1} "
        f"extended_union={extended_union} loss={overlap_loss} "
        f"pool_upper={upper_bound} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
