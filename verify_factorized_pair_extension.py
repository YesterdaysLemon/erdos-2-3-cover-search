#!/usr/bin/env python3
"""Independently verify a factor-plus-lifted-row block certificate."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def load_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def verification_primes(payload: dict) -> tuple[int, ...]:
    values = payload.get("anchor_primes")
    if values is None:
        values = [
            *payload.get("base_primes", ()),
            payload.get("extra_prime"),
        ]
    return tuple(int(value) for value in values if value is not None)


def verification_period(payload: dict) -> int:
    if payload.get("period") is not None:
        return int(payload["period"])
    if payload.get("enumerated_period") is not None:
        return int(payload["enumerated_period"])
    return int(payload.get("base_period", 1)) * int(
        payload.get(
            "residual_modulus",
            payload.get("extra_modulus", 1),
        )
    )


def verification_union(payload: dict) -> Fraction:
    for key in (
        "maximum_union_density",
        "maximum_block_union_density",
        "extended_maximum_union_density",
    ):
        if payload.get(key) is not None:
            return load_fraction(payload[key])
    raise RuntimeError("base verification omits maximum union density")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("base_verification", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    cert = json.loads(args.certificate.read_text())
    base_check = json.loads(args.base_verification.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    base_primes = tuple(int(value) for value in cert["base_anchor_primes"])
    factor_record = cert["factor_row"]
    lifted_record = cert["lifted_row"]
    factor = by_prime[int(factor_record["p"])]
    lifted = by_prime[int(lifted_record["p"])]

    factor_h = int(factor["h"])
    lifted_h = int(lifted["h"])
    residual = lifted_h // factor_h if lifted_h % factor_h == 0 else 0
    base_period = int(cert["base_period"])
    det = (
        int(factor["a"]) * int(lifted["b"])
        - int(factor["b"]) * int(lifted["a"])
    )
    joint = math.gcd(det, factor_h) == 1
    component_union = (
        Fraction(1, factor_h)
        + Fraction(1, lifted_h)
        - Fraction(1, factor_h * lifted_h)
    )
    base_union = load_fraction(cert["base_maximum_union_density"])
    base_sum = load_fraction(cert["base_individual_density_sum"])
    extended_union = base_union + (1 - base_union) * component_union
    extended_sum = base_sum + Fraction(1, factor_h) + Fraction(1, lifted_h)
    loss = extended_sum - extended_union
    total = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    upper = total - loss

    base_loss = base_sum - base_union
    base_verified = (
        bool(base_check.get("verified"))
        and verification_primes(base_check) == base_primes
        and verification_period(base_check) == base_period
        and verification_union(base_check) == base_union
        and load_fraction(base_check["forced_overlap_loss"]) == base_loss
    )
    row_records_match = all(
        all(
            int(actual[key]) == int(recorded[key])
            for key in ("p", "h", "a", "b")
        )
        for actual, recorded in (
            (factor, factor_record),
            (lifted, lifted_record),
        )
    )
    targets_free = all(
        int(row.get("target_modulus", 1)) == 1
        for row in (factor, lifted)
    )
    maps_surjective = all(
        math.gcd(int(row["a"]), int(row["b"]), int(row["h"])) == 1
        for row in (factor, lifted)
    )
    anchors_present = set(base_primes) <= by_prime.keys()
    verified = (
        str(args.pool) == cert["pool"]
        and int(cert["row_count"]) == len(rows)
        and base_verified
        and anchors_present
        and row_records_match
        and factor_h > 1
        and math.gcd(factor_h, base_period) == 1
        and residual > 1
        and math.gcd(residual, base_period * factor_h) == 1
        and targets_free
        and maps_surjective
        and joint
        and int(cert["shared_modulus"]) == factor_h
        and int(cert["residual_modulus"]) == residual
        and int(cert["shared_determinant"]) == det
        and bool(cert["shared_pair_jointly_surjective"]) == joint
        and load_fraction(cert["pair_component_union_density"])
        == component_union
        and load_fraction(cert["extended_maximum_union_density"])
        == extended_union
        and load_fraction(cert["extended_individual_density_sum"])
        == extended_sum
        and load_fraction(cert["forced_overlap_loss"]) == loss
        and load_fraction(cert["total_pool_density"]) == total
        and load_fraction(cert["pool_union_density_upper_bound"]) == upper
        and bool(cert["proved_no_cover"]) == (upper < 1)
    )
    result = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "base_verification": str(args.base_verification),
        "base_verified": base_verified,
        "base_period": base_period,
        "factor_prime": int(factor["p"]),
        "lifted_prime": int(lifted["p"]),
        "shared_modulus": factor_h,
        "residual_modulus": residual,
        "shared_determinant": det,
        "shared_pair_jointly_surjective": joint,
        "anchor_primes": [
            *base_primes,
            int(factor["p"]),
            int(lifted["p"]),
        ],
        "enumerated_period": base_period * factor_h * residual,
        "maximum_block_union_density": {
            "numerator": extended_union.numerator,
            "denominator": extended_union.denominator,
        },
        "block_individual_density_sum": {
            "numerator": extended_sum.numerator,
            "denominator": extended_sum.denominator,
        },
        "forced_overlap_loss": {
            "numerator": loss.numerator,
            "denominator": loss.denominator,
        },
        "pool_union_upper_bound": {
            "numerator": upper.numerator,
            "denominator": upper.denominator,
        },
        "proved_no_cover": upper < 1,
        "verified": verified,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base_verified={base_verified} shared={factor_h} "
        f"residual={residual} joint={joint} component_union="
        f"{component_union} loss={loss} pool_upper={upper} "
        f"verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
