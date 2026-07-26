#!/usr/bin/env python3
"""Independently verify a block-plus-matching union bound."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def enumerate_pair_density(left: dict, right: dict) -> Fraction:
    h1 = int(left["h"])
    h2 = int(right["h"])
    period = math.lcm(h1, h2)
    counts = {}
    for k in range(period):
        for l in range(period):
            key = (
                (int(left["a"]) * k + int(left["b"]) * l) % h1,
                (int(right["a"]) * k + int(right["b"]) * l) % h2,
            )
            counts[key] = counts.get(key, 0) + 1
    if len(counts) != h1 * h2 or len(set(counts.values())) != 1:
        return Fraction(0)
    return Fraction(next(iter(counts.values())), period * period)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("block_verification", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    certificate = json.loads(args.certificate.read_text())
    block_verification = json.loads(args.block_verification.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    block_primes = {
        int(p) for p in certificate["block_anchor_primes"]
    }
    block_loss = read_fraction(certificate["block_overlap_loss"])
    block_verified = (
        bool(block_verification.get("verified"))
        and block_verification.get("certificate")
        == certificate["block_certificate"]
        and read_fraction(block_verification["forced_overlap_loss"])
        == block_loss
    )
    used = set(block_primes)
    pair_loss = Fraction(0)
    pair_checks = []
    pairs_valid = True
    for pair in certificate["selected_disjoint_pairs"]:
        first, second = (int(p) for p in pair["primes"])
        if (
            first == second
            or first in used
            or second in used
            or first not in by_prime
            or second not in by_prime
        ):
            pairs_valid = False
            continue
        used.update((first, second))
        density = enumerate_pair_density(
            by_prime[first], by_prime[second]
        )
        claimed = read_fraction(pair["forced_intersection_density"])
        valid = density == claimed and density > 0
        pairs_valid &= valid
        pair_loss += density
        pair_checks.append(
            {
                "primes": [first, second],
                "density": {
                    "numerator": density.numerator,
                    "denominator": density.denominator,
                },
                "valid": valid,
            }
        )
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    total_loss = block_loss + pair_loss
    upper_bound = total_density - total_loss
    verified = (
        str(args.pool) == certificate["pool"]
        and int(certificate["row_count"]) == len(rows)
        and block_primes <= by_prime.keys()
        and block_verified
        and pairs_valid
        and read_fraction(certificate["pair_overlap_sum"]) == pair_loss
        and read_fraction(certificate["total_forced_overlap_loss"])
        == total_loss
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
        "block_verification": str(args.block_verification),
        "block_verified": block_verified,
        "pairs_valid": pairs_valid,
        "pair_checks": pair_checks,
        "block_overlap_loss": {
            "numerator": block_loss.numerator,
            "denominator": block_loss.denominator,
        },
        "pair_overlap_sum": {
            "numerator": pair_loss.numerator,
            "denominator": pair_loss.denominator,
        },
        "total_forced_overlap_loss": {
            "numerator": total_loss.numerator,
            "denominator": total_loss.denominator,
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
        f"block_verified={block_verified} pairs={len(pair_checks)} "
        f"pair_loss={pair_loss} upper={upper_bound} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
