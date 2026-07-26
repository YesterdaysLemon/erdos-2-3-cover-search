#!/usr/bin/env python3
"""Exhaustively verify a forced-overlap matching certificate."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def enumerate_intersection_density(
    left: dict,
    right: dict,
    max_cells: int,
) -> tuple[Fraction, int, int]:
    h1 = int(left["h"])
    h2 = int(right["h"])
    period = math.lcm(h1, h2)
    cells = period * period
    if cells > max_cells:
        raise RuntimeError(
            f"pair grid has {cells} cells, above guard {max_cells}"
        )
    counts = {}
    for k in range(period):
        for l in range(period):
            target_pair = (
                (int(left["a"]) * k + int(left["b"]) * l) % h1,
                (int(right["a"]) * k + int(right["b"]) * l) % h2,
            )
            counts[target_pair] = counts.get(target_pair, 0) + 1
    if len(counts) != h1 * h2 or len(set(counts.values())) != 1:
        return Fraction(0), period, cells
    return Fraction(next(iter(counts.values())), cells), period, cells


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-pair-cells", type=int, default=10_000_000)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    certificate = json.loads(args.certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    used = set()
    overlap_sum = Fraction(0)
    pair_checks = []
    pairs_valid = True
    for pair in certificate["selected_disjoint_pairs"]:
        first_prime, second_prime = (int(p) for p in pair["primes"])
        if (
            first_prime == second_prime
            or first_prime in used
            or second_prime in used
            or first_prime not in by_prime
            or second_prime not in by_prime
        ):
            pairs_valid = False
            continue
        used.update((first_prime, second_prime))
        density, period, cells = enumerate_intersection_density(
            by_prime[first_prime],
            by_prime[second_prime],
            args.max_pair_cells,
        )
        claimed = read_fraction(pair["forced_intersection_density"])
        valid = density == claimed and density > 0
        pairs_valid &= valid
        overlap_sum += density
        pair_checks.append(
            {
                "primes": [first_prime, second_prime],
                "period": period,
                "cells": cells,
                "intersection_density": {
                    "numerator": density.numerator,
                    "denominator": density.denominator,
                },
                "valid": valid,
            }
        )
    upper_bound = total_density - overlap_sum
    verified = (
        str(args.pool) == certificate["pool"]
        and int(certificate["row_count"]) == len(rows)
        and pairs_valid
        and read_fraction(certificate["total_reciprocal_density"])
        == total_density
        and read_fraction(certificate["forced_overlap_sum"]) == overlap_sum
        and read_fraction(certificate["union_density_upper_bound"])
        == upper_bound
        and certificate["proved_no_cover"]
        and upper_bound < 1
    )
    result = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "row_count": len(rows),
        "pair_checks": pair_checks,
        "pairs_disjoint_and_valid": pairs_valid,
        "total_density": {
            "numerator": total_density.numerator,
            "denominator": total_density.denominator,
        },
        "forced_overlap_sum": {
            "numerator": overlap_sum.numerator,
            "denominator": overlap_sum.denominator,
        },
        "union_upper_bound": {
            "numerator": upper_bound.numerator,
            "denominator": upper_bound.denominator,
        },
        "verified": verified,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"pairs={len(pair_checks)} overlap={overlap_sum} "
        f"upper={upper_bound} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
