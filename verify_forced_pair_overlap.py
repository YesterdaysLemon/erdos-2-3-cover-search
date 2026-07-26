#!/usr/bin/env python3
"""Independently verify a forced-pair-overlap no-cover certificate."""

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
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-enumerated-cells", type=int, default=10_000_000)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    certificate = json.loads(args.certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    first_prime, second_prime = (
        int(value) for value in certificate["anchor_primes"]
    )
    left = by_prime[first_prime]
    right = by_prime[second_prime]
    h1 = int(left["h"])
    h2 = int(right["h"])
    period = math.lcm(h1, h2)
    cells = period * period
    if cells > args.max_enumerated_cells:
        raise RuntimeError(
            f"anchor grid has {cells} cells, above enumeration guard"
        )

    image_counts = {}
    for k in range(period):
        for l in range(period):
            target_pair = (
                (int(left["a"]) * k + int(left["b"]) * l) % h1,
                (int(right["a"]) * k + int(right["b"]) * l) % h2,
            )
            image_counts[target_pair] = image_counts.get(target_pair, 0) + 1
    expected_targets = h1 * h2
    jointly_surjective = len(image_counts) == expected_targets
    uniform_fibre_size = len(set(image_counts.values())) == 1
    intersection_density = (
        Fraction(next(iter(image_counts.values())), cells)
        if jointly_surjective and uniform_fibre_size
        else Fraction(0)
    )
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    upper_bound = total_density - intersection_density
    verified = (
        str(args.pool) == certificate["pool"]
        and int(certificate["row_count"]) == len(rows)
        and read_fraction(certificate["total_reciprocal_density"])
        == total_density
        and read_fraction(certificate["forced_pair_overlap_density"])
        == intersection_density
        and read_fraction(certificate["union_density_upper_bound"])
        == upper_bound
        and certificate["joint_target_map_surjective"]
        == jointly_surjective
        and certificate["proved_no_cover"]
        and jointly_surjective
        and uniform_fibre_size
        and upper_bound < 1
    )
    result = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "enumerated_period": period,
        "enumerated_cells": cells,
        "target_pairs": len(image_counts),
        "expected_target_pairs": expected_targets,
        "jointly_surjective": jointly_surjective,
        "uniform_fibre_size": uniform_fibre_size,
        "intersection_density": {
            "numerator": intersection_density.numerator,
            "denominator": intersection_density.denominator,
        },
        "total_density": {
            "numerator": total_density.numerator,
            "denominator": total_density.denominator,
        },
        "union_upper_bound": {
            "numerator": upper_bound.numerator,
            "denominator": upper_bound.denominator,
        },
        "verified": verified,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"period={period} cells={cells} targets={len(image_counts)}/"
        f"{expected_targets} intersection={intersection_density} "
        f"upper={upper_bound} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
