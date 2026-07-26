#!/usr/bin/env python3
"""Independently verify a two-anchor block-star certificate."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def recorded_anchor_rows(block: dict) -> list[dict]:
    if block.get("anchor_rows"):
        return list(block["anchor_rows"])
    if block.get("base_certificate"):
        base = json.loads(Path(block["base_certificate"]).read_text())
        return [*recorded_anchor_rows(base), block["extra_row"]]
    raise RuntimeError("block certificate does not expose its anchor rows")


def pair_cokernel_index(first: dict, second: dict) -> int:
    h1, a1, b1 = (int(first[key]) for key in ("h", "a", "b"))
    h2, a2, b2 = (int(second[key]) for key in ("h", "a", "b"))
    presentation_minors = (
        h1 * h2,
        h1 * a2,
        h1 * b2,
        h2 * a1,
        h2 * b1,
        a1 * b2 - a2 * b1,
    )
    return math.gcd(*(abs(value) for value in presentation_minors))


def explicit_determinant(
    first: tuple[int, int, int],
    second: tuple[int, int, int],
    third: tuple[int, int, int],
) -> int:
    a, b, c = first
    d, e, f = second
    g, h, i = third
    return a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h


def triple_index_explicit(rows: tuple[dict, dict, dict]) -> int:
    h0, h1, h2 = (int(row["h"]) for row in rows)
    a0, a1, a2 = (int(row["a"]) for row in rows)
    b0, b1, b2 = (int(row["b"]) for row in rows)
    columns = (
        (h0, 0, 0),
        (0, h1, 0),
        (0, 0, h2),
        (a0, a1, a2),
        (b0, b1, b2),
    )
    determinants = []
    for indices in (
        (0, 1, 2),
        (0, 1, 3),
        (0, 1, 4),
        (0, 2, 3),
        (0, 2, 4),
        (0, 3, 4),
        (1, 2, 3),
        (1, 2, 4),
        (1, 3, 4),
        (2, 3, 4),
    ):
        determinants.append(
            abs(
                explicit_determinant(
                    columns[indices[0]],
                    columns[indices[1]],
                    columns[indices[2]],
                )
            )
        )
    index = math.gcd(*determinants)
    if index < 1:
        raise RuntimeError("three-target map has infinite cokernel")
    return index


def exhaustive_triple_index(
    rows: tuple[dict, dict, dict],
    max_cells: int,
) -> tuple[int | None, int]:
    period = math.lcm(*(int(row["h"]) for row in rows))
    cells = period * period
    if cells > max_cells:
        return None, cells
    image = {
        tuple(
            (
                int(row["a"]) * k + int(row["b"]) * l
            )
            % int(row["h"])
            for row in rows
        )
        for k in range(period)
        for l in range(period)
    }
    product = math.prod(int(row["h"]) for row in rows)
    if product % len(image):
        raise AssertionError("image order does not divide target-group order")
    return product // len(image), cells


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("block_verification", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-exhaustive-cells", type=int, default=1_000_000)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    certificate = json.loads(args.certificate.read_text())
    block_verification = json.loads(args.block_verification.read_text())
    block_certificate = json.loads(
        Path(certificate["block_certificate"]).read_text()
    )
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    block_primes = {
        int(p) for p in certificate["block_anchor_primes"]
    }
    recorded_anchors = {
        int(row["p"]): row for row in recorded_anchor_rows(block_certificate)
    }
    anchors_match = (
        set(recorded_anchors) == block_primes
        and all(
            prime in by_prime
            and all(
                int(by_prime[prime][key]) == int(recorded[key])
                for key in ("h", "a", "b")
            )
            for prime, recorded in recorded_anchors.items()
        )
    )
    block_loss = read_fraction(certificate["block_overlap_loss"])
    block_verified = (
        bool(block_verification.get("verified"))
        and block_verification.get("certificate")
        == certificate["block_certificate"]
        and read_fraction(block_verification["forced_overlap_loss"])
        == block_loss
        and set(int(p) for p in block_verification["anchor_primes"])
        == block_primes
        and anchors_match
    )

    seen_outside: set[int] = set()
    star_valid = True
    star_overlap = Fraction(0)
    edge_checks = []
    for edge in certificate["selected_star_edges"]:
        outside_prime = int(edge["outside_prime"])
        witness_primes = tuple(
            int(p) for p in edge["witness_anchor_primes"]
        )
        endpoints_valid = (
            outside_prime in by_prime
            and outside_prime not in block_primes
            and outside_prime not in seen_outside
            and len(witness_primes) in (1, 2)
            and len(set(witness_primes)) == len(witness_primes)
            and set(witness_primes) <= block_primes
        )
        seen_outside.add(outside_prime)
        if endpoints_valid:
            outside = by_prime[outside_prime]
            witnesses = tuple(by_prime[p] for p in witness_primes)
            pair_indices = tuple(
                pair_cokernel_index(outside, anchor)
                for anchor in witnesses
            )
            pair_valid = all(index == 1 for index in pair_indices)
            pair_densities = tuple(
                Fraction(
                    1,
                    int(outside["h"]) * int(anchor["h"]),
                )
                for anchor in witnesses
            )
            if len(witnesses) == 2:
                triple_rows = (outside, *witnesses)
                triple_index = triple_index_explicit(triple_rows)
                exhaustive_index, exhaustive_cells = (
                    exhaustive_triple_index(
                        triple_rows,
                        args.max_exhaustive_cells,
                    )
                )
                exhaustive_agrees = (
                    exhaustive_index is None
                    or exhaustive_index == triple_index
                )
                triple_upper = Fraction(
                    triple_index,
                    math.prod(int(row["h"]) for row in triple_rows),
                )
                lower = sum(pair_densities, Fraction(0)) - triple_upper
                recorded_triple_valid = (
                    int(edge["triple_image_index"]) == triple_index
                    and read_fraction(
                        edge["triple_intersection_upper_bound"]
                    )
                    == triple_upper
                )
            else:
                triple_index = None
                exhaustive_index = None
                exhaustive_cells = 0
                exhaustive_agrees = True
                triple_upper = Fraction(0)
                lower = pair_densities[0]
                recorded_triple_valid = (
                    "triple_image_index" not in edge
                    and "triple_intersection_upper_bound" not in edge
                )
        else:
            pair_indices = ()
            pair_valid = False
            triple_index = None
            exhaustive_index = None
            exhaustive_cells = 0
            exhaustive_agrees = False
            triple_upper = Fraction(0)
            lower = Fraction(0)
            recorded_triple_valid = False
        claimed = read_fraction(edge["block_intersection_lower_bound"])
        valid = (
            endpoints_valid
            and pair_valid
            and recorded_triple_valid
            and exhaustive_agrees
            and lower == claimed
            and lower > 0
        )
        star_valid &= valid
        star_overlap += lower
        edge_checks.append(
            {
                "outside_prime": outside_prime,
                "witness_anchor_primes": list(witness_primes),
                "pair_cokernel_indices": list(pair_indices),
                "triple_image_index": triple_index,
                "exhaustive_index": exhaustive_index,
                "exhaustive_cells": exhaustive_cells,
                "lower_bound": {
                    "numerator": lower.numerator,
                    "denominator": lower.denominator,
                },
                "valid": valid,
            }
        )

    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    total_loss = block_loss + star_overlap
    upper_bound = total_density - total_loss
    verified = (
        str(args.pool) == certificate["pool"]
        and int(certificate["row_count"]) == len(rows)
        and block_verified
        and star_valid
        and read_fraction(certificate["star_overlap_sum"])
        == star_overlap
        and read_fraction(certificate["total_forced_overlap_loss"])
        == total_loss
        and read_fraction(certificate["total_pool_density"])
        == total_density
        and read_fraction(
            certificate["pool_union_density_upper_bound"]
        )
        == upper_bound
        and bool(certificate["proved_no_cover"]) == (upper_bound < 1)
        and upper_bound < 1
    )
    result = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "block_verification": str(args.block_verification),
        "anchors_match": anchors_match,
        "block_verified": block_verified,
        "star_valid": star_valid,
        "edge_checks": edge_checks,
        "block_overlap_loss": {
            "numerator": block_loss.numerator,
            "denominator": block_loss.denominator,
        },
        "star_overlap_sum": {
            "numerator": star_overlap.numerator,
            "denominator": star_overlap.denominator,
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
        f"block_verified={block_verified} edges={len(edge_checks)} "
        f"star_valid={star_valid} star_loss={star_overlap} "
        f"upper={upper_bound} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
