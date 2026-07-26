#!/usr/bin/env python3
"""Independently replay a block-plus-overlap-forest certificate."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def enumerate_pair_density(
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
    counts: dict[tuple[int, int], int] = {}
    for k in range(period):
        for l in range(period):
            key = (
                (int(left["a"]) * k + int(left["b"]) * l) % h1,
                (int(right["a"]) * k + int(right["b"]) * l) % h2,
            )
            counts[key] = counts.get(key, 0) + 1
    if len(counts) != h1 * h2 or len(set(counts.values())) != 1:
        return Fraction(0), period, cells
    density = Fraction(next(iter(counts.values())), cells)
    return density, period, cells


class CycleChecker:
    def __init__(self, vertices: set[int]) -> None:
        self.parent = {vertex: vertex for vertex in vertices}

    def find(self, vertex: int) -> int:
        while self.parent[vertex] != vertex:
            vertex = self.parent[vertex]
        return vertex

    def add_edge(self, left: int, right: int) -> bool:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return False
        self.parent[right_root] = left_root
        return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("block_verification", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-pair-cells", type=int, default=10_000_000)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    certificate = json.loads(args.certificate.read_text())
    block_verification = json.loads(args.block_verification.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    block_primes = {
        int(p) for p in certificate["block_anchor_primes"]
    }
    outside_primes = set(by_prime) - block_primes
    block_loss = read_fraction(certificate["block_overlap_loss"])
    block_verified = (
        bool(block_verification.get("verified"))
        and block_verification.get("certificate")
        == certificate["block_certificate"]
        and read_fraction(block_verification["forced_overlap_loss"])
        == block_loss
    )

    cycle_checker = CycleChecker(outside_primes)
    forest_valid = True
    forest_overlap = Fraction(0)
    edge_checks = []
    seen_edges: set[tuple[int, int]] = set()
    for edge in certificate["selected_forest_edges"]:
        first, second = (int(p) for p in edge["primes"])
        canonical = tuple(sorted((first, second)))
        endpoints_valid = (
            first != second
            and first in outside_primes
            and second in outside_primes
            and canonical not in seen_edges
        )
        seen_edges.add(canonical)
        acyclic = (
            cycle_checker.add_edge(first, second)
            if endpoints_valid
            else False
        )
        if endpoints_valid:
            density, period, cells = enumerate_pair_density(
                by_prime[first],
                by_prime[second],
                args.max_pair_cells,
            )
        else:
            density, period, cells = Fraction(0), 0, 0
        claimed = read_fraction(edge["forced_intersection_density"])
        valid = endpoints_valid and acyclic and density == claimed and density > 0
        forest_valid &= valid
        forest_overlap += density
        edge_checks.append(
            {
                "primes": [first, second],
                "period": period,
                "cells": cells,
                "intersection_density": {
                    "numerator": density.numerator,
                    "denominator": density.denominator,
                },
                "endpoints_valid": endpoints_valid,
                "acyclic": acyclic,
                "valid": valid,
            }
        )

    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    total_loss = block_loss + forest_overlap
    upper_bound = total_density - total_loss
    verified = (
        str(args.pool) == certificate["pool"]
        and int(certificate["row_count"]) == len(rows)
        and block_primes <= by_prime.keys()
        and block_verified
        and forest_valid
        and read_fraction(certificate["forest_overlap_sum"])
        == forest_overlap
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
        "block_verified": block_verified,
        "forest_valid": forest_valid,
        "edge_checks": edge_checks,
        "block_overlap_loss": {
            "numerator": block_loss.numerator,
            "denominator": block_loss.denominator,
        },
        "forest_overlap_sum": {
            "numerator": forest_overlap.numerator,
            "denominator": forest_overlap.denominator,
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
        f"forest_valid={forest_valid} forest_loss={forest_overlap} "
        f"upper={upper_bound} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
