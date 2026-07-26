#!/usr/bin/env python3
"""Independently replay a block-node overlap-forest certificate."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from verify_block_plus_overlap_forest import (
    CycleChecker,
    enumerate_pair_density,
    read_fraction,
)


BLOCK_VERTEX = 0


def recorded_anchor_rows(block: dict) -> list[dict]:
    if block.get("anchor_rows"):
        return list(block["anchor_rows"])
    if block.get("base_certificate"):
        base = json.loads(Path(block["base_certificate"]).read_text())
        return [*recorded_anchor_rows(base), block["extra_row"]]
    raise RuntimeError("block certificate does not expose its anchor rows")


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
    block_certificate = json.loads(
        Path(certificate["block_certificate"]).read_text()
    )
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    block_primes = {
        int(p) for p in certificate["block_anchor_primes"]
    }
    outside_primes = set(by_prime) - block_primes
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
        and anchors_match
    )

    cycle_checker = CycleChecker({BLOCK_VERTEX, *outside_primes})
    forest_valid = True
    forest_overlap = Fraction(0)
    edge_checks = []
    seen_graph_edges: set[tuple[int, int]] = set()
    for edge in certificate["selected_forest_edges"]:
        kind = edge["kind"]
        witness_first, witness_second = (
            int(p) for p in edge["witness_primes"]
        )
        if kind == "block":
            graph_first = BLOCK_VERTEX
            graph_second = witness_second
            endpoints_valid = (
                witness_first in block_primes
                and witness_second in outside_primes
            )
        elif kind == "outside":
            graph_first = witness_first
            graph_second = witness_second
            endpoints_valid = (
                witness_first in outside_primes
                and witness_second in outside_primes
                and witness_first != witness_second
            )
        else:
            graph_first = graph_second = BLOCK_VERTEX
            endpoints_valid = False
        canonical = tuple(sorted((graph_first, graph_second)))
        endpoints_valid &= canonical not in seen_graph_edges
        seen_graph_edges.add(canonical)
        acyclic = (
            cycle_checker.add_edge(graph_first, graph_second)
            if endpoints_valid
            else False
        )
        if endpoints_valid:
            density, period, cells = enumerate_pair_density(
                by_prime[witness_first],
                by_prime[witness_second],
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
                "kind": kind,
                "graph_vertices": [graph_first, graph_second],
                "witness_primes": [witness_first, witness_second],
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
        "anchors_match": anchors_match,
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
