#!/usr/bin/env python3
"""Certify a union bound with an exact block represented as one forest node.

The anchor rows are replaced by their union B, whose maximum density is
provided by a separately verified block certificate.  For an outside row A,
the intersection B intersect A contains the intersection of A with any one
anchor row.  A phase-independent anchor/outside pair therefore supplies a
valid lower bound for a B--A forest edge.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from certify_block_plus_overlap_forest import (
    DisjointSet,
    fraction_payload,
    read_fraction,
    recorded_anchor_rows,
)
from certify_forced_pair_overlap import joint_index


BLOCK_VERTEX = 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("block_certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    block = json.loads(args.block_certificate.read_text())
    rows = source["choices"]
    if len({int(row["p"]) for row in rows}) != len(rows):
        raise RuntimeError("pool contains a repeated fibre prime")
    if any(int(row.get("target_modulus", 1)) != 1 for row in rows):
        raise RuntimeError("certificate requires unrestricted row targets")
    by_prime = {int(row["p"]): row for row in rows}
    anchor_primes = {int(p) for p in block["anchor_primes"]}
    anchor_rows = recorded_anchor_rows(block)
    for recorded in anchor_rows:
        prime = int(recorded["p"])
        if prime not in by_prime or any(
            int(by_prime[prime][key]) != int(recorded[key])
            for key in ("h", "a", "b")
        ):
            raise RuntimeError("block anchor differs in the supplied pool")

    block_loss = read_fraction(block["forced_overlap_loss"])
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    available = [
        row for row in rows if int(row["p"]) not in anchor_primes
    ]
    candidate_edges = []
    for index, left in enumerate(available):
        for right in available[:index]:
            if joint_index(left, right) != 1:
                continue
            weight = Fraction(1, int(left["h"]) * int(right["h"]))
            candidate_edges.append(
                (
                    weight,
                    "outside",
                    min(int(left["p"]), int(right["p"])),
                    max(int(left["p"]), int(right["p"])),
                    left,
                    right,
                )
            )
    for outside in available:
        compatible = [
            anchor
            for anchor in anchor_rows
            if joint_index(anchor, outside) == 1
        ]
        if not compatible:
            continue
        anchor = min(
            compatible,
            key=lambda row: (int(row["h"]), int(row["p"])),
        )
        weight = Fraction(1, int(anchor["h"]) * int(outside["h"]))
        candidate_edges.append(
            (
                weight,
                "block",
                int(anchor["p"]),
                int(outside["p"]),
                anchor,
                outside,
            )
        )
    candidate_edges.sort(
        key=lambda edge: (
            -edge[0],
            edge[1],
            edge[2],
            edge[3],
        )
    )

    dsu = DisjointSet(
        [BLOCK_VERTEX, *(int(row["p"]) for row in available)]
    )
    selected_edges = []
    forest_overlap = Fraction(0)
    for weight, kind, first, second, left, right in candidate_edges:
        graph_left = BLOCK_VERTEX if kind == "block" else int(left["p"])
        graph_right = int(right["p"])
        if not dsu.union(graph_left, graph_right):
            continue
        forest_overlap += weight
        selected_edges.append(
            {
                "kind": kind,
                "graph_vertices": [
                    "block" if kind == "block" else int(left["p"]),
                    int(right["p"]),
                ],
                "witness_primes": [int(left["p"]), int(right["p"])],
                "witness_rows": [
                    {
                        key: int(row[key])
                        for key in ("p", "h", "a", "b")
                    }
                    for row in (left, right)
                ],
                "forced_intersection_density": fraction_payload(weight),
            }
        )
        if total_density - block_loss - forest_overlap < 1:
            break

    total_loss = block_loss + forest_overlap
    upper_bound = total_density - total_loss
    proved = upper_bound < 1
    result = {
        "pool": str(args.pool),
        "row_count": len(rows),
        "block_certificate": str(args.block_certificate),
        "block_anchor_primes": sorted(anchor_primes),
        "block_overlap_loss": fraction_payload(block_loss),
        "selected_forest_edges": selected_edges,
        "forest_overlap_sum": fraction_payload(forest_overlap),
        "total_forced_overlap_loss": fraction_payload(total_loss),
        "total_pool_density": fraction_payload(total_density),
        "pool_union_density_upper_bound": fraction_payload(upper_bound),
        "proved_no_cover": proved,
        "argument": (
            "The exact anchor union is one forest vertex. Each block edge "
            "uses a phase-independent intersection between its outside row "
            "and one recorded anchor, which is a subset of the outside "
            "row's intersection with the anchor union. All graph edges form "
            "a forest, so their lower-bound intersection densities may be "
            "subtracted pointwise from the sum of node densities."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"anchors={len(anchor_primes)} edges={len(selected_edges)} "
        f"block_loss={block_loss} forest_loss={forest_overlap} "
        f"upper={upper_bound} proved_no_cover={proved}",
        flush=True,
    )
    return 0 if proved else 1


if __name__ == "__main__":
    raise SystemExit(main())
