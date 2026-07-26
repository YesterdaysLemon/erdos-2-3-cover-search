#!/usr/bin/env python3
"""Combine an exact anchor-block loss with an overlap forest.

For sets indexed by the vertices of a forest F,

    |union A_v| <= sum_v |A_v| - sum_{uv in F} |A_u intersect A_v|.

Indeed, at any point contained in r of the sets, the induced subgraph is a
forest and therefore has at most r-1 edges.  Here the anchor block is first
unioned exactly (or upper-bounded by its exact maximum), and the forest uses
only rows outside that block.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from certify_forced_pair_overlap import joint_index


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def recorded_anchor_rows(block: dict) -> list[dict]:
    if block.get("anchor_rows"):
        return list(block["anchor_rows"])
    if block.get("base_certificate"):
        base = json.loads(Path(block["base_certificate"]).read_text())
        return [*recorded_anchor_rows(base), block["extra_row"]]
    raise RuntimeError("block certificate does not expose its anchor rows")


class DisjointSet:
    def __init__(self, vertices: list[int]) -> None:
        self.parent = {vertex: vertex for vertex in vertices}
        self.rank = {vertex: 0 for vertex in vertices}

    def find(self, vertex: int) -> int:
        root = vertex
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[vertex] != vertex:
            parent = self.parent[vertex]
            self.parent[vertex] = root
            vertex = parent
        return root

    def union(self, left: int, right: int) -> bool:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return False
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1
        return True


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
    for recorded in recorded_anchor_rows(block):
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
    edges = sorted(
        (
            (
                Fraction(1, int(left["h"]) * int(right["h"])),
                min(int(left["p"]), int(right["p"])),
                max(int(left["p"]), int(right["p"])),
                left,
                right,
            )
            for index, left in enumerate(available)
            for right in available[:index]
            if joint_index(left, right) == 1
        ),
        key=lambda edge: (-edge[0], edge[1], edge[2]),
    )
    dsu = DisjointSet([int(row["p"]) for row in available])
    selected_edges = []
    forest_overlap = Fraction(0)
    for weight, _first, _second, left, right in edges:
        left_prime = int(left["p"])
        right_prime = int(right["p"])
        if not dsu.union(left_prime, right_prime):
            continue
        forest_overlap += weight
        selected_edges.append(
            {
                "primes": [left_prime, right_prime],
                "rows": [
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
            "The block anchors are first unioned with their certified maximum "
            "density. The remaining selected phase-independent pair "
            "intersections form a forest. At any point lying in r outside "
            "rows, the induced forest has at most r-1 edges, so subtracting "
            "every selected edge intersection from the sum of row densities "
            "is a valid pointwise union upper bound."
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
