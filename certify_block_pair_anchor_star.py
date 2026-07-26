#!/usr/bin/env python3
"""Strengthen block-to-row forest edges with two anchor witnesses.

For an outside row C and two block anchors A and B,

  |C intersect (A union B)|
    >= |C intersect A| + |C intersect B| - |C intersect A intersect B|.

When each pair target map is surjective, the first two densities are fixed.
The three-target map is a homomorphism to a product of three cyclic groups.
Its cokernel index is the gcd of the maximal minors of a presentation matrix,
so every nonempty triple fibre has density index/(h_C h_A h_B).  This gives
a phase-independent lower bound for C's intersection with the full block.
Using each outside row once produces a star, hence a valid forest.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

from certify_block_plus_overlap_forest import (
    fraction_payload,
    read_fraction,
    recorded_anchor_rows,
)
from certify_forced_pair_overlap import joint_index


def determinant3(
    first: tuple[int, int, int],
    second: tuple[int, int, int],
    third: tuple[int, int, int],
) -> int:
    return (
        first[0] * (second[1] * third[2] - second[2] * third[1])
        - second[0] * (first[1] * third[2] - first[2] * third[1])
        + third[0] * (first[1] * second[2] - first[2] * second[1])
    )


def triple_image_index_minors(rows: tuple[dict, dict, dict]) -> int:
    moduli = [int(row["h"]) for row in rows]
    columns = [
        (moduli[0], 0, 0),
        (0, moduli[1], 0),
        (0, 0, moduli[2]),
        tuple(int(row["a"]) for row in rows),
        tuple(int(row["b"]) for row in rows),
    ]
    minors = [
        abs(determinant3(*selected))
        for selected in itertools.combinations(columns, 3)
    ]
    index = math.gcd(*minors)
    if index < 1:
        raise RuntimeError("three-target map has infinite cokernel")
    return index


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
    anchors = recorded_anchor_rows(block)
    for recorded in anchors:
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
    candidate_edges = []
    for outside in rows:
        outside_prime = int(outside["p"])
        if outside_prime in anchor_primes:
            continue
        compatible = [
            anchor
            for anchor in anchors
            if joint_index(outside, anchor) == 1
        ]
        best_lower = Fraction(0)
        best_anchors: tuple[dict, ...] = ()
        best_triple_index = None
        for anchor in compatible:
            lower = Fraction(
                1,
                int(outside["h"]) * int(anchor["h"]),
            )
            if lower > best_lower:
                best_lower = lower
                best_anchors = (anchor,)
                best_triple_index = None
        for first, second in itertools.combinations(compatible, 2):
            index = triple_image_index_minors(
                (outside, first, second)
            )
            triple_upper = Fraction(
                index,
                int(outside["h"])
                * int(first["h"])
                * int(second["h"]),
            )
            lower = (
                Fraction(
                    1,
                    int(outside["h"]) * int(first["h"]),
                )
                + Fraction(
                    1,
                    int(outside["h"]) * int(second["h"]),
                )
                - triple_upper
            )
            if (
                lower > best_lower
                or (
                    lower == best_lower
                    and tuple(
                        sorted((int(first["p"]), int(second["p"])))
                    )
                    < tuple(
                        sorted(int(row["p"]) for row in best_anchors)
                    )
                )
            ):
                best_lower = lower
                best_anchors = (first, second)
                best_triple_index = index
        if best_lower <= 0:
            continue
        candidate_edges.append(
            (
                best_lower,
                outside_prime,
                outside,
                best_anchors,
                best_triple_index,
            )
        )

    candidate_edges.sort(key=lambda edge: (-edge[0], edge[1]))
    selected_edges = []
    star_overlap = Fraction(0)
    for lower, outside_prime, outside, witnesses, triple_index in candidate_edges:
        star_overlap += lower
        witness_payload = [
            {
                key: int(row[key]) for key in ("p", "h", "a", "b")
            }
            for row in witnesses
        ]
        edge = {
            "outside_prime": outside_prime,
            "outside_row": {
                key: int(outside[key]) for key in ("p", "h", "a", "b")
            },
            "witness_anchor_primes": [
                int(row["p"]) for row in witnesses
            ],
            "witness_anchor_rows": witness_payload,
            "block_intersection_lower_bound": fraction_payload(lower),
        }
        if len(witnesses) == 2:
            first, second = witnesses
            edge["pair_intersection_densities"] = [
                fraction_payload(
                    Fraction(
                        1,
                        int(outside["h"]) * int(anchor["h"]),
                    )
                )
                for anchor in witnesses
            ]
            edge["triple_image_index"] = int(triple_index)
            edge["triple_intersection_upper_bound"] = fraction_payload(
                Fraction(
                    int(triple_index),
                    int(outside["h"])
                    * int(first["h"])
                    * int(second["h"]),
                )
            )
        selected_edges.append(edge)
        if total_density - block_loss - star_overlap < 1:
            break

    total_loss = block_loss + star_overlap
    upper_bound = total_density - total_loss
    proved = upper_bound < 1
    result = {
        "pool": str(args.pool),
        "row_count": len(rows),
        "block_certificate": str(args.block_certificate),
        "block_anchor_primes": sorted(anchor_primes),
        "block_overlap_loss": fraction_payload(block_loss),
        "selected_star_edges": selected_edges,
        "star_overlap_sum": fraction_payload(star_overlap),
        "total_forced_overlap_loss": fraction_payload(total_loss),
        "total_pool_density": fraction_payload(total_density),
        "pool_union_density_upper_bound": fraction_payload(upper_bound),
        "proved_no_cover": proved,
        "argument": (
            "Each selected outside row is connected once to the exact block, "
            "so the selected graph is a star forest. A one-anchor edge uses "
            "its exact phase-independent pair intersection. A two-anchor "
            "edge uses the second Bonferroni lower bound, with the maximum "
            "triple-fibre density certified by the cokernel index of the "
            "three-target homomorphism."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"anchors={len(anchor_primes)} edges={len(selected_edges)} "
        f"block_loss={block_loss} star_loss={star_overlap} "
        f"upper={upper_bound} proved_no_cover={proved}",
        flush=True,
    )
    return 0 if proved else 1


if __name__ == "__main__":
    raise SystemExit(main())
