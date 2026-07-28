#!/usr/bin/env python3
"""Independently verify a finite-group density/overlap obstruction.

For every source row descending to the declared quotient, a phase selects a
set of density 1/h.  If two affine maps are jointly surjective onto their
two cyclic codomains, every pair of phases intersects in exact density
1/(h_i h_j).  Hunter's union bound permits subtracting intersections on any
forest.  This verifier independently computes the maximum possible forest
weight.  If the resulting upper bound is below one, no assignment of phases
can cover the declared finite group.

This checks only that finite-family statement.  It makes no global claim
about the Erdos problem and does not construct an integer m.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

from build_axis_layered_pool import sha256


def transformed_row(
    row: dict,
    direction: tuple[int, int],
    transverse: tuple[int, int],
) -> tuple[int, int, int]:
    h = int(row["h"])
    if h < 1:
        raise AssertionError("row modulus must be positive")
    source_a = int(row["a"]) % h
    source_b = int(row["b"]) % h
    if math.gcd(source_a, source_b, h) != 1:
        raise AssertionError("source row is not surjective")
    a = (
        source_a * direction[0] + source_b * direction[1]
    ) % h
    b = (
        source_a * transverse[0] + source_b * transverse[1]
    ) % h
    if math.gcd(a, b, h) != 1:
        raise AssertionError("basis transform lost surjectivity")
    return h, a, b


def joint_image_index(left: dict, right: dict) -> int:
    """Compute the joint-image index from the presentation minors."""

    left_h = int(left["h"])
    right_h = int(right["h"])
    left_a = int(left["a"]) % left_h
    left_b = int(left["b"]) % left_h
    right_a = int(right["a"]) % right_h
    right_b = int(right["b"]) % right_h
    minors = (
        left_a * right_b - left_b * right_a,
        -left_h * right_a,
        -left_h * right_b,
        right_h * left_a,
        right_h * left_b,
        left_h * right_h,
    )
    result = 0
    for value in minors:
        result = math.gcd(result, value)
    return result


def maximum_overlap_weight(slots: list[dict]) -> Fraction:
    """Use Kruskal independently of the producer's Prim implementation."""

    ordered = sorted(slots, key=lambda slot: slot["p"])
    candidates = []
    for left_index, left in enumerate(ordered):
        for right in ordered[left_index + 1 :]:
            if joint_image_index(left, right) == 1:
                candidates.append(
                    (
                        Fraction(1, left["h"] * right["h"]),
                        left["p"],
                        right["p"],
                    )
                )
    candidates.sort(key=lambda edge: (-edge[0], edge[1], edge[2]))
    parent = {slot["p"]: slot["p"] for slot in ordered}

    def root(prime: int) -> int:
        while parent[prime] != prime:
            parent[prime] = parent[parent[prime]]
            prime = parent[prime]
        return prime

    total = Fraction()
    for weight, left_prime, right_prime in candidates:
        left_root = root(left_prime)
        right_root = root(right_prime)
        if left_root == right_root:
            continue
        parent[right_root] = left_root
        total += weight
    return total


def verify_emitted_forest(
    records: list[dict],
    slots: list[dict],
    cell_count: int,
) -> Fraction:
    by_prime = {slot["p"]: slot for slot in slots}
    parent = {prime: prime for prime in by_prime}

    def root(prime: int) -> int:
        while parent[prime] != prime:
            parent[prime] = parent[parent[prime]]
            prime = parent[prime]
        return prime

    total = Fraction()
    for record in records:
        left_prime = int(record["left_prime"])
        right_prime = int(record["right_prime"])
        if (
            left_prime >= right_prime
            or left_prime not in by_prime
            or right_prime not in by_prime
        ):
            raise AssertionError("forced-overlap edge endpoint mismatch")
        left = by_prime[left_prime]
        right = by_prime[right_prime]
        if (
            int(record["left_modulus"]) != left["h"]
            or int(record["right_modulus"]) != right["h"]
            or joint_image_index(left, right) != 1
            or int(record["joint_map_cokernel_index"]) != 1
        ):
            raise AssertionError("forced-overlap edge is not surjective")
        weight = Fraction(1, left["h"] * right["h"])
        if (
            int(record["intersection_density_numerator"])
            != weight.numerator
            or int(record["intersection_density_denominator"])
            != weight.denominator
            or cell_count % weight.denominator
            or int(record["intersection_cell_count"])
            != cell_count * weight.numerator // weight.denominator
        ):
            raise AssertionError("forced-overlap edge weight mismatch")
        left_root = root(left_prime)
        right_root = root(right_prime)
        if left_root == right_root:
            raise AssertionError("forced-overlap records contain a cycle")
        parent[right_root] = left_root
        total += weight
    return total


def verify(artifact: dict) -> dict:
    if artifact.get("schema") != "reverse_group_arithmetic_line_inventory_v3":
        raise AssertionError("inventory schema mismatch")
    source_path = Path(artifact["source"])
    if sha256(source_path) != artifact["source_sha256"]:
        raise AssertionError("source SHA-256 mismatch")
    source = json.loads(source_path.read_text())
    rows = source["choices"]
    direction = tuple(
        int(value) for value in artifact["basis"]["direction"]
    )
    transverse = tuple(
        int(value) for value in artifact["basis"]["transverse"]
    )
    determinant = (
        direction[0] * transverse[1]
        - direction[1] * transverse[0]
    )
    if abs(determinant) != 1:
        raise AssertionError("basis is not unimodular")
    if determinant != int(artifact["basis"]["determinant"]):
        raise AssertionError("basis determinant mismatch")
    width = int(artifact["group"]["width"])
    height = int(artifact["group"]["height"])
    cell_count = width * height
    if (
        width < 1
        or height < 1
        or cell_count != int(artifact["group"]["cell_count"])
    ):
        raise AssertionError("group size mismatch")

    slots = []
    seen_primes = set()
    excluded = 0
    density = Fraction()
    for row in rows:
        prime = int(row["p"])
        if prime in seen_primes:
            raise AssertionError("duplicate source prime")
        seen_primes.add(prime)
        h, a, b = transformed_row(row, direction, transverse)
        if (a * width) % h or (b * height) % h:
            excluded += 1
            continue
        if cell_count % h:
            raise AssertionError("descending row has fractional fibre")
        slots.append({"p": prime, "h": h, "a": a, "b": b})
        density += Fraction(1, h)

    if (
        len(rows) != int(artifact["source_row_count"])
        or len(slots) != int(artifact["descending_row_count"])
        or excluded != int(artifact["excluded_not_descending_count"])
        or density.numerator != int(artifact["raw_density_numerator"])
        or density.denominator != int(artifact["raw_density_denominator"])
        or (density >= 1) is not artifact["density_at_least_one"]
    ):
        raise AssertionError("inventory density metadata mismatch")

    forest = artifact["forced_overlap_forest"]
    overlap = verify_emitted_forest(forest, slots, cell_count)
    maximum_overlap = maximum_overlap_weight(slots)
    if overlap != maximum_overlap:
        raise AssertionError("emitted Hunter forest is not maximum-weight")
    union_bound = density - overlap
    impossible = union_bound < 1
    expected_scalars = {
        "forced_overlap_density_numerator": overlap.numerator,
        "forced_overlap_density_denominator": overlap.denominator,
        "phase_independent_union_upper_bound_numerator": (
            union_bound.numerator
        ),
        "phase_independent_union_upper_bound_denominator": (
            union_bound.denominator
        ),
    }
    if any(
        int(artifact[key]) != value
        for key, value in expected_scalars.items()
    ):
        raise AssertionError("forced-overlap bound metadata mismatch")
    if (
        artifact["finite_group_cover_impossible_by_density_overlap"]
        is not impossible
        or artifact["claim"][
            "finite_group_cover_impossible_by_density_overlap"
        ]
        is not impossible
        or artifact["claim"]["integer_m_found"] is not False
    ):
        raise AssertionError("finite-family claim mismatch")
    return {
        "verified": True,
        "source_sha256": artifact["source_sha256"],
        "descending_row_count": len(slots),
        "raw_density_numerator": density.numerator,
        "raw_density_denominator": density.denominator,
        "forced_overlap_forest_edge_count": len(forest),
        "forced_overlap_density_numerator": overlap.numerator,
        "forced_overlap_density_denominator": overlap.denominator,
        "union_upper_bound_numerator": union_bound.numerator,
        "union_upper_bound_denominator": union_bound.denominator,
        "declared_finite_group_family_noncover_proved": impossible,
        "integer_m_found": False,
        "scope": (
            "exact Hunter union-bound verification for every phase "
            "assignment in the declared finite quotient family only"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inventory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    artifact = json.loads(args.inventory.read_text())
    report = verify(artifact)
    report["inventory"] = str(args.inventory)
    report["inventory_sha256"] = sha256(args.inventory)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"verified={report['verified']} "
        f"rows={report['descending_row_count']} "
        f"union_upper={report['union_upper_bound_numerator']}/"
        f"{report['union_upper_bound_denominator']} "
        f"noncover="
        f"{report['declared_finite_group_family_noncover_proved']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
