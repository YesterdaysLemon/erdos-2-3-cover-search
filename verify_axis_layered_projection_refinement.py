#!/usr/bin/env python3
"""Independently replay exact claims in a projected layered refinement.

The refinement artifact deliberately mixes proof-grade facts with
floating-point discovery evidence.  This verifier authenticates exact
capacity witnesses, exact anchor-dual obstructions, monotone uses of those
obstructions, and the algebra behind persisted cuts.  It reports discovery
MILP outcomes separately and never turns them into a proof.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path


def row_data(row: dict, axis: str) -> tuple[int, int, int]:
    h = int(row["h"])
    a = int(row["a"]) % h
    b = int(row["b"]) % h
    if h < 1 or math.gcd(a, b, h) != 1:
        raise AssertionError("invalid affine row")
    if axis == "k":
        active_modulus = math.gcd(b, h)
        coefficient = a
        residual = h // active_modulus
    elif axis == "l":
        active_modulus = math.gcd(a, h)
        coefficient = b
        residual = h // active_modulus
    else:
        raise AssertionError("invalid layer axis")
    if math.gcd(coefficient, active_modulus) != 1:
        raise AssertionError("active-coordinate map is not surjective")
    return active_modulus, coefficient, residual


def projected_row(
    rows: list[dict],
    axis: str,
    row_index: int,
    projection: int,
) -> dict:
    _active_modulus, _coefficient, residual = row_data(
        rows[row_index],
        axis,
    )
    modulus = math.gcd(residual, projection)
    return {
        "row_index": row_index,
        "p": int(rows[row_index]["p"]),
        "projected_modulus": modulus,
        "conditional_denominator": residual // modulus,
    }


def maximum_weighted_contribution(row: dict, weights: list[int]) -> Fraction:
    modulus = int(row["projected_modulus"])
    denominator = int(row["conditional_denominator"])
    maximum = max(
        sum(weights[cell] for cell in range(residue, len(weights), modulus))
        for residue in range(modulus)
    )
    return Fraction(maximum, denominator)


def replay_dual_evidence(
    rows: list[dict],
    axis: str,
    active_indices: frozenset[int],
    projection: int,
    evidence: dict,
) -> dict:
    if not evidence.get("oracle", "").startswith(
        "exact-anchor-projection-dual"
    ):
        raise AssertionError("cache evidence is not an exact anchor dual")
    if evidence.get("status") != "infeasible":
        raise AssertionError("exact dual does not declare infeasibility")
    cut_active = frozenset(
        int(index)
        for index in evidence["infeasible_cut_active_row_indices"]
    )
    if cut_active != active_indices:
        raise AssertionError("exact-dual cache key mismatch")
    obstruction = evidence["dual_obstruction"]
    if int(obstruction["projection_modulus"]) != projection:
        raise AssertionError("dual projection mismatch")
    base_indices = tuple(
        int(index) for index in obstruction["base_anchor_row_indices"]
    )
    branch_index = int(obstruction["branch_anchor_row_index"])
    anchor_indices = {*base_indices, branch_index}
    if len(base_indices) != 2 or len(anchor_indices) != 3:
        raise AssertionError("invalid anchor indices")
    if not anchor_indices <= active_indices:
        raise AssertionError("dual anchor is not active")
    first, second = (
        projected_row(rows, axis, index, projection)
        for index in base_indices
    )
    branch_anchor = projected_row(
        rows,
        axis,
        branch_index,
        projection,
    )
    first_modulus = int(first["projected_modulus"])
    second_modulus = int(second["projected_modulus"])
    if (
        int(first["conditional_denominator"]) != 1
        or int(second["conditional_denominator"]) != 1
        or int(branch_anchor["conditional_denominator"]) != 1
        or math.gcd(first_modulus, second_modulus) != 1
        or math.lcm(first_modulus, second_modulus) != projection
        or int(branch_anchor["projected_modulus"]) != projection
    ):
        raise AssertionError("invalid full-cell anchor normalization")
    if obstruction["base_anchor_primes"] != [first["p"], second["p"]]:
        raise AssertionError("base-anchor prime mismatch")
    if int(obstruction["branch_anchor_prime"]) != branch_anchor["p"]:
        raise AssertionError("branch-anchor prime mismatch")
    tail = [
        projected_row(rows, axis, index, projection)
        for index in sorted(active_indices - anchor_indices)
    ]
    if int(obstruction["tail_row_count"]) != len(tail):
        raise AssertionError("dual tail count mismatch")
    branches = obstruction["branches"]
    if [
        int(branch["branch_anchor_target"]) for branch in branches
    ] != list(range(projection)):
        raise AssertionError("dual branches are not exhaustive")
    minimum_gap = None
    for branch in branches:
        target = int(branch["branch_anchor_target"])
        covered = {
            cell
            for cell in range(projection)
            if cell % first_modulus == 0
            or cell % second_modulus == 0
            or cell == target
        }
        if branch["covered_cells"] != sorted(covered):
            raise AssertionError("dual covered-cell list mismatch")
        weights = [int(value) for value in branch["cell_weights"]]
        if (
            len(weights) != projection
            or any(value < 0 for value in weights)
            or any(weights[cell] for cell in covered)
        ):
            raise AssertionError("invalid exact dual weights")
        total = sum(weights)
        if total != int(branch["total_uncovered_cell_weight"]) or total <= 0:
            raise AssertionError("invalid exact dual total weight")
        maximum = sum(
            (
                maximum_weighted_contribution(row, weights)
                for row in tail
            ),
            Fraction(),
        )
        stored_maximum = Fraction(
            int(branch["maximum_tail_weight_numerator"]),
            int(branch["maximum_tail_weight_denominator"]),
        )
        if maximum != stored_maximum:
            raise AssertionError("exact dual tail maximum mismatch")
        gap = Fraction(total) - maximum
        if gap != Fraction(
            int(branch["strict_gap_numerator"]),
            int(branch["strict_gap_denominator"]),
        ):
            raise AssertionError("exact dual gap mismatch")
        if gap <= 0:
            raise AssertionError("exact dual gap is not strict")
        minimum_gap = gap if minimum_gap is None else min(minimum_gap, gap)
    assert minimum_gap is not None
    return {
        "active_row_count": len(active_indices),
        "tail_row_count": len(tail),
        "minimum_gap_numerator": minimum_gap.numerator,
        "minimum_gap_denominator": minimum_gap.denominator,
    }


def exact_projection_capacities(
    rows: list[dict],
    axis: str,
    active_indices: tuple[int, ...],
    projection: int,
    placement: list[int],
) -> list[Fraction]:
    capacities = []
    for cell in range(projection):
        capacity = Fraction()
        for row_index, phase in zip(
            active_indices,
            placement,
            strict=True,
        ):
            row = projected_row(rows, axis, row_index, projection)
            modulus = int(row["projected_modulus"])
            denominator = int(row["conditional_denominator"])
            if not 0 <= phase < modulus:
                raise AssertionError("projected phase is outside its modulus")
            if cell % modulus == phase:
                capacity += Fraction(1, denominator)
        capacities.append(capacity)
    return capacities


def replay(payload: dict) -> dict:
    if payload.get("schema") != "axis_layered_projected_refinement_v1":
        raise AssertionError("refinement schema mismatch")
    rows = payload["choices"]
    if len({int(row["p"]) for row in rows}) != len(rows):
        raise AssertionError("duplicate row prime")
    axis = payload["layer_axis"]
    pattern_period = int(payload["capacity_pattern_period"])
    placement = []
    for row in rows:
        active_modulus, coefficient, _residual = row_data(row, axis)
        active_class = int(row["layer_active_class"])
        if (
            not 0 <= active_class < active_modulus
            or int(row["target_modulus"]) != active_modulus
            or int(row["target_residue"]) % active_modulus
            != coefficient * active_class % active_modulus
        ):
            raise AssertionError("materialized active class mismatch")
        placement.append(active_class)

    refinement = payload["projection_refinement"]
    projection = int(refinement["projection_modulus"])
    exact_cache = {}
    cache_payload = refinement["monotonicity_cache"]
    for record in cache_payload.get("exact_infeasible", []):
        active = frozenset(
            int(index) for index in record["active_row_indices"]
        )
        if active in exact_cache:
            raise AssertionError("duplicate exact cache entry")
        evidence = record["evidence"]
        replay_dual_evidence(rows, axis, active, projection, evidence)
        exact_cache[active] = evidence

    seen_coordinates = set()
    status_counts = Counter()
    exact_status_counts = Counter()
    for record in refinement["projection_types"]:
        active = tuple(
            int(index) for index in record["active_row_indices"]
        )
        coordinates = [int(value) for value in record["coordinates"]]
        for coordinate in coordinates:
            if not 0 <= coordinate < pattern_period:
                raise AssertionError("projection coordinate out of range")
            if coordinate in seen_coordinates:
                raise AssertionError("projection coordinate repeated")
            reconstructed = tuple(
                index
                for index, row in enumerate(rows)
                if coordinate % row_data(row, axis)[0] == placement[index]
            )
            if reconstructed != active:
                raise AssertionError("stored active set mismatch")
            seen_coordinates.add(coordinate)
        status = str(record["status"])
        status_counts[status] += 1
        oracle = str(record.get("oracle", ""))
        if status == "feasible":
            capacities = exact_projection_capacities(
                rows,
                axis,
                active,
                projection,
                [int(value) for value in record["placement"]],
            )
            minimum = min(capacities)
            if minimum < 1 or minimum != Fraction(
                int(record["minimum_numerator"]),
                int(record["minimum_denominator"]),
            ):
                raise AssertionError("feasible projection replay mismatch")
            exact_status_counts["feasible"] += 1
        elif oracle.startswith("exact-anchor-projection-dual"):
            cut_active = frozenset(
                int(index)
                for index in record[
                    "infeasible_cut_active_row_indices"
                ]
            )
            replay_dual_evidence(
                rows,
                axis,
                cut_active,
                projection,
                record,
            )
            if not frozenset(active) <= cut_active:
                raise AssertionError("direct exact cut misses active rows")
            exact_status_counts["infeasible"] += 1
        elif oracle == "monotone-exact-infeasible-superset":
            cut_active = frozenset(
                int(index)
                for index in record[
                    "infeasible_cut_active_row_indices"
                ]
            )
            if (
                cut_active not in exact_cache
                or not frozenset(active) <= cut_active
            ):
                raise AssertionError("monotone exact evidence mismatch")
            exact_status_counts["infeasible"] += 1
        elif status not in {"infeasible", "unknown"}:
            raise AssertionError("unrecognized projection status")
    if seen_coordinates != set(range(pattern_period)):
        raise AssertionError("projection types do not partition the period")

    exact_cut_count = 0
    discovery_cut_count = 0
    for record in refinement.get("projection_benders_cuts", []):
        cut_active = frozenset(
            int(index)
            for index in record["infeasible_cut_active_row_indices"]
        )
        coordinate = int(record["coordinate"])
        if not 0 <= coordinate < pattern_period:
            raise AssertionError("projection cut coordinate out of range")
        if record.get("evidence_tier") == "exact":
            if not any(cut_active <= known for known in exact_cache):
                raise AssertionError("exact cut lacks a dual superset")
            exact_cut_count += 1
        else:
            discovery_cut_count += 1

    valid_pair_cut_count = 0
    for record in refinement.get("pair_benders_cuts", []):
        coordinate = int(record["coordinate"])
        left, right = (
            int(index) for index in record["row_indices"]
        )
        if (
            not 0 <= coordinate < pattern_period
            or not 0 <= left < len(rows)
            or not 0 <= right < len(rows)
            or left == right
        ):
            raise AssertionError("invalid persisted pair cut")
        left_residual = row_data(rows[left], axis)[2]
        right_residual = row_data(rows[right], axis)[2]
        if math.gcd(left_residual, right_residual) != 1:
            raise AssertionError("persisted pair cut is not unavoidable")
        valid_pair_cut_count += 1

    declared_projection_safe = refinement.get("projection_safe") is True
    projection_safe_verified = (
        declared_projection_safe
        and status_counts == Counter({"feasible": len(refinement["projection_types"])})
    )
    return {
        "verified_exact_claims": True,
        "projection_modulus": projection,
        "projection_type_count": len(refinement["projection_types"]),
        "projection_status_counts": dict(status_counts),
        "exact_projection_status_counts": dict(exact_status_counts),
        "exact_infeasible_cache_count": len(exact_cache),
        "exact_benders_cut_count": exact_cut_count,
        "discovery_benders_cut_count": discovery_cut_count,
        "valid_pair_cut_count": valid_pair_cut_count,
        "declared_projection_safe": declared_projection_safe,
        "projection_safe_verified": projection_safe_verified,
        "scope": (
            "exact claims inside one finite projected layered refinement; "
            "MILP discovery outcomes and the original infinite problem are "
            "not certified"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = replay(json.loads(args.input.read_text()))
    report["input"] = str(args.input)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"verified_exact_claims={report['verified_exact_claims']} "
        f"types={report['projection_type_count']} "
        f"exact_cache={report['exact_infeasible_cache_count']} "
        f"projection_safe_verified={report['projection_safe_verified']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
