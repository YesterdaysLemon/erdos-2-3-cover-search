#!/usr/bin/env python3
"""Refine a layered placement with exact small-projection capacity screens.

The master chooses one active layer class per row and requires ordinary
one-dimensional capacity in every layer column.  For each resulting active
row set, a subordinate MILP asks whether the residual phases can distribute
conditional density at least one in every cell modulo a small projection
``D``.  An exact weighted obstruction yields a rigorous monotone Benders
cut: because deleting active rows cannot repair a no-cover capacity
obstruction, at least one currently inactive row must be activated in that
column.  Floating-point MILP infeasibility is retained separately as
exploratory discovery evidence and must not be promoted to a certificate.

MILP answers are discovery evidence.  Every accepted master placement and
every accepted projected placement are replayed with exact rational
arithmetic before an artifact is written.  A projected-safe placement is
therefore a certified necessary-condition survivor, not yet a lattice cover.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from fractions import Fraction
from pathlib import Path

from build_axis_layered_pool import (
    exact_capacity_replay,
    layer_data,
    materialize_rows,
    sha256,
    unavoidable_pair_violations,
)


class ProjectionMonotoneCache:
    """Reuse exact-feasible and discovery-infeasible active-set answers.

    Projected capacity is monotone in the active rows.  A feasible set makes
    every superset feasible, while an infeasible set makes every subset
    infeasible.  Feasible witnesses are always replayed exactly before they
    enter this cache.  Infeasible entries remain solver-discovery evidence;
    they accelerate candidate search but are not publication certificates.
    """

    def __init__(self) -> None:
        self.feasible: dict[frozenset[int], dict[int, int]] = {}
        self.infeasible: dict[frozenset[int], dict] = {}
        self.exact_infeasible: dict[frozenset[int], dict] = {}

    def add_feasible(
        self,
        active_indices: tuple[int, ...],
        placement: list[int],
    ) -> None:
        active = frozenset(active_indices)
        witness = dict(zip(active_indices, placement, strict=True))
        if any(known <= active for known in self.feasible):
            return
        self.feasible = {
            known: known_witness
            for known, known_witness in self.feasible.items()
            if not active <= known
        }
        self.feasible[active] = witness

    def add_infeasible(
        self,
        active_indices: tuple[int, ...],
        evidence: dict | None = None,
    ) -> None:
        active = frozenset(active_indices)
        evidence = evidence or {}
        if evidence.get("oracle", "").startswith("exact-"):
            self.exact_infeasible = self._insert_infeasible_frontier(
                self.exact_infeasible,
                active,
                evidence,
            )
        self.infeasible = self._insert_infeasible_frontier(
            self.infeasible,
            active,
            evidence,
        )

    @staticmethod
    def _insert_infeasible_frontier(
        frontier: dict[frozenset[int], dict],
        active: frozenset[int],
        evidence: dict,
    ) -> dict[frozenset[int], dict]:
        if active in frontier:
            if (
                evidence.get("oracle", "").startswith("exact-")
                and not frontier[active]
                .get("oracle", "")
                .startswith("exact-")
            ):
                frontier[active] = evidence
            return frontier
        if any(active <= known for known in frontier):
            return frontier
        return {
            **{
                known: known_evidence
                for known, known_evidence in frontier.items()
                if not known <= active
            },
            active: evidence,
        }

    def classify(
        self,
        active_indices: tuple[int, ...],
        rows: list[dict],
        axis: str,
        projection: int,
        exact_infeasible_only: bool = False,
    ) -> dict | None:
        active = frozenset(active_indices)
        feasible_subsets = [
            known for known in self.feasible if known <= active
        ]
        if feasible_subsets:
            known = min(feasible_subsets, key=len)
            witness = self.feasible[known]
            placement = [witness.get(index, 0) for index in active_indices]
            capacities = exact_projection_capacities(
                rows,
                axis,
                active_indices,
                projection,
                placement,
            )
            minimum = min(capacities)
            if minimum < 1:
                raise AssertionError(
                    "monotone feasible witness failed exact replay"
                )
            return {
                "status": "feasible",
                "oracle": "monotone-feasible-subset",
                "placement": placement,
                "minimum_numerator": minimum.numerator,
                "minimum_denominator": minimum.denominator,
                "minimum_decimal": float(minimum),
                "dominating_active_row_indices": sorted(known),
            }
        infeasible_frontier = (
            self.exact_infeasible
            if exact_infeasible_only
            else self.infeasible
        )
        infeasible_supersets = [
            known for known in infeasible_frontier if active <= known
        ]
        if infeasible_supersets:
            known = max(infeasible_supersets, key=len)
            return {
                "status": "infeasible",
                "oracle": (
                    "monotone-exact-infeasible-superset"
                    if exact_infeasible_only
                    else "monotone-infeasible-superset"
                ),
                "infeasible_cut_active_row_indices": sorted(known),
                "dominating_evidence_oracle": infeasible_frontier[
                    known
                ].get("oracle", "unspecified-discovery"),
            }
        return None

    def add_result(
        self,
        active_indices: tuple[int, ...],
        result: dict,
    ) -> None:
        if result["status"] == "feasible":
            self.add_feasible(active_indices, result["placement"])
        elif result["status"] == "infeasible":
            cut_active = tuple(
                result.get(
                    "infeasible_cut_active_row_indices",
                    active_indices,
                )
            )
            self.add_infeasible(cut_active, result)

    def export(self) -> dict:
        return {
            "feasible": [
                {
                    "active_row_indices": sorted(active),
                    "placement_by_row": {
                        str(index): witness[index]
                        for index in sorted(witness)
                    },
                }
                for active, witness in sorted(
                    self.feasible.items(),
                    key=lambda item: (len(item[0]), sorted(item[0])),
                )
            ],
            "infeasible": [
                {
                    "active_row_indices": sorted(active),
                    "evidence": evidence,
                }
                for active, evidence in sorted(
                    self.infeasible.items(),
                    key=lambda item: (-len(item[0]), sorted(item[0])),
                )
            ],
            "exact_infeasible": [
                {
                    "active_row_indices": sorted(active),
                    "evidence": evidence,
                }
                for active, evidence in sorted(
                    self.exact_infeasible.items(),
                    key=lambda item: (-len(item[0]), sorted(item[0])),
                )
            ],
        }

    def import_payload(self, payload: dict | None) -> None:
        if not payload:
            return
        for record in payload.get("feasible", []):
            active = tuple(int(value) for value in record["active_row_indices"])
            by_row = {
                int(index): int(value)
                for index, value in record["placement_by_row"].items()
            }
            self.add_feasible(active, [by_row[index] for index in active])
        for record in payload.get("infeasible", []):
            self.add_infeasible(
                tuple(int(value) for value in record["active_row_indices"]),
                record.get("evidence"),
            )
        for record in payload.get("exact_infeasible", []):
            self.add_infeasible(
                tuple(int(value) for value in record["active_row_indices"]),
                record.get("evidence"),
            )

    def strengthen_exact_infeasible(
        self,
        rows: list[dict],
        axis: str,
        projection: int,
    ) -> dict:
        before = len(self.exact_infeasible)
        expanded_count = 0
        total_added_rows = 0
        for active, evidence in list(self.exact_infeasible.items()):
            if not (
                evidence.get("oracle", "").startswith(
                    "exact-anchor-projection-dual"
                )
                and "dual_obstruction" in evidence
            ):
                continue
            expanded = expand_anchor_projection_dual(
                rows,
                axis,
                tuple(sorted(active)),
                projection,
                evidence,
            )
            expanded_active = frozenset(
                int(index)
                for index in expanded[
                    "infeasible_cut_active_row_indices"
                ]
            )
            if expanded_active == active:
                continue
            expanded_count += 1
            total_added_rows += len(expanded_active - active)
            self.add_infeasible(tuple(sorted(expanded_active)), expanded)
        return {
            "frontier_size_before": before,
            "frontier_size_after": len(self.exact_infeasible),
            "expanded_certificate_count": expanded_count,
            "total_added_rows": total_added_rows,
        }


def active_set_groups(
    rows: list[dict],
    axis: str,
    period: int,
    placement: list[int],
) -> dict[tuple[int, ...], list[int]]:
    data = [layer_data(row, axis) for row in rows]
    groups: dict[tuple[int, ...], list[int]] = {}
    for coordinate in range(period):
        active = tuple(
            index
            for index, (_order, modulus, _coefficient, _residual) in enumerate(
                data
            )
            if coordinate % modulus == placement[index]
        )
        groups.setdefault(active, []).append(coordinate)
    return groups


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
        for row_index, phase_class in zip(
            active_indices,
            placement,
            strict=True,
        ):
            residual = layer_data(rows[row_index], axis)[3]
            modulus = math.gcd(residual, projection)
            denominator = residual // modulus
            if cell % modulus == phase_class:
                capacity += Fraction(1, denominator)
        capacities.append(capacity)
    return capacities


def solve_anchor_projection_dual(
    rows: list[dict],
    axis: str,
    active_indices: tuple[int, ...],
    projection: int,
) -> dict | None:
    """Try an exact weighted obstruction after normalizing full-cell anchors."""

    from certify_axis_layered_projection_noncover import (
        discover_branch_weights,
        maximum_weighted_tail,
    )

    projected_rows = []
    for row_index in active_indices:
        residual = layer_data(rows[row_index], axis)[3]
        modulus = math.gcd(residual, projection)
        projected_rows.append(
            {
                "row_index": row_index,
                "p": int(rows[row_index]["p"]),
                "projected_modulus": modulus,
                "conditional_denominator": residual // modulus,
            }
        )
    anchors = [
        row
        for row in projected_rows
        if row["conditional_denominator"] == 1
    ]
    branch_anchors = [
        row for row in anchors if row["projected_modulus"] == projection
    ]
    base_pairs = [
        (first, second)
        for first in anchors
        for second in anchors
        if first["row_index"] < second["row_index"]
        and math.gcd(
            int(first["projected_modulus"]),
            int(second["projected_modulus"]),
        )
        == 1
        and math.lcm(
            int(first["projected_modulus"]),
            int(second["projected_modulus"]),
        )
        == projection
    ]
    for first, second in base_pairs:
        for branch_anchor in branch_anchors:
            anchor_indices = {
                first["row_index"],
                second["row_index"],
                branch_anchor["row_index"],
            }
            if len(anchor_indices) != 3:
                continue
            first_modulus = int(first["projected_modulus"])
            second_modulus = int(second["projected_modulus"])
            tail = [
                row
                for row in projected_rows
                if row["row_index"] not in anchor_indices
            ]
            branches = []
            try:
                for target in range(projection):
                    covered = {
                        cell
                        for cell in range(projection)
                        if cell % first_modulus == 0
                        or cell % second_modulus == 0
                        or cell == target
                    }
                    weights = discover_branch_weights(
                        tail,
                        covered,
                        projection,
                    )
                    total_weight = sum(weights)
                    maximum = maximum_weighted_tail(tail, weights)
                    gap = Fraction(total_weight) - maximum
                    if gap <= 0:
                        raise RuntimeError("non-strict exact dual")
                    branches.append(
                        {
                            "branch_anchor_target": target,
                            "covered_cells": sorted(covered),
                            "cell_weights": weights,
                            "total_uncovered_cell_weight": total_weight,
                            "maximum_tail_weight_numerator": (
                                maximum.numerator
                            ),
                            "maximum_tail_weight_denominator": (
                                maximum.denominator
                            ),
                            "strict_gap_numerator": gap.numerator,
                            "strict_gap_denominator": gap.denominator,
                        }
                    )
            except RuntimeError:
                continue
            result = {
                "status": "infeasible",
                "oracle": "exact-anchor-projection-dual",
                "infeasible_cut_active_row_indices": list(active_indices),
                "dual_obstruction": {
                    "projection_modulus": projection,
                    "base_anchor_row_indices": [
                        first["row_index"],
                        second["row_index"],
                    ],
                    "base_anchor_primes": [
                        first["p"],
                        second["p"],
                    ],
                    "base_anchor_moduli": [
                        first_modulus,
                        second_modulus,
                    ],
                    "branch_anchor_row_index": (
                        branch_anchor["row_index"]
                    ),
                    "branch_anchor_prime": branch_anchor["p"],
                    "branch_anchor_modulus": projection,
                    "tail_row_count": len(tail),
                    "branches": branches,
                },
            }
            return expand_anchor_projection_dual(
                rows,
                axis,
                active_indices,
                projection,
                result,
            )
    return None


def expand_anchor_projection_dual(
    rows: list[dict],
    axis: str,
    active_indices: tuple[int, ...],
    projection: int,
    result: dict,
) -> dict:
    """Greedily enlarge an exact dual obstruction without another LP."""

    from certify_axis_layered_projection_noncover import (
        maximum_weighted_tail,
    )

    obstruction = result["dual_obstruction"]
    branches = [dict(branch) for branch in obstruction["branches"]]
    maxima = [
        Fraction(
            branch["maximum_tail_weight_numerator"],
            branch["maximum_tail_weight_denominator"],
        )
        for branch in branches
    ]
    totals = [
        Fraction(branch["total_uncovered_cell_weight"])
        for branch in branches
    ]
    active = set(active_indices)
    candidate_contributions = {}
    for row_index, raw in enumerate(rows):
        if row_index in active:
            continue
        residual = layer_data(raw, axis)[3]
        modulus = math.gcd(residual, projection)
        projected_row = {
            "projected_modulus": modulus,
            "conditional_denominator": residual // modulus,
        }
        candidate_contributions[row_index] = [
            maximum_weighted_tail(
                [projected_row],
                branch["cell_weights"],
            )
            for branch in branches
        ]

    added = []
    while True:
        fitting = []
        for row_index, contributions in candidate_contributions.items():
            gaps = [
                total - maximum
                for total, maximum in zip(totals, maxima, strict=True)
            ]
            if not all(
                contribution < gap
                for contribution, gap in zip(
                    contributions,
                    gaps,
                    strict=True,
                )
            ):
                continue
            score = max(
                (
                    contribution / gap
                    for contribution, gap in zip(
                        contributions,
                        gaps,
                        strict=True,
                    )
                ),
                default=Fraction(),
            )
            fitting.append((score, row_index, contributions))
        if not fitting:
            break
        _score, row_index, contributions = min(fitting)
        added.append(row_index)
        active.add(row_index)
        maxima = [
            maximum + contribution
            for maximum, contribution in zip(
                maxima,
                contributions,
                strict=True,
            )
        ]
        del candidate_contributions[row_index]

    if not added:
        return result
    for branch, maximum, total in zip(
        branches,
        maxima,
        totals,
        strict=True,
    ):
        gap = total - maximum
        if gap <= 0:
            raise AssertionError("expanded exact dual lost its strict gap")
        branch["maximum_tail_weight_numerator"] = maximum.numerator
        branch["maximum_tail_weight_denominator"] = maximum.denominator
        branch["strict_gap_numerator"] = gap.numerator
        branch["strict_gap_denominator"] = gap.denominator
    expanded = dict(result)
    expanded["oracle"] = "exact-anchor-projection-dual-expanded"
    expanded["infeasible_cut_active_row_indices"] = sorted(active)
    expanded_obstruction = dict(obstruction)
    expanded_obstruction["expanded_from_active_row_indices"] = list(
        active_indices
    )
    expanded_obstruction["greedy_added_row_indices"] = added
    expanded_obstruction["tail_row_count"] = (
        int(obstruction["tail_row_count"]) + len(added)
    )
    expanded_obstruction["branches"] = branches
    expanded["dual_obstruction"] = expanded_obstruction
    return expanded


def solve_projection_subproblem(
    rows: list[dict],
    axis: str,
    active_indices: tuple[int, ...],
    projection: int,
    time_limit: float,
    try_anchor_dual: bool = True,
) -> dict:
    if try_anchor_dual:
        dual = solve_anchor_projection_dual(
            rows,
            axis,
            active_indices,
            projection,
        )
        if dual is not None:
            return dual

    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import lil_matrix

    projection_rows = []
    offsets = []
    choice_variable_count = 0
    for row_index in active_indices:
        residual = layer_data(rows[row_index], axis)[3]
        modulus = math.gcd(residual, projection)
        projection_rows.append((modulus, residual // modulus))
        offsets.append(choice_variable_count)
        choice_variable_count += modulus
    minimum_variable = choice_variable_count
    variable_count = choice_variable_count + 1
    constraint_count = len(active_indices) + projection
    matrix = lil_matrix((constraint_count, variable_count))
    lower = np.empty(constraint_count)
    upper = np.empty(constraint_count)
    for index, (modulus, _denominator) in enumerate(projection_rows):
        matrix[index, offsets[index] : offsets[index] + modulus] = 1.0
        lower[index] = 1.0
        upper[index] = 1.0
    for cell in range(projection):
        constraint = len(active_indices) + cell
        for index, (modulus, denominator) in enumerate(projection_rows):
            matrix[constraint, offsets[index] + cell % modulus] = (
                1.0 / denominator
            )
        matrix[constraint, minimum_variable] = -1.0
        lower[constraint] = 0.0
        upper[constraint] = np.inf
    objective = np.zeros(variable_count)
    objective[minimum_variable] = -1.0
    integrality = np.ones(variable_count, dtype=np.int8)
    integrality[minimum_variable] = 0
    bound_lower = np.zeros(variable_count)
    bound_upper = np.ones(variable_count)
    bound_upper[minimum_variable] = sum(
        1.0 / denominator
        for _modulus, denominator in projection_rows
    )
    result = milp(
        objective,
        integrality=integrality,
        bounds=Bounds(bound_lower, bound_upper),
        constraints=LinearConstraint(matrix.tocsr(), lower, upper),
        options={
            "presolve": True,
            "time_limit": time_limit,
            "mip_rel_gap": 0.01,
        },
    )
    if result.x is None:
        return {
            "status": "unknown",
            "oracle": "scipy-highs-maximin-discovery",
            "solver_status": int(result.status),
            "solver_message": str(result.message),
        }
    placement = [
        int(
            np.argmax(
                result.x[offsets[index] : offsets[index] + modulus]
            )
        )
        for index, (modulus, _denominator) in enumerate(projection_rows)
    ]
    capacities = exact_projection_capacities(
        rows,
        axis,
        active_indices,
        projection,
        placement,
    )
    minimum = min(capacities)
    if minimum < 1:
        discovered_minimum = float(result.x[minimum_variable])
        if int(result.status) == 0 and discovered_minimum < 1.0 - 1e-7:
            return {
                "status": "infeasible",
                "oracle": "scipy-highs-maximin-discovery",
                "solver_status": int(result.status),
                "solver_message": str(result.message),
                "floating_optimum": discovered_minimum,
                "decoded_minimum_numerator": minimum.numerator,
                "decoded_minimum_denominator": minimum.denominator,
                "decoded_minimum_decimal": float(minimum),
            }
        return {
            "status": "unknown",
            "solver_status": int(result.status),
            "oracle": "scipy-highs-maximin-discovery",
            "solver_message": (
                "maximin incumbent remains below one or failed exact replay"
            ),
            "floating_incumbent": discovered_minimum,
            "decoded_minimum_numerator": minimum.numerator,
            "decoded_minimum_denominator": minimum.denominator,
            "decoded_minimum_decimal": float(minimum),
        }
    return {
        "status": "feasible",
        "oracle": "scipy-highs-maximin-exact-witness-replay",
        "solver_status": int(result.status),
        "solver_message": str(result.message),
        "placement": placement,
        "minimum_numerator": minimum.numerator,
        "minimum_denominator": minimum.denominator,
        "minimum_decimal": float(minimum),
    }


def refinement_root(
    input_path: Path,
    payload: dict,
) -> tuple[str, str]:
    """Resolve the first non-refinement ancestor with hash checks."""

    declared_root = payload.get("refinement_root_source")
    declared_hash = payload.get("refinement_root_source_sha256")
    if declared_root is not None or declared_hash is not None:
        if not isinstance(declared_root, str) or not isinstance(
            declared_hash,
            str,
        ):
            raise ValueError("incomplete declared refinement root")
        root_path = Path(declared_root)
        if sha256(root_path) != declared_hash:
            raise ValueError("declared refinement root SHA-256 mismatch")
        return declared_root, declared_hash

    current_path = input_path
    current = payload
    seen = set()
    while current.get("schema") == "axis_layered_projected_refinement_v1":
        key = str(current_path.resolve())
        if key in seen:
            raise ValueError("refinement ancestry contains a cycle")
        seen.add(key)
        parent_path = Path(current["refinement_source"])
        parent_hash = str(current["refinement_source_sha256"])
        if sha256(parent_path) != parent_hash:
            raise ValueError("refinement parent SHA-256 mismatch")
        current_path = parent_path
        current = json.loads(parent_path.read_text())
    return str(current_path), sha256(current_path)


def discover_refined_placement(
    rows: list[dict],
    axis: str,
    period: int,
    projection: int,
    time_limit: float,
    subproblem_time_limit: float,
    cut_round_limit: int,
    seed_cache_payload: dict | None = None,
    seed_projection_cuts: list[dict] | None = None,
    seed_pair_cuts: list[dict] | None = None,
) -> tuple[list[int], dict]:
    dependency_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dependency_path))
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import lil_matrix, vstack

    data = [layer_data(row, axis) for row in rows]
    offsets = []
    variable_count = 0
    for _order, modulus, _coefficient, _residual in data:
        offsets.append(variable_count)
        variable_count += modulus
    objective = np.zeros(variable_count)
    for index, (_order, modulus, _coefficient, _residual) in enumerate(data):
        active_class = int(rows[index].get("layer_active_class", -1))
        if 0 <= active_class < modulus:
            objective[offsets[index] + active_class] = -1.0
    constraint_count = len(rows) + period
    matrix = lil_matrix((constraint_count, variable_count))
    lower = np.empty(constraint_count)
    upper = np.empty(constraint_count)
    for index, (_order, modulus, _coefficient, _residual) in enumerate(data):
        matrix[index, offsets[index] : offsets[index] + modulus] = 1.0
        lower[index] = 1.0
        upper[index] = 1.0
    for coordinate in range(period):
        constraint = len(rows) + coordinate
        for index, (_order, modulus, _coefficient, residual) in enumerate(
            data
        ):
            matrix[constraint, offsets[index] + coordinate % modulus] = (
                1.0 / residual
            )
        lower[constraint] = 1.0
        upper[constraint] = np.inf
    constraint_matrix = matrix.tocsr()
    started = time.monotonic()
    added_cut_keys: set[tuple[int, ...]] = set()
    projection_cut_records: dict[tuple[int, ...], dict] = {}
    pair_cut_records: dict[tuple[int, ...], dict] = {}
    records = []
    placement = None
    final_projection_records = []
    final_pair_violations = []
    monotone_cache = ProjectionMonotoneCache()
    monotone_cache.import_payload(seed_cache_payload)
    exact_cache_strengthening = monotone_cache.strengthen_exact_infeasible(
        rows,
        axis,
        projection,
    )
    initial_feasible_cache_size = len(monotone_cache.feasible)
    initial_infeasible_cache_size = len(monotone_cache.infeasible)
    initial_exact_infeasible_cache_size = len(
        monotone_cache.exact_infeasible
    )

    seed_rows = []
    seed_lowers = []
    strengthened_seed_projection_cut_count = 0
    for record in seed_projection_cuts or []:
        coordinate = int(record["coordinate"]) % period
        cut_active = {
            int(index)
            for index in record["infeasible_cut_active_row_indices"]
        }
        if any(index < 0 or index >= len(rows) for index in cut_active):
            raise ValueError("seed projection cut has an invalid row index")
        exact_supersets = [
            known
            for known in monotone_cache.exact_infeasible
            if cut_active <= known
        ]
        evidence_tier = str(
            record.get("evidence_tier", "discovery")
        )
        if exact_supersets:
            strongest = max(exact_supersets, key=len)
            if strongest != cut_active:
                strengthened_seed_projection_cut_count += 1
            cut_active = set(strongest)
            evidence_tier = "exact"
        variables = tuple(
            offsets[index] + coordinate % modulus
            for index, (
                _order,
                modulus,
                _coefficient,
                _residual,
            ) in enumerate(data)
            if index not in cut_active
        )
        if not variables:
            raise RuntimeError(
                "seed projection-infeasible superset contains every row"
            )
        if variables in added_cut_keys:
            continue
        added_cut_keys.add(variables)
        projection_cut_records[variables] = {
            "coordinate": coordinate,
            "infeasible_cut_active_row_indices": sorted(cut_active),
            "evidence_tier": evidence_tier,
        }
        seed_rows.append(variables)
        seed_lowers.append(1.0)
    for record in seed_pair_cuts or []:
        coordinate = int(record["coordinate"]) % period
        left, right = (
            int(index) for index in record["row_indices"]
        )
        if (
            left < 0
            or right < 0
            or left >= len(rows)
            or right >= len(rows)
            or left == right
        ):
            raise ValueError("seed pair cut has invalid row indices")
        left_period = data[left][3]
        right_period = data[right][3]
        if math.gcd(left_period, right_period) != 1:
            raise ValueError("seed pair cut periods are not coprime")
        overlap = Fraction(1, left_period * right_period)
        variables = tuple(
            offsets[index] + coordinate % modulus
            for index, (
                _order,
                modulus,
                _coefficient,
                _residual,
            ) in enumerate(data)
        )
        key = (-1, coordinate, min(left, right), max(left, right))
        if key in added_cut_keys:
            continue
        added_cut_keys.add(key)
        coefficients = {
            variable: 1.0 / data[index][3]
            for index, variable in enumerate(variables)
        }
        coefficients[variables[left]] -= float(overlap)
        coefficients[variables[right]] -= float(overlap)
        pair_cut_records[key] = {
            "coordinate": coordinate,
            "row_indices": [left, right],
        }
        seed_rows.append(coefficients)
        seed_lowers.append(1.0 - float(overlap))
    if seed_rows:
        additions = lil_matrix((len(seed_rows), variable_count))
        for row_index, specification in enumerate(seed_rows):
            if isinstance(specification, dict):
                for variable, coefficient in specification.items():
                    additions[row_index, variable] = coefficient
            else:
                for variable in specification:
                    additions[row_index, variable] = 1.0
        constraint_matrix = vstack(
            (constraint_matrix, additions.tocsr()),
            format="csr",
        )
        lower = np.concatenate((lower, np.asarray(seed_lowers)))
        upper = np.concatenate(
            (upper, np.full(len(seed_rows), np.inf))
        )
    initial_projection_cut_count = len(projection_cut_records)
    initial_pair_cut_count = len(pair_cut_records)

    for cut_round in range(cut_round_limit + 1):
        remaining_time = max(0.001, time_limit - (time.monotonic() - started))
        result = milp(
            objective,
            integrality=np.ones(variable_count, dtype=np.int8),
            bounds=Bounds(
                np.zeros(variable_count),
                np.ones(variable_count),
            ),
            constraints=LinearConstraint(
                constraint_matrix,
                lower,
                upper,
            ),
            options={
                "presolve": True,
                "time_limit": remaining_time,
                "mip_rel_gap": 0.05,
            },
        )
        if result.x is None:
            raise RuntimeError(
                "refined placement master failed: "
                f"status={result.status} message={result.message}"
            )
        placement = [
            int(
                np.argmax(
                    result.x[offsets[index] : offsets[index] + modulus]
                )
            )
            for index, (_order, modulus, _coefficient, _residual) in enumerate(
                data
            )
        ]
        exact_capacity_replay(rows, axis, period, placement)
        pair_violations = unavoidable_pair_violations(
            rows,
            axis,
            period,
            placement,
        )
        groups = active_set_groups(rows, axis, period, placement)
        projection_cache = {}
        monotone_hit_count = 0
        exact_dual_count = 0
        direct_oracle_count = 0
        for active in groups:
            answer = monotone_cache.classify(
                active,
                rows,
                axis,
                projection,
                exact_infeasible_only=True,
            )
            if answer is None:
                answer = solve_anchor_projection_dual(
                    rows,
                    axis,
                    active,
                    projection,
                )
                if answer is not None:
                    exact_dual_count += 1
                    direct_oracle_count += 1
                    monotone_cache.add_result(active, answer)
            if answer is None:
                answer = monotone_cache.classify(
                    active,
                    rows,
                    axis,
                    projection,
                )
            if answer is None:
                answer = solve_projection_subproblem(
                    rows,
                    axis,
                    active,
                    projection,
                    subproblem_time_limit,
                    try_anchor_dual=False,
                )
                direct_oracle_count += 1
                monotone_cache.add_result(active, answer)
            else:
                if answer["oracle"].startswith("monotone-"):
                    monotone_hit_count += 1
            projection_cache[active] = answer
        projection_records = [
            {
                "active_row_indices": list(active),
                "coordinates": coordinates,
                **projection_cache[active],
            }
            for active, coordinates in groups.items()
        ]
        infeasible_coordinates = [
            coordinate
            for active, coordinates in groups.items()
            if projection_cache[active]["status"] == "infeasible"
            for coordinate in coordinates
        ]
        unknown_type_count = sum(
            record["status"] == "unknown" for record in projection_records
        )
        records.append(
            {
                "cut_round": cut_round,
                "master_status": int(result.status),
                "master_message": str(result.message),
                "layer_active_set_type_count": len(groups),
                "projection_infeasible_coordinate_count": len(
                    infeasible_coordinates
                ),
                "projection_unknown_type_count": unknown_type_count,
                "monotone_cache_hit_count": monotone_hit_count,
                "exact_dual_count": exact_dual_count,
                "direct_oracle_count": direct_oracle_count,
                "feasible_cache_size": len(monotone_cache.feasible),
                "infeasible_cache_size": len(monotone_cache.infeasible),
                "exact_infeasible_cache_size": len(
                    monotone_cache.exact_infeasible
                ),
                "pair_violation_count": len(pair_violations),
                "seconds": time.monotonic() - started,
            }
        )
        print(
            f"round={cut_round} types={len(groups)} "
            f"projection_bad={len(infeasible_coordinates)} "
            f"projection_unknown={unknown_type_count} "
            f"cache_hits={monotone_hit_count} "
            f"pair_bad={len(pair_violations)} "
            f"elapsed_s={time.monotonic() - started:.3f}",
            flush=True,
        )
        final_projection_records = projection_records
        final_pair_violations = pair_violations
        if (
            not infeasible_coordinates
            and not pair_violations
            and unknown_type_count == 0
        ):
            break
        if (
            cut_round == cut_round_limit
            or time.monotonic() - started >= time_limit
        ):
            break

        cut_rows = []
        cut_lowers = []
        for active, coordinates in groups.items():
            answer = projection_cache[active]
            if answer["status"] != "infeasible":
                continue
            cut_active = set(
                answer.get(
                    "infeasible_cut_active_row_indices",
                    active,
                )
            )
            for coordinate in coordinates:
                variables = tuple(
                    offsets[index] + coordinate % modulus
                    for index, (
                        _order,
                        modulus,
                        _coefficient,
                        _residual,
                    ) in enumerate(data)
                    if index not in cut_active
                )
                if not variables:
                    raise RuntimeError(
                        "projection-infeasible superset contains every row"
                    )
                if variables in added_cut_keys:
                    continue
                added_cut_keys.add(variables)
                projection_cut_records[variables] = {
                    "coordinate": coordinate,
                    "infeasible_cut_active_row_indices": sorted(
                        cut_active
                    ),
                    "evidence_tier": (
                        "exact"
                        if (
                            answer["oracle"].startswith("exact-")
                            or answer["oracle"].startswith(
                                "monotone-exact-"
                            )
                        )
                        else "discovery"
                    ),
                }
                cut_rows.append(variables)
                cut_lowers.append(1.0)
        for violation in pair_violations:
            coordinate = int(violation["coordinate"])
            left, right = map(int, violation["row_indices"])
            overlap = Fraction(
                int(violation["overlap_numerator"]),
                int(violation["overlap_denominator"]),
            )
            variables = tuple(
                offsets[index] + coordinate % modulus
                for index, (_order, modulus, _coefficient, _residual) in enumerate(
                    data
                )
            )
            key = (-1, coordinate, left, right)
            if key in added_cut_keys:
                continue
            added_cut_keys.add(key)
            coefficients = {
                variable: 1.0 / data[index][3]
                for index, variable in enumerate(variables)
            }
            coefficients[variables[left]] -= float(overlap)
            coefficients[variables[right]] -= float(overlap)
            pair_cut_records[key] = {
                "coordinate": coordinate,
                "row_indices": [left, right],
            }
            cut_rows.append(coefficients)
            cut_lowers.append(1.0 - float(overlap))
        if not cut_rows:
            break
        additions = lil_matrix((len(cut_rows), variable_count))
        for row_index, specification in enumerate(cut_rows):
            if isinstance(specification, dict):
                for variable, coefficient in specification.items():
                    additions[row_index, variable] = coefficient
            else:
                for variable in specification:
                    additions[row_index, variable] = 1.0
        constraint_matrix = vstack(
            (constraint_matrix, additions.tocsr()),
            format="csr",
        )
        lower = np.concatenate((lower, np.asarray(cut_lowers)))
        upper = np.concatenate(
            (upper, np.full(len(cut_rows), np.inf))
        )
    assert placement is not None
    projection_safe = all(
        record["status"] == "feasible"
        for record in final_projection_records
    )
    return placement, {
        "engine": "scipy-highs-layer-capacity-with-projection-benders-cuts",
        "projection_modulus": projection,
        "base_constraints": constraint_count,
        "variables": variable_count,
        "cut_round_limit": cut_round_limit,
        "cuts_added": len(added_cut_keys),
        "initial_projection_cut_count": initial_projection_cut_count,
        "strengthened_seed_projection_cut_count": (
            strengthened_seed_projection_cut_count
        ),
        "new_projection_cut_count": (
            len(projection_cut_records) - initial_projection_cut_count
        ),
        "projection_benders_cuts": list(
            projection_cut_records.values()
        ),
        "initial_pair_cut_count": initial_pair_cut_count,
        "new_pair_cut_count": len(pair_cut_records) - initial_pair_cut_count,
        "pair_benders_cuts": list(pair_cut_records.values()),
        "initial_feasible_cache_size": initial_feasible_cache_size,
        "initial_infeasible_cache_size": initial_infeasible_cache_size,
        "initial_exact_infeasible_cache_size": (
            initial_exact_infeasible_cache_size
        ),
        "exact_cache_strengthening": exact_cache_strengthening,
        "monotonicity_cache": monotone_cache.export(),
        "rounds": records,
        "projection_safe": projection_safe,
        "pair_safe": not final_pair_violations,
        "projection_types": final_projection_records,
        "seconds": time.monotonic() - started,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--projection", type=int, required=True)
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--subproblem-time-limit", type=float, default=1.0)
    parser.add_argument("--cut-rounds", type=int, default=20)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (
        args.projection < 1
        or args.time_limit <= 0
        or args.subproblem_time_limit <= 0
        or args.cut_rounds < 0
    ):
        raise SystemExit("invalid positive limit or projection")
    source = json.loads(args.input.read_text())
    if source.get("schema") not in {
        "axis_layered_pool_v1",
        "axis_layered_projected_refinement_v1",
    }:
        raise SystemExit("input is not an axis-layered pool artifact")
    rows = [dict(row) for row in source["choices"]]
    axis = source["layer_axis"]
    pattern_period = int(source["capacity_pattern_period"])
    source_refinement = source.get("projection_refinement", {})
    seed_cache_payload = None
    seed_projection_cuts = None
    seed_pair_cuts = None
    if (
        int(source_refinement.get("projection_modulus", -1))
        == args.projection
    ):
        seed_cache_payload = source_refinement.get("monotonicity_cache")
        if seed_cache_payload is None:
            seed_cache_payload = {"feasible": [], "infeasible": []}
            for record in source_refinement.get("projection_types", []):
                active = [
                    int(value)
                    for value in record["active_row_indices"]
                ]
                if record.get("status") == "feasible":
                    seed_cache_payload["feasible"].append(
                        {
                            "active_row_indices": active,
                            "placement_by_row": {
                                str(index): int(value)
                                for index, value in zip(
                                    active,
                                    record["placement"],
                                    strict=True,
                                )
                            },
                        }
                    )
                elif record.get("status") == "infeasible":
                    seed_cache_payload["infeasible"].append(
                        {
                            "active_row_indices": record.get(
                                "infeasible_cut_active_row_indices",
                                active,
                            )
                        }
                    )
        seed_projection_cuts = source_refinement.get(
            "projection_benders_cuts"
        )
        if seed_projection_cuts is None:
            seed_projection_cuts = []
            for record in source_refinement.get("projection_types", []):
                if record.get("status") != "infeasible":
                    continue
                cut_active = record.get(
                    "infeasible_cut_active_row_indices",
                    record["active_row_indices"],
                )
                seed_projection_cuts.extend(
                    {
                        "coordinate": int(coordinate),
                        "infeasible_cut_active_row_indices": cut_active,
                    }
                    for coordinate in record["coordinates"]
                )
        seed_pair_cuts = source_refinement.get("pair_benders_cuts", [])
    placement, discovery = discover_refined_placement(
        rows,
        axis,
        pattern_period,
        args.projection,
        args.time_limit,
        args.subproblem_time_limit,
        args.cut_rounds,
        seed_cache_payload,
        seed_projection_cuts,
        seed_pair_cuts,
    )
    full_replay = exact_capacity_replay(
        rows,
        axis,
        int(source["layer_period"]),
        placement,
    )
    materialized = materialize_rows(rows, axis, placement)
    root_source, root_source_hash = refinement_root(args.input, source)
    result = {
        **source,
        "schema": "axis_layered_projected_refinement_v1",
        "refinement_source": str(args.input),
        "refinement_source_sha256": sha256(args.input),
        "refinement_root_source": root_source,
        "refinement_root_source_sha256": root_source_hash,
        "projection_refinement": discovery,
        "exact_capacity_replay": full_replay,
        "choices": materialized,
        "scope": (
            "discovery placement with exact ordinary and small-projection "
            "capacity replay; projected phases remain independently relaxed "
            "between layer columns, so this is not a lattice cover"
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"projection_safe={discovery['projection_safe']} "
        f"pair_safe={discovery['pair_safe']} "
        f"rows={len(rows)} output={args.output}",
        flush=True,
    )
    return 0 if discovery["projection_safe"] and discovery["pair_safe"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
