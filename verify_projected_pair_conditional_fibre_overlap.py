#!/usr/bin/env python3
"""Independent replay of a structured conditional-fibre certificate."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

try:
    import numpy as np
except ModuleNotFoundError:  # pragma: no cover - explicit runtime guard
    np = None

from verify_projected_pair_start_leaf_block import brute_path_leaf_density


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def recorded_row(row: dict) -> dict:
    return {key: int(row[key]) for key in ("p", "h", "a", "b")}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if np is None:
        raise RuntimeError("NumPy is required for exact replay")

    source = json.loads(args.pool.read_text())
    cert = json.loads(args.certificate.read_text())
    block = json.loads(Path(cert["block_certificate"]).read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    structured = block
    structured_path = Path(cert["block_certificate"])
    structured_schemas = {
        "projected_pair_normalized_start_leaf_v2",
        "projected_pair_normalized_start_leaf_v3",
    }
    seen_paths = set()
    while structured.get("schema") not in structured_schemas:
        base_name = structured.get("base_certificate")
        if not base_name:
            raise RuntimeError("block has no projected-pair structured base")
        next_path = Path(base_name)
        if not next_path.is_absolute() and not next_path.exists():
            next_path = structured_path.parent / next_path
        resolved = next_path.resolve()
        if resolved in seen_paths:
            raise RuntimeError("cycle in block certificate chain")
        seen_paths.add(resolved)
        structured_path = next_path
        structured = json.loads(structured_path.read_text())
    structured_schema = structured["schema"]
    has_paired_projection = (
        structured_schema == "projected_pair_normalized_start_leaf_v3"
    )
    base_primes = tuple(int(value) for value in structured["base_primes"])
    projected_prime = int(structured["projected_prime"])
    lifted_prime = int(structured["lifted_prime"])
    core_independent_primes = tuple(
        int(value) for value in structured["independent_projected_primes"]
    )
    paired_projected_prime = (
        int(structured["paired_projected_prime"])
        if has_paired_projection else None
    )
    paired_shared_prime = (
        int(structured["paired_shared_prime"])
        if has_paired_projection else None
    )
    leaf_prime = int(structured["normalized_start_shared_prime"])
    outside_prime = int(cert["outside_prime"])
    structured_anchor_primes = (
        *base_primes,
        projected_prime,
        lifted_prime,
        *core_independent_primes,
        *((paired_projected_prime, paired_shared_prime)
          if has_paired_projection else ()),
        leaf_prime,
    )
    block_anchor_primes = tuple(int(value) for value in block["anchor_primes"])
    extra_independent_primes = tuple(
        prime
        for prime in block_anchor_primes
        if prime not in structured_anchor_primes
    )
    independent_primes = (
        *core_independent_primes,
        *extra_independent_primes,
    )
    include_leaf = bool(cert.get("start_leaf_included", True))
    anchor_primes = tuple(
        prime
        for prime in block_anchor_primes
        if include_leaf or prime != leaf_prime
    )
    missing = {outside_prime, *block_anchor_primes} - by_prime.keys()
    if missing:
        raise RuntimeError(f"certificate rows absent from pool: {sorted(missing)}")

    bases = [by_prime[prime] for prime in base_primes]
    projected = by_prime[projected_prime]
    lifted = by_prime[lifted_prime]
    independents = [by_prime[prime] for prime in independent_primes]
    paired_projected = (
        by_prime[paired_projected_prime]
        if has_paired_projection else None
    )
    paired_shared = (
        by_prime[paired_shared_prime]
        if has_paired_projection else None
    )
    leaf = by_prime[leaf_prime]
    outside = by_prime[outside_prime]
    anchors = [by_prime[prime] for prime in anchor_primes]

    base_period = math.lcm(*(int(row["h"]) for row in bases))
    base_cells = base_period**2
    projection = math.gcd(base_period, int(projected["h"]))
    shared = int(projected["h"]) // projection
    lifted_residual = int(lifted["h"]) // shared
    independent_projections = tuple(
        math.gcd(base_period, int(row["h"])) for row in independents
    )
    independent_residuals = tuple(
        int(row["h"]) // projection_value
        for row, projection_value in zip(
            independents,
            independent_projections,
        )
    )
    paired_projection = (
        math.gcd(base_period, int(paired_projected["h"]))
        if has_paired_projection else None
    )
    paired_residual = (
        int(paired_projected["h"]) // paired_projection
        if has_paired_projection else None
    )
    paired_determinant = (
        int(paired_projected["a"]) * int(paired_shared["b"])
        - int(paired_projected["b"]) * int(paired_shared["a"])
        if has_paired_projection else None
    )
    leaf_projection = math.gcd(base_period, int(leaf["h"]))
    leaf_residual = int(leaf["h"]) // leaf_projection
    determinants = (
        int(projected["a"]) * int(lifted["b"])
        - int(projected["b"]) * int(lifted["a"]),
        int(projected["a"]) * int(leaf["b"])
        - int(projected["b"]) * int(leaf["a"]),
        int(lifted["a"]) * int(leaf["b"])
        - int(lifted["b"]) * int(leaf["a"]),
    )
    residuals = (
        shared,
        lifted_residual,
        *independent_residuals,
        *((paired_residual,) if has_paired_projection else ()),
    )
    structural = (
        structured.get("schema") in structured_schemas
        and len(by_prime) == len(rows)
        and len(set(anchor_primes)) == len(anchor_primes)
        and block_anchor_primes
        == (*structured_anchor_primes, *extra_independent_primes)
        and outside_prime not in anchor_primes
        and base_period % int(outside["h"]) == 0
        and all(int(row.get("target_modulus", 1)) == 1 for row in anchors)
        and min(
            projection,
            leaf_projection,
            *independent_projections[:len(core_independent_primes)],
            *((paired_projection,) if has_paired_projection else ()),
        )
        > 1
        and all(
            value == 1
            for value in independent_projections[
                len(core_independent_primes):
            ]
        )
        and min(residuals) > 1
        and len(set(residuals)) == len(residuals)
        and all(
            math.gcd(left, right) == 1
            for left, right in itertools.combinations(residuals, 2)
        )
        and all(math.gcd(value, base_period) == 1 for value in residuals)
        and int(projected["h"]) == projection * shared
        and int(lifted["h"]) == shared * lifted_residual
        and all(
            int(row["h"]) == projection_value * residual_value
            for row, projection_value, residual_value in zip(
                independents,
                independent_projections,
                independent_residuals,
            )
        )
        and (
            not has_paired_projection
            or (
                int(paired_projected["h"])
                == paired_projection * paired_residual
                and int(paired_shared["h"]) == paired_residual
                and math.gcd(paired_determinant, paired_residual) == 1
            )
        )
        and int(leaf["h"]) == leaf_projection * shared
        and leaf_residual == shared
        and all(math.gcd(value, shared) == 1 for value in determinants)
    )

    normalization_primes = tuple(
        int(value) for value in cert["normalization_primes"]
    )
    normalized_indices = {
        base_primes.index(prime) for prime in normalization_primes
    }
    normalizers = [
        outside,
        *[by_prime[prime] for prime in normalization_primes],
        *([leaf] if include_leaf else []),
    ]
    normalization_period = math.lcm(
        *(int(row["h"]) for row in normalizers)
    )
    # The verifier reverses the plane-loop order and reconstructs the image
    # directly, independently of the certificate's recorded image size.
    normalization_image = {
        tuple(
            (int(row["a"]) * first + int(row["b"]) * second) % int(row["h"])
            for row in normalizers
        )
        for second in range(normalization_period)
        for first in range(normalization_period)
    }
    normalization_target_count = math.prod(
        int(row["h"]) for row in normalizers
    )
    normalization_surjective = (
        len(normalization_image) == normalization_target_count
    )
    target_ranges = [
        range(1) if index in normalized_indices else range(int(row["h"]))
        for index, row in enumerate(bases)
    ]
    target_tuples = math.prod(len(values) for values in target_ranges)

    path_values = {
        "inactive": brute_path_leaf_density(
            projected,
            lifted,
            leaf,
            shared,
            lifted_residual,
            projected_active=False,
            leaf_active=False,
            compatible=False,
        ),
        "one": brute_path_leaf_density(
            projected,
            lifted,
            leaf,
            shared,
            lifted_residual,
            projected_active=True,
            leaf_active=False,
            compatible=False,
        ),
        "both_compatible": brute_path_leaf_density(
            projected,
            lifted,
            leaf,
            shared,
            lifted_residual,
            projected_active=True,
            leaf_active=True,
            compatible=True,
        ),
        "both_incompatible": brute_path_leaf_density(
            projected,
            lifted,
            leaf,
            shared,
            lifted_residual,
            projected_active=True,
            leaf_active=True,
            compatible=False,
        ),
    }
    if path_values["both_incompatible"] > path_values["both_compatible"]:
        raise RuntimeError("incompatible residual targets are not minimal")

    projection_moduli = (
        projection,
        *independent_projections,
        *((paired_projection,) if has_paired_projection else ()),
    )
    weights = {}
    for state in itertools.product(
        (False, True),
        repeat=len(projection_moduli) + int(include_leaf),
    ):
        main_active = state[0]
        independent_activity = state[1:1 + len(independents)]
        paired_active = (
            state[1 + len(independents)]
            if has_paired_projection else False
        )
        leaf_active = state[-1] if include_leaf else False
        if main_active and leaf_active:
            path_density = path_values["both_incompatible"]
        elif main_active or leaf_active:
            path_density = path_values["one"]
        else:
            path_density = path_values["inactive"]
        uncovered_weight = 1 - path_density
        for active, residual_value in zip(
            independent_activity,
            independent_residuals,
        ):
            if active:
                uncovered_weight *= Fraction(
                    residual_value - 1,
                    residual_value,
                )
        if has_paired_projection:
            paired_density = (
                Fraction(2, paired_residual)
                - Fraction(1, paired_residual * paired_residual)
                if paired_active
                else Fraction(1, paired_residual)
            )
            uncovered_weight *= 1 - paired_density
        weights[state] = 1 - uncovered_weight
    denominator = math.lcm(
        *(value.denominator for value in weights.values())
    )

    def transposed_masks(row: dict, modulus: int) -> list[int]:
        masks = [0] * modulus
        for second in range(base_period):
            for first in range(base_period):
                value = (
                    int(row["a"]) * first + int(row["b"]) * second
                ) % modulus
                masks[value] |= 1 << (second * base_period + first)
        return masks

    base_masks = [
        transposed_masks(row, int(row["h"])) for row in bases
    ]
    outside_mask = transposed_masks(outside, int(outside["h"]))[0]
    outside_base_points = outside_mask.bit_count()
    projection_rows = [
        projected,
        *independents,
        *((paired_projected,) if has_paired_projection else ()),
    ]
    projection_masks = [
        transposed_masks(row, modulus)
        for row, modulus in zip(projection_rows, projection_moduli)
    ]
    leaf_active_mask = (
        transposed_masks(leaf, leaf_projection)[0]
        if include_leaf
        else None
    )
    full_grid_mask = (1 << base_cells) - 1
    observed_codes = tuple(
        itertools.product(
            *(range(modulus) for modulus in projection_moduli),
            *((False, True),) if include_leaf else (),
        )
    )
    observed_masks = []
    for code in observed_codes:
        cells = full_grid_mask
        for masks, value in zip(
            projection_masks,
            code[:-1] if include_leaf else code,
        ):
            cells &= masks[value]
        if include_leaf:
            cells &= (
                leaf_active_mask
                if code[-1]
                else full_grid_mask ^ leaf_active_mask
            )
        observed_masks.append(cells)
    partition_valid = (
        sum(mask.bit_count() for mask in observed_masks) == base_cells
    )
    choices = tuple(
        itertools.product(
            *(range(modulus) for modulus in projection_moduli)
        )
    )
    # The verifier contracts in the transpose orientation: category rows by
    # projection choices, whereas the certifier stores choices by categories.
    cell_choice_scores = np.array(
        [
            [
                int(
                    weights[
                        (
                            *(
                                observed == target
                                for observed, target in zip(
                                    code[:-1] if include_leaf else code,
                                    choice,
                                )
                            ),
                            *((code[-1],) if include_leaf else ()),
                        )
                    ]
                    * denominator
                )
                for choice in choices
            ]
            for code in observed_codes
        ],
        dtype=np.int64,
    )

    minimum_score = None
    minimizing_targets = None
    minimizing_projections = None
    minimizing_base_covered = None
    minimizing_counts = None
    checked = 0
    for targets in itertools.product(*target_ranges):
        base_union = 0
        for masks, target in zip(base_masks, targets):
            base_union |= masks[target]
        base_union &= outside_mask
        uncovered = outside_mask ^ base_union
        histogram = np.fromiter(
            (
                (uncovered & category_mask).bit_count()
                for category_mask in observed_masks
            ),
            dtype=np.int64,
            count=len(observed_masks),
        )
        scores = histogram @ cell_choice_scores
        scores += base_union.bit_count() * denominator
        choice_index = int(np.argmin(scores))
        score = int(scores[choice_index])
        checked += len(choices)
        if minimum_score is None or score < minimum_score:
            minimum_score = score
            minimizing_targets = targets
            minimizing_projections = choices[choice_index]
            minimizing_base_covered = base_union.bit_count()
            minimizing_counts = {
                "".join("1" if flag else "0" for flag in state): 0
                for state in itertools.product(
                    (False, True),
                    repeat=len(projection_moduli) + int(include_leaf),
                )
            }
            for code, count in zip(observed_codes, histogram.tolist()):
                state = (
                    *(
                        observed == target
                        for observed, target in zip(
                            code[:-1] if include_leaf else code,
                            minimizing_projections,
                        )
                    ),
                    *((code[-1],) if include_leaf else ()),
                )
                key = "".join("1" if flag else "0" for flag in state)
                minimizing_counts[key] += int(count)

    intersection = Fraction(
        int(minimum_score),
        base_cells * denominator,
    )
    anchors_match = (
        list(cert["anchor_primes"]) == list(anchor_primes)
        and len(cert["anchor_rows"]) == len(anchors)
        and all(
            int(record["p"]) in by_prime
            and recorded_row(by_prime[int(record["p"])])
            == recorded_row(record)
            for record in cert["anchor_rows"]
        )
        and {
            int(record["p"]): recorded_row(record)
            for record in block["anchor_rows"]
        }
        == {
            prime: recorded_row(by_prime[prime])
            for prime in block_anchor_primes
        }
    )
    recorded_weights = {
        key: read_fraction(value)
        for key, value in cert["residual_union_weights"].items()
    }
    expected_weights = {
        "".join("1" if flag else "0" for flag in state): value
        for state, value in weights.items()
    }
    verified = (
        cert.get("schema")
        == "projected_pair_conditional_fibre_overlap_v1"
        and cert.get("structured_base_schema", structured_schema)
        == structured_schema
        and str(args.pool) == cert["pool"]
        and int(cert["row_count"]) == len(rows)
        and recorded_row(cert["outside_row"]) == recorded_row(outside)
        and anchors_match
        and structural
        and normalization_surjective
        and int(cert["normalization_period"]) == normalization_period
        and int(cert["normalization_target_count"])
        == normalization_target_count
        and bool(cert["normalization_jointly_surjective"])
        and int(cert["base_period"]) == base_period
        and int(cert["base_cells"]) == base_cells
        and int(cert["outside_base_points"]) == outside_base_points
        and outside_base_points == base_cells // int(outside["h"])
        and list(cert["projection_moduli"]) == list(projection_moduli)
        and int(cert["fixed_leaf_projection_modulus"])
        == leaf_projection
        and cert["fixed_leaf_projection_target"]
        == (0 if include_leaf else None)
        and list(cert["residual_moduli"]) == list(residuals)
        and cert["residual_target_relation"]
        == ("incompatible" if include_leaf else "leaf_omitted")
        and all(
            read_fraction(cert["path_leaf_densities"][key]) == value
            for key, value in path_values.items()
        )
        and recorded_weights == expected_weights
        and partition_valid
        and int(cert["target_tuples"]) == target_tuples
        and int(cert["projection_target_combinations_per_base"])
        == len(choices)
        and int(cert["target_combinations_checked"]) == checked
        and list(cert["minimizing_base_targets"])
        == list(minimizing_targets)
        and list(cert["minimizing_projection_targets"])
        == list(minimizing_projections)
        and int(cert["minimizing_base_covered_cells"])
        == minimizing_base_covered
        and cert["minimizing_uncovered_category_counts"]
        == minimizing_counts
        and int(cert["minimum_intersection_score"]) == minimum_score
        and int(cert["score_denominator"]) == denominator
        and read_fraction(cert["forced_intersection_density"])
        == intersection
    )
    result = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "block_certificate": cert["block_certificate"],
        "outside_prime": outside_prime,
        "anchor_count": len(anchors),
        "base_period": base_period,
        "normalization_primes": list(normalization_primes),
        "target_combinations_checked": checked,
        "forced_intersection_density": {
            "numerator": intersection.numerator,
            "denominator": intersection.denominator,
        },
        "verified": verified,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"outside={outside_prime} anchors={len(anchors)} "
        f"checked={checked} intersection={intersection} "
        f"verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
