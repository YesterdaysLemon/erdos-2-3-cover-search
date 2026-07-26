#!/usr/bin/env python3
"""Independent replay of a projected pair with a normalized start leaf."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

try:
    import numpy as np
except ModuleNotFoundError:
    np = None


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def brute_path_leaf_density(
    projected: dict,
    lifted: dict,
    leaf: dict,
    shared: int,
    lifted_residual: int,
    *,
    projected_active: bool,
    leaf_active: bool,
    compatible: bool,
) -> Fraction:
    """Count the residual union directly on the transposed product grid."""
    covered = 0
    cells = shared**2 * lifted_residual**2
    lifted_shared_target = 0 if compatible else 1
    for shared_second in range(shared):
        for shared_first in range(shared):
            projected_event = (
                (
                    int(projected["a"]) * shared_first
                    + int(projected["b"]) * shared_second
                )
                % shared
                == 0
            )
            leaf_event = (
                (
                    int(leaf["a"]) * shared_first
                    + int(leaf["b"]) * shared_second
                )
                % shared
                == 0
            )
            lifted_shared_event = (
                (
                    int(lifted["a"]) * shared_first
                    + int(lifted["b"]) * shared_second
                )
                % shared
                == lifted_shared_target
            )
            for lifted_second in range(lifted_residual):
                for lifted_first in range(lifted_residual):
                    lifted_event = (
                        lifted_shared_event
                        and (
                            int(lifted["a"]) * lifted_first
                            + int(lifted["b"]) * lifted_second
                        )
                        % lifted_residual
                        == 0
                    )
                    covered += int(
                        lifted_event
                        or (projected_active and projected_event)
                        or (leaf_active and leaf_event)
                    )
    return Fraction(covered, cells)


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
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    base_primes = tuple(int(value) for value in cert["base_primes"])
    projected_prime = int(cert["projected_prime"])
    lifted_prime = int(cert["lifted_prime"])
    independent_primes = tuple(
        int(value)
        for value in cert.get(
            "independent_projected_primes",
            [cert["independent_projected_prime"]],
        )
    )
    paired_projected_prime = cert.get("paired_projected_prime")
    paired_shared_prime = cert.get("paired_shared_prime")
    paired_primes = (
        (
            int(paired_projected_prime),
            int(paired_shared_prime),
        )
        if paired_projected_prime is not None
        and paired_shared_prime is not None
        else ()
    )
    leaf_prime = int(cert["normalized_start_shared_prime"])
    anchor_primes = (
        *base_primes,
        projected_prime,
        lifted_prime,
        *independent_primes,
        *paired_primes,
        leaf_prime,
    )
    missing = set(anchor_primes) - by_prime.keys()
    if missing:
        raise RuntimeError(f"certificate anchors absent: {sorted(missing)}")
    bases = [by_prime[prime] for prime in base_primes]
    projected = by_prime[projected_prime]
    lifted = by_prime[lifted_prime]
    independents = [by_prime[prime] for prime in independent_primes]
    paired_projected = (
        by_prime[int(paired_projected_prime)]
        if paired_projected_prime is not None
        else None
    )
    paired_shared = (
        by_prime[int(paired_shared_prime)]
        if paired_shared_prime is not None
        else None
    )
    leaf = by_prime[leaf_prime]
    anchors = [
        *bases,
        projected,
        lifted,
        *independents,
        *(
            [paired_projected, paired_shared]
            if paired_projected is not None
            else []
        ),
        leaf,
    ]

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
    paired_projection = None
    paired_residual = None
    paired_determinant = None
    if paired_projected is not None:
        paired_projection = math.gcd(
            base_period,
            int(paired_projected["h"]),
        )
        paired_residual = (
            int(paired_projected["h"]) // paired_projection
        )
        paired_determinant = (
            int(paired_projected["a"]) * int(paired_shared["b"])
            - int(paired_projected["b"]) * int(paired_shared["a"])
        )
    leaf_projection = math.gcd(base_period, int(leaf["h"]))
    leaf_residual = int(leaf["h"]) // leaf_projection
    projected_lifted_determinant = (
        int(projected["a"]) * int(lifted["b"])
        - int(projected["b"]) * int(lifted["a"])
    )
    projected_leaf_determinant = (
        int(projected["a"]) * int(leaf["b"])
        - int(projected["b"]) * int(leaf["a"])
    )
    lifted_leaf_determinant = (
        int(lifted["a"]) * int(leaf["b"])
        - int(lifted["b"]) * int(leaf["a"])
    )
    residuals = (
        shared,
        lifted_residual,
        *independent_residuals,
        *([paired_residual] if paired_residual is not None else []),
    )
    structural = (
        len(set(anchor_primes)) == len(anchor_primes)
        and all(int(row.get("target_modulus", 1)) == 1 for row in anchors)
        and min(
            projection,
            leaf_projection,
            *independent_projections,
            *([paired_projection] if paired_projection is not None else []),
        )
        > 1
        and min(residuals) > 1
        and len(set(residuals)) == len(residuals)
        and all(
            math.gcd(first, second) == 1
            for first, second in itertools.combinations(residuals, 2)
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
        and int(leaf["h"]) == leaf_projection * shared
        and leaf_residual == shared
        and (
            paired_projected is None
            or (
                int(paired_projected["h"])
                == paired_projection * paired_residual
                and int(paired_shared["h"]) == paired_residual
                and math.gcd(paired_determinant, paired_residual) == 1
            )
        )
        and all(
            math.gcd(determinant, shared) == 1
            for determinant in (
                projected_lifted_determinant,
                projected_leaf_determinant,
                lifted_leaf_determinant,
            )
        )
    )

    normalization_primes = tuple(
        int(value) for value in cert["normalization_primes"]
    )
    normalized_indices = {
        base_primes.index(prime) for prime in normalization_primes
    }
    normalizers = [
        *[by_prime[prime] for prime in normalization_primes],
        leaf,
    ]
    normalization_period = math.lcm(
        *(int(row["h"]) for row in normalizers)
    )
    normalization_image = {
        tuple(
            (int(row["a"]) * k + int(row["b"]) * l) % int(row["h"])
            for row in normalizers
        )
        for l in range(normalization_period)
        for k in range(normalization_period)
    }
    normalization_surjective = len(normalization_image) == math.prod(
        int(row["h"]) for row in normalizers
    )
    target_ranges = [
        range(1) if index in normalized_indices else range(int(row["h"]))
        for index, row in enumerate(bases)
    ]
    target_tuples = math.prod(len(values) for values in target_ranges)

    path_densities = {}
    for label, compatible in (
        ("compatible", True),
        ("incompatible", False),
    ):
        for projected_active, leaf_active in itertools.product(
            (False, True),
            repeat=2,
        ):
            path_densities[
                (label, projected_active, leaf_active)
            ] = brute_path_leaf_density(
                projected,
                lifted,
                leaf,
                shared,
                lifted_residual,
                projected_active=projected_active,
                leaf_active=leaf_active,
                compatible=compatible,
            )
    if any(
        path_densities[("compatible", first, second)]
        < path_densities[("incompatible", first, second)]
        for first, second in itertools.product((False, True), repeat=2)
    ):
        raise RuntimeError("compatible residual targets are not maximal")
    weights = {}
    projection_moduli = (
        projection,
        *independent_projections,
        *([paired_projection] if paired_projection is not None else []),
    )
    for state in itertools.product(
        (False, True),
        repeat=len(projection_moduli) + 1,
    ):
        main_active = state[0]
        independent_activity = state[1:1 + len(independents)]
        paired_active = (
            state[1 + len(independents)]
            if paired_projected is not None
            else None
        )
        leaf_active = state[-1]
        uncovered = (
            1
            - path_densities[
                ("compatible", main_active, leaf_active)
            ]
        )
        for active, residual_value in zip(
            independent_activity,
            independent_residuals,
        ):
            if active:
                uncovered *= Fraction(
                    residual_value - 1,
                    residual_value,
                )
        if paired_projected is not None:
            paired_density = (
                Fraction(2, paired_residual)
                - Fraction(1, paired_residual * paired_residual)
                if paired_active
                else Fraction(1, paired_residual)
            )
            uncovered *= 1 - paired_density
        weights[state] = 1 - uncovered
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
    projection_rows = [
        projected,
        *independents,
        *(
            [paired_projected]
            if paired_projected is not None
            else []
        ),
    ]
    projection_masks = [
        transposed_masks(row, modulus)
        for row, modulus in zip(projection_rows, projection_moduli)
    ]
    leaf_active_mask = transposed_masks(leaf, leaf_projection)[0]
    full_grid_mask = (1 << base_cells) - 1
    observed_codes = tuple(
        itertools.product(
            *(range(modulus) for modulus in projection_moduli),
            (False, True),
        )
    )
    observed_masks = []
    for code in observed_codes:
        cells = full_grid_mask
        for masks, value in zip(projection_masks, code[:-1]):
            cells &= masks[value]
        cells &= (
            leaf_active_mask
            if code[-1]
            else full_grid_mask ^ leaf_active_mask
        )
        observed_masks.append(cells)
    choices = tuple(
        itertools.product(
            *(range(modulus) for modulus in projection_moduli),
        )
    )
    cell_choice_scores = np.array(
        [
            [
                int(
                    weights[
                        (
                            *(
                                observed == target
                                for observed, target in zip(
                                    code[:-1],
                                    choice,
                                )
                            ),
                            code[-1],
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

    maximum_score = -1
    maximizing_targets = None
    maximizing_projections = None
    maximizing_base_covered = None
    maximizing_counts = None
    checked = 0
    for targets in itertools.product(*target_ranges):
        base_union = 0
        for masks, target in zip(base_masks, targets):
            base_union |= masks[target]
        uncovered = full_grid_mask ^ base_union
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
        choice_index = int(np.argmax(scores))
        score = int(scores[choice_index])
        checked += len(choices)
        if score > maximum_score:
            maximum_score = score
            maximizing_targets = targets
            maximizing_projections = choices[choice_index]
            maximizing_base_covered = base_union.bit_count()
            maximizing_counts = {
                "".join("1" if flag else "0" for flag in state): 0
                for state in itertools.product(
                    (False, True),
                    repeat=len(projection_moduli) + 1,
                )
            }
            for code, count in zip(observed_codes, histogram.tolist()):
                state = (
                    *(
                        observed == target
                        for observed, target in zip(
                            code[:-1],
                            maximizing_projections,
                        )
                    ),
                    code[-1],
                )
                key = "".join("1" if flag else "0" for flag in state)
                maximizing_counts[key] += int(count)

    maximum = Fraction(maximum_score, base_cells * denominator)
    individual_sum = sum(
        (Fraction(1, int(row["h"])) for row in anchors),
        Fraction(0),
    )
    loss = individual_sum - maximum
    total = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    upper = total - loss
    recorded_path = {
        "inactive": path_densities[("compatible", False, False)],
        "one": path_densities[("compatible", True, False)],
        "both_compatible": path_densities[
            ("compatible", True, True)
        ],
        "both_incompatible": path_densities[
            ("incompatible", True, True)
        ],
    }
    anchors_match = (
        list(cert["anchor_primes"]) == list(anchor_primes)
        and len(cert["anchor_rows"]) == len(anchors)
        and all(
            int(record["p"]) in by_prime
            and all(
                int(by_prime[int(record["p"])][key]) == int(record[key])
                for key in ("h", "a", "b")
            )
            for record in cert["anchor_rows"]
        )
    )
    verified = (
        cert.get("schema")
        == (
            "projected_pair_normalized_start_leaf_v3"
            if paired_projected is not None
            else "projected_pair_normalized_start_leaf_v1"
            if len(independent_primes) == 1
            else "projected_pair_normalized_start_leaf_v2"
        )
        and str(args.pool) == cert["pool"]
        and int(cert["row_count"]) == len(rows)
        and anchors_match
        and structural
        and normalization_surjective
        and int(cert["normalization_period"]) == normalization_period
        and bool(cert["normalization_jointly_surjective"])
        == normalization_surjective
        and int(cert["base_period"]) == base_period
        and int(cert["base_cells"]) == base_cells
        and list(cert["projection_moduli"]) == list(projection_moduli)
        and int(cert["fixed_leaf_projection_modulus"])
        == leaf_projection
        and int(cert["fixed_leaf_projection_target"]) == 0
        and int(cert["shared_residual_modulus"]) == shared
        and int(cert["lifted_residual_modulus"]) == lifted_residual
        and int(cert["independent_residual_modulus"])
        == independent_residuals[0]
        and list(
            cert.get(
                "independent_residual_moduli",
                [cert["independent_residual_modulus"]],
            )
        )
        == list(independent_residuals)
        and cert.get("paired_projection_modulus") == paired_projection
        and cert.get("paired_residual_modulus") == paired_residual
        and cert.get("paired_determinant") == paired_determinant
        and cert.get("paired_maps_jointly_surjective")
        == (True if paired_projected is not None else None)
        and int(cert["projected_lifted_determinant"])
        == projected_lifted_determinant
        and int(cert["projected_leaf_determinant"])
        == projected_leaf_determinant
        and int(cert["lifted_leaf_determinant"])
        == lifted_leaf_determinant
        and all(
            read_fraction(cert["path_leaf_densities"][key]) == value
            for key, value in recorded_path.items()
        )
        and bool(cert["compatible_residual_targets_dominate"])
        and {
            key: read_fraction(value)
            for key, value in cert["residual_union_weights"].items()
        }
        == {
            "".join("1" if flag else "0" for flag in state): value
            for state, value in weights.items()
        }
        and int(cert["enumerated_period"])
        == base_period * math.prod(residuals)
        and int(cert["target_tuples"]) == target_tuples
        and int(cert["projection_target_combinations_per_base"])
        == len(choices)
        and int(cert["target_combinations_checked"]) == checked
        and list(cert["maximizing_base_targets"]) == list(maximizing_targets)
        and list(cert["maximizing_projection_targets"])
        == list(maximizing_projections)
        and int(cert["maximizing_base_covered_cells"])
        == maximizing_base_covered
        and cert["maximizing_uncovered_category_counts"]
        == maximizing_counts
        and read_fraction(cert["maximum_block_union_density"]) == maximum
        and read_fraction(cert["block_individual_density_sum"])
        == individual_sum
        and read_fraction(cert["forced_overlap_loss"]) == loss
        and read_fraction(cert["total_pool_density"]) == total
        and read_fraction(cert["pool_union_density_upper_bound"]) == upper
        and bool(cert["proved_no_cover"]) == (upper < 1)
    )
    result = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "anchor_primes": list(anchor_primes),
        "base_period": base_period,
        "projection_moduli": list(projection_moduli),
        "fixed_leaf_projection_modulus": leaf_projection,
        "residual_moduli": list(residuals),
        "enumerated_period": base_period * math.prod(residuals),
        "target_tuples": target_tuples,
        "target_combinations_checked": checked,
        "maximum_block_union_density": {
            "numerator": maximum.numerator,
            "denominator": maximum.denominator,
        },
        "forced_overlap_loss": {
            "numerator": loss.numerator,
            "denominator": loss.denominator,
        },
        "pool_union_density_upper_bound": {
            "numerator": upper.numerator,
            "denominator": upper.denominator,
        },
        "verified": verified,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base_targets={target_tuples} projection_targets={len(choices)} "
        f"checked={checked} maximum={maximum} loss={loss} "
        f"pool_upper={upper} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
