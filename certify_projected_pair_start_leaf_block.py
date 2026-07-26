#!/usr/bin/env python3
"""Certify a projected pair with a normalized leaf on its shared residual."""

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


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def path_leaf_densities(
    shared: int,
    lifted_residual: int,
) -> dict[str, Fraction]:
    inactive = Fraction(1, shared * lifted_residual)
    one = (
        Fraction(1, shared)
        + inactive
        - Fraction(1, shared * shared * lifted_residual)
    )
    both_compatible = (
        Fraction(2, shared)
        - Fraction(1, shared * shared)
        + inactive
        - Fraction(1, shared * shared * lifted_residual)
    )
    both_incompatible = (
        both_compatible
        - Fraction(1, shared * shared * lifted_residual)
    )
    return {
        "inactive": inactive,
        "one": one,
        "both_compatible": both_compatible,
        "both_incompatible": both_incompatible,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--base-primes", required=True)
    parser.add_argument("--projected-prime", type=int, required=True)
    parser.add_argument("--lifted-prime", type=int, required=True)
    parser.add_argument(
        "--independent-projected-prime",
        type=int,
        required=True,
    )
    parser.add_argument("--second-independent-projected-prime", type=int)
    parser.add_argument("--third-independent-projected-prime", type=int)
    parser.add_argument("--paired-projected-prime", type=int)
    parser.add_argument("--paired-shared-prime", type=int)
    parser.add_argument(
        "--normalized-start-shared-prime",
        type=int,
        required=True,
    )
    parser.add_argument("--normalize-primes", required=True)
    parser.add_argument("--max-base-cells", type=int, default=1_000_000)
    parser.add_argument("--max-target-tuples", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if np is None:
        raise RuntimeError("NumPy is required for exact score replay")

    base_primes = tuple(
        int(value) for value in args.base_primes.split(",") if value
    )
    normalization_primes = tuple(
        int(value) for value in args.normalize_primes.split(",") if value
    )
    independent_primes = (
        args.independent_projected_prime,
        *(
            (args.second_independent_projected_prime,)
            if args.second_independent_projected_prime is not None
            else ()
        ),
        *(
            (args.third_independent_projected_prime,)
            if args.third_independent_projected_prime is not None
            else ()
        ),
    )
    paired_primes = (
        (
            args.paired_projected_prime,
            args.paired_shared_prime,
        )
        if args.paired_projected_prime is not None
        and args.paired_shared_prime is not None
        else ()
    )
    anchor_primes = (
        *base_primes,
        args.projected_prime,
        args.lifted_prime,
        *independent_primes,
        *paired_primes,
        args.normalized_start_shared_prime,
    )
    if (
        len(base_primes) < 2
        or len(set(anchor_primes)) != len(anchor_primes)
        or len(normalization_primes) != 2
        or not set(normalization_primes) <= set(base_primes)
        or (
            (args.paired_projected_prime is None)
            != (args.paired_shared_prime is None)
        )
    ):
        raise SystemExit("invalid anchors or normalization")

    source = json.loads(args.pool.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    missing = set(anchor_primes) - by_prime.keys()
    if missing:
        raise RuntimeError(f"anchor primes are absent: {sorted(missing)}")
    bases = [by_prime[prime] for prime in base_primes]
    projected = by_prime[args.projected_prime]
    lifted = by_prime[args.lifted_prime]
    independents = [by_prime[prime] for prime in independent_primes]
    paired_projected = (
        by_prime[args.paired_projected_prime]
        if args.paired_projected_prime is not None
        else None
    )
    paired_shared = (
        by_prime[args.paired_shared_prime]
        if args.paired_shared_prime is not None
        else None
    )
    start_leaf = by_prime[args.normalized_start_shared_prime]
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
        start_leaf,
    ]
    if any(int(row.get("target_modulus", 1)) != 1 for row in anchors):
        raise RuntimeError("anchor targets must be unrestricted")
    if any(
        math.gcd(int(row["a"]), int(row["b"]), int(row["h"])) != 1
        for row in anchors[len(bases):]
    ):
        raise RuntimeError("an extra row target map is not surjective")

    base_period = math.lcm(*(int(row["h"]) for row in bases))
    base_cells = base_period * base_period
    if base_cells > args.max_base_cells:
        raise RuntimeError("base grid exceeds guard")
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
    leaf_projection = math.gcd(base_period, int(start_leaf["h"]))
    leaf_residual = int(start_leaf["h"]) // leaf_projection
    projected_lifted_determinant = (
        int(projected["a"]) * int(lifted["b"])
        - int(projected["b"]) * int(lifted["a"])
    )
    projected_leaf_determinant = (
        int(projected["a"]) * int(start_leaf["b"])
        - int(projected["b"]) * int(start_leaf["a"])
    )
    lifted_leaf_determinant = (
        int(lifted["a"]) * int(start_leaf["b"])
        - int(lifted["b"]) * int(start_leaf["a"])
    )
    residuals = (
        shared,
        lifted_residual,
        *independent_residuals,
        *([paired_residual] if paired_residual is not None else []),
    )
    if (
        min(
            projection,
            leaf_projection,
            *independent_projections,
            *([paired_projection] if paired_projection is not None else []),
        )
        <= 1
        or min(residuals) <= 1
        or len(set(residuals)) != len(residuals)
        or any(
            math.gcd(first, second) != 1
            for first, second in itertools.combinations(residuals, 2)
        )
        or any(math.gcd(value, base_period) != 1 for value in residuals)
        or int(projected["h"]) != projection * shared
        or int(lifted["h"]) != shared * lifted_residual
        or any(
            int(row["h"]) != projection_value * residual_value
            for row, projection_value, residual_value in zip(
                independents,
                independent_projections,
                independent_residuals,
            )
        )
        or int(start_leaf["h"]) != leaf_projection * shared
        or leaf_residual != shared
        or (
            paired_projected is not None
            and (
                int(paired_projected["h"])
                != paired_projection * paired_residual
                or int(paired_shared["h"]) != paired_residual
                or math.gcd(paired_determinant, paired_residual) != 1
            )
        )
        or any(
            math.gcd(determinant, shared) != 1
            for determinant in (
                projected_lifted_determinant,
                projected_leaf_determinant,
                lifted_leaf_determinant,
            )
        )
    ):
        raise RuntimeError("rows do not form the required residual triangle")

    normalizer_indices = {
        base_primes.index(prime) for prime in normalization_primes
    }
    normalizers = [
        *[by_prime[prime] for prime in normalization_primes],
        start_leaf,
    ]
    normalization_period = math.lcm(
        *(int(row["h"]) for row in normalizers)
    )
    normalization_image = {
        tuple(
            (int(row["a"]) * k + int(row["b"]) * l) % int(row["h"])
            for row in normalizers
        )
        for k in range(normalization_period)
        for l in range(normalization_period)
    }
    normalization_surjective = len(normalization_image) == math.prod(
        int(row["h"]) for row in normalizers
    )
    if not normalization_surjective:
        raise RuntimeError("normalizer target map is not jointly surjective")
    target_ranges = [
        range(1) if index in normalizer_indices else range(int(row["h"]))
        for index, row in enumerate(bases)
    ]
    target_tuples = math.prod(len(values) for values in target_ranges)
    if target_tuples > args.max_target_tuples:
        raise RuntimeError("base target space exceeds guard")

    def line_masks(row: dict, modulus: int) -> list[int]:
        masks = [0] * modulus
        for k in range(base_period):
            for l in range(base_period):
                value = (
                    int(row["a"]) * k + int(row["b"]) * l
                ) % modulus
                masks[value] |= 1 << (k * base_period + l)
        return masks

    base_masks = [line_masks(row, int(row["h"])) for row in bases]
    projection_rows = [
        projected,
        *independents,
        *(
            [paired_projected]
            if paired_projected is not None
            else []
        ),
    ]
    projection_moduli = (
        projection,
        *independent_projections,
        *([paired_projection] if paired_projection is not None else []),
    )
    projection_masks = [
        line_masks(row, modulus)
        for row, modulus in zip(projection_rows, projection_moduli)
    ]
    leaf_active_mask = line_masks(start_leaf, leaf_projection)[0]
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
    if sum(mask.bit_count() for mask in observed_masks) != base_cells:
        raise RuntimeError("projection categories do not partition base grid")

    path_densities = path_leaf_densities(shared, lifted_residual)
    if (
        path_densities["both_compatible"]
        < path_densities["both_incompatible"]
    ):
        raise RuntimeError("compatible residual triangle is not maximal")
    choices = tuple(
        itertools.product(
            *(range(modulus) for modulus in projection_moduli),
        )
    )
    weights = {}
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
        if main_active and leaf_active:
            path_density = path_densities["both_compatible"]
        elif main_active or leaf_active:
            path_density = path_densities["one"]
        else:
            path_density = path_densities["inactive"]
        uncovered = 1 - path_density
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
                - Fraction(1, paired_residual**2)
                if paired_active
                else Fraction(1, paired_residual)
            )
            uncovered *= 1 - paired_density
        weights[state] = 1 - uncovered
    denominator = math.lcm(
        *(value.denominator for value in weights.values())
    )
    if base_cells * denominator >= int(np.iinfo(np.int64).max):
        raise RuntimeError("exact scores exceed int64")
    score_matrix = np.array(
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
                for code in observed_codes
            ]
            for choice in choices
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
        scores = score_matrix @ histogram
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

    maximum = Fraction(
        maximum_score,
        base_cells * denominator,
    )
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
    result = {
        "schema": (
            "projected_pair_normalized_start_leaf_v3"
            if paired_projected is not None
            else "projected_pair_normalized_start_leaf_v1"
            if len(independent_primes) == 1
            else "projected_pair_normalized_start_leaf_v2"
        ),
        "pool": str(args.pool),
        "row_count": len(rows),
        "anchor_primes": list(anchor_primes),
        "anchor_rows": [
            {key: int(row[key]) for key in ("p", "h", "a", "b")}
            for row in anchors
        ],
        "base_primes": list(base_primes),
        "projected_prime": args.projected_prime,
        "lifted_prime": args.lifted_prime,
        "independent_projected_prime":
        args.independent_projected_prime,
        "independent_projected_primes": list(independent_primes),
        "paired_projected_prime": args.paired_projected_prime,
        "paired_shared_prime": args.paired_shared_prime,
        "normalized_start_shared_prime":
        args.normalized_start_shared_prime,
        "normalization_primes": list(normalization_primes),
        "normalization_period": normalization_period,
        "normalization_jointly_surjective": normalization_surjective,
        "base_period": base_period,
        "base_cells": base_cells,
        "projection_moduli": list(projection_moduli),
        "fixed_leaf_projection_modulus": leaf_projection,
        "fixed_leaf_projection_target": 0,
        "shared_residual_modulus": shared,
        "lifted_residual_modulus": lifted_residual,
        "independent_residual_modulus": independent_residuals[0],
        "independent_residual_moduli": list(independent_residuals),
        "paired_projection_modulus": paired_projection,
        "paired_residual_modulus": paired_residual,
        "paired_determinant": paired_determinant,
        "paired_maps_jointly_surjective": (
            True if paired_projected is not None else None
        ),
        "projected_lifted_determinant":
        projected_lifted_determinant,
        "projected_leaf_determinant": projected_leaf_determinant,
        "lifted_leaf_determinant": lifted_leaf_determinant,
        "all_shared_pair_maps_jointly_surjective": True,
        "path_leaf_densities": {
            key: fraction_payload(value)
            for key, value in path_densities.items()
        },
        "compatible_residual_targets_dominate": True,
        "residual_union_weights": {
            "".join("1" if flag else "0" for flag in state):
            fraction_payload(value)
            for state, value in weights.items()
        },
        "enumerated_period": base_period * math.prod(residuals),
        "target_tuples": target_tuples,
        "projection_target_combinations_per_base": len(choices),
        "target_combinations_checked": checked,
        "maximizing_base_targets": list(maximizing_targets),
        "maximizing_projection_targets": list(maximizing_projections),
        "maximizing_base_covered_cells": maximizing_base_covered,
        "maximizing_uncovered_category_counts": maximizing_counts,
        "maximum_block_union_density": fraction_payload(maximum),
        "block_individual_density_sum": fraction_payload(individual_sum),
        "forced_overlap_loss": fraction_payload(loss),
        "total_pool_density": fraction_payload(total),
        "pool_union_density_upper_bound": fraction_payload(upper),
        "proved_no_cover": upper < 1,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base_targets={target_tuples} projection_targets={len(choices)} "
        f"checked={checked} maximum={maximum} loss={loss} "
        f"pool_upper={upper}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
