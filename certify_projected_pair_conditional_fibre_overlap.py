#!/usr/bin/env python3
"""Certify a conditional fibre intersection with a structured anchor block.

The anchor block is one previously certified projected-pair/start-leaf block.
The outside modulus must divide the block's base period.  Translation
normalizes the outside target, one or more named base-anchor targets, and the
start-leaf target to zero.  The remaining base and projection targets are
enumerated exactly.  On the shared residual component the incompatible target
relation is the pointwise minimum, so the resulting score is the exact minimum
intersection of the outside fibre with the anchor union.
"""

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

from certify_projected_pair_start_leaf_block import (
    fraction_payload,
    path_leaf_densities,
)


def parse_primes(value: str) -> tuple[int, ...]:
    return tuple(int(item) for item in value.split(",") if item)


def recorded_row(row: dict) -> dict:
    return {key: int(row[key]) for key in ("p", "h", "a", "b")}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("block_certificate", type=Path)
    parser.add_argument("--outside-prime", type=int, required=True)
    parser.add_argument(
        "--normalize-base-primes",
        required=True,
        help="comma-separated base anchors normalized with the outside/leaf",
    )
    parser.add_argument(
        "--omit-start-leaf",
        action="store_true",
        help="certify intersection with the block subset excluding the leaf",
    )
    parser.add_argument("--max-base-cells", type=int, default=1_000_000)
    parser.add_argument("--max-target-tuples", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if np is None:
        raise RuntimeError("NumPy is required for exact score replay")

    source = json.loads(args.pool.read_text())
    block = json.loads(args.block_certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    if len(by_prime) != len(rows):
        raise RuntimeError("pool contains repeated fibre primes")
    structured = block
    structured_path = args.block_certificate
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
    if block_anchor_primes != (
        *structured_anchor_primes,
        *extra_independent_primes,
    ):
        raise RuntimeError("factorized extension does not append its extra rows")
    independent_primes = (
        *core_independent_primes,
        *extra_independent_primes,
    )
    include_leaf = not args.omit_start_leaf
    anchor_primes = tuple(
        prime
        for prime in block_anchor_primes
        if include_leaf or prime != leaf_prime
    )
    listed = {args.outside_prime, *block_anchor_primes}
    if (
        len(listed) != len(block_anchor_primes) + 1
        or not listed <= by_prime.keys()
    ):
        raise RuntimeError("outside and anchor primes must be present/distinct")

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
    outside = by_prime[args.outside_prime]
    anchors = [by_prime[prime] for prime in anchor_primes]
    block_rows = {
        int(row["p"]): recorded_row(row)
        for row in block["anchor_rows"]
    }
    if (
        set(block_rows) != set(block_anchor_primes)
        or any(
            block_rows[prime] != recorded_row(by_prime[prime])
            for prime in block_anchor_primes
        )
    ):
        raise RuntimeError("block anchors differ from the pool")
    if any(int(row.get("target_modulus", 1)) != 1 for row in anchors):
        raise RuntimeError("anchor targets must be unrestricted")

    base_period = math.lcm(*(int(row["h"]) for row in bases))
    base_cells = base_period**2
    if base_cells > args.max_base_cells:
        raise RuntimeError("base grid exceeds guard")
    if base_period % int(outside["h"]):
        raise RuntimeError("outside modulus does not divide the base period")

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
        min(
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
    if not structural:
        raise RuntimeError("block rows fail the projected-pair structure")

    normalization_primes = parse_primes(args.normalize_base_primes)
    if (
        not normalization_primes
        or len(set(normalization_primes)) != len(normalization_primes)
        or not set(normalization_primes) <= set(base_primes)
    ):
        raise RuntimeError("invalid base normalization primes")
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
    normalization_image = {
        tuple(
            (int(row["a"]) * k + int(row["b"]) * l) % int(row["h"])
            for row in normalizers
        )
        for k in range(normalization_period)
        for l in range(normalization_period)
    }
    normalization_target_count = math.prod(
        int(row["h"]) for row in normalizers
    )
    if len(normalization_image) != normalization_target_count:
        raise RuntimeError("outside/anchor normalization is not surjective")

    target_ranges = [
        range(1) if index in normalized_indices else range(int(row["h"]))
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
    outside_mask = line_masks(outside, int(outside["h"]))[0]
    outside_base_points = outside_mask.bit_count()
    if outside_base_points != base_cells // int(outside["h"]):
        raise RuntimeError("outside fibre has the wrong base-grid size")

    projection_rows = [
        projected,
        *independents,
        *((paired_projected,) if has_paired_projection else ()),
    ]
    projection_moduli = (
        projection,
        *independent_projections,
        *((paired_projection,) if has_paired_projection else ()),
    )
    projection_masks = [
        line_masks(row, modulus)
        for row, modulus in zip(projection_rows, projection_moduli)
    ]
    leaf_active_mask = (
        line_masks(leaf, leaf_projection)[0]
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
    if sum(mask.bit_count() for mask in observed_masks) != base_cells:
        raise RuntimeError("projection categories do not partition base grid")

    path_densities = path_leaf_densities(shared, lifted_residual)
    if (
        path_densities["both_incompatible"]
        > path_densities["both_compatible"]
    ):
        raise RuntimeError("incompatible shared targets are not minimal")
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
            path_density = path_densities["both_incompatible"]
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
                uncovered *= Fraction(residual_value - 1, residual_value)
        if has_paired_projection:
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

    choices = tuple(
        itertools.product(
            *(range(modulus) for modulus in projection_moduli)
        )
    )
    score_matrix = np.array(
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
                for code in observed_codes
            ]
            for choice in choices
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
        scores = score_matrix @ histogram
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
    result = {
        "schema": "projected_pair_conditional_fibre_overlap_v1",
        "pool": str(args.pool),
        "block_certificate": str(args.block_certificate),
        "structured_base_certificate": str(structured_path),
        "structured_base_schema": structured_schema,
        "row_count": len(rows),
        "outside_prime": args.outside_prime,
        "outside_row": recorded_row(outside),
        "anchor_primes": list(anchor_primes),
        "anchor_rows": [recorded_row(row) for row in anchors],
        "base_primes": list(base_primes),
        "projected_prime": projected_prime,
        "lifted_prime": lifted_prime,
        "independent_projected_primes": list(independent_primes),
        "core_independent_projected_primes":
        list(core_independent_primes),
        "extra_independent_primes": list(extra_independent_primes),
        "paired_projected_prime": paired_projected_prime,
        "paired_shared_prime": paired_shared_prime,
        "normalized_start_shared_prime": leaf_prime,
        "start_leaf_included": include_leaf,
        "normalization_primes": list(normalization_primes),
        "normalization_period": normalization_period,
        "normalization_target_count": normalization_target_count,
        "normalization_jointly_surjective": True,
        "base_period": base_period,
        "base_cells": base_cells,
        "outside_base_points": outside_base_points,
        "projection_moduli": list(projection_moduli),
        "fixed_leaf_projection_modulus": leaf_projection,
        "fixed_leaf_projection_target": 0 if include_leaf else None,
        "residual_moduli": list(residuals),
        "residual_target_relation": (
            "incompatible" if include_leaf else "leaf_omitted"
        ),
        "path_leaf_densities": {
            key: fraction_payload(value)
            for key, value in path_densities.items()
        },
        "residual_union_weights": {
            "".join("1" if flag else "0" for flag in state):
            fraction_payload(value)
            for state, value in weights.items()
        },
        "target_tuples": target_tuples,
        "projection_target_combinations_per_base": len(choices),
        "target_combinations_checked": checked,
        "minimizing_base_targets": list(minimizing_targets),
        "minimizing_projection_targets": list(minimizing_projections),
        "minimizing_base_covered_cells": minimizing_base_covered,
        "minimizing_uncovered_category_counts": minimizing_counts,
        "minimum_intersection_score": int(minimum_score),
        "score_denominator": denominator,
        "forced_intersection_density": fraction_payload(intersection),
        "argument": (
            "Translate the outside fibre, the named base anchors, and the "
            "start leaf to target zero; the recorded joint image is "
            "surjective. Restrict the exact base-plane partition to the "
            "outside fibre and enumerate every remaining anchor target. "
            "Residual weights count the exact projected union. Incompatible "
            "targets minimize the shared residual triangle pointwise. The "
            "smallest score is therefore the exact minimum intersection of "
            "the outside fibre with the full anchor-block union."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"outside={args.outside_prime} anchors={len(anchors)} "
        f"base_targets={target_tuples} projections={len(choices)} "
        f"checked={checked} intersection={intersection}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
