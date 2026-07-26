#!/usr/bin/env python3
"""Independent replay for a projected triple-independent block certificate."""

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


def main_pair_densities_explicit(
    shared: int,
    lifted_residual: int,
) -> tuple[Fraction, Fraction]:
    projected = Fraction(1, shared)
    lifted = Fraction(1, shared * lifted_residual)
    intersection = Fraction(1, shared * shared * lifted_residual)
    return projected + lifted - intersection, lifted


def main_chain_densities_explicit(
    shared: int,
    lifted_residual: int,
    tail_residual: int,
) -> tuple[Fraction, Fraction]:
    densities = (
        Fraction(1, shared),
        Fraction(1, shared * lifted_residual),
        Fraction(1, lifted_residual * tail_residual),
    )
    pair_intersections = (
        Fraction(1, shared * shared * lifted_residual),
        Fraction(1, shared * lifted_residual * tail_residual),
        Fraction(
            1,
            shared * lifted_residual * lifted_residual * tail_residual,
        ),
    )
    triple = Fraction(
        1,
        shared * shared * lifted_residual * lifted_residual * tail_residual,
    )
    active = sum(densities, Fraction(0))
    active -= sum(pair_intersections, Fraction(0))
    active += triple
    inactive = densities[1] + densities[2] - pair_intersections[2]
    return active, inactive


def explicit_weights(
    active_main: Fraction,
    inactive_main: Fraction,
    independent_residuals: tuple[int, ...],
    event_active_densities: tuple[Fraction, ...] | None = None,
    event_inactive_densities: tuple[Fraction, ...] | None = None,
) -> dict[tuple[bool, ...], Fraction]:
    active_densities = event_active_densities or tuple(
        Fraction(1, residual) for residual in independent_residuals
    )
    inactive_densities = event_inactive_densities or tuple(
        Fraction(0) for _ in independent_residuals
    )
    answer = {}
    for bitmask in range(1 << (1 + len(independent_residuals))):
        state = tuple(
            bool(bitmask & (1 << index))
            for index in range(1 + len(independent_residuals))
        )
        complement = 1 - (
            active_main if state[0] else inactive_main
        )
        for index in range(len(independent_residuals)):
            density = (
                active_densities[index]
                if state[index + 1]
                else inactive_densities[index]
            )
            complement *= 1 - density
        answer[state] = 1 - complement
    return answer


def explicit_weights_with_first_projected_pair(
    active_main: Fraction,
    inactive_main: Fraction,
    independent_residuals: tuple[int, ...],
    event_active_densities: tuple[Fraction, ...],
    event_inactive_densities: tuple[Fraction, ...],
) -> dict[tuple[bool, ...], Fraction]:
    answer = {}
    state_count = 2 + len(independent_residuals)
    first_residual = independent_residuals[0]
    for bitmask in range(1 << state_count):
        state = tuple(
            bool(bitmask & (1 << index))
            for index in range(state_count)
        )
        complement = 1 - (
            active_main if state[0] else inactive_main
        )
        enabled_count = int(state[1]) + int(state[-1])
        if enabled_count == 1:
            first_density = Fraction(1, first_residual)
        elif enabled_count == 2:
            first_density = (
                Fraction(2, first_residual)
                - Fraction(1, first_residual * first_residual)
            )
        else:
            first_density = Fraction(0)
        complement *= 1 - first_density
        for index in range(1, len(independent_residuals)):
            density = (
                event_active_densities[index]
                if state[index + 1]
                else event_inactive_densities[index]
            )
            complement *= 1 - density
        answer[state] = 1 - complement
    return answer


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    cert = json.loads(args.certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    base_primes = tuple(int(value) for value in cert["base_primes"])
    independent_primes = tuple(
        int(value)
        for value in cert["independent_projected_primes"]
    )
    factorized_independent_indices = tuple(
        int(value)
        for value in cert.get("factorized_independent_indices", [])
    )
    projected_prime = int(cert["projected_prime"])
    lifted_prime = int(cert["lifted_prime"])
    tail_prime = cert.get("tail_prime")
    third_shared_lifted_prime = cert.get(
        "third_shared_lifted_prime"
    )
    first_shared_projected_prime = cert.get(
        "first_shared_projected_prime"
    )
    shared_paired_index = (
        int(cert.get("shared_paired_index", 2))
        if third_shared_lifted_prime is not None
        else None
    )
    all_primes = (
        *base_primes,
        projected_prime,
        lifted_prime,
        *([int(tail_prime)] if tail_prime is not None else []),
        *independent_primes,
        *(
            [int(third_shared_lifted_prime)]
            if third_shared_lifted_prime is not None
            else []
        ),
        *(
            [int(first_shared_projected_prime)]
            if first_shared_projected_prime is not None
            else []
        ),
    )
    missing = set(all_primes) - by_prime.keys()
    if missing:
        raise RuntimeError(f"certificate rows absent: {sorted(missing)}")
    bases = [by_prime[prime] for prime in base_primes]
    projected = by_prime[projected_prime]
    lifted = by_prime[lifted_prime]
    tail = by_prime[int(tail_prime)] if tail_prime is not None else None
    independents = [by_prime[prime] for prime in independent_primes]
    third_shared_lifted = (
        by_prime[int(third_shared_lifted_prime)]
        if third_shared_lifted_prime is not None
        else None
    )
    first_shared_projected = (
        by_prime[int(first_shared_projected_prime)]
        if first_shared_projected_prime is not None
        else None
    )
    if (
        third_shared_lifted is not None
        and not 0 <= shared_paired_index < len(independents)
    ):
        raise RuntimeError(
            "shared lifted row has an invalid paired projected index"
        )
    if first_shared_projected is not None and shared_paired_index == 0:
        raise RuntimeError(
            "the first shared projection conflicts with paired index zero"
        )
    anchors = [
        *bases,
        projected,
        lifted,
        *([tail] if tail else []),
        *independents,
        *([third_shared_lifted] if third_shared_lifted else []),
        *([first_shared_projected] if first_shared_projected else []),
    ]

    base_period = math.lcm(*(int(row["h"]) for row in bases))
    base_cells = base_period**2
    projection = math.gcd(base_period, int(projected["h"]))
    shared = int(projected["h"]) // projection
    lifted_residual = (
        int(lifted["h"]) // shared
        if int(lifted["h"]) % shared == 0
        else 0
    )
    determinant = (
        int(projected["a"]) * int(lifted["b"])
        - int(projected["b"]) * int(lifted["a"])
    )
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
    projection_moduli = (projection, *independent_projections)
    tail_residual = None
    tail_determinant = None
    if tail is not None:
        tail_residual = (
            int(tail["h"]) // lifted_residual
            if int(tail["h"]) % lifted_residual == 0
            else 0
        )
        tail_determinant = (
            int(lifted["a"]) * int(tail["b"])
            - int(lifted["b"]) * int(tail["a"])
        )
    residuals = (
        shared,
        lifted_residual,
        *([tail_residual] if tail_residual is not None else []),
        *independent_residuals,
    )
    if tail_residual is None:
        active_main, inactive_main = main_pair_densities_explicit(
            shared,
            lifted_residual,
        )
    else:
        active_main, inactive_main = main_chain_densities_explicit(
            shared,
            lifted_residual,
            tail_residual,
        )
    event_active_densities = tuple(
        Fraction(1, residual) for residual in independent_residuals
    )
    event_inactive_densities = tuple(
        Fraction(0) for _ in independent_residuals
    )
    third_shared_determinant = None
    if third_shared_lifted is not None:
        third_residual = independent_residuals[shared_paired_index]
        third_shared_determinant = (
            int(independents[shared_paired_index]["a"])
            * int(third_shared_lifted["b"])
            - int(independents[shared_paired_index]["b"])
            * int(third_shared_lifted["a"])
        )
        active_densities = list(event_active_densities)
        inactive_densities = list(event_inactive_densities)
        active_densities[shared_paired_index] = (
            Fraction(2, third_residual)
            - Fraction(1, third_residual * third_residual)
        )
        inactive_densities[shared_paired_index] = Fraction(
            1,
            third_residual,
        )
        event_active_densities = tuple(active_densities)
        event_inactive_densities = tuple(inactive_densities)
    first_shared_projection = None
    first_shared_residual = None
    first_shared_determinant = None
    if first_shared_projected is not None:
        first_shared_projection = math.gcd(
            base_period,
            int(first_shared_projected["h"]),
        )
        first_shared_residual = (
            int(first_shared_projected["h"]) // first_shared_projection
        )
        first_shared_determinant = (
            int(independents[0]["a"])
            * int(first_shared_projected["b"])
            - int(independents[0]["b"])
            * int(first_shared_projected["a"])
        )
        projection_moduli = (
            *projection_moduli,
            first_shared_projection,
        )
        weights = explicit_weights_with_first_projected_pair(
            active_main,
            inactive_main,
            independent_residuals,
            event_active_densities,
            event_inactive_densities,
        )
    else:
        weights = explicit_weights(
            active_main,
            inactive_main,
            independent_residuals,
            event_active_densities,
            event_inactive_densities,
        )
    if np is None:
        raise RuntimeError("NumPy is required for exact histogram replay")

    normalization_primes = tuple(
        int(value) for value in cert["normalization_primes"]
    )
    normalized_indices = {
        base_primes.index(prime) for prime in normalization_primes
    }
    normalizers = [by_prime[prime] for prime in normalization_primes]
    normalization_period = (
        math.lcm(*(int(row["h"]) for row in normalizers))
        if normalizers
        else 1
    )
    normalization_image = {
        tuple(
            (
                int(row["a"]) * k + int(row["b"]) * l
            )
            % int(row["h"])
            for row in normalizers
        )
        for k in range(normalization_period)
        for l in range(normalization_period)
    }
    normalization_surjective = (
        not normalizers
        or len(normalization_image)
        == math.prod(int(row["h"]) for row in normalizers)
    )
    target_ranges = [
        range(1) if index in normalized_indices else range(int(row["h"]))
        for index, row in enumerate(bases)
    ]

    # This implementation deliberately transposes the bit index relative to
    # the certifier and reconstructs each exact category directly.
    def transposed_line_masks(row: dict, modulus: int) -> list[int]:
        masks = [0] * modulus
        for l in range(base_period):
            for k in range(base_period):
                target = (
                    int(row["a"]) * k + int(row["b"]) * l
                ) % modulus
                masks[target] |= 1 << (l * base_period + k)
        return masks

    base_masks = [
        transposed_line_masks(row, int(row["h"])) for row in bases
    ]
    full_grid_mask = (1 << base_cells) - 1
    projection_rows = [
        projected,
        *independents,
        *(
            [first_shared_projected]
            if first_shared_projected is not None
            else []
        ),
    ]
    projection_masks = [
        transposed_line_masks(row, modulus)
        for row, modulus in zip(projection_rows, projection_moduli)
    ]
    maximum_score = -1
    maximizing_targets = None
    maximizing_projections = None
    maximizing_base_covered = None
    maximizing_counts = None
    checked = 0
    target_tuples = math.prod(len(values) for values in target_ranges)
    factorized_projection_axes = tuple(
        1 + index for index in factorized_independent_indices
    )
    core_projection_axes = tuple(
        index
        for index in range(len(projection_moduli))
        if index not in factorized_projection_axes
    )
    core_projection_moduli = tuple(
        projection_moduli[index] for index in core_projection_axes
    )
    factorized_projection_moduli = tuple(
        projection_moduli[index] for index in factorized_projection_axes
    )
    if factorized_independent_indices:
        core_weights = {}
        for core_state in itertools.product(
            (False, True),
            repeat=len(core_projection_axes),
        ):
            full_state = [False] * len(projection_moduli)
            for axis, enabled in zip(core_projection_axes, core_state):
                full_state[axis] = enabled
            core_weights[core_state] = weights[tuple(full_state)]
        factorized_residuals = tuple(
            independent_residuals[index]
            for index in factorized_independent_indices
        )
        extra_weights = {}
        for extra_state in itertools.product(
            (False, True),
            repeat=len(factorized_projection_axes),
        ):
            complement = Fraction(1)
            for enabled, residual in zip(
                extra_state,
                factorized_residuals,
            ):
                if enabled:
                    complement *= Fraction(residual - 1, residual)
            extra_weights[extra_state] = 1 - complement
        core_weight_denominator = math.lcm(
            *(value.denominator for value in core_weights.values())
        )
        extra_weight_denominator = math.lcm(
            *(value.denominator for value in extra_weights.values())
        )
        score_denominator = (
            core_weight_denominator * extra_weight_denominator
        )
        if (
            base_cells * score_denominator
            >= int(np.iinfo(np.int64).max)
        ):
            raise RuntimeError("factorized verifier score exceeds int64")
        core_codes = tuple(
            itertools.product(
                *(range(modulus) for modulus in core_projection_moduli)
            )
        )
        extra_codes = tuple(
            itertools.product(
                *(
                    range(modulus)
                    for modulus in factorized_projection_moduli
                )
            )
        )

        def transposed_masks_for_axes(
            codes: tuple[tuple[int, ...], ...],
            axes: tuple[int, ...],
        ) -> list[int]:
            answer = []
            for code in codes:
                cells = full_grid_mask
                for axis, target in zip(axes, code):
                    cells &= projection_masks[axis][target]
                answer.append(cells)
            return answer

        core_code_masks = transposed_masks_for_axes(
            core_codes,
            core_projection_axes,
        )
        extra_code_masks = transposed_masks_for_axes(
            extra_codes,
            factorized_projection_axes,
        )
        if (
            sum(mask.bit_count() for mask in core_code_masks) != base_cells
            or sum(mask.bit_count() for mask in extra_code_masks)
            != base_cells
        ):
            raise RuntimeError(
                "verifier factorized values do not partition the grid"
            )
        core_cell_choice_weights = np.empty(
            (len(core_codes), len(core_codes)),
            dtype=np.int64,
        )
        for cell_index, code in enumerate(core_codes):
            for choice_index, choice in enumerate(core_codes):
                state = tuple(
                    observed == target
                    for observed, target in zip(code, choice)
                )
                value = core_weights[state]
                core_cell_choice_weights[cell_index, choice_index] = (
                    value.numerator
                    * (core_weight_denominator // value.denominator)
                )
        extra_cell_choice_weights = np.empty(
            (len(extra_codes), len(extra_codes)),
            dtype=np.int64,
        )
        for cell_index, code in enumerate(extra_codes):
            for choice_index, choice in enumerate(extra_codes):
                state = tuple(
                    observed == target
                    for observed, target in zip(code, choice)
                )
                value = extra_weights[state]
                extra_cell_choice_weights[cell_index, choice_index] = (
                    value.numerator
                    * (extra_weight_denominator // value.denominator)
                )
        for base_index, targets in enumerate(
            itertools.product(*target_ranges),
            start=1,
        ):
            base_union = 0
            for masks, target in zip(base_masks, targets):
                base_union |= masks[target]
            base_covered = base_union.bit_count()
            uncovered = full_grid_mask ^ base_union
            histogram = np.empty(
                (len(core_codes), len(extra_codes)),
                dtype=np.int64,
            )
            for core_index, core_mask in enumerate(core_code_masks):
                core_uncovered = uncovered & core_mask
                histogram[core_index, :] = [
                    (core_uncovered & extra_mask).bit_count()
                    for extra_mask in extra_code_masks
                ]
            total_by_core = histogram.sum(axis=1)
            core_scores = (
                total_by_core @ core_cell_choice_weights
            )
            core_complements = (
                core_weight_denominator - core_cell_choice_weights
            )
            extension_scores = (
                (core_complements.T @ histogram)
                @ extra_cell_choice_weights
            )
            scores = (
                core_scores[:, None] * extra_weight_denominator
                + extension_scores
            )
            scores += base_covered * score_denominator
            flat_choice = int(np.argmax(scores))
            core_choice_index, extra_choice_index = np.unravel_index(
                flat_choice,
                scores.shape,
            )
            score = int(scores[core_choice_index, extra_choice_index])
            projections_list = [0] * len(projection_moduli)
            for axis, target in zip(
                core_projection_axes,
                core_codes[core_choice_index],
            ):
                projections_list[axis] = target
            for axis, target in zip(
                factorized_projection_axes,
                extra_codes[extra_choice_index],
            ):
                projections_list[axis] = target
            projections = tuple(projections_list)
            checked += len(core_codes) * len(extra_codes)
            if score > maximum_score:
                maximum_score = score
                maximizing_targets = targets
                maximizing_projections = projections
                maximizing_base_covered = base_covered
                maximizing_counts = {
                    "".join("1" if flag else "0" for flag in state): 0
                    for state in itertools.product(
                        (False, True),
                        repeat=len(projection_moduli),
                    )
                }
                for core_index, core_code in enumerate(core_codes):
                    for extra_index, extra_code in enumerate(extra_codes):
                        full_code = [0] * len(projection_moduli)
                        for axis, observed in zip(
                            core_projection_axes,
                            core_code,
                        ):
                            full_code[axis] = observed
                        for axis, observed in zip(
                            factorized_projection_axes,
                            extra_code,
                        ):
                            full_code[axis] = observed
                        state = tuple(
                            observed == target
                            for observed, target in zip(
                                full_code,
                                projections,
                            )
                        )
                        key = "".join(
                            "1" if flag else "0" for flag in state
                        )
                        maximizing_counts[key] += int(
                            histogram[core_index, extra_index]
                        )
            if base_index % 500 == 0:
                print(
                    f"verified_base_assignments={base_index}/"
                    f"{target_tuples} projection_combinations={checked} "
                    f"best_score={maximum_score}",
                    flush=True,
                )
    else:
        weight_denominator = math.lcm(
            *(value.denominator for value in weights.values())
        )
        score_denominator = weight_denominator
        weight_numerators = {
            state: value.numerator
            * (weight_denominator // value.denominator)
            for state, value in weights.items()
        }
        if (
            base_cells * weight_denominator
            >= int(np.iinfo(np.int64).max)
        ):
            raise RuntimeError("exact verifier score exceeds int64")
        projection_codes = tuple(
            itertools.product(
                *(range(modulus) for modulus in projection_moduli)
            )
        )
        code_masks = []
        for code in projection_codes:
            cells = full_grid_mask
            for index, target in enumerate(code):
                cells &= projection_masks[index][target]
            code_masks.append(cells)
        if sum(mask.bit_count() for mask in code_masks) != base_cells:
            raise RuntimeError(
                "verifier projection cells do not partition grid"
            )
        cell_choice_scores = np.empty(
            (len(projection_codes), len(projection_codes)),
            dtype=np.int64,
        )
        for cell_index, code in enumerate(projection_codes):
            for choice_index, choice in enumerate(projection_codes):
                state = tuple(
                    observed == target
                    for observed, target in zip(code, choice)
                )
                cell_choice_scores[cell_index, choice_index] = (
                    weight_numerators[state]
                )
        for base_index, targets in enumerate(
            itertools.product(*target_ranges),
            start=1,
        ):
            base_union = 0
            for masks, target in zip(base_masks, targets):
                base_union |= masks[target]
            base_covered = base_union.bit_count()
            uncovered = full_grid_mask ^ base_union
            base_score = base_covered * weight_denominator
            histogram = np.asarray(
                [
                    (uncovered & mask).bit_count()
                    for mask in code_masks
                ],
                dtype=np.int64,
            )
            scores = histogram @ cell_choice_scores
            scores += base_score
            choice_index = int(np.argmax(scores))
            score = int(scores[choice_index])
            projections = projection_codes[choice_index]
            checked += len(projection_codes)
            if score > maximum_score:
                maximum_score = score
                maximizing_targets = targets
                maximizing_projections = projections
                maximizing_base_covered = base_covered
                maximizing_counts = {
                    "".join("1" if flag else "0" for flag in state): 0
                    for state in itertools.product(
                        (False, True),
                        repeat=len(projection_moduli),
                    )
                }
                for code, count in zip(
                    projection_codes,
                    histogram.tolist(),
                ):
                    state = tuple(
                        observed == target
                        for observed, target in zip(code, projections)
                    )
                    key = "".join(
                        "1" if flag else "0" for flag in state
                    )
                    maximizing_counts[key] += int(count)
            if base_index % 500 == 0:
                print(
                    f"verified_base_assignments={base_index}/"
                    f"{target_tuples} projection_combinations={checked} "
                    f"best_score={maximum_score}",
                    flush=True,
                )

    maximum = Fraction(
        maximum_score,
        base_cells * score_denominator,
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
    anchors_match = (
        list(cert["anchor_primes"]) == list(all_primes)
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
    structural = (
        len(independents) in (2, 3, 4)
        and len(set(all_primes)) == len(all_primes)
        and len(set(factorized_independent_indices))
        == len(factorized_independent_indices)
        and all(
            0 <= index < len(independents)
            for index in factorized_independent_indices
        )
        and (
            third_shared_lifted is None
            or shared_paired_index
            not in factorized_independent_indices
        )
        and (
            first_shared_projected is None
            or 0 not in factorized_independent_indices
        )
        and all(int(row.get("target_modulus", 1)) == 1 for row in anchors)
        and all(
            math.gcd(int(row["a"]), int(row["b"]), int(row["h"])) == 1
            for row in [
                projected,
                lifted,
                *([tail] if tail else []),
                *independents,
                *([third_shared_lifted] if third_shared_lifted else []),
                *(
                    [first_shared_projected]
                    if first_shared_projected
                    else []
                ),
            ]
        )
        and shared > 1
        and lifted_residual > 1
        and math.gcd(shared, base_period) == 1
        and math.gcd(int(lifted["h"]), base_period) == 1
        and math.gcd(determinant, shared) == 1
        and all(
            residual > 1
            and math.gcd(residual, base_period) == 1
            and int(row["h"]) == projection_value * residual
            for row, projection_value, residual in zip(
                independents,
                independent_projections,
                independent_residuals,
            )
        )
        and all(
            math.gcd(first, second) == 1
            for first, second in itertools.combinations(residuals, 2)
        )
        and (
            tail is None
            or (
                tail_residual > 1
                and int(tail["h"]) == lifted_residual * tail_residual
                and math.gcd(
                    int(tail["h"]),
                    base_period * shared,
                )
                == 1
                and math.gcd(lifted_residual, tail_residual) == 1
                and math.gcd(tail_determinant, lifted_residual) == 1
            )
        )
        and (
            third_shared_lifted is None
            or (
                int(third_shared_lifted["h"])
                == independent_residuals[shared_paired_index]
                and math.gcd(
                    third_shared_determinant,
                    independent_residuals[shared_paired_index],
                )
                == 1
            )
        )
        and (
            first_shared_projected is None
            or (
                first_shared_residual == independent_residuals[0]
                and int(first_shared_projected["h"])
                == first_shared_projection * first_shared_residual
                and math.gcd(
                    first_shared_determinant,
                    first_shared_residual,
                )
                == 1
            )
        )
    )
    recorded_weights = {
        key: read_fraction(value)
        for key, value in cert["residual_union_weights"].items()
    }
    expected_weights = {
        "".join("1" if flag else "0" for flag in state): value
        for state, value in weights.items()
    }
    expected_schema = (
        "projected_triple_event_block_v4"
        if factorized_independent_indices
        else (
            "projected_triple_event_block_v3"
            if first_shared_projected is not None
            else (
                "projected_triple_event_block_v2"
                if tail is not None or third_shared_lifted is not None
                else "projected_triple_independent_block_v1"
            )
        )
    )
    recorded_active_densities = cert.get(
        "event_active_residual_densities"
    )
    recorded_inactive_densities = cert.get(
        "event_inactive_residual_densities"
    )
    event_densities_match = (
        recorded_active_densities is None
        or [
            read_fraction(value) for value in recorded_active_densities
        ]
        == list(event_active_densities)
    ) and (
        recorded_inactive_densities is None
        or [
            read_fraction(value) for value in recorded_inactive_densities
        ]
        == list(event_inactive_densities)
    )
    verified = (
        cert.get("schema") == expected_schema
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
        and list(cert.get("factorized_independent_indices", []))
        == list(factorized_independent_indices)
        and list(
            cert.get(
                "core_projection_axes",
                list(range(len(projection_moduli))),
            )
        )
        == list(core_projection_axes)
        and list(cert.get("factorized_projection_axes", []))
        == list(factorized_projection_axes)
        and list(
            cert.get(
                "core_projection_moduli",
                list(projection_moduli),
            )
        )
        == list(core_projection_moduli)
        and list(cert.get("factorized_projection_moduli", []))
        == list(factorized_projection_moduli)
        and int(cert["shared_residual_modulus"]) == shared
        and int(cert["lifted_residual_modulus"]) == lifted_residual
        and cert.get("tail_residual_modulus") == tail_residual
        and cert.get("shared_paired_index", shared_paired_index)
        == shared_paired_index
        and list(cert["independent_residual_moduli"])
        == list(independent_residuals)
        and int(cert["shared_determinant"]) == determinant
        and bool(cert["shared_pair_jointly_surjective"])
        and cert.get("tail_determinant") == tail_determinant
        and bool(cert.get("tail_pair_jointly_surjective"))
        == (tail is not None)
        and cert.get("third_shared_determinant")
        == third_shared_determinant
        and bool(cert.get("third_shared_pair_jointly_surjective"))
        == (third_shared_lifted is not None)
        and cert.get("first_shared_projection_modulus")
        == first_shared_projection
        and cert.get("first_shared_residual_modulus")
        == first_shared_residual
        and cert.get("first_shared_determinant")
        == first_shared_determinant
        and bool(cert.get("first_shared_pair_jointly_surjective"))
        == (first_shared_projected is not None)
        and event_densities_match
        and read_fraction(cert["active_main_residual_union_density"])
        == active_main
        and read_fraction(cert["inactive_main_residual_union_density"])
        == inactive_main
        and recorded_weights == expected_weights
        and int(cert["enumerated_period"])
        == base_period * math.prod(residuals)
        and int(cert["target_tuples"]) == target_tuples
        and int(cert["projection_target_combinations_per_base"])
        == math.prod(projection_moduli)
        and int(cert["target_combinations_checked"]) == checked
        and list(maximizing_targets) == cert["maximizing_base_targets"]
        and list(maximizing_projections)
        == cert["maximizing_projection_targets"]
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
        "anchor_primes": list(all_primes),
        "base_period": base_period,
        "base_cells": base_cells,
        "projection_moduli": list(projection_moduli),
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
        f"base_targets={target_tuples} "
        f"projection_targets={math.prod(projection_moduli)} checked={checked} "
        f"maximum={maximum} loss={loss} pool_upper={upper} "
        f"verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
