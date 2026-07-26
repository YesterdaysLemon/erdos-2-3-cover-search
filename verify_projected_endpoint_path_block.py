#!/usr/bin/env python3
"""Independent transposed-grid replay of an endpoint-path block."""

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


def explicit_endpoint_weights(
    densities: tuple[Fraction, Fraction, Fraction, Fraction, Fraction],
    independent_density: Fraction,
    paired_active_density: Fraction | None = None,
    paired_inactive_density: Fraction | None = None,
    first_shared_pair_density: Fraction | None = None,
    second_independent_density: Fraction | None = None,
) -> dict[tuple[bool, ...], Fraction]:
    paired = (
        paired_active_density is not None
        and paired_inactive_density is not None
    )
    answer = {}
    first_shared = first_shared_pair_density is not None
    second_independent = second_independent_density is not None
    state_count = (
        3
        + int(paired)
        + int(first_shared)
        + int(second_independent)
    )
    for bitmask in range(1 << state_count):
        state = tuple(
            bool(bitmask & (1 << index))
            for index in range(state_count)
        )
        state_index = 3
        paired_active = state[state_index] if paired else False
        state_index += int(paired)
        first_shared_active = (
            state[state_index] if first_shared else False
        )
        state_index += int(first_shared)
        second_independent_active = (
            state[state_index] if second_independent else False
        )
        factors = [densities[2], densities[3]]
        factors.append(densities[1])
        if state[0]:
            factors.append(densities[0])
        if first_shared:
            enabled_count = (
                int(state[1]) + int(first_shared_active)
            )
            if enabled_count == 1:
                factors.append(independent_density)
            elif enabled_count == 2:
                factors.append(first_shared_pair_density)
        elif state[1]:
            factors.append(independent_density)
        if state[2]:
            factors.append(densities[4])
        if paired:
            factors.append(
                paired_active_density
                if paired_active
                else paired_inactive_density
            )
        if second_independent_active:
            factors.append(second_independent_density)
        uncovered = math.prod((1 - value for value in factors), start=Fraction(1))
        answer[state] = 1 - uncovered
    return answer


def brute_equation_intersection_density(
    prime: int,
    equations: tuple[tuple[int, int, int], ...],
) -> Fraction:
    """Count an affine system directly, with coordinates transposed."""
    satisfying = 0
    for second in range(prime):
        for first in range(prime):
            if all(
                (a * first + b * second - target) % prime == 0
                for a, b, target in equations
            ):
                satisfying += 1
    return Fraction(satisfying, prime * prime)


def brute_endpoint_path_with_start_leaf_density(
    path_rows: tuple[dict, dict, dict, dict, dict],
    start_leaf: dict,
    residuals: tuple[int, int, int, int],
    *,
    projected_enabled: bool,
    endpoint_enabled: bool,
    start_leaf_enabled: bool,
    compatible: bool,
) -> Fraction:
    """Replay the six-event residual union by direct component counts."""
    rows = (*path_rows, start_leaf)
    incidence = (
        (0, 1, 5),
        (1, 2),
        (2, 3),
        (3, 4),
    )
    enabled = [1, 2, 3]
    if projected_enabled:
        enabled.append(0)
    if endpoint_enabled:
        enabled.append(4)
    if start_leaf_enabled:
        enabled.append(5)
    enabled.sort()
    union = Fraction(0)
    for bitmask in range(1, 1 << len(enabled)):
        selected = {
            enabled[index]
            for index in range(len(enabled))
            if bitmask & (1 << index)
        }
        intersection = Fraction(1)
        for component_index in range(4):
            equations = []
            for event in incidence[component_index]:
                if event not in selected:
                    continue
                target = (
                    1
                    if (
                        component_index == 0
                        and event == 1
                        and not compatible
                    )
                    else 0
                )
                equations.append(
                    (
                        int(rows[event]["a"]),
                        int(rows[event]["b"]),
                        target,
                    )
                )
            intersection *= brute_equation_intersection_density(
                residuals[component_index],
                tuple(equations),
            )
        union += (
            intersection
            if bitmask.bit_count() % 2
            else -intersection
        )
    return union


def explicit_endpoint_weights_with_start_leaf(
    path_rows: tuple[dict, dict, dict, dict, dict],
    start_leaf: dict,
    residuals: tuple[int, int, int, int],
    independent_residual: int,
    paired_residual: int | None,
    *,
    first_shared: bool,
) -> tuple[
    dict[tuple[bool, ...], Fraction],
    dict[str, dict[str, Fraction]],
]:
    """Independent residual replay for a normalized start-component leaf."""
    tables = {"compatible": {}, "incompatible": {}}
    for label, compatible in (
        ("compatible", True),
        ("incompatible", False),
    ):
        for main_active in (False, True):
            for endpoint_active in (False, True):
                for leaf_active in (False, True):
                    state_key = (
                        f"{int(main_active)}"
                        f"{int(endpoint_active)}"
                        f"{int(leaf_active)}"
                    )
                    tables[label][state_key] = (
                        brute_endpoint_path_with_start_leaf_density(
                            path_rows,
                            start_leaf,
                            residuals,
                            projected_enabled=main_active,
                            endpoint_enabled=endpoint_active,
                            start_leaf_enabled=leaf_active,
                            compatible=compatible,
                        )
                    )
    if any(
        tables["compatible"][state] < tables["incompatible"][state]
        for state in tables["compatible"]
    ):
        raise RuntimeError("compatible leaf replay is not pointwise maximal")

    weights = {}
    for bitmask in range(1 << 6):
        state = tuple(bool(bitmask & (1 << index)) for index in range(6))
        (
            main_active,
            independent_active,
            endpoint_active,
            paired_active,
            first_shared_active,
            leaf_active,
        ) = state
        path_state = (
            f"{int(main_active)}"
            f"{int(endpoint_active)}"
            f"{int(leaf_active)}"
        )
        path_density = tables["compatible"][path_state]
        independent_count = (
            int(independent_active) + int(first_shared_active)
        )
        if first_shared and independent_count == 2:
            independent_density = (
                Fraction(2, independent_residual)
                - Fraction(1, independent_residual * independent_residual)
            )
        elif independent_count:
            independent_density = Fraction(1, independent_residual)
        else:
            independent_density = Fraction(0)
        if paired_residual is None:
            paired_density = Fraction(0)
        elif paired_active:
            paired_density = (
                Fraction(2, paired_residual)
                - Fraction(1, paired_residual * paired_residual)
            )
        else:
            paired_density = Fraction(1, paired_residual)
        weights[state] = 1 - (
            (1 - path_density)
            * (1 - independent_density)
            * (1 - paired_density)
        )
    return weights, tables


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if np is None:
        raise RuntimeError("NumPy is required for exact histogram replay")

    source = json.loads(args.pool.read_text())
    cert = json.loads(args.certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    base_primes = tuple(int(value) for value in cert["base_primes"])
    path_primes = tuple(int(value) for value in cert["path_primes"])
    independent_prime = int(cert["independent_projected_prime"])
    paired_projected_prime = cert.get("paired_projected_prime")
    paired_shared_prime = cert.get("paired_shared_prime")
    first_shared_projected_prime = cert.get(
        "first_shared_projected_prime"
    )
    second_independent_projected_prime = cert.get(
        "second_independent_projected_prime"
    )
    normalized_start_shared_prime = cert.get(
        "normalized_start_shared_prime"
    )
    if (paired_projected_prime is None) != (paired_shared_prime is None):
        raise RuntimeError("certificate has an incomplete paired event")
    all_primes = (
        *base_primes,
        *path_primes[:4],
        independent_prime,
        path_primes[4],
        *(
            (
                int(paired_projected_prime),
                int(paired_shared_prime),
            )
            if paired_projected_prime is not None
            else ()
        ),
        *(
            (int(first_shared_projected_prime),)
            if first_shared_projected_prime is not None
            else ()
        ),
        *(
            (int(second_independent_projected_prime),)
            if second_independent_projected_prime is not None
            else ()
        ),
        *(
            (int(normalized_start_shared_prime),)
            if normalized_start_shared_prime is not None
            else ()
        ),
    )
    missing = set(all_primes) - by_prime.keys()
    if missing:
        raise RuntimeError(f"certificate rows absent: {sorted(missing)}")
    bases = [by_prime[prime] for prime in base_primes]
    path_rows = [by_prime[prime] for prime in path_primes]
    independent = by_prime[independent_prime]
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
    first_shared_projected = (
        by_prime[int(first_shared_projected_prime)]
        if first_shared_projected_prime is not None
        else None
    )
    second_independent = (
        by_prime[int(second_independent_projected_prime)]
        if second_independent_projected_prime is not None
        else None
    )
    normalized_start_shared = (
        by_prime[int(normalized_start_shared_prime)]
        if normalized_start_shared_prime is not None
        else None
    )
    anchors = [
        *bases,
        *path_rows[:4],
        independent,
        path_rows[4],
        *(
            [paired_projected, paired_shared]
            if paired_projected is not None
            else []
        ),
        *(
            [first_shared_projected]
            if first_shared_projected is not None
            else []
        ),
        *(
            [second_independent]
            if second_independent is not None
            else []
        ),
        *(
            [normalized_start_shared]
            if normalized_start_shared is not None
            else []
        ),
    ]
    base_period = math.lcm(*(int(row["h"]) for row in bases))
    base_cells = base_period**2
    projected, lifted, tail, second_tail, endpoint = path_rows
    projection = math.gcd(base_period, int(projected["h"]))
    residual0 = int(projected["h"]) // projection
    residual1 = (
        int(lifted["h"]) // residual0
        if int(lifted["h"]) % residual0 == 0
        else 0
    )
    residual2 = (
        int(tail["h"]) // residual1
        if residual1 and int(tail["h"]) % residual1 == 0
        else 0
    )
    residual3 = (
        int(second_tail["h"]) // residual2
        if residual2 and int(second_tail["h"]) % residual2 == 0
        else 0
    )
    independent_projection = math.gcd(
        base_period,
        int(independent["h"]),
    )
    independent_residual = int(independent["h"]) // independent_projection
    endpoint_projection = math.gcd(base_period, int(endpoint["h"]))
    endpoint_residual = int(endpoint["h"]) // endpoint_projection
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
    first_shared_projection = None
    first_shared_residual = None
    first_shared_determinant = None
    if first_shared_projected is not None:
        first_shared_projection = math.gcd(
            base_period,
            int(first_shared_projected["h"]),
        )
        first_shared_residual = (
            int(first_shared_projected["h"])
            // first_shared_projection
        )
        first_shared_determinant = (
            int(independent["a"]) * int(first_shared_projected["b"])
            - int(independent["b"]) * int(first_shared_projected["a"])
        )
    second_independent_projection = None
    second_independent_residual = None
    if second_independent is not None:
        second_independent_projection = math.gcd(
            base_period,
            int(second_independent["h"]),
        )
        second_independent_residual = (
            int(second_independent["h"])
            // second_independent_projection
        )
    normalized_start_projection = None
    normalized_start_residual = None
    normalized_start_projected_determinant = None
    normalized_start_lifted_determinant = None
    if normalized_start_shared is not None:
        normalized_start_projection = math.gcd(
            base_period,
            int(normalized_start_shared["h"]),
        )
        normalized_start_residual = (
            int(normalized_start_shared["h"])
            // normalized_start_projection
        )
        normalized_start_projected_determinant = (
            int(projected["a"]) * int(normalized_start_shared["b"])
            - int(projected["b"]) * int(normalized_start_shared["a"])
        )
        normalized_start_lifted_determinant = (
            int(lifted["a"]) * int(normalized_start_shared["b"])
            - int(lifted["b"]) * int(normalized_start_shared["a"])
        )
    unique_residuals = (
        residual0,
        residual1,
        residual2,
        residual3,
        independent_residual,
        *([paired_residual] if paired_residual is not None else []),
        *(
            [second_independent_residual]
            if second_independent_residual is not None
            else []
        ),
    )
    shared_components = unique_residuals[:4]
    determinants = tuple(
        int(first["a"]) * int(second["b"])
        - int(first["b"]) * int(second["a"])
        for first, second in zip(path_rows, path_rows[1:])
    )
    path_densities = (
        Fraction(1, residual0),
        Fraction(1, residual0 * residual1),
        Fraction(1, residual1 * residual2),
        Fraction(1, residual2 * residual3),
        Fraction(1, residual3),
    )
    start_leaf_path_tables = None
    if normalized_start_shared is not None:
        start_leaf_core_weights, start_leaf_path_tables = (
            explicit_endpoint_weights_with_start_leaf(
                tuple(path_rows),
                normalized_start_shared,
                tuple(shared_components),
                independent_residual,
                paired_residual,
                first_shared=first_shared_projected is not None,
            )
        )
        weights = {
            (*state[:5], second_active, state[5]): (
                1
                - (1 - value)
                * (
                    1
                    - (
                        Fraction(1, second_independent_residual)
                        if second_active
                        else Fraction(0)
                    )
                )
            )
            for state, value in start_leaf_core_weights.items()
            for second_active in (False, True)
        }
    else:
        weights = explicit_endpoint_weights(
            path_densities,
            Fraction(1, independent_residual),
            (
                Fraction(2, paired_residual)
                - Fraction(1, paired_residual * paired_residual)
                if paired_residual is not None
                else None
            ),
            (
                Fraction(1, paired_residual)
                if paired_residual is not None
                else None
            ),
            (
                Fraction(2, independent_residual)
                - Fraction(
                    1,
                    independent_residual * independent_residual,
                )
                if first_shared_projected is not None
                else None
            ),
            (
                Fraction(1, second_independent_residual)
                if second_independent_residual is not None
                else None
            ),
        )
    projection_moduli = (
        projection,
        independent_projection,
        endpoint_projection,
        *([paired_projection] if paired_projection is not None else []),
        *(
            [first_shared_projection]
            if first_shared_projection is not None
            else []
        ),
        *(
            [second_independent_projection]
            if second_independent_projection is not None
            else []
        ),
    )
    normalization_primes = tuple(
        int(value) for value in cert["normalization_primes"]
    )
    normalized_indices = {
        base_primes.index(prime) for prime in normalization_primes
    }
    normalizers = [
        *[by_prime[prime] for prime in normalization_primes],
        *(
            [normalized_start_shared]
            if normalized_start_shared is not None
            else []
        ),
    ]
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
    target_tuples = math.prod(len(values) for values in target_ranges)

    def transposed_masks(row: dict, modulus: int) -> list[int]:
        masks = [0] * modulus
        for l in range(base_period):
            for k in range(base_period):
                target = (
                    int(row["a"]) * k + int(row["b"]) * l
                ) % modulus
                masks[target] |= 1 << (l * base_period + k)
        return masks

    base_masks = [
        transposed_masks(row, int(row["h"])) for row in bases
    ]
    projection_rows = (
        projected,
        independent,
        endpoint,
        *([paired_projected] if paired_projected is not None else []),
        *(
            [first_shared_projected]
            if first_shared_projected is not None
            else []
        ),
        *(
            [second_independent]
            if second_independent is not None
            else []
        ),
    )
    projection_masks = [
        transposed_masks(row, modulus)
        for row, modulus in zip(projection_rows, projection_moduli)
    ]
    full_grid_mask = (1 << base_cells) - 1
    maximum_score = -1
    maximizing_targets = None
    maximizing_projections = None
    maximizing_base_covered = None
    maximizing_counts = None
    checked = 0
    projection_target_count = math.prod(projection_moduli)
    if second_independent is not None:
        core_projection_moduli = projection_moduli[:-1]
        factorized_projection_moduli = projection_moduli[-1:]
        if normalized_start_shared is not None:
            core_weights = {
                state: weights[(*state[:5], False, state[5])]
                for state in itertools.product(
                    (False, True),
                    repeat=len(core_projection_moduli) + 1,
                )
            }
        else:
            core_weights = {
                state: weights[(*state, False)]
                for state in itertools.product(
                    (False, True),
                    repeat=len(core_projection_moduli),
                )
            }
        extra_weights = {
            (False,): Fraction(0),
            (True,): Fraction(1, second_independent_residual),
        }
        core_denominator = math.lcm(
            *(value.denominator for value in core_weights.values())
        )
        extra_denominator = second_independent_residual
        score_denominator = core_denominator * extra_denominator
        score_ceiling = base_cells * score_denominator
        if score_ceiling >= int(np.iinfo(np.uint64).max):
            raise RuntimeError("factorized verifier scores exceed uint64")
        score_dtype = (
            np.int64
            if score_ceiling < int(np.iinfo(np.int64).max)
            else np.uint64
        )
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
        core_masks = []
        for code in core_codes:
            cells = full_grid_mask
            for masks, target in zip(projection_masks[:-1], code):
                cells &= masks[target]
            core_masks.append(cells)
        if normalized_start_shared is not None:
            normalized_start_active_mask = transposed_masks(
                normalized_start_shared,
                normalized_start_projection,
            )[0]
            core_observed_codes = []
            core_observed_masks = []
            for code, cells in zip(core_codes, core_masks):
                core_observed_codes.extend(
                    ((*code, False), (*code, True))
                )
                core_observed_masks.extend(
                    (
                        cells & (full_grid_mask ^ normalized_start_active_mask),
                        cells & normalized_start_active_mask,
                    )
                )
        else:
            core_observed_codes = list(core_codes)
            core_observed_masks = core_masks
        extra_masks = []
        for code in extra_codes:
            cells = full_grid_mask
            for masks, target in zip(projection_masks[-1:], code):
                cells &= masks[target]
            extra_masks.append(cells)
        if (
            sum(mask.bit_count() for mask in core_observed_masks)
            != base_cells
            or sum(mask.bit_count() for mask in extra_masks) != base_cells
        ):
            raise RuntimeError(
                "verifier factorized cells do not partition the grid"
            )
        core_cell_choice_weights = np.empty(
            (len(core_observed_codes), len(core_codes)),
            dtype=score_dtype,
        )
        for cell_index, code in enumerate(core_observed_codes):
            for choice_index, choice in enumerate(core_codes):
                state = tuple(
                    observed == target
                    for observed, target in zip(
                        code[:len(core_projection_moduli)],
                        choice,
                    )
                )
                if normalized_start_shared is not None:
                    state = (*state, bool(code[-1]))
                value = core_weights[state]
                core_cell_choice_weights[cell_index, choice_index] = (
                    value.numerator
                    * (core_denominator // value.denominator)
                )
        extra_cell_choice_weights = np.empty(
            (len(extra_codes), len(extra_codes)),
            dtype=score_dtype,
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
                    * (extra_denominator // value.denominator)
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
                (len(core_observed_codes), len(extra_codes)),
                dtype=score_dtype,
            )
            for core_index, core_mask in enumerate(core_observed_masks):
                core_uncovered = uncovered & core_mask
                histogram[core_index, :] = [
                    (core_uncovered & extra_mask).bit_count()
                    for extra_mask in extra_masks
                ]
            total_by_core = histogram.sum(axis=1)
            core_scores = total_by_core @ core_cell_choice_weights
            core_complements = (
                core_denominator - core_cell_choice_weights
            )
            extension_scores = (
                (core_complements.T @ histogram)
                @ extra_cell_choice_weights
            )
            scores = (
                core_scores[:, None] * extra_denominator
                + extension_scores
            )
            scores += base_covered * score_denominator
            flat_choice = int(np.argmax(scores))
            core_choice_index, extra_choice_index = np.unravel_index(
                flat_choice,
                scores.shape,
            )
            score = int(scores[core_choice_index, extra_choice_index])
            projections = (
                *core_codes[core_choice_index],
                *extra_codes[extra_choice_index],
            )
            checked += projection_target_count
            if score > maximum_score:
                maximum_score = score
                maximizing_targets = targets
                maximizing_projections = projections
                maximizing_base_covered = base_covered
                maximizing_counts = {
                    "".join("1" if flag else "0" for flag in state): 0
                    for state in itertools.product(
                        (False, True),
                        repeat=(
                            len(projection_moduli)
                            + int(normalized_start_shared is not None)
                        ),
                    )
                }
                for core_index, core_code in enumerate(
                    core_observed_codes
                ):
                    for extra_index, extra_code in enumerate(extra_codes):
                        core_state = tuple(
                            observed == target
                            for observed, target in zip(
                                core_code[
                                    :len(core_projection_moduli)
                                ],
                                projections[
                                    :len(core_projection_moduli)
                                ],
                            )
                        )
                        extra_state = tuple(
                            observed == target
                            for observed, target in zip(
                                extra_code,
                                projections[
                                    len(core_projection_moduli):
                                ],
                            )
                        )
                        state = (*core_state, *extra_state)
                        if normalized_start_shared is not None:
                            state = (*state, bool(core_code[-1]))
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
        denominator = math.lcm(
            *(value.denominator for value in weights.values())
        )
        score_denominator = denominator
        if base_cells * denominator >= int(np.iinfo(np.int64).max):
            raise RuntimeError("verifier scores exceed int64")
        numerators = {
            state: int(value * denominator)
            for state, value in weights.items()
        }
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
                    numerators[state]
                )
        for base_index, targets in enumerate(
            itertools.product(*target_ranges),
            start=1,
        ):
            base_union = 0
            for masks, target in zip(base_masks, targets):
                base_union |= masks[target]
            uncovered = full_grid_mask ^ base_union
            histogram = np.asarray(
                [
                    (uncovered & cells).bit_count()
                    for cells in code_masks
                ],
                dtype=np.int64,
            )
            scores = histogram @ cell_choice_scores
            scores += base_union.bit_count() * denominator
            choice_index = int(np.argmax(scores))
            score = int(scores[choice_index])
            checked += len(projection_codes)
            if score > maximum_score:
                maximum_score = score
                maximizing_targets = targets
                maximizing_projections = projection_codes[choice_index]
                maximizing_base_covered = base_union.bit_count()
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
                        for observed, target in zip(
                            code,
                            maximizing_projections,
                        )
                    )
                    key = "".join(
                        "1" if flag else "0" for flag in state
                    )
                    maximizing_counts[key] += int(count)
            if base_index % 1000 == 0:
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
        len(set(all_primes)) == len(all_primes)
        and all(int(row.get("target_modulus", 1)) == 1 for row in anchors)
        and all(
            math.gcd(int(row["a"]), int(row["b"]), int(row["h"])) == 1
            for row in anchors[len(bases):]
        )
        and endpoint_residual == residual3
        and all(residual > 1 for residual in unique_residuals)
        and all(
            math.gcd(first, second) == 1
            for first, second in itertools.combinations(unique_residuals, 2)
        )
        and all(
            math.gcd(residual, base_period) == 1
            for residual in unique_residuals
        )
        and int(lifted["h"]) == residual0 * residual1
        and int(tail["h"]) == residual1 * residual2
        and int(second_tail["h"]) == residual2 * residual3
        and int(independent["h"])
        == independent_projection * independent_residual
        and int(endpoint["h"]) == endpoint_projection * residual3
        and all(
            math.gcd(determinant, shared) == 1
            for determinant, shared in zip(
                determinants,
                shared_components,
            )
        )
        and (
            paired_projected is None
            or (
                paired_residual > 1
                and int(paired_projected["h"])
                == paired_projection * paired_residual
                and int(paired_shared["h"]) == paired_residual
                and math.gcd(paired_determinant, paired_residual) == 1
            )
        )
        and (
            first_shared_projected is None
            or (
                first_shared_residual == independent_residual
                and int(first_shared_projected["h"])
                == first_shared_projection * first_shared_residual
                and math.gcd(
                    first_shared_determinant,
                    first_shared_residual,
                )
                == 1
            )
        )
        and (
            second_independent is None
            or (
                second_independent_residual > 1
                and int(second_independent["h"])
                == second_independent_projection
                * second_independent_residual
            )
        )
        and (
            normalized_start_shared is None
            or (
                second_independent is not None
                and bool(normalization_primes)
                and normalized_start_projection > 1
                and normalized_start_residual == residual0
                and int(normalized_start_shared["h"])
                == normalized_start_projection * residual0
                and math.gcd(
                    normalized_start_projected_determinant,
                    residual0,
                )
                == 1
                and math.gcd(
                    normalized_start_lifted_determinant,
                    residual0,
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
    recorded_start_leaf_path_tables = (
        {
            compatibility: {
                state: read_fraction(value)
                for state, value in table.items()
            }
            for compatibility, table in cert[
                "normalized_start_path_union_densities"
            ].items()
        }
        if cert.get("normalized_start_path_union_densities") is not None
        else None
    )
    verified = (
        cert.get("schema")
        == (
            "projected_endpoint_path_block_v4"
            if normalized_start_shared is not None
            else (
                "projected_endpoint_path_block_v3"
                if first_shared_projected is not None
                or second_independent is not None
                else (
                    "projected_endpoint_path_block_v2"
                    if paired_projected is not None
                    else "projected_endpoint_path_block_v1"
                )
            )
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
        and list(cert["path_residual_moduli"]) == list(shared_components)
        and int(cert["independent_residual_modulus"])
        == independent_residual
        and cert.get("paired_residual_modulus") == paired_residual
        and cert.get("paired_determinant") == paired_determinant
        and cert.get("paired_maps_jointly_surjective")
        == (True if paired_projected is not None else None)
        and cert.get("first_shared_projection_modulus")
        == first_shared_projection
        and cert.get("first_shared_residual_modulus")
        == first_shared_residual
        and cert.get("first_shared_determinant")
        == first_shared_determinant
        and cert.get("first_shared_maps_jointly_surjective")
        == (True if first_shared_projected is not None else None)
        and cert.get("second_independent_projection_modulus")
        == second_independent_projection
        and cert.get("second_independent_residual_modulus")
        == second_independent_residual
        and cert.get("normalized_start_projection_modulus")
        == normalized_start_projection
        and cert.get("normalized_start_projection_target")
        == (0 if normalized_start_shared is not None else None)
        and cert.get("normalized_start_residual_modulus")
        == normalized_start_residual
        and cert.get("normalized_start_projected_determinant")
        == normalized_start_projected_determinant
        and cert.get("normalized_start_lifted_determinant")
        == normalized_start_lifted_determinant
        and cert.get("normalized_start_pair_maps_jointly_surjective")
        == (True if normalized_start_shared is not None else None)
        and cert.get("normalized_start_compatible_targets_dominate")
        == (True if normalized_start_shared is not None else None)
        and recorded_start_leaf_path_tables == start_leaf_path_tables
        and list(cert["path_determinants"]) == list(determinants)
        and bool(cert["path_adjacent_maps_jointly_surjective"])
        and [
            read_fraction(value) for value in cert["path_event_densities"]
        ]
        == list(path_densities)
        and recorded_weights == expected_weights
        and int(cert["enumerated_period"])
        == base_period * math.prod(unique_residuals)
        and int(cert["target_tuples"]) == target_tuples
        and int(cert["projection_target_combinations_per_base"])
        == projection_target_count
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
        "projection_moduli": list(projection_moduli),
        "residual_moduli": list(unique_residuals),
        "enumerated_period": base_period * math.prod(unique_residuals),
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
        f"base_targets={target_tuples} projection_targets="
        f"{projection_target_count} checked={checked} maximum={maximum} "
        f"loss={loss} pool_upper={upper} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
