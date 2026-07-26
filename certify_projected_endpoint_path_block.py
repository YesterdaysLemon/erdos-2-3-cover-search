#!/usr/bin/env python3
"""Certify a projected path with conditional events at both endpoints.

The residual event graph is

  projected -- lifted -- tail -- second_tail -- endpoint

on four pairwise-coprime CRT components.  Adjacent row maps must be jointly
surjective on their shared component.  The five event indicators are then
mutually independent for every enabled subset.  The first and last events
are enabled by separate base-grid projections.  One further projected row
uses a new independent residual component.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

from certify_projected_triple_independent_block import (
    factorized_projection_scores,
)

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


def endpoint_path_weights(
    path_densities: tuple[Fraction, Fraction, Fraction, Fraction, Fraction],
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
    if (paired_active_density is None) != (paired_inactive_density is None):
        raise ValueError("paired densities must be supplied together")
    first_shared = first_shared_pair_density is not None
    second_independent = second_independent_density is not None
    answer = {}
    for state in itertools.product(
        (False, True),
        repeat=(
            3
            + int(paired)
            + int(first_shared)
            + int(second_independent)
        ),
    ):
        main_active, independent_active, endpoint_active = state[:3]
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
        enabled = (
            main_active,
            True,
            True,
            True,
            endpoint_active,
        )
        uncovered = Fraction(1)
        for active, density in zip(enabled, path_densities):
            if active:
                uncovered *= 1 - density
        if first_shared:
            enabled_count = (
                int(independent_active) + int(first_shared_active)
            )
            if enabled_count == 1:
                uncovered *= 1 - independent_density
            elif enabled_count == 2:
                uncovered *= 1 - first_shared_pair_density
        elif independent_active:
            uncovered *= 1 - independent_density
        if paired:
            paired_density = (
                paired_active_density
                if paired_active
                else paired_inactive_density
            )
            uncovered *= 1 - paired_density
        if second_independent_active:
            uncovered *= 1 - second_independent_density
        answer[state] = 1 - uncovered
    return answer


def linear_equation_intersection_density(
    prime: int,
    equations: tuple[tuple[int, int, int], ...],
) -> Fraction:
    """Solve affine equations in two variables over a prime field."""
    if not equations:
        return Fraction(1)
    matrix = [
        [a % prime, b % prime, target % prime]
        for a, b, target in equations
    ]
    pivot_row = 0
    for column in range(2):
        pivot = next(
            (
                row
                for row in range(pivot_row, len(matrix))
                if matrix[row][column] % prime
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = (
            matrix[pivot],
            matrix[pivot_row],
        )
        inverse = pow(matrix[pivot_row][column], -1, prime)
        matrix[pivot_row] = [
            value * inverse % prime for value in matrix[pivot_row]
        ]
        for row in range(len(matrix)):
            if row == pivot_row:
                continue
            multiplier = matrix[row][column] % prime
            if multiplier:
                matrix[row] = [
                    (value - multiplier * pivot_value) % prime
                    for value, pivot_value in zip(
                        matrix[row],
                        matrix[pivot_row],
                    )
                ]
        pivot_row += 1
    if any(
        row[0] % prime == 0
        and row[1] % prime == 0
        and row[2] % prime != 0
        for row in matrix
    ):
        return Fraction(0)
    return Fraction(1, prime**pivot_row)


def endpoint_path_with_start_leaf_density(
    path_rows: tuple[dict, dict, dict, dict, dict],
    start_leaf: dict,
    residuals: tuple[int, int, int, int],
    *,
    projected_enabled: bool,
    endpoint_enabled: bool,
    start_leaf_enabled: bool,
    compatible: bool,
) -> Fraction:
    """Exact six-event path/leaf union by componentwise elimination."""
    event_rows = (*path_rows, start_leaf)
    component_events = (
        (0, 1, 5),
        (1, 2),
        (2, 3),
        (3, 4),
    )
    enabled = {1, 2, 3}
    if projected_enabled:
        enabled.add(0)
    if endpoint_enabled:
        enabled.add(4)
    if start_leaf_enabled:
        enabled.add(5)
    answer = Fraction(0)
    enabled_tuple = tuple(sorted(enabled))
    for size in range(1, len(enabled_tuple) + 1):
        sign = 1 if size % 2 else -1
        for selected_tuple in itertools.combinations(enabled_tuple, size):
            selected = set(selected_tuple)
            intersection = Fraction(1)
            for component_index, (prime, touching) in enumerate(
                zip(residuals, component_events)
            ):
                equations = []
                for event in touching:
                    if event not in selected:
                        continue
                    target = 0
                    if (
                        component_index == 0
                        and event == 1
                        and not compatible
                    ):
                        target = 1
                    row = event_rows[event]
                    equations.append(
                        (int(row["a"]), int(row["b"]), target)
                    )
                intersection *= linear_equation_intersection_density(
                    prime,
                    tuple(equations),
                )
                if not intersection:
                    break
            answer += sign * intersection
    return answer


def endpoint_path_weights_with_start_leaf(
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
    """Add a fixed-target projected leaf on the path's first component."""
    path_tables = {}
    for compatible in (False, True):
        label = "compatible" if compatible else "incompatible"
        path_tables[label] = {}
        for projected_enabled, endpoint_enabled, leaf_enabled in (
            itertools.product((False, True), repeat=3)
        ):
            key = (
                f"{int(projected_enabled)}"
                f"{int(endpoint_enabled)}"
                f"{int(leaf_enabled)}"
            )
            path_tables[label][key] = (
                endpoint_path_with_start_leaf_density(
                    path_rows,
                    start_leaf,
                    residuals,
                    projected_enabled=projected_enabled,
                    endpoint_enabled=endpoint_enabled,
                    start_leaf_enabled=leaf_enabled,
                    compatible=compatible,
                )
            )
    if any(
        path_tables["compatible"][key]
        < path_tables["incompatible"][key]
        for key in path_tables["compatible"]
    ):
        raise RuntimeError(
            "compatible start-leaf targets do not dominate pointwise"
        )

    weights = {}
    for state in itertools.product((False, True), repeat=6):
        (
            projected_active,
            independent_active,
            endpoint_active,
            paired_active,
            first_shared_active,
            leaf_active,
        ) = state
        path_key = (
            f"{int(projected_active)}"
            f"{int(endpoint_active)}"
            f"{int(leaf_active)}"
        )
        path_density = path_tables["compatible"][path_key]
        independent_enabled = (
            int(independent_active) + int(first_shared_active)
        )
        if first_shared and independent_enabled == 2:
            independent_density = (
                Fraction(2, independent_residual)
                - Fraction(1, independent_residual**2)
            )
        elif independent_enabled:
            independent_density = Fraction(1, independent_residual)
        else:
            independent_density = Fraction(0)
        if paired_residual is None:
            paired_density = Fraction(0)
        elif paired_active:
            paired_density = (
                Fraction(2, paired_residual)
                - Fraction(1, paired_residual**2)
            )
        else:
            paired_density = Fraction(1, paired_residual)
        weights[state] = 1 - (
            (1 - path_density)
            * (1 - independent_density)
            * (1 - paired_density)
        )
    return weights, path_tables


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--base-primes", required=True)
    parser.add_argument("--projected-prime", type=int, required=True)
    parser.add_argument("--lifted-prime", type=int, required=True)
    parser.add_argument("--tail-prime", type=int, required=True)
    parser.add_argument("--second-tail-prime", type=int, required=True)
    parser.add_argument(
        "--independent-projected-prime",
        type=int,
        required=True,
    )
    parser.add_argument(
        "--endpoint-projected-prime",
        type=int,
        required=True,
    )
    parser.add_argument("--paired-projected-prime", type=int)
    parser.add_argument("--paired-shared-prime", type=int)
    parser.add_argument("--first-shared-projected-prime", type=int)
    parser.add_argument("--second-independent-projected-prime", type=int)
    parser.add_argument("--normalized-start-shared-prime", type=int)
    parser.add_argument("--normalize-primes", default="")
    parser.add_argument("--max-base-cells", type=int, default=1_000_000)
    parser.add_argument("--max-target-tuples", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    base_primes = tuple(
        int(value) for value in args.base_primes.split(",") if value
    )
    normalization_primes = tuple(
        int(value)
        for value in args.normalize_primes.split(",")
        if value
    )
    extra_primes = (
        args.projected_prime,
        args.lifted_prime,
        args.tail_prime,
        args.second_tail_prime,
        args.independent_projected_prime,
        args.endpoint_projected_prime,
        *(
            (
                args.paired_projected_prime,
                args.paired_shared_prime,
            )
            if args.paired_projected_prime is not None
            and args.paired_shared_prime is not None
            else ()
        ),
        *(
            (args.first_shared_projected_prime,)
            if args.first_shared_projected_prime is not None
            else ()
        ),
        *(
            (args.second_independent_projected_prime,)
            if args.second_independent_projected_prime is not None
            else ()
        ),
        *(
            (args.normalized_start_shared_prime,)
            if args.normalized_start_shared_prime is not None
            else ()
        ),
    )
    all_primes = (*base_primes, *extra_primes)
    if (
        len(base_primes) < 2
        or len(set(all_primes)) != len(all_primes)
        or (
            (args.paired_projected_prime is None)
            != (args.paired_shared_prime is None)
        )
        or (
            normalization_primes
            and (
                len(normalization_primes) != 2
                or not set(normalization_primes) <= set(base_primes)
            )
        )
    ):
        raise SystemExit("invalid base, normalization, or extra primes")
    if np is None:
        raise RuntimeError("NumPy is required for exact histogram scoring")

    source = json.loads(args.pool.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    missing = set(all_primes) - by_prime.keys()
    if missing:
        raise RuntimeError(f"anchor primes are absent: {sorted(missing)}")
    bases = [by_prime[prime] for prime in base_primes]
    (
        projected,
        lifted,
        tail,
        second_tail,
        independent,
        endpoint,
    ) = (by_prime[prime] for prime in extra_primes[:6])
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
    first_shared_projected = (
        by_prime[args.first_shared_projected_prime]
        if args.first_shared_projected_prime is not None
        else None
    )
    second_independent = (
        by_prime[args.second_independent_projected_prime]
        if args.second_independent_projected_prime is not None
        else None
    )
    normalized_start_shared = (
        by_prime[args.normalized_start_shared_prime]
        if args.normalized_start_shared_prime is not None
        else None
    )
    anchors = [
        *bases,
        projected,
        lifted,
        tail,
        second_tail,
        independent,
        endpoint,
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
    if any(int(row.get("target_modulus", 1)) != 1 for row in anchors):
        raise RuntimeError("anchor targets must be unrestricted")
    if any(
        math.gcd(int(row["a"]), int(row["b"]), int(row["h"])) != 1
        for row in anchors[len(bases):]
    ):
        raise RuntimeError("an extra row is not target-surjective")

    base_period = math.lcm(*(int(row["h"]) for row in bases))
    base_cells = base_period**2
    if base_cells > args.max_base_cells:
        raise RuntimeError("base grid exceeds guard")
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
    independent_residual = (
        int(independent["h"]) // independent_projection
    )
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
    residuals = (
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
    path_rows = (projected, lifted, tail, second_tail, endpoint)
    path_shared = (residual0, residual1, residual2, residual3)
    determinants = tuple(
        int(first["a"]) * int(second["b"])
        - int(first["b"]) * int(second["a"])
        for first, second in zip(path_rows, path_rows[1:])
    )
    if (
        any(residual <= 1 for residual in residuals)
        or endpoint_residual != residual3
        or any(
            math.gcd(first, second) != 1
            for first, second in itertools.combinations(residuals, 2)
        )
        or any(math.gcd(residual, base_period) != 1 for residual in residuals)
        or int(lifted["h"]) != residual0 * residual1
        or int(tail["h"]) != residual1 * residual2
        or int(second_tail["h"]) != residual2 * residual3
        or int(independent["h"])
        != independent_projection * independent_residual
        or int(endpoint["h"]) != endpoint_projection * residual3
        or any(
            math.gcd(determinant, shared) != 1
            for determinant, shared in zip(determinants, path_shared)
        )
        or (
            paired_projected is not None
            and (
                paired_residual <= 1
                or int(paired_projected["h"])
                != paired_projection * paired_residual
                or int(paired_shared["h"]) != paired_residual
                or math.gcd(paired_determinant, paired_residual) != 1
            )
        )
        or (
            first_shared_projected is not None
            and (
                first_shared_residual != independent_residual
                or math.gcd(
                    first_shared_determinant,
                    independent_residual,
                )
                != 1
            )
        )
        or (
            second_independent is not None
            and (
                second_independent_residual <= 1
                or int(second_independent["h"])
                != second_independent_projection
                * second_independent_residual
            )
        )
        or (
            normalized_start_shared is not None
            and (
                second_independent is None
                or not normalization_primes
                or normalized_start_projection <= 1
                or normalized_start_residual != residual0
                or int(normalized_start_shared["h"])
                != normalized_start_projection * residual0
                or math.gcd(
                    normalized_start_projected_determinant,
                    residual0,
                )
                != 1
                or math.gcd(
                    normalized_start_lifted_determinant,
                    residual0,
                )
                != 1
            )
        )
    ):
        raise RuntimeError("rows do not form the required CRT endpoint path")

    normalizer_indices = {
        base_primes.index(prime) for prime in normalization_primes
    }
    normalization_period = 1
    normalization_surjective = not normalization_primes
    if normalization_primes:
        normalizers = [
            *[by_prime[prime] for prime in normalization_primes],
            *(
                [normalized_start_shared]
                if normalized_start_shared is not None
                else []
            ),
        ]
        normalization_period = math.lcm(
            *(int(row["h"]) for row in normalizers)
        )
        image = {
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
        normalization_surjective = len(image) == math.prod(
            int(row["h"]) for row in normalizers
        )
        if not normalization_surjective:
            raise RuntimeError("normalization map is not jointly surjective")
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
                target = (
                    int(row["a"]) * k + int(row["b"]) * l
                ) % modulus
                masks[target] |= 1 << (k * base_period + l)
        return masks

    base_masks = [line_masks(row, int(row["h"])) for row in bases]
    full_grid_mask = (1 << base_cells) - 1
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
    projection_masks = [
        line_masks(row, modulus)
        for row, modulus in zip(projection_rows, projection_moduli)
    ]
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
            endpoint_path_weights_with_start_leaf(
                path_rows,
                normalized_start_shared,
                path_shared,
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
        weights = endpoint_path_weights(
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
            raise RuntimeError("factorized scores exceed uint64")
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
            normalized_start_active_mask = line_masks(
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
                "factorized projection cells do not partition the grid"
            )
        core_choice_weights = np.empty(
            (len(core_codes), len(core_observed_codes)),
            dtype=score_dtype,
        )
        for choice_index, choice in enumerate(core_codes):
            for cell_index, code in enumerate(core_observed_codes):
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
                core_choice_weights[choice_index, cell_index] = (
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
            scores = factorized_projection_scores(
                core_choice_weights,
                histogram,
                extra_cell_choice_weights,
                core_denominator,
                extra_denominator,
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
                    f"base_assignments={base_index}/{target_tuples} "
                    f"projection_combinations={checked} "
                    f"best_score={maximum_score}",
                    flush=True,
                )
    else:
        denominator = math.lcm(
            *(value.denominator for value in weights.values())
        )
        score_denominator = denominator
        if base_cells * denominator >= int(np.iinfo(np.int64).max):
            raise RuntimeError("exact scores exceed int64")
        weight_numerators = {
            state: int(value * denominator)
            for state, value in weights.items()
        }
        projection_codes = tuple(
            itertools.product(
                *(range(modulus) for modulus in projection_moduli)
            )
        )
        code_masks = []
        for code in projection_codes:
            cells = full_grid_mask
            for masks, target in zip(projection_masks, code):
                cells &= masks[target]
            code_masks.append(cells)
        if sum(mask.bit_count() for mask in code_masks) != base_cells:
            raise RuntimeError(
                "projection cells do not partition the base grid"
            )
        score_matrix = np.empty(
            (len(projection_codes), len(projection_codes)),
            dtype=np.int64,
        )
        for choice_index, choice in enumerate(projection_codes):
            for cell_index, code in enumerate(projection_codes):
                state = tuple(
                    observed == target
                    for observed, target in zip(code, choice)
                )
                score_matrix[choice_index, cell_index] = (
                    weight_numerators[state]
                )
        for base_index, targets in enumerate(
            itertools.product(*target_ranges),
            start=1,
        ):
            base_union = 0
            for masks, target in zip(base_masks, targets):
                base_union |= masks[target]
            uncovered = full_grid_mask ^ base_union
            histogram = np.fromiter(
                (
                    (uncovered & cells).bit_count()
                    for cells in code_masks
                ),
                dtype=np.int64,
                count=len(code_masks),
            )
            scores = score_matrix @ histogram
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
                    f"base_assignments={base_index}/{target_tuples} "
                    f"projection_combinations={checked} "
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
    result = {
        "schema": (
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
        ),
        "pool": str(args.pool),
        "row_count": len(rows),
        "anchor_primes": list(all_primes),
        "anchor_rows": [
            {key: int(row[key]) for key in ("p", "h", "a", "b")}
            for row in anchors
        ],
        "base_primes": list(base_primes),
        "path_primes": [
            args.projected_prime,
            args.lifted_prime,
            args.tail_prime,
            args.second_tail_prime,
            args.endpoint_projected_prime,
        ],
        "independent_projected_prime": args.independent_projected_prime,
        "paired_projected_prime": args.paired_projected_prime,
        "paired_shared_prime": args.paired_shared_prime,
        "first_shared_projected_prime":
        args.first_shared_projected_prime,
        "second_independent_projected_prime":
        args.second_independent_projected_prime,
        "normalized_start_shared_prime":
        args.normalized_start_shared_prime,
        "normalization_primes": list(normalization_primes),
        "normalization_period": normalization_period,
        "normalization_jointly_surjective": normalization_surjective,
        "base_period": base_period,
        "base_cells": base_cells,
        "projection_moduli": list(projection_moduli),
        "path_residual_moduli": list(path_shared),
        "independent_residual_modulus": independent_residual,
        "paired_residual_modulus": paired_residual,
        "paired_determinant": paired_determinant,
        "paired_maps_jointly_surjective": (
            True if paired_projected is not None else None
        ),
        "first_shared_projection_modulus": first_shared_projection,
        "first_shared_residual_modulus": first_shared_residual,
        "first_shared_determinant": first_shared_determinant,
        "first_shared_maps_jointly_surjective": (
            True if first_shared_projected is not None else None
        ),
        "second_independent_projection_modulus":
        second_independent_projection,
        "second_independent_residual_modulus":
        second_independent_residual,
        "normalized_start_projection_modulus":
        normalized_start_projection,
        "normalized_start_projection_target": (
            0 if normalized_start_shared is not None else None
        ),
        "normalized_start_residual_modulus":
        normalized_start_residual,
        "normalized_start_projected_determinant":
        normalized_start_projected_determinant,
        "normalized_start_lifted_determinant":
        normalized_start_lifted_determinant,
        "normalized_start_pair_maps_jointly_surjective": (
            True if normalized_start_shared is not None else None
        ),
        "normalized_start_compatible_targets_dominate": (
            True if normalized_start_shared is not None else None
        ),
        "normalized_start_path_union_densities": (
            {
                compatibility: {
                    state: fraction_payload(value)
                    for state, value in table.items()
                }
                for compatibility, table in start_leaf_path_tables.items()
            }
            if start_leaf_path_tables is not None
            else None
        ),
        "path_determinants": list(determinants),
        "path_adjacent_maps_jointly_surjective": True,
        "path_event_densities": [
            fraction_payload(value) for value in path_densities
        ],
        "residual_union_weights": {
            "".join("1" if flag else "0" for flag in state):
            fraction_payload(value)
            for state, value in weights.items()
        },
        "enumerated_period": base_period * math.prod(residuals),
        "target_tuples": target_tuples,
        "projection_target_combinations_per_base":
        projection_target_count,
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
        f"base_targets={target_tuples} projection_targets="
        f"{projection_target_count} checked={checked} maximum={maximum} "
        f"loss={loss} pool_upper={upper}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
