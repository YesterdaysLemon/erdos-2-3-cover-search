#!/usr/bin/env python3
"""Certify a base block with one residual pair and projected fibre rows.

The main projected row has modulus d*s and the lifted row has modulus s*r.
Their determinant is required to be a unit modulo s.  Each of three further
rows has modulus d_i*r_i, where d_i divides the enumerated base period and
the new residual r_i is coprime to the base and every other residual.

All base targets and all four projection targets are enumerated exactly.
Residual coverage is then a CRT product, avoiding enumeration of the much
larger full period.
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
except ModuleNotFoundError:  # Unit helpers remain usable without NumPy.
    np = None


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def main_pair_densities(shared: int, lifted_residual: int) -> tuple[Fraction, Fraction]:
    active = (
        Fraction(1, shared)
        + Fraction(1, shared * lifted_residual)
        - Fraction(1, shared * shared * lifted_residual)
    )
    inactive = Fraction(1, shared * lifted_residual)
    return active, inactive


def main_chain_densities(
    shared: int,
    lifted_residual: int,
    tail_residual: int,
) -> tuple[Fraction, Fraction]:
    first = Fraction(1, shared)
    middle = Fraction(1, shared * lifted_residual)
    tail = Fraction(1, lifted_residual * tail_residual)
    first_middle = Fraction(1, shared * shared * lifted_residual)
    first_tail = Fraction(1, shared * lifted_residual * tail_residual)
    middle_tail = Fraction(
        1,
        shared * lifted_residual * lifted_residual * tail_residual,
    )
    triple = Fraction(
        1,
        shared * shared * lifted_residual * lifted_residual * tail_residual,
    )
    return (
        first
        + middle
        + tail
        - first_middle
        - first_tail
        - middle_tail
        + triple,
        middle + tail - middle_tail,
    )


def residual_union_weights(
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
    weights = {}
    for state in itertools.product(
        (False, True),
        repeat=1 + len(independent_residuals),
    ):
        uncovered = 1 - (active_main if state[0] else inactive_main)
        for enabled, active_density, inactive_density in zip(
            state[1:],
            active_densities,
            inactive_densities,
        ):
            uncovered *= 1 - (
                active_density if enabled else inactive_density
            )
        weights[state] = 1 - uncovered
    return weights


def residual_union_weights_with_first_projected_pair(
    active_main: Fraction,
    inactive_main: Fraction,
    independent_residuals: tuple[int, ...],
    event_active_densities: tuple[Fraction, ...],
    event_inactive_densities: tuple[Fraction, ...],
) -> dict[tuple[bool, ...], Fraction]:
    """Merge two conditional projected lines on the first residual."""
    first_residual = independent_residuals[0]
    weights = {}
    for state in itertools.product(
        (False, True),
        repeat=2 + len(independent_residuals),
    ):
        uncovered = 1 - (active_main if state[0] else inactive_main)
        first_enabled = state[1]
        shared_enabled = state[-1]
        enabled_count = int(first_enabled) + int(shared_enabled)
        if enabled_count == 1:
            first_density = Fraction(1, first_residual)
        elif enabled_count == 2:
            first_density = (
                Fraction(2, first_residual)
                - Fraction(1, first_residual * first_residual)
            )
        else:
            first_density = Fraction(0)
        uncovered *= 1 - first_density
        for enabled, active_density, inactive_density in zip(
            state[2:-1],
            event_active_densities[1:],
            event_inactive_densities[1:],
        ):
            uncovered *= 1 - (
                active_density if enabled else inactive_density
            )
        weights[state] = 1 - uncovered
    return weights


def split_three_events(
    universe: int,
    events: tuple[int, int, int],
) -> dict[tuple[bool, bool, bool], int]:
    """Partition a finite nonnegative bitset by three events."""
    if universe < 0 or any(event < 0 for event in events):
        raise ValueError("bitset partitions require finite nonnegative masks")
    pieces: dict[tuple[bool, ...], int] = {(): universe}
    for event in events:
        next_pieces = {}
        for state, cells in pieces.items():
            next_pieces[(*state, False)] = cells & ~event
            next_pieces[(*state, True)] = cells & event
        pieces = next_pieces
    return pieces


def factorized_projection_scores(
    core_choice_weights,
    joint_histogram,
    extra_cell_choice_weights,
    core_weight_denominator: int,
    extra_weight_denominator: int,
):
    """Score every core/extra target pair without a full dense tensor.

    ``core_choice_weights`` has core target choices on rows and observed
    core codes on columns. ``joint_histogram`` counts uncovered cells by
    observed core and extra codes. ``extra_cell_choice_weights`` has
    observed extra codes on rows and extra target choices on columns.
    """
    total_by_core = joint_histogram.sum(axis=1)
    core_score = core_choice_weights @ total_by_core
    core_complements = core_weight_denominator - core_choice_weights
    extension_score = (
        (core_complements @ joint_histogram)
        @ extra_cell_choice_weights
    )
    return (
        core_score[:, None] * extra_weight_denominator
        + extension_score
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--base-primes", required=True)
    parser.add_argument("--projected-prime", type=int, required=True)
    parser.add_argument("--lifted-prime", type=int, required=True)
    parser.add_argument("--tail-prime", type=int)
    parser.add_argument(
        "--independent-projected-primes",
        required=True,
        help="two to four comma-separated fibre primes",
    )
    parser.add_argument(
        "--third-shared-lifted-prime",
        type=int,
        help=(
            "a row on the third projected residual; when present, the "
            "third event is a jointly-surjective two-line pair"
        ),
    )
    parser.add_argument(
        "--first-shared-projected-prime",
        type=int,
        help=(
            "a second projected row sharing the first independent "
            "residual component"
        ),
    )
    parser.add_argument(
        "--factorized-independent-indices",
        default="",
        help=(
            "zero-based indices of plain independent projected rows to "
            "score as a separate CRT factor"
        ),
    )
    parser.add_argument("--normalize-primes", default="")
    parser.add_argument("--max-base-cells", type=int, default=1_000_000)
    parser.add_argument("--max-target-tuples", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    base_primes = tuple(
        int(value) for value in args.base_primes.split(",") if value
    )
    independent_primes = tuple(
        int(value)
        for value in args.independent_projected_primes.split(",")
        if value
    )
    factorized_independent_indices = tuple(
        int(value)
        for value in args.factorized_independent_indices.split(",")
        if value
    )
    normalization_primes = tuple(
        int(value)
        for value in args.normalize_primes.split(",")
        if value
    )
    shared_paired_index = (
        min(2, len(independent_primes) - 1)
        if args.third_shared_lifted_prime is not None
        else None
    )
    all_primes = (
        *base_primes,
        args.projected_prime,
        args.lifted_prime,
        *([args.tail_prime] if args.tail_prime is not None else []),
        *independent_primes,
        *(
            [args.third_shared_lifted_prime]
            if args.third_shared_lifted_prime is not None
            else []
        ),
        *(
            [args.first_shared_projected_prime]
            if args.first_shared_projected_prime is not None
            else []
        ),
    )
    if (
        len(base_primes) < 2
        or len(independent_primes) not in (2, 3, 4)
        or len(set(all_primes)) != len(all_primes)
        or (
            args.third_shared_lifted_prime is not None
            and len(independent_primes) < 2
        )
        or (
            args.first_shared_projected_prime is not None
            and shared_paired_index == 0
        )
        or len(set(factorized_independent_indices))
        != len(factorized_independent_indices)
        or any(
            not 0 <= index < len(independent_primes)
            for index in factorized_independent_indices
        )
        or (
            args.third_shared_lifted_prime is not None
            and shared_paired_index in factorized_independent_indices
        )
        or (
            args.first_shared_projected_prime is not None
            and 0 in factorized_independent_indices
        )
        or (
            normalization_primes
            and (
                len(normalization_primes) != 2
                or not set(normalization_primes) <= set(base_primes)
            )
        )
    ):
        raise SystemExit("invalid base, normalization, or anchor primes")

    source = json.loads(args.pool.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    missing = set(all_primes) - by_prime.keys()
    if missing:
        raise RuntimeError(f"anchor primes are absent: {sorted(missing)}")
    bases = [by_prime[prime] for prime in base_primes]
    projected = by_prime[args.projected_prime]
    lifted = by_prime[args.lifted_prime]
    tail = by_prime[args.tail_prime] if args.tail_prime is not None else None
    independents = [by_prime[prime] for prime in independent_primes]
    third_shared_lifted = (
        by_prime[args.third_shared_lifted_prime]
        if args.third_shared_lifted_prime is not None
        else None
    )
    first_shared_projected = (
        by_prime[args.first_shared_projected_prime]
        if args.first_shared_projected_prime is not None
        else None
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
    if any(int(row.get("target_modulus", 1)) != 1 for row in anchors):
        raise RuntimeError("anchor targets must be unrestricted")
    if any(
        math.gcd(int(row["a"]), int(row["b"]), int(row["h"])) != 1
        for row in [
            projected,
            lifted,
            *([tail] if tail else []),
            *independents,
            *([third_shared_lifted] if third_shared_lifted else []),
            *([first_shared_projected] if first_shared_projected else []),
        ]
    ):
        raise RuntimeError("an extra-row target map is not surjective")

    base_period = math.lcm(*(int(row["h"]) for row in bases))
    base_cells = base_period**2
    if base_cells > args.max_base_cells:
        raise RuntimeError("base grid exceeds guard")

    projected_h = int(projected["h"])
    projection = math.gcd(base_period, projected_h)
    shared = projected_h // projection
    lifted_h = int(lifted["h"])
    lifted_residual = (
        lifted_h // shared if lifted_h % shared == 0 else 0
    )
    determinant = (
        int(projected["a"]) * int(lifted["b"])
        - int(projected["b"]) * int(lifted["a"])
    )
    if (
        shared <= 1
        or math.gcd(shared, base_period) != 1
        or lifted_residual <= 1
        or math.gcd(lifted_h, base_period) != 1
        or math.gcd(shared, lifted_residual) != 1
        or math.gcd(determinant, shared) != 1
    ):
        raise RuntimeError("main projected/lifted pair is not a CRT pair")
    tail_residual = None
    tail_determinant = None
    if tail is not None:
        tail_h = int(tail["h"])
        tail_residual = (
            tail_h // lifted_residual
            if tail_h % lifted_residual == 0
            else 0
        )
        tail_determinant = (
            int(lifted["a"]) * int(tail["b"])
            - int(lifted["b"]) * int(tail["a"])
        )
        if (
            tail_residual <= 1
            or math.gcd(tail_h, base_period * shared) != 1
            or math.gcd(lifted_residual, tail_residual) != 1
            or math.gcd(tail_determinant, lifted_residual) != 1
        ):
            raise RuntimeError("tail row is not a valid residual extension")

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
    residuals = (
        shared,
        lifted_residual,
        *([tail_residual] if tail_residual is not None else []),
        *independent_residuals,
    )
    if any(
        residual <= 1
        or math.gcd(residual, base_period) != 1
        or int(row["h"]) != projection_value * residual
        for row, projection_value, residual in zip(
            independents,
            independent_projections,
            independent_residuals,
        )
    ) or any(
        math.gcd(first, second) != 1
        for first, second in itertools.combinations(residuals, 2)
    ):
        raise RuntimeError("independent rows do not have new CRT residuals")
    third_shared_determinant = None
    if third_shared_lifted is not None:
        third_residual = independent_residuals[shared_paired_index]
        third_shared_determinant = (
            int(independents[shared_paired_index]["a"])
            * int(third_shared_lifted["b"])
            - int(independents[shared_paired_index]["b"])
            * int(third_shared_lifted["a"])
        )
        if (
            int(third_shared_lifted["h"]) != third_residual
            or math.gcd(
                int(third_shared_lifted["a"]),
                int(third_shared_lifted["b"]),
                third_residual,
            )
            != 1
            or math.gcd(third_shared_determinant, third_residual) != 1
        ):
            raise RuntimeError(
                "third projected and lifted rows are not a shared pair"
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
            int(first_shared_projected["h"]) // first_shared_projection
        )
        first_shared_determinant = (
            int(independents[0]["a"])
            * int(first_shared_projected["b"])
            - int(independents[0]["b"])
            * int(first_shared_projected["a"])
        )
        if (
            first_shared_residual != independent_residuals[0]
            or math.gcd(
                first_shared_determinant,
                first_shared_residual,
            )
            != 1
        ):
            raise RuntimeError(
                "first projected rows are not a shared residual pair"
            )

    normalizer_indices = {
        base_primes.index(prime) for prime in normalization_primes
    }
    normalization_period = 1
    normalization_surjective = not normalization_primes
    if normalization_primes:
        normalizers = [by_prime[prime] for prime in normalization_primes]
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
    projection_rows = [
        projected,
        *independents,
        *(
            [first_shared_projected]
            if first_shared_projected is not None
            else []
        ),
    ]
    projection_moduli = (
        projection,
        *independent_projections,
        *(
            [first_shared_projection]
            if first_shared_projection is not None
            else []
        ),
    )
    projection_masks = [
        line_masks(row, modulus)
        for row, modulus in zip(projection_rows, projection_moduli)
    ]
    if tail_residual is None:
        active_main, inactive_main = main_pair_densities(
            shared,
            lifted_residual,
        )
    else:
        active_main, inactive_main = main_chain_densities(
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
    if third_shared_lifted is not None:
        third_residual = independent_residuals[shared_paired_index]
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
    if first_shared_projected is None:
        weights = residual_union_weights(
            active_main,
            inactive_main,
            independent_residuals,
            event_active_densities,
            event_inactive_densities,
        )
    else:
        weights = residual_union_weights_with_first_projected_pair(
            active_main,
            inactive_main,
            independent_residuals,
            event_active_densities,
            event_inactive_densities,
        )
    if np is None:
        raise RuntimeError(
            "NumPy is required for the exact projection-histogram engine"
        )
    maximum_score = -1
    maximizing_targets = None
    maximizing_projections = None
    maximizing_base_covered = None
    maximizing_counts = None
    target_combinations_checked = 0
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
            uncovered_residual = Fraction(1)
            for enabled, residual in zip(
                extra_state,
                factorized_residuals,
            ):
                if enabled:
                    uncovered_residual *= 1 - Fraction(1, residual)
            extra_weights[extra_state] = 1 - uncovered_residual
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
            raise RuntimeError("factorized projection score exceeds int64")
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

        def code_masks_for_axes(
            codes: tuple[tuple[int, ...], ...],
            axes: tuple[int, ...],
        ) -> list[int]:
            answer = []
            for code in codes:
                mask = full_grid_mask
                for axis, target in zip(axes, code):
                    mask &= projection_masks[axis][target]
                answer.append(mask)
            return answer

        core_code_masks = code_masks_for_axes(
            core_codes,
            core_projection_axes,
        )
        extra_code_masks = code_masks_for_axes(
            extra_codes,
            factorized_projection_axes,
        )
        if (
            sum(mask.bit_count() for mask in core_code_masks) != base_cells
            or sum(mask.bit_count() for mask in extra_code_masks)
            != base_cells
        ):
            raise RuntimeError(
                "factorized projection values do not partition the grid"
            )
        core_choice_weights = np.empty(
            (len(core_codes), len(core_codes)),
            dtype=np.int64,
        )
        for choice_index, choice in enumerate(core_codes):
            for cell_index, code in enumerate(core_codes):
                state = tuple(
                    observed == target
                    for observed, target in zip(code, choice)
                )
                value = core_weights[state]
                core_choice_weights[choice_index, cell_index] = (
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
            scores = factorized_projection_scores(
                core_choice_weights,
                histogram,
                extra_cell_choice_weights,
                core_weight_denominator,
                extra_weight_denominator,
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
            target_combinations_checked += (
                len(core_codes) * len(extra_codes)
            )
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
                    f"base_assignments={base_index}/{target_tuples} "
                    f"projection_combinations="
                    f"{target_combinations_checked} "
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
            raise RuntimeError("exact projection score exceeds int64")
        projection_codes = tuple(
            itertools.product(
                *(range(modulus) for modulus in projection_moduli)
            )
        )
        projection_code_masks = []
        for code in projection_codes:
            mask = full_grid_mask
            for masks, target in zip(projection_masks, code):
                mask &= masks[target]
            projection_code_masks.append(mask)
        if (
            sum(mask.bit_count() for mask in projection_code_masks)
            != base_cells
        ):
            raise RuntimeError(
                "projection-value cells do not partition the grid"
            )
        weight_matrix = np.empty(
            (len(projection_codes), len(projection_codes)),
            dtype=np.int64,
        )
        for choice_index, choice in enumerate(projection_codes):
            for cell_index, code in enumerate(projection_codes):
                state = tuple(
                    cell_target == chosen_target
                    for cell_target, chosen_target in zip(code, choice)
                )
                weight_matrix[choice_index, cell_index] = (
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
            histogram = np.fromiter(
                (
                    (uncovered & mask).bit_count()
                    for mask in projection_code_masks
                ),
                dtype=np.int64,
                count=len(projection_code_masks),
            )
            scores = weight_matrix @ histogram
            scores += base_score
            choice_index = int(np.argmax(scores))
            score = int(scores[choice_index])
            projections = projection_codes[choice_index]
            target_combinations_checked += len(projection_codes)
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
                        cell_target == chosen_target
                        for cell_target, chosen_target in zip(
                            code,
                            projections,
                        )
                    )
                    key = "".join(
                        "1" if flag else "0" for flag in state
                    )
                    maximizing_counts[key] += int(count)
            if base_index % 500 == 0:
                print(
                    f"base_assignments={base_index}/{target_tuples} "
                    f"projection_combinations="
                    f"{target_combinations_checked} "
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
            "projected_triple_event_block_v4"
            if factorized_independent_indices
            else (
                "projected_triple_event_block_v3"
                if first_shared_projected is not None
                else (
                    "projected_triple_event_block_v2"
                    if tail is not None
                    or third_shared_lifted is not None
                    else "projected_triple_independent_block_v1"
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
        "projected_prime": args.projected_prime,
        "lifted_prime": args.lifted_prime,
        "tail_prime": args.tail_prime,
        "independent_projected_primes": list(independent_primes),
        "factorized_independent_indices": list(
            factorized_independent_indices
        ),
        "core_projection_axes": list(core_projection_axes),
        "factorized_projection_axes": list(
            factorized_projection_axes
        ),
        "core_projection_moduli": list(core_projection_moduli),
        "factorized_projection_moduli": list(
            factorized_projection_moduli
        ),
        "third_shared_lifted_prime": args.third_shared_lifted_prime,
        "shared_paired_index": shared_paired_index,
        "first_shared_projected_prime":
        args.first_shared_projected_prime,
        "normalization_primes": list(normalization_primes),
        "normalization_period": normalization_period,
        "normalization_jointly_surjective": normalization_surjective,
        "base_period": base_period,
        "base_cells": base_cells,
        "projection_moduli": list(projection_moduli),
        "shared_residual_modulus": shared,
        "lifted_residual_modulus": lifted_residual,
        "tail_residual_modulus": tail_residual,
        "independent_residual_moduli": list(independent_residuals),
        "shared_determinant": determinant,
        "shared_pair_jointly_surjective": True,
        "tail_determinant": tail_determinant,
        "tail_pair_jointly_surjective": tail is not None,
        "third_shared_determinant": third_shared_determinant,
        "third_shared_pair_jointly_surjective": (
            third_shared_lifted is not None
        ),
        "first_shared_projection_modulus": first_shared_projection,
        "first_shared_residual_modulus": first_shared_residual,
        "first_shared_determinant": first_shared_determinant,
        "first_shared_pair_jointly_surjective": (
            first_shared_projected is not None
        ),
        "event_active_residual_densities": [
            fraction_payload(value) for value in event_active_densities
        ],
        "event_inactive_residual_densities": [
            fraction_payload(value) for value in event_inactive_densities
        ],
        "active_main_residual_union_density": fraction_payload(active_main),
        "inactive_main_residual_union_density": fraction_payload(inactive_main),
        "residual_union_weights": {
            "".join("1" if flag else "0" for flag in state): fraction_payload(value)
            for state, value in weights.items()
        },
        "enumerated_period": base_period * math.prod(residuals),
        "target_tuples": target_tuples,
        "projection_target_combinations_per_base": math.prod(
            projection_moduli
        ),
        "target_combinations_checked": target_combinations_checked,
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
        f"{math.prod(projection_moduli)} checked="
        f"{target_combinations_checked} maximum={maximum} loss={loss} "
        f"pool_upper={upper}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
