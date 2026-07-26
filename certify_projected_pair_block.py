#!/usr/bin/env python3
"""Extend a base block by a projected row and a lifted residual row.

The first extra row has modulus d*s, where d is its projection to the base
period and s is a new coprime residual component.  The second extra row has
modulus s*r and is coprime to the base period.  When their reductions modulo
s are jointly surjective, their exact union can be counted above each base
cell while the base target space is exhaustively enumerated.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def pair_union_density(shared: int, lifted_residual: int) -> Fraction:
    return (
        Fraction(1, shared)
        + Fraction(1, shared * lifted_residual)
        - Fraction(1, shared * shared * lifted_residual)
    )


def projected_pair_density(
    base_cells: int,
    base_covered: int,
    projected_compatible_uncovered: int,
    shared: int,
    lifted_residual: int,
) -> Fraction:
    incompatible_uncovered = (
        base_cells - base_covered - projected_compatible_uncovered
    )
    return (
        Fraction(base_covered, base_cells)
        + Fraction(projected_compatible_uncovered, base_cells)
        * pair_union_density(shared, lifted_residual)
        + Fraction(incompatible_uncovered, base_cells)
        * Fraction(1, shared * lifted_residual)
    )


def residual_chain_densities(
    first: int,
    middle: int,
    last: int,
) -> tuple[Fraction, Fraction]:
    """Return residual union densities with the first line active/inactive."""
    first_row = Fraction(1, first)
    middle_row = Fraction(1, first * middle)
    last_row = Fraction(1, middle * last)
    first_middle = Fraction(1, first * first * middle)
    first_last = Fraction(1, first * middle * last)
    middle_last = Fraction(1, first * middle * middle * last)
    triple = Fraction(
        1,
        first * first * middle * middle * last,
    )
    active = (
        first_row
        + middle_row
        + last_row
        - first_middle
        - first_last
        - middle_last
        + triple
    )
    inactive = middle_row + last_row - middle_last
    return active, inactive


def projected_chain_density(
    base_cells: int,
    base_covered: int,
    projected_compatible_uncovered: int,
    active_density: Fraction,
    inactive_density: Fraction,
) -> Fraction:
    incompatible_uncovered = (
        base_cells - base_covered - projected_compatible_uncovered
    )
    return (
        Fraction(base_covered, base_cells)
        + Fraction(projected_compatible_uncovered, base_cells)
        * active_density
        + Fraction(incompatible_uncovered, base_cells)
        * inactive_density
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--base-primes", required=True)
    parser.add_argument("--projected-prime", type=int, required=True)
    parser.add_argument("--lifted-prime", type=int, required=True)
    parser.add_argument("--tail-prime", type=int)
    parser.add_argument("--second-tail-prime", type=int)
    parser.add_argument("--independent-projected-prime", type=int)
    parser.add_argument("--second-independent-projected-prime", type=int)
    parser.add_argument("--normalize-primes", default="")
    parser.add_argument("--max-base-cells", type=int, default=10_000_000)
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
    if (
        len(base_primes) < 2
        or len(set(base_primes)) != len(base_primes)
        or args.projected_prime in base_primes
        or args.lifted_prime in base_primes
        or args.projected_prime == args.lifted_prime
        or (
            args.tail_prime is not None
            and (
                args.tail_prime in base_primes
                or args.tail_prime
                in (args.projected_prime, args.lifted_prime)
            )
        )
        or (
            args.independent_projected_prime is not None
            and (
                args.independent_projected_prime in base_primes
                or args.independent_projected_prime
                in (
                    args.projected_prime,
                    args.lifted_prime,
                    args.tail_prime,
                )
            )
        )
        or (
            args.second_tail_prime is not None
            and (
                args.tail_prime is None
                or args.second_tail_prime in base_primes
                or args.second_tail_prime
                in (
                    args.projected_prime,
                    args.lifted_prime,
                    args.tail_prime,
                    args.independent_projected_prime,
                )
            )
        )
        or (
            args.second_independent_projected_prime is not None
            and (
                args.independent_projected_prime is None
                or args.second_independent_projected_prime in base_primes
                or args.second_independent_projected_prime
                in (
                    args.projected_prime,
                    args.lifted_prime,
                    args.tail_prime,
                    args.second_tail_prime,
                    args.independent_projected_prime,
                )
            )
        )
    ):
        raise SystemExit("invalid base or residual-chain primes")
    if normalization_primes and (
        len(normalization_primes) != 2
        or not set(normalization_primes) <= set(base_primes)
    ):
        raise SystemExit("normalization must be two base primes")

    source = json.loads(args.pool.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    required = set(base_primes) | {
        args.projected_prime,
        args.lifted_prime,
    }
    if args.tail_prime is not None:
        required.add(args.tail_prime)
    if args.independent_projected_prime is not None:
        required.add(args.independent_projected_prime)
    if args.second_tail_prime is not None:
        required.add(args.second_tail_prime)
    if args.second_independent_projected_prime is not None:
        required.add(args.second_independent_projected_prime)
    missing = required - by_prime.keys()
    if missing:
        raise RuntimeError(f"anchor primes are absent: {sorted(missing)}")
    bases = [by_prime[prime] for prime in base_primes]
    projected = by_prime[args.projected_prime]
    lifted = by_prime[args.lifted_prime]
    tail = (
        by_prime[args.tail_prime]
        if args.tail_prime is not None
        else None
    )
    independent_projected = (
        by_prime[args.independent_projected_prime]
        if args.independent_projected_prime is not None
        else None
    )
    second_tail = (
        by_prime[args.second_tail_prime]
        if args.second_tail_prime is not None
        else None
    )
    second_independent_projected = (
        by_prime[args.second_independent_projected_prime]
        if args.second_independent_projected_prime is not None
        else None
    )
    if any(
        int(row.get("target_modulus", 1)) != 1
        for row in [
            *bases,
            projected,
            lifted,
            *([tail] if tail else []),
            *(
                [independent_projected]
                if independent_projected
                else []
            ),
            *([second_tail] if second_tail else []),
            *(
                [second_independent_projected]
                if second_independent_projected
                else []
            ),
        ]
    ):
        raise RuntimeError("anchor targets must be unrestricted")
    if any(
        math.gcd(int(row["a"]), int(row["b"]), int(row["h"])) != 1
        for row in [
            projected,
            lifted,
            *([tail] if tail else []),
            *(
                [independent_projected]
                if independent_projected
                else []
            ),
            *([second_tail] if second_tail else []),
            *(
                [second_independent_projected]
                if second_independent_projected
                else []
            ),
        ]
    ):
        raise RuntimeError("extra-row target map is not surjective")

    base_period = math.lcm(*(int(row["h"]) for row in bases))
    base_cells = base_period**2
    if base_cells > args.max_base_cells:
        raise RuntimeError("base grid exceeds guard")
    projected_h = int(projected["h"])
    projection = math.gcd(base_period, projected_h)
    shared = projected_h // projection
    if (
        projected_h != projection * shared
        or shared <= 1
        or math.gcd(shared, base_period) != 1
    ):
        raise RuntimeError("projected row lacks one new coprime residual")
    lifted_h = int(lifted["h"])
    if lifted_h % shared:
        raise RuntimeError("lifted row does not contain shared residual")
    lifted_residual = lifted_h // shared
    if (
        lifted_residual <= 1
        or math.gcd(lifted_h, base_period) != 1
        or math.gcd(shared, lifted_residual) != 1
    ):
        raise RuntimeError("lifted row is not fully independent of base")
    determinant = (
        int(projected["a"]) * int(lifted["b"])
        - int(projected["b"]) * int(lifted["a"])
    )
    if math.gcd(determinant, shared) != 1:
        raise RuntimeError(
            "extra rows are not jointly surjective on shared residual"
        )
    tail_residual = None
    tail_determinant = None
    if tail is not None:
        tail_h = int(tail["h"])
        if tail_h % lifted_residual:
            raise RuntimeError(
                "tail row does not contain lifted residual"
            )
        tail_residual = tail_h // lifted_residual
        if (
            tail_residual <= 1
            or math.gcd(tail_h, base_period * shared) != 1
            or math.gcd(lifted_residual, tail_residual) != 1
        ):
            raise RuntimeError(
                "tail row is not a single new residual extension"
            )
        tail_determinant = (
            int(lifted["a"]) * int(tail["b"])
            - int(lifted["b"]) * int(tail["a"])
        )
        if math.gcd(tail_determinant, lifted_residual) != 1:
            raise RuntimeError(
                "lifted and tail rows are not jointly surjective "
                "on their shared residual"
            )
    second_tail_residual = None
    second_tail_determinant = None
    if second_tail is not None:
        second_tail_h = int(second_tail["h"])
        if second_tail_h % int(tail_residual):
            raise RuntimeError(
                "second tail does not contain prior tail residual"
            )
        second_tail_residual = (
            second_tail_h // int(tail_residual)
        )
        if (
            second_tail_residual <= 1
            or math.gcd(
                second_tail_h,
                base_period * shared * lifted_residual,
            )
            != 1
            or math.gcd(
                int(tail_residual),
                second_tail_residual,
            )
            != 1
        ):
            raise RuntimeError(
                "second tail is not one new residual extension"
            )
        second_tail_determinant = (
            int(tail["a"]) * int(second_tail["b"])
            - int(tail["b"]) * int(second_tail["a"])
        )
        if (
            math.gcd(
                second_tail_determinant,
                int(tail_residual),
            )
            != 1
        ):
            raise RuntimeError(
                "tail rows are not jointly surjective on shared residual"
            )
    independent_projection = None
    independent_residual = None
    if independent_projected is not None:
        independent_h = int(independent_projected["h"])
        independent_projection = math.gcd(base_period, independent_h)
        independent_residual = independent_h // independent_projection
        residual_product = (
            shared
            * lifted_residual
            * (int(tail_residual) if tail_residual else 1)
            * (
                int(second_tail_residual)
                if second_tail_residual
                else 1
            )
        )
        if (
            independent_residual <= 1
            or math.gcd(independent_residual, base_period) != 1
            or math.gcd(independent_residual, residual_product) != 1
            or independent_h
            != independent_projection * independent_residual
        ):
            raise RuntimeError(
                "independent projected row lacks one new CRT residual"
            )
    second_independent_projection = None
    second_independent_residual = None
    if second_independent_projected is not None:
        second_independent_h = int(
            second_independent_projected["h"]
        )
        second_independent_projection = math.gcd(
            base_period,
            second_independent_h,
        )
        second_independent_residual = (
            second_independent_h // second_independent_projection
        )
        prior_residual_product = (
            shared
            * lifted_residual
            * (int(tail_residual) if tail_residual else 1)
            * (
                int(second_tail_residual)
                if second_tail_residual
                else 1
            )
            * int(independent_residual)
        )
        if (
            second_independent_residual <= 1
            or math.gcd(
                second_independent_residual,
                base_period,
            )
            != 1
            or math.gcd(
                second_independent_residual,
                prior_residual_product,
            )
            != 1
            or second_independent_h
            != (
                second_independent_projection
                * second_independent_residual
            )
        ):
            raise RuntimeError(
                "second independent projected row lacks one new residual"
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

    base_masks = [
        line_masks(row, int(row["h"])) for row in bases
    ]
    projected_masks = line_masks(projected, projection)
    independent_masks = (
        line_masks(independent_projected, independent_projection)
        if independent_projected is not None
        else None
    )
    second_independent_masks = (
        line_masks(
            second_independent_projected,
            second_independent_projection,
        )
        if second_independent_projected is not None
        else None
    )
    maximum = Fraction(-1)
    maximizing_targets = None
    maximizing_projection = None
    maximizing_base_covered = None
    maximizing_compatible = None
    maximizing_independent_projection = None
    maximizing_second_independent_projection = None
    maximizing_category_counts = None
    if tail is None:
        active_residual_density = pair_union_density(
            shared,
            lifted_residual,
        )
        inactive_residual_density = Fraction(
            1,
            shared * lifted_residual,
        )
    else:
        active_residual_density, inactive_residual_density = (
            residual_chain_densities(
                shared,
                lifted_residual,
                int(tail_residual),
            )
        )
    if second_tail is not None:
        second_tail_density = Fraction(
            1,
            int(tail_residual) * int(second_tail_residual),
        )
        active_residual_density = 1 - (
            (1 - active_residual_density)
            * (1 - second_tail_density)
        )
        inactive_residual_density = 1 - (
            (1 - inactive_residual_density)
            * (1 - second_tail_density)
        )
    triple_projection_weights = None
    if second_independent_masks is not None:
        first_independent_density = Fraction(
            1,
            int(independent_residual),
        )
        second_independent_density = Fraction(
            1,
            int(second_independent_residual),
        )
        triple_projection_weights = {}
        for main_active in (False, True):
            for first_active in (False, True):
                for second_active in (False, True):
                    uncovered = 1 - (
                        active_residual_density
                        if main_active
                        else inactive_residual_density
                    )
                    if first_active:
                        uncovered *= 1 - first_independent_density
                    if second_active:
                        uncovered *= 1 - second_independent_density
                    triple_projection_weights[
                        (main_active, first_active, second_active)
                    ] = 1 - uncovered
    for targets in itertools.product(*target_ranges):
        union_mask = 0
        for masks, target in zip(base_masks, targets):
            union_mask |= masks[target]
        base_covered = union_mask.bit_count()
        if independent_masks is None:
            projection_target, compatible = max(
                (
                    (target, (mask & ~union_mask).bit_count())
                    for target, mask in enumerate(projected_masks)
                ),
                key=lambda item: (item[1], -item[0]),
            )
            density = projected_chain_density(
                base_cells,
                base_covered,
                compatible,
                active_residual_density,
                inactive_residual_density,
            )
            candidates = [
                (
                    density,
                    projection_target,
                    None,
                    None,
                    compatible,
                    None,
                )
            ]
        elif second_independent_masks is None:
            independent_line_density = Fraction(
                1,
                int(independent_residual),
            )
            active_with_independent = 1 - (
                (1 - active_residual_density)
                * (1 - independent_line_density)
            )
            inactive_with_independent = 1 - (
                (1 - inactive_residual_density)
                * (1 - independent_line_density)
            )
            candidates = []
            uncovered_mask = ~union_mask
            for projection_target, projection_mask in enumerate(
                projected_masks
            ):
                projected_uncovered = (
                    projection_mask & uncovered_mask
                )
                projected_count = projected_uncovered.bit_count()
                for (
                    independent_target,
                    independent_mask,
                ) in enumerate(independent_masks):
                    independent_uncovered = (
                        independent_mask & uncovered_mask
                    )
                    both = (
                        projected_uncovered & independent_uncovered
                    ).bit_count()
                    projected_only = projected_count - both
                    independent_only = (
                        independent_uncovered.bit_count() - both
                    )
                    neither = (
                        base_cells
                        - base_covered
                        - both
                        - projected_only
                        - independent_only
                    )
                    density = (
                        Fraction(base_covered, base_cells)
                        + Fraction(both, base_cells)
                        * active_with_independent
                        + Fraction(projected_only, base_cells)
                        * active_residual_density
                        + Fraction(independent_only, base_cells)
                        * inactive_with_independent
                        + Fraction(neither, base_cells)
                        * inactive_residual_density
                    )
                    candidates.append(
                        (
                            density,
                            projection_target,
                            independent_target,
                            None,
                            projected_count,
                            {
                                "both": both,
                                "projected_only": projected_only,
                                "independent_only": independent_only,
                                "neither": neither,
                            },
                        )
                    )
        else:
            candidates = []
            uncovered_mask = ~union_mask
            for projection_target, projection_mask in enumerate(
                projected_masks
            ):
                main_uncovered = projection_mask & uncovered_mask
                main_count = main_uncovered.bit_count()
                for (
                    independent_target,
                    independent_mask,
                ) in enumerate(independent_masks):
                    first_uncovered = (
                        independent_mask & uncovered_mask
                    )
                    first_count = first_uncovered.bit_count()
                    main_first = (
                        main_uncovered & first_uncovered
                    )
                    main_first_count = main_first.bit_count()
                    for (
                        second_independent_target,
                        second_independent_mask,
                    ) in enumerate(second_independent_masks):
                        second_uncovered = (
                            second_independent_mask & uncovered_mask
                        )
                        second_count = second_uncovered.bit_count()
                        main_second_count = (
                            main_uncovered & second_uncovered
                        ).bit_count()
                        first_second_count = (
                            first_uncovered & second_uncovered
                        ).bit_count()
                        all_three = (
                            main_first & second_uncovered
                        ).bit_count()
                        main_first_only = (
                            main_first_count - all_three
                        )
                        main_second_only = (
                            main_second_count - all_three
                        )
                        first_second_only = (
                            first_second_count - all_three
                        )
                        main_only = (
                            main_count
                            - main_first_count
                            - main_second_count
                            + all_three
                        )
                        first_only = (
                            first_count
                            - main_first_count
                            - first_second_count
                            + all_three
                        )
                        second_only = (
                            second_count
                            - main_second_count
                            - first_second_count
                            + all_three
                        )
                        none = (
                            base_cells
                            - base_covered
                            - all_three
                            - main_first_only
                            - main_second_only
                            - first_second_only
                            - main_only
                            - first_only
                            - second_only
                        )
                        counts = {
                            (True, True, True): all_three,
                            (True, True, False): main_first_only,
                            (True, False, True): main_second_only,
                            (False, True, True): first_second_only,
                            (True, False, False): main_only,
                            (False, True, False): first_only,
                            (False, False, True): second_only,
                            (False, False, False): none,
                        }
                        density = Fraction(base_covered, base_cells)
                        density += sum(
                            (
                                Fraction(count, base_cells)
                                * triple_projection_weights[state]
                                for state, count in counts.items()
                            ),
                            Fraction(0),
                        )
                        candidates.append(
                            (
                                density,
                                projection_target,
                                independent_target,
                                second_independent_target,
                                main_count,
                                {
                                    "".join(
                                        "1" if flag else "0"
                                        for flag in state
                                    ): count
                                    for state, count in counts.items()
                                },
                            )
                        )
        best = max(
            candidates,
            key=lambda item: (
                item[0],
                -item[1],
                -(item[2] if item[2] is not None else 0),
                -(item[3] if item[3] is not None else 0),
            ),
        )
        (
            density,
            projection_target,
            independent_target,
            second_independent_target,
            compatible,
            category_counts,
        ) = best
        if density > maximum:
            maximum = density
            maximizing_targets = targets
            maximizing_projection = projection_target
            maximizing_independent_projection = independent_target
            maximizing_second_independent_projection = (
                second_independent_target
            )
            maximizing_base_covered = base_covered
            maximizing_compatible = compatible
            maximizing_category_counts = category_counts

    individual_sum = sum(
        (Fraction(1, int(row["h"])) for row in bases),
        Fraction(1, projected_h)
        + Fraction(1, lifted_h)
        + (Fraction(1, int(tail["h"])) if tail else Fraction(0)),
        # This start value is extended below when a second projected row
        # is present.
    )
    if independent_projected is not None:
        individual_sum += Fraction(
            1,
            int(independent_projected["h"]),
        )
    if second_independent_projected is not None:
        individual_sum += Fraction(
            1,
            int(second_independent_projected["h"]),
        )
    if second_tail is not None:
        individual_sum += Fraction(1, int(second_tail["h"]))
    loss = individual_sum - maximum
    total = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    upper = total - loss
    anchor_rows = [
        {
            key: int(row[key]) for key in ("p", "h", "a", "b")
        }
        for row in [
            *bases,
            projected,
            lifted,
            *([tail] if tail else []),
            *(
                [independent_projected]
                if independent_projected
                else []
            ),
            *(
                [second_independent_projected]
                if second_independent_projected
                else []
            ),
            *([second_tail] if second_tail else []),
        ]
    ]
    result = {
        "pool": str(args.pool),
        "row_count": len(rows),
        "anchor_primes": [
            *base_primes,
            args.projected_prime,
            args.lifted_prime,
            *([args.tail_prime] if args.tail_prime is not None else []),
            *(
                [args.independent_projected_prime]
                if args.independent_projected_prime is not None
                else []
            ),
            *(
                [args.second_independent_projected_prime]
                if args.second_independent_projected_prime is not None
                else []
            ),
            *(
                [args.second_tail_prime]
                if args.second_tail_prime is not None
                else []
            ),
        ],
        "anchor_rows": anchor_rows,
        "base_primes": list(base_primes),
        "projected_prime": args.projected_prime,
        "lifted_prime": args.lifted_prime,
        "tail_prime": args.tail_prime,
        "second_tail_prime": args.second_tail_prime,
        "independent_projected_prime": (
            args.independent_projected_prime
        ),
        "second_independent_projected_prime": (
            args.second_independent_projected_prime
        ),
        "normalization_primes": list(normalization_primes),
        "normalization_period": normalization_period,
        "normalization_jointly_surjective": normalization_surjective,
        "base_period": base_period,
        "base_cells": base_cells,
        "projection_modulus": projection,
        "shared_residual_modulus": shared,
        "lifted_residual_modulus": lifted_residual,
        "shared_determinant": determinant,
        "shared_pair_jointly_surjective": True,
        "pair_component_union_density": fraction_payload(
            pair_union_density(shared, lifted_residual)
        ),
        "tail_residual_modulus": tail_residual,
        "tail_determinant": tail_determinant,
        "tail_pair_jointly_surjective": tail is not None,
        "second_tail_residual_modulus": second_tail_residual,
        "second_tail_determinant": second_tail_determinant,
        "second_tail_pair_jointly_surjective": (
            second_tail is not None
        ),
        "independent_projection_modulus": independent_projection,
        "independent_residual_modulus": independent_residual,
        "second_independent_projection_modulus": (
            second_independent_projection
        ),
        "second_independent_residual_modulus": (
            second_independent_residual
        ),
        "active_residual_union_density": fraction_payload(
            active_residual_density
        ),
        "inactive_residual_union_density": fraction_payload(
            inactive_residual_density
        ),
        "enumerated_period": (
            base_period
            * shared
            * lifted_residual
            * (int(tail_residual) if tail_residual else 1)
            * (
                int(second_tail_residual)
                if second_tail_residual
                else 1
            )
            * (
                int(independent_residual)
                if independent_residual
                else 1
            )
            * (
                int(second_independent_residual)
                if second_independent_residual
                else 1
            )
        ),
        "target_tuples": target_tuples,
        "maximizing_base_targets": list(maximizing_targets),
        "maximizing_projection_target": maximizing_projection,
        "maximizing_independent_projection_target": (
            maximizing_independent_projection
        ),
        "maximizing_second_independent_projection_target": (
            maximizing_second_independent_projection
        ),
        "maximizing_base_covered_cells": maximizing_base_covered,
        "maximizing_projected_compatible_uncovered_cells": (
            maximizing_compatible
        ),
        "maximizing_uncovered_category_counts": (
            maximizing_category_counts
        ),
        "maximum_block_union_density": fraction_payload(maximum),
        "block_individual_density_sum": fraction_payload(individual_sum),
        "forced_overlap_loss": fraction_payload(loss),
        "total_pool_density": fraction_payload(total),
        "pool_union_density_upper_bound": fraction_payload(upper),
        "proved_no_cover": upper < 1,
        "argument": (
            "The base targets and projected lower target are exhaustively "
            "enumerated. Above every uncovered base cell, ordinary "
            "inclusion-exclusion gives the exact residual-chain union. "
            "The adjacent two-target maps are jointly surjective on their "
            "shared components."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base_period={base_period} projection={projection} "
        f"shared={shared} lifted_residual={lifted_residual} "
        f"tail_residual={tail_residual} "
        f"second_tail_residual={second_tail_residual} "
        f"independent_residual={independent_residual} "
        f"second_independent_residual={second_independent_residual} "
        f"targets={target_tuples} max_union={maximum} loss={loss} "
        f"pool_upper={upper} proved_no_cover={upper < 1}",
        flush=True,
    )
    return 0 if upper < 1 else 1


if __name__ == "__main__":
    raise SystemExit(main())
