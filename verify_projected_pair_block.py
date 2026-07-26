#!/usr/bin/env python3
"""Independently verify a projected-plus-lifted-row block certificate."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


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
    bases = [by_prime[prime] for prime in base_primes]
    projected = by_prime[int(cert["projected_prime"])]
    lifted = by_prime[int(cert["lifted_prime"])]
    tail_prime = cert.get("tail_prime")
    tail = (
        by_prime[int(tail_prime)]
        if tail_prime is not None
        else None
    )
    independent_projected_prime = cert.get(
        "independent_projected_prime"
    )
    independent_projected = (
        by_prime[int(independent_projected_prime)]
        if independent_projected_prime is not None
        else None
    )
    second_tail_prime = cert.get("second_tail_prime")
    second_tail = (
        by_prime[int(second_tail_prime)]
        if second_tail_prime is not None
        else None
    )
    second_independent_projected_prime = cert.get(
        "second_independent_projected_prime"
    )
    second_independent_projected = (
        by_prime[int(second_independent_projected_prime)]
        if second_independent_projected_prime is not None
        else None
    )
    if (
        second_independent_projected is not None
        and independent_projected is None
    ):
        raise RuntimeError(
            "second independent projection requires the first"
        )
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
    joint = math.gcd(determinant, shared) == 1
    component_union = (
        Fraction(1, shared)
        + Fraction(1, int(lifted["h"]))
        - Fraction(1, shared * int(lifted["h"]))
    )
    tail_residual = None
    tail_determinant = None
    tail_joint = False
    if tail is None:
        active_residual_density = component_union
        inactive_residual_density = Fraction(1, int(lifted["h"]))
    else:
        tail_residual = (
            int(tail["h"]) // lifted_residual
            if int(tail["h"]) % lifted_residual == 0
            else 0
        )
        tail_determinant = (
            int(lifted["a"]) * int(tail["b"])
            - int(lifted["b"]) * int(tail["a"])
        )
        tail_joint = math.gcd(
            tail_determinant,
            lifted_residual,
        ) == 1
        first_row = Fraction(1, shared)
        middle_row = Fraction(1, shared * lifted_residual)
        last_row = Fraction(1, lifted_residual * tail_residual)
        first_middle = Fraction(
            1,
            shared * shared * lifted_residual,
        )
        first_last = Fraction(
            1,
            shared * lifted_residual * tail_residual,
        )
        middle_last = Fraction(
            1,
            shared
            * lifted_residual
            * lifted_residual
            * tail_residual,
        )
        triple = Fraction(
            1,
            shared
            * shared
            * lifted_residual
            * lifted_residual
            * tail_residual,
        )
        active_residual_density = (
            first_row
            + middle_row
            + last_row
            - first_middle
            - first_last
            - middle_last
            + triple
        )
        inactive_residual_density = (
            middle_row + last_row - middle_last
        )
    second_tail_residual = None
    second_tail_determinant = None
    second_tail_joint = False
    if second_tail is not None:
        second_tail_residual = (
            int(second_tail["h"]) // tail_residual
            if tail_residual
            and int(second_tail["h"]) % tail_residual == 0
            else 0
        )
        second_tail_determinant = (
            int(tail["a"]) * int(second_tail["b"])
            - int(tail["b"]) * int(second_tail["a"])
        )
        second_tail_joint = math.gcd(
            second_tail_determinant,
            tail_residual,
        ) == 1
        second_tail_density = Fraction(
            1,
            tail_residual * second_tail_residual,
        )
        active_residual_density = 1 - (
            (1 - active_residual_density)
            * (1 - second_tail_density)
        )
        inactive_residual_density = 1 - (
            (1 - inactive_residual_density)
            * (1 - second_tail_density)
        )
    independent_projection = None
    independent_residual = None
    if independent_projected is not None:
        independent_projection = math.gcd(
            base_period,
            int(independent_projected["h"]),
        )
        independent_residual = (
            int(independent_projected["h"]) // independent_projection
        )
        independent_line_density = Fraction(
            1,
            independent_residual,
        )
        active_with_independent = 1 - (
            (1 - active_residual_density)
            * (1 - independent_line_density)
        )
        inactive_with_independent = 1 - (
            (1 - inactive_residual_density)
            * (1 - independent_line_density)
        )
    else:
        active_with_independent = None
        inactive_with_independent = None
    second_independent_projection = None
    second_independent_residual = None
    triple_projection_weights = None
    if second_independent_projected is not None:
        second_independent_projection = math.gcd(
            base_period,
            int(second_independent_projected["h"]),
        )
        second_independent_residual = (
            int(second_independent_projected["h"])
            // second_independent_projection
        )
        first_independent_density = Fraction(1, independent_residual)
        second_independent_density = Fraction(
            1,
            second_independent_residual,
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

    normalization_primes = tuple(
        int(value) for value in cert["normalization_primes"]
    )
    normalized_indices = {
        base_primes.index(prime) for prime in normalization_primes
    }
    normalizers = [by_prime[prime] for prime in normalization_primes]
    normalization_period = math.lcm(
        *(int(row["h"]) for row in normalizers)
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
    normalization_surjective = len(normalization_image) == math.prod(
        int(row["h"]) for row in normalizers
    )

    subset_indices = {
        mask: tuple(
            index
            for index in range(len(bases))
            if mask & (1 << index)
        )
        for mask in range(1, 1 << len(bases))
    }
    base_tables = {mask: {} for mask in subset_indices}
    projected_tables = {mask: {} for mask in subset_indices}
    projected_sizes = [0] * projection
    independent_tables = (
        {mask: {} for mask in subset_indices}
        if independent_projected is not None
        else None
    )
    joint_tables = (
        {mask: {} for mask in subset_indices}
        if independent_projected is not None
        else None
    )
    independent_sizes = (
        [0] * independent_projection
        if independent_projected is not None
        else None
    )
    joint_sizes = (
        [
            [0] * independent_projection
            for _ in range(projection)
        ]
        if independent_projected is not None
        else None
    )
    second_independent_tables = (
        {mask: {} for mask in subset_indices}
        if second_independent_projected is not None
        else None
    )
    projected_second_tables = (
        {mask: {} for mask in subset_indices}
        if second_independent_projected is not None
        else None
    )
    independent_second_tables = (
        {mask: {} for mask in subset_indices}
        if second_independent_projected is not None
        else None
    )
    triple_projection_tables = (
        {mask: {} for mask in subset_indices}
        if second_independent_projected is not None
        else None
    )
    second_independent_sizes = (
        [0] * second_independent_projection
        if second_independent_projected is not None
        else None
    )
    projected_second_sizes = (
        [
            [0] * second_independent_projection
            for _ in range(projection)
        ]
        if second_independent_projected is not None
        else None
    )
    independent_second_sizes = (
        [
            [0] * second_independent_projection
            for _ in range(independent_projection)
        ]
        if second_independent_projected is not None
        else None
    )
    triple_projection_sizes = (
        [
            [
                [0] * second_independent_projection
                for _ in range(independent_projection)
            ]
            for _ in range(projection)
        ]
        if second_independent_projected is not None
        else None
    )
    for k in range(base_period):
        for l in range(base_period):
            values = tuple(
                (
                    int(row["a"]) * k + int(row["b"]) * l
                )
                % int(row["h"])
                for row in bases
            )
            projected_target = (
                int(projected["a"]) * k
                + int(projected["b"]) * l
            ) % projection
            projected_sizes[projected_target] += 1
            if independent_projected is not None:
                independent_target = (
                    int(independent_projected["a"]) * k
                    + int(independent_projected["b"]) * l
                ) % independent_projection
                independent_sizes[independent_target] += 1
                joint_sizes[projected_target][independent_target] += 1
            if second_independent_projected is not None:
                second_independent_target = (
                    int(second_independent_projected["a"]) * k
                    + int(second_independent_projected["b"]) * l
                ) % second_independent_projection
                second_independent_sizes[
                    second_independent_target
                ] += 1
                projected_second_sizes[projected_target][
                    second_independent_target
                ] += 1
                independent_second_sizes[independent_target][
                    second_independent_target
                ] += 1
                triple_projection_sizes[projected_target][
                    independent_target
                ][second_independent_target] += 1
            for mask, indices in subset_indices.items():
                key = tuple(values[index] for index in indices)
                table = base_tables[mask]
                table[key] = table.get(key, 0) + 1
                projected_key = (projected_target, *key)
                projected_table = projected_tables[mask]
                projected_table[projected_key] = (
                    projected_table.get(projected_key, 0) + 1
                )
                if independent_projected is not None:
                    independent_key = (independent_target, *key)
                    independent_table = independent_tables[mask]
                    independent_table[independent_key] = (
                        independent_table.get(independent_key, 0) + 1
                    )
                    joint_key = (
                        projected_target,
                        independent_target,
                        *key,
                    )
                    joint_table = joint_tables[mask]
                    joint_table[joint_key] = (
                        joint_table.get(joint_key, 0) + 1
                    )
                if second_independent_projected is not None:
                    second_key = (second_independent_target, *key)
                    second_table = second_independent_tables[mask]
                    second_table[second_key] = (
                        second_table.get(second_key, 0) + 1
                    )
                    projected_second_key = (
                        projected_target,
                        second_independent_target,
                        *key,
                    )
                    projected_second_table = projected_second_tables[
                        mask
                    ]
                    projected_second_table[projected_second_key] = (
                        projected_second_table.get(
                            projected_second_key,
                            0,
                        )
                        + 1
                    )
                    independent_second_key = (
                        independent_target,
                        second_independent_target,
                        *key,
                    )
                    independent_second_table = (
                        independent_second_tables[mask]
                    )
                    independent_second_table[independent_second_key] = (
                        independent_second_table.get(
                            independent_second_key,
                            0,
                        )
                        + 1
                    )
                    triple_key = (
                        projected_target,
                        independent_target,
                        second_independent_target,
                        *key,
                    )
                    triple_table = triple_projection_tables[mask]
                    triple_table[triple_key] = (
                        triple_table.get(triple_key, 0) + 1
                    )

    target_ranges = [
        range(1) if index in normalized_indices else range(int(row["h"]))
        for index, row in enumerate(bases)
    ]
    maximum = Fraction(-1)
    maximizing_targets = None
    maximizing_projection = None
    maximizing_base_covered = None
    maximizing_compatible = None
    maximizing_independent_projection = None
    maximizing_second_independent_projection = None
    maximizing_category_counts = None
    for targets in itertools.product(*target_ranges):
        base_covered = 0
        projected_union = [0] * projection
        independent_union = (
            [0] * independent_projection
            if independent_projected is not None
            else None
        )
        joint_union = (
            [
                [0] * independent_projection
                for _ in range(projection)
            ]
            if independent_projected is not None
            else None
        )
        second_independent_union = (
            [0] * second_independent_projection
            if second_independent_projected is not None
            else None
        )
        projected_second_union = (
            [
                [0] * second_independent_projection
                for _ in range(projection)
            ]
            if second_independent_projected is not None
            else None
        )
        independent_second_union = (
            [
                [0] * second_independent_projection
                for _ in range(independent_projection)
            ]
            if second_independent_projected is not None
            else None
        )
        triple_projection_union = (
            [
                [
                    [0] * second_independent_projection
                    for _ in range(independent_projection)
                ]
                for _ in range(projection)
            ]
            if second_independent_projected is not None
            else None
        )
        for mask, indices in subset_indices.items():
            key = tuple(targets[index] for index in indices)
            sign = 1 if len(indices) % 2 else -1
            base_covered += sign * base_tables[mask].get(key, 0)
            table = projected_tables[mask]
            for target in range(projection):
                projected_union[target] += sign * table.get(
                    (target, *key),
                    0,
                )
            if independent_projected is not None:
                independent_table = independent_tables[mask]
                for target in range(independent_projection):
                    independent_union[target] += sign * (
                        independent_table.get((target, *key), 0)
                    )
                joint_table = joint_tables[mask]
                for projected_target in range(projection):
                    for independent_target in range(
                        independent_projection
                    ):
                        joint_union[projected_target][
                            independent_target
                        ] += sign * joint_table.get(
                            (
                                projected_target,
                                independent_target,
                                *key,
                            ),
                            0,
                        )
            if second_independent_projected is not None:
                second_table = second_independent_tables[mask]
                for target in range(second_independent_projection):
                    second_independent_union[target] += sign * (
                        second_table.get((target, *key), 0)
                    )
                projected_second_table = projected_second_tables[mask]
                for projected_target in range(projection):
                    for second_target in range(
                        second_independent_projection
                    ):
                        projected_second_union[projected_target][
                            second_target
                        ] += sign * projected_second_table.get(
                            (projected_target, second_target, *key),
                            0,
                        )
                independent_second_table = independent_second_tables[
                    mask
                ]
                for independent_target in range(
                    independent_projection
                ):
                    for second_target in range(
                        second_independent_projection
                    ):
                        independent_second_union[independent_target][
                            second_target
                        ] += sign * independent_second_table.get(
                            (independent_target, second_target, *key),
                            0,
                        )
                triple_table = triple_projection_tables[mask]
                for projected_target in range(projection):
                    for independent_target in range(
                        independent_projection
                    ):
                        for second_target in range(
                            second_independent_projection
                        ):
                            triple_projection_union[
                                projected_target
                            ][independent_target][
                                second_target
                            ] += sign * triple_table.get(
                                (
                                    projected_target,
                                    independent_target,
                                    second_target,
                                    *key,
                                ),
                                0,
                            )
        if independent_projected is None:
            projection_target = max(
                range(projection),
                key=lambda target: (
                    projected_sizes[target] - projected_union[target],
                    -target,
                ),
            )
            independent_target = None
            second_independent_target = None
            compatible = (
                projected_sizes[projection_target]
                - projected_union[projection_target]
            )
            incompatible = base_cells - base_covered - compatible
            category_counts = None
            density = (
                Fraction(base_covered, base_cells)
                + Fraction(compatible, base_cells)
                * active_residual_density
                + Fraction(incompatible, base_cells)
                * inactive_residual_density
            )
        elif second_independent_projected is None:
            candidates = []
            for projection_target in range(projection):
                projected_count = (
                    projected_sizes[projection_target]
                    - projected_union[projection_target]
                )
                for independent_target in range(
                    independent_projection
                ):
                    independent_count = (
                        independent_sizes[independent_target]
                        - independent_union[independent_target]
                    )
                    both = (
                        joint_sizes[projection_target][independent_target]
                        - joint_union[projection_target][independent_target]
                    )
                    projected_only = projected_count - both
                    independent_only = independent_count - both
                    neither = (
                        base_cells
                        - base_covered
                        - both
                        - projected_only
                        - independent_only
                    )
                    candidate_density = (
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
                            candidate_density,
                            projection_target,
                            independent_target,
                            projected_count,
                            {
                                "both": both,
                                "projected_only": projected_only,
                                "independent_only": independent_only,
                                "neither": neither,
                            },
                        )
                    )
            (
                density,
                projection_target,
                independent_target,
                compatible,
                category_counts,
            ) = max(
                candidates,
                key=lambda item: (
                    item[0],
                    -item[1],
                    -item[2],
                ),
            )
            second_independent_target = None
        else:
            candidates = []
            for projection_target in range(projection):
                main_count = (
                    projected_sizes[projection_target]
                    - projected_union[projection_target]
                )
                for independent_target in range(
                    independent_projection
                ):
                    first_count = (
                        independent_sizes[independent_target]
                        - independent_union[independent_target]
                    )
                    main_first_count = (
                        joint_sizes[projection_target][independent_target]
                        - joint_union[projection_target][independent_target]
                    )
                    for second_independent_target in range(
                        second_independent_projection
                    ):
                        second_count = (
                            second_independent_sizes[
                                second_independent_target
                            ]
                            - second_independent_union[
                                second_independent_target
                            ]
                        )
                        main_second_count = (
                            projected_second_sizes[projection_target][
                                second_independent_target
                            ]
                            - projected_second_union[projection_target][
                                second_independent_target
                            ]
                        )
                        first_second_count = (
                            independent_second_sizes[independent_target][
                                second_independent_target
                            ]
                            - independent_second_union[independent_target][
                                second_independent_target
                            ]
                        )
                        all_three = (
                            triple_projection_sizes[projection_target][
                                independent_target
                            ][second_independent_target]
                            - triple_projection_union[projection_target][
                                independent_target
                            ][second_independent_target]
                        )
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
                        candidate_density = Fraction(
                            base_covered,
                            base_cells,
                        )
                        candidate_density += sum(
                            (
                                Fraction(count, base_cells)
                                * triple_projection_weights[state]
                                for state, count in counts.items()
                            ),
                            Fraction(0),
                        )
                        candidates.append(
                            (
                                candidate_density,
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
            (
                density,
                projection_target,
                independent_target,
                second_independent_target,
                compatible,
                category_counts,
            ) = max(
                candidates,
                key=lambda item: (
                    item[0],
                    -item[1],
                    -item[2],
                    -item[3],
                ),
            )
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
        Fraction(1, int(projected["h"]))
        + Fraction(1, int(lifted["h"]))
        + (Fraction(1, int(tail["h"])) if tail else Fraction(0)),
    )
    if second_tail is not None:
        individual_sum += Fraction(1, int(second_tail["h"]))
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
    loss = individual_sum - maximum
    total = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    upper = total - loss
    anchors_match = all(
        all(
            int(by_prime[int(record["p"])][key]) == int(record[key])
            for key in ("h", "a", "b")
        )
        for record in cert["anchor_rows"]
    )
    structural = (
        shared > 1
        and math.gcd(shared, base_period) == 1
        and lifted_residual > 1
        and math.gcd(int(lifted["h"]), base_period) == 1
        and math.gcd(shared, lifted_residual) == 1
        and all(
            int(row.get("target_modulus", 1)) == 1
            for row in [
                *bases,
                projected,
                lifted,
                *([tail] if tail else []),
                *([second_tail] if second_tail else []),
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
            ]
        )
        and all(
            math.gcd(int(row["a"]), int(row["b"]), int(row["h"])) == 1
            for row in (
                projected,
                lifted,
                *([tail] if tail else []),
                *([second_tail] if second_tail else []),
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
            )
        )
        and joint
        and (
            tail is None
            or (
                tail_residual > 1
                and math.gcd(
                    int(tail["h"]),
                    base_period * shared,
                )
                == 1
                and math.gcd(lifted_residual, tail_residual) == 1
                and tail_joint
            )
        )
        and (
            independent_projected is None
            or (
                independent_residual > 1
                and int(independent_projected["h"])
                == independent_projection * independent_residual
                and math.gcd(independent_residual, base_period) == 1
                and math.gcd(
                    independent_residual,
                    shared
                    * lifted_residual
                    * (tail_residual if tail_residual else 1)
                    * (
                        second_tail_residual
                        if second_tail_residual
                        else 1
                    ),
                )
                == 1
            )
        )
        and (
            second_tail is None
            or (
                tail is not None
                and second_tail_residual > 1
                and int(second_tail["h"])
                == tail_residual * second_tail_residual
                and math.gcd(
                    int(second_tail["h"]),
                    base_period * shared * lifted_residual,
                )
                == 1
                and math.gcd(
                    tail_residual,
                    second_tail_residual,
                )
                == 1
                and second_tail_joint
            )
        )
        and (
            second_independent_projected is None
            or (
                independent_projected is not None
                and second_independent_residual > 1
                and int(second_independent_projected["h"])
                == (
                    second_independent_projection
                    * second_independent_residual
                )
                and math.gcd(
                    second_independent_residual,
                    base_period,
                )
                == 1
                and math.gcd(
                    second_independent_residual,
                    shared
                    * lifted_residual
                    * (tail_residual if tail_residual else 1)
                    * (
                        second_tail_residual
                        if second_tail_residual
                        else 1
                    )
                    * independent_residual,
                )
                == 1
            )
        )
    )
    verified = (
        str(args.pool) == cert["pool"]
        and int(cert["row_count"]) == len(rows)
        and anchors_match
        and structural
        and normalization_surjective
        and int(cert["base_period"]) == base_period
        and int(cert["base_cells"]) == base_cells
        and int(cert["projection_modulus"]) == projection
        and int(cert["shared_residual_modulus"]) == shared
        and int(cert["lifted_residual_modulus"]) == lifted_residual
        and int(cert["shared_determinant"]) == determinant
        and bool(cert["shared_pair_jointly_surjective"]) == joint
        and read_fraction(cert["pair_component_union_density"])
        == component_union
        and cert.get("tail_residual_modulus") == tail_residual
        and cert.get("tail_determinant") == tail_determinant
        and bool(cert.get("tail_pair_jointly_surjective"))
        == (tail is not None and tail_joint)
        and cert.get("second_tail_residual_modulus")
        == second_tail_residual
        and cert.get("second_tail_determinant")
        == second_tail_determinant
        and bool(cert.get("second_tail_pair_jointly_surjective"))
        == (second_tail is not None and second_tail_joint)
        and cert.get("independent_projection_modulus")
        == independent_projection
        and cert.get("independent_residual_modulus")
        == independent_residual
        and cert.get("second_independent_projection_modulus")
        == second_independent_projection
        and cert.get("second_independent_residual_modulus")
        == second_independent_residual
        and read_fraction(
            cert.get(
                "active_residual_union_density",
                cert["pair_component_union_density"],
            )
        )
        == active_residual_density
        and read_fraction(
            cert.get(
                "inactive_residual_union_density",
                {
                    "numerator": 1,
                    "denominator": int(lifted["h"]),
                },
            )
        )
        == inactive_residual_density
        and int(cert["enumerated_period"])
        == (
            base_period
            * shared
            * lifted_residual
            * (tail_residual if tail_residual else 1)
            * (
                second_tail_residual
                if second_tail_residual
                else 1
            )
            * (independent_residual if independent_residual else 1)
            * (
                second_independent_residual
                if second_independent_residual
                else 1
            )
        )
        and int(cert["target_tuples"])
        == math.prod(len(values) for values in target_ranges)
        and list(maximizing_targets)
        == cert["maximizing_base_targets"]
        and int(cert["maximizing_projection_target"])
        == maximizing_projection
        and cert.get("maximizing_independent_projection_target")
        == maximizing_independent_projection
        and cert.get(
            "maximizing_second_independent_projection_target"
        )
        == maximizing_second_independent_projection
        and int(cert["maximizing_base_covered_cells"])
        == maximizing_base_covered
        and int(
            cert["maximizing_projected_compatible_uncovered_cells"]
        )
        == maximizing_compatible
        and cert.get("maximizing_uncovered_category_counts")
        == maximizing_category_counts
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
        "anchor_primes": [
            *base_primes,
            int(projected["p"]),
            int(lifted["p"]),
            *([int(tail["p"])] if tail else []),
            *(
                [int(independent_projected["p"])]
                if independent_projected
                else []
            ),
            *(
                [int(second_independent_projected["p"])]
                if second_independent_projected
                else []
            ),
            *([int(second_tail["p"])] if second_tail else []),
        ],
        "base_primes": list(base_primes),
        "projected_prime": int(projected["p"]),
        "lifted_prime": int(lifted["p"]),
        "tail_prime": int(tail["p"]) if tail else None,
        "second_tail_prime": (
            int(second_tail["p"]) if second_tail else None
        ),
        "independent_projected_prime": (
            int(independent_projected["p"])
            if independent_projected
            else None
        ),
        "second_independent_projected_prime": (
            int(second_independent_projected["p"])
            if second_independent_projected
            else None
        ),
        "normalization_jointly_surjective": normalization_surjective,
        "base_period": base_period,
        "base_cells": base_cells,
        "projection_modulus": projection,
        "shared_residual_modulus": shared,
        "lifted_residual_modulus": lifted_residual,
        "shared_determinant": determinant,
        "shared_pair_jointly_surjective": joint,
        "pair_component_union_density": {
            "numerator": component_union.numerator,
            "denominator": component_union.denominator,
        },
        "tail_residual_modulus": tail_residual,
        "tail_determinant": tail_determinant,
        "tail_pair_jointly_surjective": tail is not None and tail_joint,
        "second_tail_residual_modulus": second_tail_residual,
        "second_tail_determinant": second_tail_determinant,
        "second_tail_pair_jointly_surjective": (
            second_tail is not None and second_tail_joint
        ),
        "independent_projection_modulus": independent_projection,
        "independent_residual_modulus": independent_residual,
        "second_independent_projection_modulus": (
            second_independent_projection
        ),
        "second_independent_residual_modulus": (
            second_independent_residual
        ),
        "active_residual_union_density": {
            "numerator": active_residual_density.numerator,
            "denominator": active_residual_density.denominator,
        },
        "inactive_residual_union_density": {
            "numerator": inactive_residual_density.numerator,
            "denominator": inactive_residual_density.denominator,
        },
        "enumerated_period": (
            base_period
            * shared
            * lifted_residual
            * (tail_residual if tail_residual else 1)
            * (
                second_tail_residual
                if second_tail_residual
                else 1
            )
            * (independent_residual if independent_residual else 1)
            * (
                second_independent_residual
                if second_independent_residual
                else 1
            )
        ),
        "target_tuples": math.prod(
            len(values) for values in target_ranges
        ),
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
        "maximum_block_union_density": {
            "numerator": maximum.numerator,
            "denominator": maximum.denominator,
        },
        "block_individual_density_sum": {
            "numerator": individual_sum.numerator,
            "denominator": individual_sum.denominator,
        },
        "forced_overlap_loss": {
            "numerator": loss.numerator,
            "denominator": loss.denominator,
        },
        "pool_union_upper_bound": {
            "numerator": upper.numerator,
            "denominator": upper.denominator,
        },
        "proved_no_cover": upper < 1,
        "verified": verified,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base_period={base_period} projection={projection} "
        f"shared={shared} lifted_residual={lifted_residual} "
        f"tail_residual={tail_residual} "
        f"second_tail_residual={second_tail_residual} "
        f"independent_residual={independent_residual} "
        f"second_independent_residual={second_independent_residual} "
        f"targets={result['target_tuples']} max_union={maximum} "
        f"loss={loss} pool_upper={upper} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
