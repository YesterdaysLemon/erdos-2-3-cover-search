#!/usr/bin/env python3
"""Extend an anchor block by one row through its projected lower component."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


def best_projected_target(
    projected_masks: list[int],
    union_mask: int,
) -> tuple[int, int]:
    """Return the target adding the most new cells, then the least target."""
    return max(
        (
            (
                target,
                (mask & ~union_mask).bit_count(),
            )
            for target, mask in enumerate(projected_masks)
        ),
        key=lambda item: (item[1], -item[0]),
    )


def normalized_shared_projection_counts(
    uncovered_mask: int,
    projected_mask: int,
    normalized_mask: int,
) -> tuple[int, int]:
    """Count cells with exactly one or both residual lines enabled."""
    single = (
        uncovered_mask & (projected_mask ^ normalized_mask)
    ).bit_count()
    both = (
        uncovered_mask & projected_mask & normalized_mask
    ).bit_count()
    return single, both


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--base-primes", required=True)
    parser.add_argument("--extra-prime", type=int, required=True)
    parser.add_argument("--normalized-shared-prime", type=int)
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
    all_primes = (
        *base_primes,
        args.extra_prime,
        *(
            (args.normalized_shared_prime,)
            if args.normalized_shared_prime is not None
            else ()
        ),
    )
    if (
        len(base_primes) < 2
        or len(set(all_primes)) != len(all_primes)
    ):
        raise SystemExit("invalid base or extra primes")
    if normalization_primes and (
        len(normalization_primes) != 2
        or not set(normalization_primes) <= set(base_primes)
    ):
        raise SystemExit("normalization must be two base primes")

    source = json.loads(args.pool.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    missing = set(all_primes) - by_prime.keys()
    if missing:
        raise RuntimeError(f"anchor primes are absent: {sorted(missing)}")
    bases = [by_prime[prime] for prime in base_primes]
    extra = by_prime[args.extra_prime]
    normalized_shared = (
        by_prime[args.normalized_shared_prime]
        if args.normalized_shared_prime is not None
        else None
    )
    if any(
        int(row.get("target_modulus", 1)) != 1
        for row in [
            *bases,
            extra,
            *([normalized_shared] if normalized_shared is not None else []),
        ]
    ):
        raise RuntimeError("anchor targets must be unrestricted")
    base_period = math.lcm(*(int(row["h"]) for row in bases))
    base_cells = base_period**2
    if base_cells > args.max_base_cells:
        raise RuntimeError("base grid exceeds guard")
    extra_h = int(extra["h"])
    projection_modulus = math.gcd(base_period, extra_h)
    full_period = math.lcm(base_period, extra_h)
    residual_modulus = full_period // base_period
    if extra_h != projection_modulus * residual_modulus:
        raise AssertionError("unexpected period decomposition")
    if math.gcd(int(extra["a"]), int(extra["b"]), extra_h) != 1:
        raise RuntimeError("extra row is not surjective")
    normalized_projection = None
    normalized_residual = None
    shared_determinant = None
    if normalized_shared is not None:
        normalized_projection = math.gcd(
            base_period,
            int(normalized_shared["h"]),
        )
        normalized_residual = (
            int(normalized_shared["h"]) // normalized_projection
        )
        shared_determinant = (
            int(extra["a"]) * int(normalized_shared["b"])
            - int(extra["b"]) * int(normalized_shared["a"])
        )
        if (
            not normalization_primes
            or normalized_projection <= 1
            or normalized_residual != residual_modulus
            or int(normalized_shared["h"])
            != normalized_projection * residual_modulus
            or math.gcd(shared_determinant, residual_modulus) != 1
        ):
            raise RuntimeError(
                "normalized row does not share the projected residual"
            )

    normalizer_indices = {
        base_primes.index(prime) for prime in normalization_primes
    }
    normalization_surjective = not normalization_primes
    normalization_period = 1
    if normalization_primes:
        normalizers = [
            *[by_prime[prime] for prime in normalization_primes],
            *(
                [normalized_shared]
                if normalized_shared is not None
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
        result = [0] * modulus
        for k in range(base_period):
            for l in range(base_period):
                target = (
                    int(row["a"]) * k + int(row["b"]) * l
                ) % modulus
                result[target] |= 1 << (k * base_period + l)
        return result

    base_masks = [
        line_masks(row, int(row["h"])) for row in bases
    ]
    full_grid_mask = (1 << base_cells) - 1
    projected_masks = line_masks(extra, projection_modulus)
    normalized_active_mask = (
        line_masks(normalized_shared, normalized_projection)[0]
        if normalized_shared is not None
        else None
    )
    maximum_union = Fraction(-1)
    maximizing_targets = None
    maximizing_projection_target = None
    maximizing_base_covered = None
    maximizing_projected_new = None
    for targets in itertools.product(*target_ranges):
        union_mask = 0
        for masks, target in zip(base_masks, targets):
            union_mask |= masks[target]
        base_covered = union_mask.bit_count()
        if normalized_shared is None:
            projected_target, projected_new = best_projected_target(
                projected_masks,
                union_mask,
            )
            projected_single = projected_new
            projected_both = 0
            union_density = (
                Fraction(base_covered, base_cells)
                + Fraction(
                    projected_new,
                    base_cells * residual_modulus,
                )
            )
        else:
            uncovered = full_grid_mask ^ union_mask
            choices = []
            for projected_target, projected_mask in enumerate(
                projected_masks
            ):
                projected_single, projected_both = (
                    normalized_shared_projection_counts(
                        uncovered,
                        projected_mask,
                        normalized_active_mask,
                    )
                )
                union_density = (
                    Fraction(base_covered, base_cells)
                    + Fraction(
                        projected_single,
                        base_cells * residual_modulus,
                    )
                    + Fraction(
                        projected_both
                        * (2 * residual_modulus - 1),
                        base_cells * residual_modulus**2,
                    )
                )
                choices.append(
                    (
                        union_density,
                        -projected_target,
                        projected_target,
                        projected_single,
                        projected_both,
                    )
                )
            (
                union_density,
                _,
                projected_target,
                projected_single,
                projected_both,
            ) = max(choices)
            projected_new = projected_single + projected_both
        if union_density > maximum_union:
            maximum_union = union_density
            maximizing_targets = targets
            maximizing_projection_target = projected_target
            maximizing_base_covered = base_covered
            maximizing_projected_new = projected_new
            maximizing_projected_single = projected_single
            maximizing_projected_both = projected_both

    individual_sum = sum(
        (Fraction(1, int(row["h"])) for row in bases),
        Fraction(1, extra_h)
        + (
            Fraction(1, int(normalized_shared["h"]))
            if normalized_shared is not None
            else Fraction(0)
        ),
    )
    overlap_loss = individual_sum - maximum_union
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    upper_bound = total_density - overlap_loss
    proved = upper_bound < 1
    anchor_rows = [
        {
            key: int(row[key]) for key in ("p", "h", "a", "b")
        }
        for row in [
            *bases,
            extra,
            *([normalized_shared] if normalized_shared is not None else []),
        ]
    ]
    result = {
        "schema": (
            "projected_row_normalized_shared_v2"
            if normalized_shared is not None
            else "projected_row_block_v1"
        ),
        "pool": str(args.pool),
        "row_count": len(rows),
        "anchor_primes": list(all_primes),
        "anchor_rows": anchor_rows,
        "base_primes": list(base_primes),
        "extra_prime": args.extra_prime,
        "normalized_shared_prime": args.normalized_shared_prime,
        "normalization_primes": list(normalization_primes),
        "normalization_period": normalization_period,
        "normalization_jointly_surjective": normalization_surjective,
        "base_period": base_period,
        "base_cells": base_cells,
        "projection_modulus": projection_modulus,
        "residual_modulus": residual_modulus,
        "normalized_projection_modulus": normalized_projection,
        "normalized_projection_target": (
            0 if normalized_shared is not None else None
        ),
        "normalized_residual_modulus": normalized_residual,
        "shared_determinant": shared_determinant,
        "shared_maps_jointly_surjective": (
            True if normalized_shared is not None else None
        ),
        "enumerated_period": full_period,
        "target_tuples": target_tuples,
        "maximizing_base_targets": list(maximizing_targets),
        "maximizing_projection_target": maximizing_projection_target,
        "maximizing_base_covered_cells": maximizing_base_covered,
        "maximizing_projected_new_base_cells": maximizing_projected_new,
        "maximizing_projected_single_base_cells": (
            maximizing_projected_single
            if normalized_shared is not None
            else None
        ),
        "maximizing_projected_both_base_cells": (
            maximizing_projected_both
            if normalized_shared is not None
            else None
        ),
        "maximum_block_union_density": fraction_payload(maximum_union),
        "block_individual_density_sum": fraction_payload(individual_sum),
        "forced_overlap_loss": fraction_payload(overlap_loss),
        "total_pool_density": fraction_payload(total_density),
        "pool_union_density_upper_bound": fraction_payload(upper_bound),
        "proved_no_cover": proved,
        "argument": (
            "The extra row projects to one affine line modulo the gcd of "
            "its modulus and the base period. Above every compatible base "
            "cell, its primitive top-component equation covers exactly "
            "1/residual_modulus of the lifts. The listed maximum exhausts "
            "all normalized base targets and all projected extra targets."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base_period={base_period} projection={projection_modulus} "
        f"residual={residual_modulus} targets={target_tuples} "
        f"max_union={maximum_union} loss={overlap_loss} "
        f"pool_upper={upper_bound} proved_no_cover={proved}",
        flush=True,
    )
    return 0 if proved else 1


if __name__ == "__main__":
    raise SystemExit(main())
