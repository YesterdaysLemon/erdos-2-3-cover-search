#!/usr/bin/env python3
"""Independently verify a projected-row anchor-block certificate."""

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
    certificate = json.loads(args.certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    base_primes = tuple(int(p) for p in certificate["base_primes"])
    bases = [by_prime[prime] for prime in base_primes]
    extra = by_prime[int(certificate["extra_prime"])]
    normalized_shared_prime = certificate.get("normalized_shared_prime")
    normalized_shared = (
        by_prime[int(normalized_shared_prime)]
        if normalized_shared_prime is not None
        else None
    )
    base_period = math.lcm(*(int(row["h"]) for row in bases))
    base_cells = base_period**2
    projection = math.gcd(base_period, int(extra["h"]))
    full_period = math.lcm(base_period, int(extra["h"]))
    residual = full_period // base_period
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
    normalization_primes = tuple(
        int(p) for p in certificate["normalization_primes"]
    )
    normalized_indices = {
        base_primes.index(prime) for prime in normalization_primes
    }
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

    # Independent inclusion-exclusion tables for the base union and for its
    # intersection with every projected extra-row target.
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
    projected_sizes = (
        {
            (target, active): 0
            for target in range(projection)
            for active in (False, True)
        }
        if normalized_shared is not None
        else [0] * projection
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
                int(extra["a"]) * k + int(extra["b"]) * l
            ) % projection
            normalized_active = (
                (
                    int(normalized_shared["a"]) * k
                    + int(normalized_shared["b"]) * l
                )
                % normalized_projection
                == 0
                if normalized_shared is not None
                else None
            )
            if normalized_shared is not None:
                projected_sizes[
                    (projected_target, normalized_active)
                ] += 1
            else:
                projected_sizes[projected_target] += 1
            for mask, indices in subset_indices.items():
                key = tuple(values[index] for index in indices)
                base_table = base_tables[mask]
                base_table[key] = base_table.get(key, 0) + 1
                projected_key = (
                    (
                        projected_target,
                        normalized_active,
                        *key,
                    )
                    if normalized_shared is not None
                    else (projected_target, *key)
                )
                projected_table = projected_tables[mask]
                projected_table[projected_key] = (
                    projected_table.get(projected_key, 0) + 1
                )

    target_ranges = [
        range(1) if index in normalized_indices else range(int(row["h"]))
        for index, row in enumerate(bases)
    ]
    maximum_union = Fraction(-1)
    maximizing_targets = None
    maximizing_projection = None
    for targets in itertools.product(*target_ranges):
        base_covered = 0
        projected_union = (
            {
                (target, active): 0
                for target in range(projection)
                for active in (False, True)
            }
            if normalized_shared is not None
            else [0] * projection
        )
        for mask, indices in subset_indices.items():
            key = tuple(targets[index] for index in indices)
            sign = 1 if len(indices) % 2 else -1
            base_covered += sign * base_tables[mask].get(key, 0)
            table = projected_tables[mask]
            for projected_target in range(projection):
                if normalized_shared is not None:
                    for normalized_active in (False, True):
                        projected_union[
                            (projected_target, normalized_active)
                        ] += sign * table.get(
                            (
                                projected_target,
                                normalized_active,
                                *key,
                            ),
                            0,
                        )
                else:
                    projected_union[projected_target] += sign * table.get(
                        (projected_target, *key),
                        0,
                    )
        if normalized_shared is None:
            projected_new = max(
                projected_sizes[target] - projected_union[target]
                for target in range(projection)
            )
            projected_target = max(
                range(projection),
                key=lambda target: (
                    projected_sizes[target] - projected_union[target],
                    -target,
                ),
            )
            projected_single = projected_new
            projected_both = 0
            union = (
                Fraction(base_covered, base_cells)
                + Fraction(projected_new, base_cells * residual)
            )
        else:
            choices = []
            for projected_target in range(projection):
                uncovered = {
                    (target, active): (
                        projected_sizes[(target, active)]
                        - projected_union[(target, active)]
                    )
                    for target in range(projection)
                    for active in (False, True)
                }
                projected_both = uncovered[(projected_target, True)]
                projected_single = (
                    uncovered[(projected_target, False)]
                    + sum(
                        uncovered[(other, True)]
                        for other in range(projection)
                        if other != projected_target
                    )
                )
                union = (
                    Fraction(base_covered, base_cells)
                    + Fraction(
                        projected_single,
                        base_cells * residual,
                    )
                    + Fraction(
                        projected_both * (2 * residual - 1),
                        base_cells * residual**2,
                    )
                )
                choices.append(
                    (
                        union,
                        -projected_target,
                        projected_target,
                        projected_single,
                        projected_both,
                    )
                )
            (
                union,
                _,
                projected_target,
                projected_single,
                projected_both,
            ) = max(choices)
        if union > maximum_union:
            maximum_union = union
            maximizing_targets = targets
            maximizing_projection = projected_target
            maximizing_projected_single = projected_single
            maximizing_projected_both = projected_both

    individual_sum = sum(
        (Fraction(1, int(row["h"])) for row in bases),
        Fraction(1, int(extra["h"]))
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
    anchors_match = all(
        all(
            int(by_prime[int(record["p"])][key]) == int(record[key])
            for key in ("h", "a", "b")
        )
        for record in certificate["anchor_rows"]
    ) and list(certificate["anchor_primes"]) == [
        *base_primes,
        int(extra["p"]),
        *(
            (int(normalized_shared["p"]),)
            if normalized_shared is not None
            else ()
        ),
    ]
    structural = (
        normalized_shared is None
        or (
            bool(normalization_primes)
            and normalized_projection > 1
            and normalized_residual == residual
            and int(normalized_shared["h"])
            == normalized_projection * residual
            and math.gcd(shared_determinant, residual) == 1
        )
    )
    verified = (
        certificate.get("schema", "projected_row_block_v1")
        == (
            "projected_row_normalized_shared_v2"
            if normalized_shared is not None
            else "projected_row_block_v1"
        )
        and
        str(args.pool) == certificate["pool"]
        and int(certificate["row_count"]) == len(rows)
        and anchors_match
        and structural
        and normalization_surjective
        and int(certificate["base_period"]) == base_period
        and int(certificate["base_cells"]) == base_cells
        and int(certificate["projection_modulus"]) == projection
        and int(certificate["residual_modulus"]) == residual
        and certificate.get("normalized_projection_modulus")
        == normalized_projection
        and certificate.get("normalized_projection_target")
        == (0 if normalized_shared is not None else None)
        and certificate.get("normalized_residual_modulus")
        == normalized_residual
        and certificate.get("shared_determinant") == shared_determinant
        and certificate.get("shared_maps_jointly_surjective")
        == (True if normalized_shared is not None else None)
        and int(certificate["enumerated_period"]) == full_period
        and int(certificate["target_tuples"])
        == math.prod(len(values) for values in target_ranges)
        and read_fraction(certificate["maximum_block_union_density"])
        == maximum_union
        and certificate.get("maximizing_projected_single_base_cells")
        == (
            maximizing_projected_single
            if normalized_shared is not None
            else None
        )
        and certificate.get("maximizing_projected_both_base_cells")
        == (
            maximizing_projected_both
            if normalized_shared is not None
            else None
        )
        and read_fraction(certificate["block_individual_density_sum"])
        == individual_sum
        and read_fraction(certificate["forced_overlap_loss"])
        == overlap_loss
        and read_fraction(certificate["total_pool_density"])
        == total_density
        and read_fraction(certificate["pool_union_density_upper_bound"])
        == upper_bound
        and bool(certificate["proved_no_cover"]) == (upper_bound < 1)
    )
    result = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "anchor_primes": [
            *base_primes,
            int(extra["p"]),
            *(
                (int(normalized_shared["p"]),)
                if normalized_shared is not None
                else ()
            ),
        ],
        "base_primes": list(base_primes),
        "extra_prime": int(extra["p"]),
        "normalized_shared_prime": normalized_shared_prime,
        "normalization_jointly_surjective": normalization_surjective,
        "base_period": base_period,
        "base_cells": base_cells,
        "projection_modulus": projection,
        "residual_modulus": residual,
        "normalized_projection_modulus": normalized_projection,
        "normalized_residual_modulus": normalized_residual,
        "shared_determinant": shared_determinant,
        "enumerated_period": full_period,
        "target_tuples": math.prod(len(values) for values in target_ranges),
        "maximizing_base_targets": list(maximizing_targets),
        "maximizing_projection_target": maximizing_projection,
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
        "maximum_block_union_density": {
            "numerator": maximum_union.numerator,
            "denominator": maximum_union.denominator,
        },
        "forced_overlap_loss": {
            "numerator": overlap_loss.numerator,
            "denominator": overlap_loss.denominator,
        },
        "pool_union_upper_bound": {
            "numerator": upper_bound.numerator,
            "denominator": upper_bound.denominator,
        },
        "proved_no_cover": upper_bound < 1,
        "verified": verified,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base_period={base_period} projection={projection} "
        f"residual={residual} targets={result['target_tuples']} "
        f"max_union={maximum_union} loss={overlap_loss} "
        f"pool_upper={upper_bound} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
