#!/usr/bin/env python3
"""Certify a phase-independent upper bound for a small anchor block."""

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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--primes", required=True)
    parser.add_argument(
        "--normalize-primes",
        default="",
        help=(
            "two jointly-surjective anchor primes whose targets are fixed "
            "to zero by translation symmetry"
        ),
    )
    parser.add_argument("--max-cells", type=int, default=10_000_000)
    parser.add_argument("--max-target-tuples", type=int, default=10_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    anchor_primes = tuple(
        int(value) for value in args.primes.split(",") if value
    )
    if len(anchor_primes) < 2 or len(set(anchor_primes)) != len(anchor_primes):
        raise SystemExit("--primes must contain at least two distinct values")
    normalization_primes = tuple(
        int(value)
        for value in args.normalize_primes.split(",")
        if value
    )
    if normalization_primes and (
        len(normalization_primes) != 2
        or len(set(normalization_primes)) != 2
        or not set(normalization_primes) <= set(anchor_primes)
    ):
        raise SystemExit(
            "--normalize-primes must be two distinct listed anchors"
        )

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    missing = set(anchor_primes) - by_prime.keys()
    if missing:
        raise RuntimeError(f"anchor primes are absent: {sorted(missing)}")
    if any(int(row.get("target_modulus", 1)) != 1 for row in rows):
        raise RuntimeError("certificate requires unrestricted row targets")
    anchors = [by_prime[prime] for prime in anchor_primes]
    normalizer_indices = {
        anchor_primes.index(prime) for prime in normalization_primes
    }
    normalization_surjective = not normalization_primes
    normalization_period = 1
    normalization_cells = 1
    if normalization_primes:
        normalizers = [by_prime[prime] for prime in normalization_primes]
        normalization_period = math.lcm(
            *(int(row["h"]) for row in normalizers)
        )
        normalization_cells = normalization_period**2
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
            raise RuntimeError(
                "declared normalization target map is not jointly surjective"
            )
    period = math.lcm(*(int(row["h"]) for row in anchors))
    cells = period * period
    target_tuples = math.prod(
        1 if index in normalizer_indices else int(row["h"])
        for index, row in enumerate(anchors)
    )
    if cells > args.max_cells:
        raise RuntimeError(f"anchor grid has {cells} cells, above guard")
    if target_tuples > args.max_target_tuples:
        raise RuntimeError(
            f"anchor target space has {target_tuples} tuples, above guard"
        )

    masks = []
    for row in anchors:
        h = int(row["h"])
        row_masks = [0] * h
        for k in range(period):
            for l in range(period):
                target = (
                    int(row["a"]) * k + int(row["b"]) * l
                ) % h
                row_masks[target] |= 1 << (k * period + l)
        masks.append(row_masks)

    maximum_covered = -1
    maximizing_targets = None
    target_ranges = [
        range(1) if index in normalizer_indices else range(int(row["h"]))
        for index, row in enumerate(anchors)
    ]
    for targets in itertools.product(*target_ranges):
        union_mask = 0
        for row_masks, target in zip(masks, targets):
            union_mask |= row_masks[target]
        covered = union_mask.bit_count()
        if covered > maximum_covered:
            maximum_covered = covered
            maximizing_targets = targets

    maximum_block_union = Fraction(maximum_covered, cells)
    block_density_sum = sum(
        (Fraction(1, int(row["h"])) for row in anchors),
        Fraction(0),
    )
    forced_overlap_loss = block_density_sum - maximum_block_union
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    union_upper_bound = total_density - forced_overlap_loss
    proved = union_upper_bound < 1
    result = {
        "pool": str(args.pool),
        "row_count": len(rows),
        "anchor_primes": list(anchor_primes),
        "anchor_rows": [
            {
                key: int(row[key])
                for key in ("p", "h", "a", "b")
            }
            for row in anchors
        ],
        "normalization_primes": list(normalization_primes),
        "normalization_period": normalization_period,
        "normalization_cells": normalization_cells,
        "normalization_jointly_surjective": normalization_surjective,
        "enumerated_period": period,
        "enumerated_cells": cells,
        "target_tuples": target_tuples,
        "maximum_covered_cells": maximum_covered,
        "maximizing_targets": list(maximizing_targets),
        "maximum_block_union_density": fraction_payload(
            maximum_block_union
        ),
        "block_individual_density_sum": fraction_payload(
            block_density_sum
        ),
        "forced_overlap_loss": fraction_payload(forced_overlap_loss),
        "total_pool_density": fraction_payload(total_density),
        "pool_union_density_upper_bound": fraction_payload(
            union_upper_bound
        ),
        "proved_no_cover": proved,
        "argument": (
            "After a proof-safe translation normalization, every remaining "
            "anchor-target tuple was enumerated on the full anchor period "
            "plane. Grouping those anchor rows first loses at least the "
            "stated density relative to their individual density sum; adding "
            "every other row by the union bound gives the final pool upper "
            "bound."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"anchors={anchor_primes} period={period} cells={cells} "
        f"target_tuples={target_tuples} max_union={maximum_block_union} "
        f"loss={forced_overlap_loss} pool_upper={union_upper_bound} "
        f"proved_no_cover={proved}",
        flush=True,
    )
    return 0 if proved else 1


if __name__ == "__main__":
    raise SystemExit(main())
