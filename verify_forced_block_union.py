#!/usr/bin/env python3
"""Independently replay a small forced-block union certificate."""

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
    anchor_primes = tuple(
        int(prime) for prime in certificate["anchor_primes"]
    )
    anchors = [by_prime[prime] for prime in anchor_primes]
    normalization_primes = tuple(
        int(prime)
        for prime in certificate.get("normalization_primes", ())
    )
    normalizer_indices = {
        anchor_primes.index(prime) for prime in normalization_primes
    }
    normalization_period = 1
    normalization_cells = 1
    normalization_surjective = not normalization_primes
    if normalization_primes:
        if len(normalization_primes) != 2:
            raise RuntimeError("certificate normalization is not a pair")
        normalizers = [by_prime[prime] for prime in normalization_primes]
        normalization_period = math.lcm(
            *(int(row["h"]) for row in normalizers)
        )
        normalization_cells = normalization_period**2
        image = set()
        for k in range(normalization_period):
            for l in range(normalization_period):
                image.add(
                    tuple(
                        (
                            int(row["a"]) * k + int(row["b"]) * l
                        )
                        % int(row["h"])
                        for row in normalizers
                    )
                )
        normalization_surjective = len(image) == math.prod(
            int(row["h"]) for row in normalizers
        )
    period = math.lcm(*(int(row["h"]) for row in anchors))
    cells = period * period

    # This verifier deliberately uses inclusion-exclusion tables rather than
    # the certificate generator's bitset unions.  For every nonempty subset
    # of anchors, count the simultaneous target fibres once on the full grid.
    # Any proposed target tuple is then evaluated by the exact alternating
    # sum of those independently tabulated intersection counts.
    subset_indices = {
        mask: tuple(
            index
            for index in range(len(anchors))
            if mask & (1 << index)
        )
        for mask in range(1, 1 << len(anchors))
    }
    intersection_tables = {
        mask: {} for mask in subset_indices
    }
    for k in range(period):
        for l in range(period):
            values = tuple(
                (
                    int(row["a"]) * k + int(row["b"]) * l
                )
                % int(row["h"])
                for row in anchors
            )
            for mask, indices in subset_indices.items():
                key = tuple(values[index] for index in indices)
                table = intersection_tables[mask]
                table[key] = table.get(key, 0) + 1

    maximum_covered = -1
    maximizing_targets = None
    target_ranges = [
        range(1) if index in normalizer_indices else range(int(row["h"]))
        for index, row in enumerate(anchors)
    ]
    for targets in itertools.product(*target_ranges):
        covered = 0
        for mask, indices in subset_indices.items():
            key = tuple(targets[index] for index in indices)
            intersection = intersection_tables[mask].get(key, 0)
            covered += (
                intersection
                if len(indices) % 2
                else -intersection
            )
        if covered > maximum_covered:
            maximum_covered = covered
            maximizing_targets = targets

    maximum_union = Fraction(maximum_covered, cells)
    block_sum = sum(
        (Fraction(1, int(row["h"])) for row in anchors),
        Fraction(0),
    )
    overlap_loss = block_sum - maximum_union
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    upper_bound = total_density - overlap_loss
    verified = (
        str(args.pool) == certificate["pool"]
        and int(certificate["row_count"]) == len(rows)
        and int(certificate["enumerated_period"]) == period
        and int(certificate["enumerated_cells"]) == cells
        and int(certificate["target_tuples"])
        == math.prod(len(values) for values in target_ranges)
        and normalization_surjective
        and certificate.get("normalization_jointly_surjective", True)
        == normalization_surjective
        and int(certificate.get("normalization_period", 1))
        == normalization_period
        and int(certificate.get("normalization_cells", 1))
        == normalization_cells
        and int(certificate["maximum_covered_cells"]) == maximum_covered
        and read_fraction(certificate["maximum_block_union_density"])
        == maximum_union
        and read_fraction(certificate["block_individual_density_sum"])
        == block_sum
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
        "anchor_primes": list(anchor_primes),
        "normalization_primes": list(normalization_primes),
        "normalization_period": normalization_period,
        "normalization_cells": normalization_cells,
        "normalization_jointly_surjective": normalization_surjective,
        "period": period,
        "cells": cells,
        "target_tuples": math.prod(len(values) for values in target_ranges),
        "maximum_covered_cells": maximum_covered,
        "maximizing_targets": list(maximizing_targets),
        "maximum_union_density": {
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
        f"anchors={anchor_primes} max_union={maximum_union} "
        f"loss={overlap_loss} pool_upper={upper_bound} "
        f"verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
