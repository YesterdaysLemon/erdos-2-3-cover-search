#!/usr/bin/env python3
"""Independent replay of a small-subgroup conditional fibre certificate."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

try:
    import numpy as np
except ImportError:  # pragma: no cover - explicit runtime guard
    np = None


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def pair_index(first: dict, second: dict) -> int:
    h1, a1, b1 = (int(first[key]) for key in ("h", "a", "b"))
    h2, a2, b2 = (int(second[key]) for key in ("h", "a", "b"))
    return math.gcd(
        h1 * h2,
        h1 * a2,
        h1 * b2,
        h2 * a1,
        h2 * b1,
        a1 * b2 - a2 * b1,
    )


def generator_order(
    generator: tuple[int, ...],
    moduli: tuple[int, ...],
) -> int:
    return math.lcm(
        *(
            modulus // math.gcd(value, modulus)
            for value, modulus in zip(generator, moduli)
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if np is None:
        raise RuntimeError("NumPy is required for independent replay")

    source = json.loads(args.pool.read_text())
    cert = json.loads(args.certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    outside = by_prime[int(cert["outside_prime"])]
    anchor_primes = tuple(int(p) for p in cert["anchor_primes"])
    anchors = [by_prime[prime] for prime in anchor_primes]
    normalizer_index = anchor_primes.index(
        int(cert["normalizer_anchor_prime"])
    )
    moduli = tuple(int(row["h"]) for row in anchors)
    generators = tuple(
        tuple(int(value) for value in generator)
        for generator in cert["kernel_image_generators"]
    )
    # Independently enumerate all pairs of generator multiples.  This differs
    # from the generator's breadth-first subgroup traversal.
    first_order = generator_order(generators[0], moduli)
    first_subgroup = {
        tuple(
            (first * generators[0][axis]) % moduli[axis]
            for axis in range(len(moduli))
        )
        for first in range(first_order)
    }
    # Add disjoint cosets of the first cyclic subgroup until the second
    # generator returns to an existing coset. This visits the subgroup once
    # instead of scanning a potentially very redundant order1*order2 box.
    image = set(first_subgroup)
    shift = (0,) * len(moduli)
    while True:
        shift = tuple(
            (shift[axis] + generators[1][axis]) % moduli[axis]
            for axis in range(len(moduli))
        )
        if shift in image:
            break
        image.update(
            tuple(
                (value[axis] + shift[axis]) % moduli[axis]
                for axis in range(len(moduli))
            )
            for value in first_subgroup
        )
    shape = moduli[:normalizer_index] + moduli[normalizer_index + 1:]
    counts = np.zeros(shape, dtype=np.int64)
    normalizer_covered = 0
    for values in image:
        if values[normalizer_index] == 0:
            normalizer_covered += 1
            continue
        index = values[:normalizer_index] + values[normalizer_index + 1:]
        counts[index] += 1
    transformed = counts
    # Use matrix contractions, rather than the generator's sum-minus-cell
    # transform, to replay every "observed value differs from target" query.
    for axis, modulus in enumerate(shape):
        matrix = np.ones((modulus, modulus), dtype=np.int64)
        np.fill_diagonal(matrix, 0)
        transformed = np.tensordot(
            matrix,
            transformed,
            axes=(1, axis),
        )
        transformed = np.moveaxis(transformed, 0, axis)
    flat_index = int(np.argmax(transformed))
    maximizing_targets = tuple(
        int(value)
        for value in np.unravel_index(flat_index, transformed.shape)
    )
    maximum_uncovered = int(transformed[maximizing_targets])
    minimum_covered = len(image) - maximum_uncovered
    intersection = Fraction(
        minimum_covered,
        int(outside["h"]) * len(image),
    )
    rows_match = (
        cert["outside_row"]
        == {key: int(outside[key]) for key in ("p", "h", "a", "b")}
        and cert["anchor_rows"]
        == [
            {key: int(row[key]) for key in ("p", "h", "a", "b")}
            for row in anchors
        ]
    )
    verified = (
        cert.get("schema") == "conditional_fibre_overlap_v1"
        and str(args.pool) == cert["pool"]
        and len(by_prime) == len(rows)
        and rows_match
        and pair_index(outside, anchors[normalizer_index]) == 1
        and int(cert["kernel_image_size"]) == len(image)
        and int(cert["normalizer_covered_image_points"])
        == normalizer_covered
        and list(cert["target_tensor_shape"]) == list(shape)
        and int(cert["target_combinations"]) == math.prod(shape)
        and list(cert["maximizing_targets"]) == list(maximizing_targets)
        and int(cert["maximum_uncovered_image_points"])
        == maximum_uncovered
        and int(cert["minimum_covered_image_points"]) == minimum_covered
        and read_fraction(cert["forced_intersection_density"])
        == intersection
    )
    result = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "outside_prime": int(outside["p"]),
        "anchor_primes": list(anchor_primes),
        "kernel_image_size": len(image),
        "target_combinations": math.prod(shape),
        "forced_intersection_density": {
            "numerator": intersection.numerator,
            "denominator": intersection.denominator,
        },
        "verified": verified,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"outside={outside['p']} image={len(image)} "
        f"targets={math.prod(shape)} intersection={intersection} "
        f"verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
