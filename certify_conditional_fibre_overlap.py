#!/usr/bin/env python3
"""Certify an exact small-subgroup conditional fibre intersection."""

from __future__ import annotations

import argparse
import collections
import json
import math
from fractions import Fraction
from pathlib import Path

try:
    import numpy as np
except ImportError:  # pragma: no cover - explicit runtime guard
    np = None

from certify_forced_pair_overlap import joint_index


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def extended_gcd(left: int, right: int) -> tuple[int, int, int]:
    if right == 0:
        return abs(left), 1 if left >= 0 else -1, 0
    gcd, x, y = extended_gcd(right, left % right)
    return gcd, y, x - (left // right) * y


def kernel_image(
    outside: dict,
    anchors: list[dict],
) -> tuple[set[tuple[int, ...]], tuple[tuple[int, ...], ...]]:
    a = int(outside["a"])
    b = int(outside["b"])
    h = int(outside["h"])
    gcd, u, v = extended_gcd(a, b)
    if math.gcd(gcd, h) != 1:
        raise RuntimeError("outside row is not primitive modulo h")
    # These two integer vectors are a determinant-h basis of the target-zero
    # kernel of a*k+b*l modulo h.
    basis = ((b // gcd, -a // gcd), (h * u, h * v))
    moduli = tuple(int(row["h"]) for row in anchors)
    generators = tuple(
        tuple(
            (
                int(row["a"]) * k + int(row["b"]) * l
            ) % modulus
            for row, modulus in zip(anchors, moduli)
        )
        for k, l in basis
    )
    zero = (0,) * len(anchors)
    image = {zero}
    queue = collections.deque([zero])
    while queue:
        current = queue.popleft()
        for generator in generators:
            successor = tuple(
                (value + step) % modulus
                for value, step, modulus
                in zip(current, generator, moduli)
            )
            if successor not in image:
                image.add(successor)
                queue.append(successor)
    return image, generators


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--outside-prime", type=int, required=True)
    parser.add_argument("--anchor-primes", required=True)
    parser.add_argument(
        "--normalizer-anchor-prime",
        type=int,
        required=True,
    )
    parser.add_argument(
        "--max-target-combinations",
        type=int,
        default=5_000_000,
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if np is None:
        raise RuntimeError("NumPy is required for exact complement replay")

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    anchor_primes = tuple(
        int(value) for value in args.anchor_primes.split(",") if value
    )
    if (
        len(anchor_primes) < 2
        or len(set(anchor_primes)) != len(anchor_primes)
        or args.normalizer_anchor_prime not in anchor_primes
        or args.outside_prime in anchor_primes
    ):
        raise SystemExit("invalid outside/anchor prime selection")
    missing = (
        {args.outside_prime, *anchor_primes} - by_prime.keys()
    )
    if missing:
        raise RuntimeError(f"primes are absent from the pool: {sorted(missing)}")
    outside = by_prime[args.outside_prime]
    anchors = [by_prime[prime] for prime in anchor_primes]
    normalizer_index = anchor_primes.index(
        args.normalizer_anchor_prime
    )
    if joint_index(outside, anchors[normalizer_index]) != 1:
        raise RuntimeError(
            "outside/normalizer target map is not jointly surjective"
        )

    image, generators = kernel_image(outside, anchors)
    moduli = tuple(int(row["h"]) for row in anchors)
    shape = moduli[:normalizer_index] + moduli[normalizer_index + 1:]
    target_combinations = math.prod(shape)
    if target_combinations > args.max_target_combinations:
        raise RuntimeError(
            f"target tensor has {target_combinations} cells, above guard"
        )
    counts = np.zeros(shape, dtype=np.int64)
    normalizer_covered = 0
    for values in image:
        if values[normalizer_index] == 0:
            normalizer_covered += 1
            continue
        index = values[:normalizer_index] + values[normalizer_index + 1:]
        counts[index] += 1
    transformed = counts
    for axis in range(transformed.ndim):
        transformed = (
            np.expand_dims(transformed.sum(axis=axis), axis)
            - transformed
        )
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
    result = {
        "schema": "conditional_fibre_overlap_v1",
        "pool": str(args.pool),
        "row_count": len(rows),
        "outside_prime": args.outside_prime,
        "outside_row": {
            key: int(outside[key]) for key in ("p", "h", "a", "b")
        },
        "anchor_primes": list(anchor_primes),
        "anchor_rows": [
            {key: int(row[key]) for key in ("p", "h", "a", "b")}
            for row in anchors
        ],
        "normalizer_anchor_prime": args.normalizer_anchor_prime,
        "kernel_image_generators": [list(value) for value in generators],
        "kernel_image_size": len(image),
        "normalizer_covered_image_points": normalizer_covered,
        "target_tensor_shape": list(shape),
        "target_combinations": target_combinations,
        "maximizing_targets": list(maximizing_targets),
        "maximum_uncovered_image_points": maximum_uncovered,
        "minimum_covered_image_points": minimum_covered,
        "forced_intersection_density": fraction_payload(intersection),
        "argument": (
            "After normalizing the outside and one anchor target to zero, "
            "the two displayed kernel-image generators enumerate exactly the "
            "anchor values on the outside fibre. Successive complement "
            "transforms evaluate every remaining anchor-target tuple. The "
            "largest uncovered image count therefore gives the exact minimum "
            "intersection with the selected anchor subunion."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"outside={args.outside_prime} anchors={len(anchors)} "
        f"image={len(image)} targets={target_combinations} "
        f"intersection={intersection}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
