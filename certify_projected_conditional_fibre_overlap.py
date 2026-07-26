#!/usr/bin/env python3
"""Certify a conditional fibre-to-anchor-union intersection lower bound.

The outside fibre and one anchor target are normalized to zero.  Base-anchor
targets are handled exactly on a finite base plane.  Each projected anchor
has one residual congruence, while the optional projected/shared pair has
two independent residual lines on the same prime component.
"""

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

from certify_forced_pair_overlap import joint_index


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def parse_primes(value: str) -> tuple[int, ...]:
    return tuple(int(item) for item in value.split(",") if item)


def parse_projected(value: str) -> tuple[tuple[int, int, int], ...]:
    if not value:
        return ()
    parsed = []
    for item in value.split(","):
        prime, projection, residual = (
            int(part) for part in item.split(":")
        )
        parsed.append((prime, projection, residual))
    return tuple(parsed)


def primitive_fibre_points(
    row: dict,
    period: int,
) -> object:
    """Enumerate one target-zero fibre without scanning the full plane."""
    h = int(row["h"])
    a = int(row["a"])
    b = int(row["b"])
    if period % h:
        raise RuntimeError("outside modulus does not divide the base period")
    if math.gcd(b, h) == 1:
        inverse = pow(b, -1, h)
        for k in range(period):
            first_l = (-a * k * inverse) % h
            for l in range(first_l, period, h):
                yield k, l
    elif math.gcd(a, h) == 1:
        inverse = pow(a, -1, h)
        for l in range(period):
            first_k = (-b * l * inverse) % h
            for k in range(first_k, period, h):
                yield k, l
    else:
        if math.gcd(a, b, h) != 1:
            raise RuntimeError("outside row is not primitive modulo h")
        residues = [
            (k, l)
            for k in range(h)
            for l in range(h)
            if (a * k + b * l) % h == 0
        ]
        if len(residues) != h:
            raise RuntimeError("primitive outside fibre has wrong residue size")
        lift_count = period // h
        for residue_k, residue_l in residues:
            for k_lift in range(lift_count):
                k = residue_k + h * k_lift
                for l_lift in range(lift_count):
                    yield k, residue_l + h * l_lift


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--outside-prime", type=int, required=True)
    parser.add_argument(
        "--normalizer-anchor-prime",
        type=int,
        required=True,
    )
    parser.add_argument("--base-anchor-primes", required=True)
    parser.add_argument(
        "--projected-anchors",
        default="",
        help="comma-separated prime:projection:residual triples",
    )
    parser.add_argument("--paired-projected-prime", type=int)
    parser.add_argument("--paired-shared-prime", type=int)
    parser.add_argument("--paired-projection-modulus", type=int)
    parser.add_argument("--paired-residual-modulus", type=int)
    parser.add_argument("--base-period", type=int, required=True)
    parser.add_argument(
        "--max-tensor-cells",
        type=int,
        default=20_000_000,
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if np is None:
        raise RuntimeError("NumPy is required for exact tensor replay")

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    if len(by_prime) != len(rows):
        raise RuntimeError("pool contains repeated fibre primes")
    base_primes = parse_primes(args.base_anchor_primes)
    projected_specs = parse_projected(args.projected_anchors)
    paired_values = (
        args.paired_projected_prime,
        args.paired_shared_prime,
        args.paired_projection_modulus,
        args.paired_residual_modulus,
    )
    if any(value is None for value in paired_values) != all(
        value is None for value in paired_values
    ):
        raise SystemExit("all paired-prime options must be supplied together")
    paired = paired_values[0] is not None
    listed_primes = (
        args.outside_prime,
        args.normalizer_anchor_prime,
        *base_primes,
        *(prime for prime, _projection, _residual in projected_specs),
        *(
            (
                args.paired_projected_prime,
                args.paired_shared_prime,
            )
            if paired
            else ()
        ),
    )
    if len(set(listed_primes)) != len(listed_primes):
        raise RuntimeError("all outside and anchor primes must be distinct")
    missing = set(listed_primes) - by_prime.keys()
    if missing:
        raise RuntimeError(f"primes are absent from the pool: {sorted(missing)}")

    outside = by_prime[args.outside_prime]
    normalizer = by_prime[args.normalizer_anchor_prime]
    base_rows = [by_prime[prime] for prime in base_primes]
    projected_rows = [
        (by_prime[prime], projection, residual)
        for prime, projection, residual in projected_specs
    ]
    period = args.base_period
    base_all = [outside, normalizer, *base_rows]
    if any(period % int(row["h"]) for row in base_all):
        raise RuntimeError("outside and base anchors must divide base period")
    if joint_index(outside, normalizer) != 1:
        raise RuntimeError(
            "outside/normalizer target map is not jointly surjective"
        )

    residual_moduli = []
    for row, projection, residual in projected_rows:
        h = int(row["h"])
        if h != projection * residual:
            raise RuntimeError("projected anchor h is not projection*residual")
        if math.gcd(h, period) != projection:
            raise RuntimeError(
                "projected modulus is not gcd(anchor h, base period)"
            )
        if math.gcd(period // projection, residual) != 1:
            raise RuntimeError("residual lift is not independent of base")
        if math.gcd(int(row["a"]), int(row["b"]), residual) != 1:
            raise RuntimeError("projected residual equation is not primitive")
        residual_moduli.append(residual)

    paired_rows = None
    if paired:
        paired_projected = by_prime[args.paired_projected_prime]
        paired_shared = by_prime[args.paired_shared_prime]
        projection = args.paired_projection_modulus
        residual = args.paired_residual_modulus
        if (
            int(paired_projected["h"]) != projection * residual
            or math.gcd(int(paired_projected["h"]), period) != projection
            or int(paired_shared["h"]) != residual
            or math.gcd(residual, period) != 1
        ):
            raise RuntimeError("paired projection/residual factorization fails")
        determinant = (
            int(paired_projected["a"]) * int(paired_shared["b"])
            - int(paired_shared["a"]) * int(paired_projected["b"])
        )
        if math.gcd(determinant, residual) != 1:
            raise RuntimeError("paired residual maps are not independent")
        residual_moduli.append(residual)
        paired_rows = (
            paired_projected,
            paired_shared,
            projection,
            residual,
            determinant,
        )
    if any(
        math.gcd(left, right) != 1
        for index, left in enumerate(residual_moduli)
        for right in residual_moduli[:index]
    ):
        raise RuntimeError("distinct residual components are not coprime")

    projection_rows = [
        (row, projection, residual, "single")
        for row, projection, residual in projected_rows
    ]
    if paired_rows:
        projection_rows.append(
            (
                paired_rows[0],
                paired_rows[2],
                paired_rows[3],
                "paired",
            )
        )
    shape = (
        *(int(row["h"]) for row in base_rows),
        *(projection for _row, projection, _residual, _kind
          in projection_rows),
    )
    tensor_cells = math.prod(shape)
    if tensor_cells > args.max_tensor_cells:
        raise RuntimeError(
            f"tensor has {tensor_cells} cells, above guard"
        )
    counts = np.zeros(shape, dtype=np.int64)
    outside_base_points = period * period // int(outside["h"])
    normalizer_covered = 0
    enumerated_outside_points = 0
    for k, l in primitive_fibre_points(outside, period):
        enumerated_outside_points += 1
        if (
            int(normalizer["a"]) * k + int(normalizer["b"]) * l
        ) % int(normalizer["h"]) == 0:
            normalizer_covered += 1
            continue
        index = (
            *(
                (
                    int(row["a"]) * k + int(row["b"]) * l
                ) % int(row["h"])
                for row in base_rows
            ),
            *(
                (
                    int(row["a"]) * k + int(row["b"]) * l
                ) % projection
                for row, projection, _residual, _kind
                in projection_rows
            ),
        )
        counts[index] += 1
    if enumerated_outside_points != outside_base_points:
        raise RuntimeError("outside fibre enumeration has wrong size")

    transformed = counts
    residual_denominator = 1
    for axis, (_row, projection, residual, kind) in reversed(
        list(enumerate(projection_rows, start=len(base_rows)))
    ):
        if kind == "single":
            matrix = np.full(
                (projection, projection),
                residual,
                dtype=np.int64,
            )
            np.fill_diagonal(matrix, residual - 1)
            residual_denominator *= residual
        else:
            matrix = np.full(
                (projection, projection),
                residual * (residual - 1),
                dtype=np.int64,
            )
            np.fill_diagonal(matrix, (residual - 1) ** 2)
            residual_denominator *= residual**2
        transformed = np.tensordot(
            matrix,
            transformed,
            axes=(1, axis),
        )
        transformed = np.moveaxis(transformed, 0, axis)
    for axis in range(len(base_rows)):
        transformed = (
            np.expand_dims(transformed.sum(axis=axis), axis)
            - transformed
        )
    flat_index = int(np.argmax(transformed))
    maximizing_targets = tuple(
        int(value)
        for value in np.unravel_index(flat_index, transformed.shape)
    )
    maximum_uncovered_score = int(transformed[maximizing_targets])
    total_outside_score = (
        period * period // int(outside["h"])
    ) * residual_denominator
    intersection = Fraction(
        total_outside_score - maximum_uncovered_score,
        residual_denominator * period * period,
    )
    anchors = [
        normalizer,
        *base_rows,
        *(row for row, _projection, _residual in projected_rows),
        *(
            [paired_rows[0], paired_rows[1]]
            if paired_rows
            else []
        ),
    ]
    result = {
        "schema": "projected_conditional_fibre_overlap_v1",
        "pool": str(args.pool),
        "row_count": len(rows),
        "outside_prime": args.outside_prime,
        "outside_row": {
            key: int(outside[key]) for key in ("p", "h", "a", "b")
        },
        "anchor_primes": [int(row["p"]) for row in anchors],
        "anchor_rows": [
            {key: int(row[key]) for key in ("p", "h", "a", "b")}
            for row in anchors
        ],
        "normalizer_anchor_prime": args.normalizer_anchor_prime,
        "base_anchor_primes": list(base_primes),
        "projected_anchors": [
            {
                "prime": int(row["p"]),
                "projection_modulus": projection,
                "residual_modulus": residual,
            }
            for row, projection, residual in projected_rows
        ],
        "paired_projected_prime": (
            args.paired_projected_prime if paired else None
        ),
        "paired_shared_prime": (
            args.paired_shared_prime if paired else None
        ),
        "paired_projection_modulus": (
            args.paired_projection_modulus if paired else None
        ),
        "paired_residual_modulus": (
            args.paired_residual_modulus if paired else None
        ),
        "paired_determinant": (
            paired_rows[4] if paired_rows else None
        ),
        "base_period": period,
        "base_cells": period * period,
        "outside_base_points": outside_base_points,
        "normalizer_covered_base_points": normalizer_covered,
        "uncovered_count_tensor_shape": list(shape),
        "uncovered_count_tensor_sum": int(counts.sum()),
        "uncovered_count_tensor_nonzero": int(np.count_nonzero(counts)),
        "residual_denominator": residual_denominator,
        "target_combinations": tensor_cells,
        "maximizing_targets": list(maximizing_targets),
        "maximum_uncovered_score": maximum_uncovered_score,
        "forced_intersection_density": fraction_payload(intersection),
        "argument": (
            "Normalize the outside fibre and the named anchor to target zero. "
            "On the base plane, points covered by that anchor are removed. "
            "Complement transforms enumerate every remaining base-anchor "
            "target, while residual matrices give the exact uncovered count "
            "for each projected target. The shared pair uses two independent "
            "residual lines. Therefore the maximum uncovered score gives the "
            "minimum intersection of the outside fibre with this anchor "
            "subunion."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"outside={args.outside_prime} anchors={len(anchors)} "
        f"base_points={outside_base_points} tensor_cells={tensor_cells} "
        f"maximum_uncovered_score={maximum_uncovered_score} "
        f"intersection={intersection}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
