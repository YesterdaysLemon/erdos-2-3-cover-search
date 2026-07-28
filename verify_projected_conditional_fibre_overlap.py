#!/usr/bin/env python3
"""Independent replay of a projected conditional-fibre certificate."""

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


def primitive_kernel_residues(row: dict) -> object:
    """Enumerate the target-zero kernel modulo h as one cyclic subgroup."""
    h = int(row["h"])
    a = int(row["a"])
    b = int(row["b"])
    if math.gcd(a, b, h) != 1:
        raise RuntimeError("outside row is not primitive modulo h")
    order = math.lcm(
        h // math.gcd(b, h),
        h // math.gcd(a, h),
    )
    if order != h:
        raise RuntimeError("primitive kernel generator has wrong order")
    for multiplier in range(h):
        yield (
            (multiplier * b) % h,
            (-multiplier * a) % h,
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
    normalizer = by_prime[int(cert["normalizer_anchor_prime"])]
    base_rows = [
        by_prime[int(prime)] for prime in cert["base_anchor_primes"]
    ]
    projected_rows = [
        (
            by_prime[int(record["prime"])],
            int(record["projection_modulus"]),
            int(record["residual_modulus"]),
        )
        for record in cert["projected_anchors"]
    ]
    paired = cert.get("paired_projected_prime") is not None
    paired_rows = None
    if paired:
        paired_rows = (
            by_prime[int(cert["paired_projected_prime"])],
            by_prime[int(cert["paired_shared_prime"])],
            int(cert["paired_projection_modulus"]),
            int(cert["paired_residual_modulus"]),
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
    period = int(cert["base_period"])
    shape = (
        *(int(row["h"]) for row in base_rows),
        *(projection for _row, projection, _residual in projected_rows),
        *([paired_rows[2]] if paired_rows else []),
    )
    counts = np.zeros(shape, dtype=np.int64)
    outside_points = 0
    normalizer_covered = 0
    def record_outside_point(k: int, l: int) -> None:
        nonlocal outside_points, normalizer_covered
        outside_points += 1
        if (
            int(normalizer["a"]) * k + int(normalizer["b"]) * l
        ) % int(normalizer["h"]) == 0:
            normalizer_covered += 1
            return
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
                for row, projection, _residual in projected_rows
            ),
            *(
                [
                    (
                        int(paired_rows[0]["a"]) * k
                        + int(paired_rows[0]["b"]) * l
                    ) % paired_rows[2]
                ]
                if paired_rows
                else []
            ),
        )
        counts[index] += 1

    if period * period <= 10_000_000:
        enumeration_mode = "complete_base_plane"
        # Deliberately scan the complete base plane, independently of the
        # certificate generator's primitive-fibre parametrization.
        for k in range(period):
            for l in range(period):
                if (
                    int(outside["a"]) * k + int(outside["b"]) * l
                ) % int(outside["h"]):
                    continue
                record_outside_point(k, l)
    else:
        enumeration_mode = "cyclic_kernel_lifts"
        # The vector (b,-a) generates the target-zero kernel modulo h when
        # gcd(a,b,h)=1.  Enumerate that h-cycle and lift it through the base
        # period.  This is independent of the generator's choice of an
        # invertible coefficient and avoids an unnecessary h-by-h scan.
        outside_h = int(outside["h"])
        lift_count = period // outside_h
        for residue_k, residue_l in primitive_kernel_residues(outside):
            for k_lift in range(lift_count):
                k = residue_k + outside_h * k_lift
                for l_lift in range(lift_count):
                    record_outside_point(
                        k,
                        residue_l + outside_h * l_lift,
                    )

    matrices = []
    residual_denominator = 1
    for _row, projection, residual in projected_rows:
        matrix = np.full(
            (projection, projection),
            residual,
            dtype=np.int64,
        )
        np.fill_diagonal(matrix, residual - 1)
        matrices.append(matrix)
        residual_denominator *= residual
    if paired_rows:
        projection = paired_rows[2]
        residual = paired_rows[3]
        matrix = np.full(
            (projection, projection),
            residual * (residual - 1),
            dtype=np.int64,
        )
        np.fill_diagonal(matrix, (residual - 1) ** 2)
        matrices.append(matrix)
        residual_denominator *= residual**2

    transformed = counts
    # Replay the separable map with explicit einsum axis labels, rather than
    # the generator's tensordot/moveaxis implementation.
    for offset, matrix in reversed(list(enumerate(matrices))):
        axis = len(base_rows) + offset
        transformed = np.einsum(
            matrix,
            [axis, transformed.ndim],
            transformed,
            list(range(transformed.ndim)),
            list(range(axis))
            + [transformed.ndim]
            + list(range(axis + 1, transformed.ndim)),
            optimize=True,
        )
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
    total_outside_score = outside_points * residual_denominator
    intersection = Fraction(
        total_outside_score - maximum_uncovered_score,
        residual_denominator * period * period,
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
    structural = (
        pair_index(outside, normalizer) == 1
        and all(period % int(row["h"]) == 0
                for row in [outside, normalizer, *base_rows])
        and all(
            int(row["h"]) == projection * residual
            and math.gcd(int(row["h"]), period) == projection
            and math.gcd(period // projection, residual) == 1
            for row, projection, residual in projected_rows
        )
    )
    if paired_rows:
        determinant = (
            int(paired_rows[0]["a"]) * int(paired_rows[1]["b"])
            - int(paired_rows[1]["a"]) * int(paired_rows[0]["b"])
        )
        structural = (
            structural
            and int(paired_rows[0]["h"])
            == paired_rows[2] * paired_rows[3]
            and int(paired_rows[1]["h"]) == paired_rows[3]
            and math.gcd(determinant, paired_rows[3]) == 1
            and int(cert["paired_determinant"]) == determinant
        )
    verified = (
        cert.get("schema") == "projected_conditional_fibre_overlap_v1"
        and str(args.pool) == cert["pool"]
        and len(by_prime) == len(rows)
        and rows_match
        and structural
        and int(cert["base_cells"]) == period * period
        and int(cert["outside_base_points"]) == outside_points
        and int(cert["normalizer_covered_base_points"])
        == normalizer_covered
        and list(cert["uncovered_count_tensor_shape"]) == list(shape)
        and int(cert["uncovered_count_tensor_sum"]) == int(counts.sum())
        and int(cert["uncovered_count_tensor_nonzero"])
        == int(np.count_nonzero(counts))
        and int(cert["residual_denominator"]) == residual_denominator
        and int(cert["target_combinations"]) == math.prod(shape)
        and list(cert["maximizing_targets"]) == list(maximizing_targets)
        and int(cert["maximum_uncovered_score"])
        == maximum_uncovered_score
        and read_fraction(cert["forced_intersection_density"])
        == intersection
    )
    result = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "outside_prime": int(outside["p"]),
        "anchor_primes": [int(row["p"]) for row in anchors],
        "base_period": period,
        "base_enumeration_mode": enumeration_mode,
        "outside_base_points": outside_points,
        "target_combinations": math.prod(shape),
        "maximum_uncovered_score": maximum_uncovered_score,
        "forced_intersection_density": {
            "numerator": intersection.numerator,
            "denominator": intersection.denominator,
        },
        "verified": verified,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"outside={outside['p']} base_points={outside_points} "
        f"targets={math.prod(shape)} intersection={intersection} "
        f"verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
