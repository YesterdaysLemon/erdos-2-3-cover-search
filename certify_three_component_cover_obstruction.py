#!/usr/bin/env python3
"""Certify a phase-independent obstruction for a three-component line cover.

Rows are grouped by the exact squarefree moduli x, y, z, xy, xz, yz, xyz.
The proof combines:

* the Alon--Furedi lower bound for holes left by the x-only lines;
* the Jamison--Brouwer--Schrijver affine-plane cover threshold; and
* a pair-incidence bound on points meeting many projected cross-lines.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path


def direction(row: dict, prime: int) -> tuple[int, int]:
    a = int(row["a"]) % prime
    b = int(row["b"]) % prime
    if a:
        inverse = pow(a, -1, prime)
        return 1, b * inverse % prime
    if b:
        return 0, 1
    raise RuntimeError(
        f"row p={row['p']} is degenerate modulo {prime}"
    )


def alon_furedi_grid_bound(prime: int, degree: int) -> int:
    required_sum = 2 * prime - degree
    products = [
        first * (required_sum - first)
        for first in range(1, prime + 1)
        if 1 <= required_sum - first <= prime
    ]
    if not products:
        return 0
    return min(products)


def high_incidence_bound(
    rows: list[dict],
    prime: int,
    threshold: int,
) -> dict:
    counts = Counter(direction(row, prime) for row in rows)
    parallel_pairs = sum(
        count * (count - 1) // 2 for count in counts.values()
    )
    all_pairs = len(rows) * (len(rows) - 1) // 2
    nonparallel_pairs = all_pairs - parallel_pairs
    pair_incidence_upper = (
        nonparallel_pairs + prime * parallel_pairs
    )
    threshold_pairs = threshold * (threshold - 1) // 2
    point_upper = pair_incidence_upper // threshold_pairs
    return {
        "row_count": len(rows),
        "threshold": threshold,
        "direction_counts": {
            f"{first},{second}": count
            for (first, second), count in sorted(counts.items())
        },
        "parallel_pairs": parallel_pairs,
        "nonparallel_pairs": nonparallel_pairs,
        "pair_incidence_upper": pair_incidence_upper,
        "pairs_per_high_incidence_point": threshold_pairs,
        "high_incidence_point_upper_bound": point_upper,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--x-prime", type=int, required=True)
    parser.add_argument("--y-prime", type=int, required=True)
    parser.add_argument("--z-prime", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    x, y, z = args.x_prime, args.y_prime, args.z_prime
    if len({x, y, z}) != 3:
        raise SystemExit("component primes must be distinct")
    payload = json.loads(args.pool.read_text())
    allowed = (x, y, z, x * y, x * z, y * z, x * y * z)
    groups = {
        modulus: [
            row
            for row in payload["choices"]
            if int(row["h"]) == modulus
        ]
        for modulus in allowed
    }
    actual = {int(row["h"]) for row in payload["choices"]}
    unsupported = sorted(actual - set(allowed))
    if unsupported:
        raise RuntimeError(f"unsupported moduli: {unsupported}")
    if any(
        int(row["target_modulus"]) != 1
        for row in payload["choices"]
    ):
        raise RuntimeError("all phases must be freely selectable")

    x_rows = groups[x]
    y_rows = groups[y]
    z_rows = groups[z]
    xy_rows = groups[x * y]
    xz_rows = groups[x * z]
    yz_rows = groups[y * z]
    xyz_rows = groups[x * y * z]
    y_direction_counts = Counter(
        direction(row, y) for row in [*y_rows, *xy_rows]
    )
    z_direction_counts = Counter(
        direction(row, z)
        for row in [*z_rows, *xz_rows, *yz_rows, *xyz_rows]
    )
    x_direction_counts = Counter(
        direction(row, x) for row in x_rows
    )

    x_cover_threshold = 2 * x - 1
    y_cover_threshold = 2 * y - 1
    z_cover_threshold = 2 * z - 1
    if len(x_rows) >= x_cover_threshold:
        raise RuntimeError("x-only row count is outside this proof regime")
    if max(x_direction_counts.values(), default=0) >= x:
        raise RuntimeError("x-only rows may contain a parallel class")
    if max(y_direction_counts.values(), default=0) >= y:
        raise RuntimeError("y rows may contain a parallel class")
    if max(z_direction_counts.values(), default=0) >= z:
        raise RuntimeError("z rows may contain a parallel class")

    x_hole_lower = alon_furedi_grid_bound(x, len(x_rows))
    xz_threshold = (
        z_cover_threshold
        - len(z_rows)
        - len(yz_rows)
        - len(xyz_rows)
    )
    xy_threshold = y_cover_threshold - len(y_rows)
    if xz_threshold < 2 or xy_threshold < 2:
        raise RuntimeError("high-incidence thresholds are too small")
    xz_bound = high_incidence_bound(
        xz_rows,
        x,
        xz_threshold,
    )
    xy_bound = high_incidence_bound(
        xy_rows,
        x,
        xy_threshold,
    )
    x_hole_upper_if_cover = (
        xz_bound["high_incidence_point_upper_bound"]
        + xy_bound["high_incidence_point_upper_bound"]
    )
    proved = x_hole_lower > x_hole_upper_if_cover

    result = {
        "pool": str(args.pool),
        "components": {"x": x, "y": y, "z": z},
        "row_counts": {
            str(modulus): len(groups[modulus])
            for modulus in allowed
        },
        "affine_cover_thresholds": {
            str(x): x_cover_threshold,
            str(y): y_cover_threshold,
            str(z): z_cover_threshold,
        },
        "maximum_direction_multiplicities": {
            "x_only": max(x_direction_counts.values(), default=0),
            "y_and_xy": max(y_direction_counts.values(), default=0),
            "z_incident": max(z_direction_counts.values(), default=0),
        },
        "alon_furedi": {
            "grid": [x, x],
            "degree": len(x_rows),
            "hole_lower_bound": x_hole_lower,
        },
        "xz_high_incidence": xz_bound,
        "xy_high_incidence": xy_bound,
        "x_hole_upper_bound_if_cover_exists": (
            x_hole_upper_if_cover
        ),
        "proved_no_cover": proved,
        "scope": "all phase assignments for the supplied squarefree pool",
        "argument": (
            "Every x-hole must either meet enough xz projections to reach "
            "the nonparallel z-cover threshold even after all yz and xyz "
            "rows activate, or meet enough xy projections to let the y-lines "
            "cover. Pair-incidence bounds limit these two exceptional x-point "
            "sets. Their union is smaller than the Alon--Furedi lower bound "
            "for x-holes."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"x_holes_lb={x_hole_lower} "
        f"xz_threshold={xz_threshold} "
        f"xz_high_points_ub="
        f"{xz_bound['high_incidence_point_upper_bound']} "
        f"xy_threshold={xy_threshold} "
        f"xy_high_points_ub="
        f"{xy_bound['high_incidence_point_upper_bound']} "
        f"x_holes_ub_if_cover={x_hole_upper_if_cover}",
        flush=True,
    )
    print(
        "PROVED no three-component cover"
        if proved
        else "NO obstruction from these bounds",
        flush=True,
    )
    return 0 if proved else 1


if __name__ == "__main__":
    raise SystemExit(main())
