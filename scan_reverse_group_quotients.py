#!/usr/bin/env python3
"""Rank small arithmetic-matched quotient families without cell enumeration.

For each declared unimodular basis and pair of quotient periods, this scans
the authenticated source rows, keeps precisely the affine predicates that
descend to the quotient, and computes:

* exact reciprocal density;
* a maximum-weight forced-intersection Hunter-forest bound; and
* the exact independent-random-phase first moment for unrestricted rows.

The scan is a finite-family preflight.  A surviving record is not a cover.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

from build_axis_layered_pool import sha256
from inventory_reverse_group_lines import (
    maximum_forced_overlap_forest,
    transformed_coefficients,
)


DEFAULT_BASES = [
    ("id", (1, 0), (0, 1)),
    ("d11", (1, 1), (0, 1)),
    ("d1m1", (1, -1), (0, 1)),
    ("d12", (1, 2), (0, 1)),
    ("d1m2", (1, -2), (0, 1)),
    ("d13", (1, 3), (0, 1)),
    ("d1m3", (1, -3), (0, 1)),
    ("d21", (2, 1), (1, 1)),
    ("d2m1", (2, -1), (1, 0)),
    ("d31", (3, 1), (-1, 0)),
    ("d3m1", (3, -1), (1, 0)),
    ("d32", (3, 2), (1, 1)),
    ("d3m2", (3, -2), (-1, 1)),
    ("d41", (4, 1), (-1, 0)),
    ("d4m1", (4, -1), (1, 0)),
    ("d43", (4, 3), (1, 1)),
]

DEFAULT_PERIODS = [
    420,
    630,
    720,
    840,
    1260,
    1680,
    2520,
    3360,
    3780,
    3960,
    4620,
    5040,
    7560,
    10080,
    15120,
]


def exact_fields(prefix: str, value: Fraction) -> dict:
    return {
        f"{prefix}_numerator": value.numerator,
        f"{prefix}_denominator": value.denominator,
        f"{prefix}_decimal": float(value),
    }


def scan(
    payload: dict,
    source_path: Path,
    bases: list[tuple[str, tuple[int, int], tuple[int, int]]],
    periods: list[int],
    cell_limit: int,
    minimum_density: Fraction,
) -> dict:
    if json.loads(source_path.read_text()) != payload:
        raise ValueError("source payload differs from authenticated file")
    if cell_limit < 1:
        raise ValueError("cell limit must be positive")
    if not bases or len({basis[0] for basis in bases}) != len(bases):
        raise ValueError("basis names must be nonempty and unique")
    periods = sorted(set(int(period) for period in periods))
    if not periods or periods[0] < 1:
        raise ValueError("periods must be positive")
    rows = payload["choices"]
    primes = [int(row["p"]) for row in rows]
    if len(set(primes)) != len(primes):
        raise ValueError("duplicate source prime")

    records = []
    configuration_count = 0
    density_at_least_one_count = 0
    finite_noncover_count = 0
    for name, direction, transverse in bases:
        determinant = (
            direction[0] * transverse[1]
            - direction[1] * transverse[0]
        )
        if abs(determinant) != 1:
            raise ValueError(f"basis {name} is not unimodular")
        transformed = []
        for row in rows:
            h, a, b = transformed_coefficients(
                row,
                direction,
                transverse,
            )
            target_modulus = int(row.get("target_modulus", 1))
            if target_modulus < 1 or h % target_modulus:
                raise ValueError(f"invalid target restriction p={row['p']}")
            transformed.append(
                {
                    "p": int(row["p"]),
                    "h": h,
                    "a": a,
                    "b": b,
                    "ord_a": h // math.gcd(a, h),
                    "ord_b": h // math.gcd(b, h),
                    "unrestricted": target_modulus == 1,
                }
            )
        for width in periods:
            for height in periods:
                cell_count = width * height
                if cell_count > cell_limit:
                    continue
                configuration_count += 1
                descending = [
                    row
                    for row in transformed
                    if width % row["ord_a"] == 0
                    and height % row["ord_b"] == 0
                ]
                density = sum(
                    (Fraction(1, row["h"]) for row in descending),
                    Fraction(),
                )
                if density >= 1:
                    density_at_least_one_count += 1
                line_types = [
                    {
                        "h": row["h"],
                        "a": row["a"],
                        "b": row["b"],
                        "primes": [row["p"]],
                    }
                    for row in descending
                ]
                forest, overlap = maximum_forced_overlap_forest(
                    line_types,
                    cell_count,
                )
                union_upper = density - overlap
                noncover = union_upper < 1
                if noncover:
                    finite_noncover_count += 1
                if density < minimum_density:
                    continue
                unrestricted = [
                    row for row in descending if row["unrestricted"]
                ]
                avoidance = Fraction(1)
                for row in unrestricted:
                    avoidance *= Fraction(row["h"] - 1, row["h"])
                expected = cell_count * avoidance
                if expected < 1 and noncover:
                    raise AssertionError(
                        "first-moment existence conflicts with "
                        "Hunter noncover bound"
                    )
                records.append(
                    {
                        "basis": {
                            "name": name,
                            "direction": list(direction),
                            "transverse": list(transverse),
                            "determinant": determinant,
                        },
                        "group": {
                            "width": width,
                            "height": height,
                            "cell_count": cell_count,
                        },
                        "descending_row_count": len(descending),
                        "descending_rows": descending,
                        **exact_fields("raw_density", density),
                        "forced_overlap_forest": forest,
                        **exact_fields(
                            "forced_overlap_density",
                            overlap,
                        ),
                        **exact_fields(
                            "phase_independent_union_upper_bound",
                            union_upper,
                        ),
                        "finite_group_cover_impossible_by_density_overlap": (
                            noncover
                        ),
                        "unrestricted_row_count": len(unrestricted),
                        **exact_fields(
                            "expected_uncovered_cell_count",
                            expected,
                        ),
                        "first_moment_cover_exists": expected < 1,
                    }
                )
    records.sort(
        key=lambda record: (
            -Fraction(
                record["phase_independent_union_upper_bound_numerator"],
                record["phase_independent_union_upper_bound_denominator"],
            ),
            record["group"]["cell_count"],
            record["basis"]["name"],
            record["group"]["width"],
            record["group"]["height"],
        )
    )
    return {
        "schema": "reverse_group_quotient_scan_v2",
        "source": str(source_path),
        "source_sha256": sha256(source_path),
        "periods": periods,
        "cell_limit": cell_limit,
        **exact_fields("minimum_emitted_density", minimum_density),
        "basis_count": len(bases),
        "configuration_count": configuration_count,
        "density_below_one_count": (
            configuration_count - density_at_least_one_count
        ),
        "density_at_least_one_count": density_at_least_one_count,
        "finite_noncover_count": finite_noncover_count,
        "emitted_record_count": len(records),
        "records": records,
        "claim": {
            "scan_arithmetic_exact": True,
            "surviving_record_is_a_cover": False,
            "integer_m_found": False,
        },
        "scope": (
            "exact no-cell-enumeration preflight over only the declared "
            "bases, periods, source rows, and cell limit"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--period", action="append", type=int)
    parser.add_argument(
        "--basis",
        action="append",
        nargs=5,
        metavar=("NAME", "D0", "D1", "T0", "T1"),
    )
    parser.add_argument("--cell-limit", type=int, default=100_000_000)
    parser.add_argument("--minimum-density", default="1")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    bases = (
        [
            (
                raw[0],
                (int(raw[1]), int(raw[2])),
                (int(raw[3]), int(raw[4])),
            )
            for raw in args.basis
        ]
        if args.basis
        else DEFAULT_BASES
    )
    payload = json.loads(args.pool.read_text())
    result = scan(
        payload,
        args.pool,
        bases,
        args.period or DEFAULT_PERIODS,
        args.cell_limit,
        Fraction(args.minimum_density),
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"configurations={result['configuration_count']} "
        f"density_survivors={result['density_at_least_one_count']} "
        f"finite_noncovers={result['finite_noncover_count']} "
        f"emitted={result['emitted_record_count']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
