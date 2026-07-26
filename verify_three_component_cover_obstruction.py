#!/usr/bin/env python3
"""Independent arithmetic replay of a three-component obstruction."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path


def normalized_direction(row: dict, prime: int) -> tuple[int, int]:
    a = int(row["a"]) % prime
    b = int(row["b"]) % prime
    if a:
        return 1, b * pow(a, -1, prime) % prime
    if b:
        return 0, 1
    raise RuntimeError("degenerate component line")


def incidence(rows: list[dict], prime: int, threshold: int) -> int:
    directions = Counter(
        normalized_direction(row, prime) for row in rows
    )
    parallel = sum(
        count * (count - 1) // 2
        for count in directions.values()
    )
    all_pairs = len(rows) * (len(rows) - 1) // 2
    pair_upper = all_pairs - parallel + prime * parallel
    return pair_upper // (threshold * (threshold - 1) // 2)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    pool = json.loads(Path(certificate["pool"]).read_text())
    x = int(certificate["components"]["x"])
    y = int(certificate["components"]["y"])
    z = int(certificate["components"]["z"])
    moduli = (x, y, z, x * y, x * z, y * z, x * y * z)
    groups = {
        modulus: [
            row
            for row in pool["choices"]
            if int(row["h"]) == modulus
        ]
        for modulus in moduli
    }
    x_rows = groups[x]
    y_rows = groups[y]
    z_rows = groups[z]
    xy_rows = groups[x * y]
    xz_rows = groups[x * z]
    yz_rows = groups[y * z]
    xyz_rows = groups[x * y * z]

    required_sum = 2 * x - len(x_rows)
    af_candidates = [
        first * (required_sum - first)
        for first in range(1, x + 1)
        if 1 <= required_sum - first <= x
    ]
    af_bound = min(af_candidates)
    xz_threshold = (
        (2 * z - 1)
        - len(z_rows)
        - len(yz_rows)
        - len(xyz_rows)
    )
    xy_threshold = (2 * y - 1) - len(y_rows)
    xz_points = incidence(xz_rows, x, xz_threshold)
    xy_points = incidence(xy_rows, x, xy_threshold)
    recomputed_upper = xz_points + xy_points
    y_max_direction = max(
        Counter(
            normalized_direction(row, y)
            for row in [*y_rows, *xy_rows]
        ).values()
    )
    z_max_direction = max(
        Counter(
            normalized_direction(row, z)
            for row in [
                *z_rows,
                *xz_rows,
                *yz_rows,
                *xyz_rows,
            ]
        ).values()
    )
    checks = {
        "all_moduli_supported": all(
            int(row["h"]) in moduli for row in pool["choices"]
        ),
        "row_counts_match": all(
            len(groups[modulus])
            == int(certificate["row_counts"][str(modulus)])
            for modulus in moduli
        ),
        "y_has_no_parallel_class": y_max_direction < y,
        "z_has_no_parallel_class": z_max_direction < z,
        "alon_furedi_bound_matches": (
            af_bound
            == int(
                certificate["alon_furedi"]["hole_lower_bound"]
            )
        ),
        "xz_threshold_matches": (
            xz_threshold
            == int(certificate["xz_high_incidence"]["threshold"])
        ),
        "xy_threshold_matches": (
            xy_threshold
            == int(certificate["xy_high_incidence"]["threshold"])
        ),
        "xz_point_bound_matches": (
            xz_points
            == int(
                certificate["xz_high_incidence"][
                    "high_incidence_point_upper_bound"
                ]
            )
        ),
        "xy_point_bound_matches": (
            xy_points
            == int(
                certificate["xy_high_incidence"][
                    "high_incidence_point_upper_bound"
                ]
            )
        ),
        "hole_bounds_contradict": af_bound > recomputed_upper,
        "certificate_claims_no_cover": bool(
            certificate["proved_no_cover"]
        ),
    }
    passed = all(checks.values())
    result = {
        "certificate": str(args.certificate),
        "checks": checks,
        "recomputed": {
            "alon_furedi_hole_lower_bound": af_bound,
            "xz_threshold": xz_threshold,
            "xz_high_incidence_points": xz_points,
            "xy_threshold": xy_threshold,
            "xy_high_incidence_points": xy_points,
            "x_hole_upper_if_cover": recomputed_upper,
        },
        "verified": passed,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS" if passed else "FAIL", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
