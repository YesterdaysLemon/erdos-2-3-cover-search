#!/usr/bin/env python3
"""Independently verify a lifted component-cell obstruction."""

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
    raise RuntimeError("degenerate row")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    parent = json.loads(Path(certificate["parent_pool"]).read_text())
    child = json.loads(
        Path(certificate["conditioned_pool"]).read_text()
    )
    obstruction = json.loads(
        Path(certificate["conditioned_obstruction"]).read_text()
    )
    q = int(certificate["component_prime"])
    pure = [row for row in parent["choices"] if int(row["h"]) == q]
    residual = [
        row for row in parent["choices"] if int(row["h"]) != q
    ]
    child_rows = {
        int(row["p"]): row for row in child["choices"]
    }
    structural_match = (
        sorted(child_rows)
        == sorted(int(row["p"]) for row in residual)
    )
    if structural_match:
        for row in residual:
            h = int(row["h"])
            a = int(row["a"]) % h
            b = int(row["b"]) % h
            common = math.gcd(q * a, q * b, h)
            modulus = h // common
            expected = (
                modulus,
                q * a // common % modulus,
                q * b // common % modulus,
            )
            actual_row = child_rows[int(row["p"])]
            actual = (
                int(actual_row["h"]),
                int(actual_row["a"]),
                int(actual_row["b"]),
            )
            if expected != actual:
                structural_match = False
                break
    directions = Counter(
        normalized_direction(row, q) for row in pure
    )
    maximum_direction = max(directions.values(), default=0)
    threshold = 2 * q - 1
    checks = {
        "child_period_matches": (
            int(child["cell_condition"]["period"]) == q
        ),
        "conditioned_obstruction_proved": bool(
            obstruction["proved_no_cover"]
        ),
        "conditioned_geometry_matches": structural_match,
        "pure_count_matches": (
            len(pure) == int(certificate["pure_component_rows"])
        ),
        "residual_count_matches": (
            len(residual) == int(certificate["residual_rows"])
        ),
        "cover_threshold_matches": (
            threshold == int(certificate["affine_cover_threshold"])
        ),
        "too_few_for_nonparallel_cover": len(pure) < threshold,
        "no_parallel_class": maximum_direction < q,
        "maximum_direction_matches": (
            maximum_direction
            == int(certificate["maximum_direction_count"])
        ),
        "certificate_claims_no_cover": bool(
            certificate["proved_no_cover"]
        ),
    }
    passed = all(checks.values())
    result = {
        "certificate": str(args.certificate),
        "checks": checks,
        "recomputed": {
            "pure_rows": len(pure),
            "residual_rows": len(residual),
            "cover_threshold": threshold,
            "maximum_direction_count": maximum_direction,
        },
        "verified": passed,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS" if passed else "FAIL", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
