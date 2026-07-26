#!/usr/bin/env python3
"""Independently replay a squarefree coarse-reduction certificate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    source = json.loads(Path(certificate["pool"]).read_text())
    coarse = json.loads(Path(certificate["coarse_pool"]).read_text())
    period = int(certificate["coarse_period"])
    first, second = map(int, certificate["high_digit_primes"])
    groups = {1: [], first: [], second: []}
    unsupported = []
    for row in source["choices"]:
        h = int(row["h"])
        refinement = h // math.gcd(h, period)
        if refinement in groups:
            groups[refinement].append(row)
        else:
            unsupported.append(int(row["p"]))

    first_holes = first * first - len(groups[first]) * first
    second_holes = second * second - len(groups[second]) * second
    source_coarse_primes = sorted(
        int(row["p"]) for row in groups[1]
    )
    materialized_primes = sorted(
        int(row["p"]) for row in coarse["choices"]
    )
    checks = {
        "no_unsupported_rows": not unsupported,
        "coarse_rows_divide_period": all(
            period % int(row["h"]) == 0 for row in groups[1]
        ),
        "coarse_pool_exactly_matches": (
            source_coarse_primes == materialized_primes
        ),
        "first_factor_has_holes": first_holes > 0,
        "second_factor_has_holes": second_holes > 0,
        "counts_match": (
            len(groups[1])
            == int(certificate["row_counts"]["coarse"])
            and len(groups[first])
            == int(certificate["row_counts"][str(first)])
            and len(groups[second])
            == int(certificate["row_counts"][str(second)])
        ),
        "hole_bounds_match": (
            first_holes
            == int(
                certificate["refinement_hole_lower_bounds"][
                    str(first)
                ]
            )
            and second_holes
            == int(
                certificate["refinement_hole_lower_bounds"][
                    str(second)
                ]
            )
        ),
        "product_bound_matches": (
            first_holes * second_holes
            == int(certificate["product_hole_lower_bound"])
        ),
        "certificate_claims_equivalence": bool(
            certificate["proved_equivalence"]
        ),
    }
    passed = all(checks.values())
    result = {
        "certificate": str(args.certificate),
        "checks": checks,
        "recomputed_counts": {
            "coarse": len(groups[1]),
            str(first): len(groups[first]),
            str(second): len(groups[second]),
        },
        "recomputed_hole_lower_bounds": {
            str(first): first_holes,
            str(second): second_holes,
        },
        "verified": passed,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS" if passed else "FAIL", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
