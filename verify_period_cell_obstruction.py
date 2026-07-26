#!/usr/bin/env python3
"""Independently replay a period-cell obstruction certificate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def geometry(row: dict, period: int) -> tuple[int, int, int]:
    h = int(row["h"])
    a = int(row["a"]) % h
    b = int(row["b"]) % h
    divisor = math.gcd(period * a, period * b, h)
    reduced = h // divisor
    return (
        reduced,
        period * a // divisor % reduced,
        period * b // divisor % reduced,
    )


def references_pool(obstruction: dict, pool: Path) -> bool:
    wanted = pool.name
    return any(
        value and Path(value).name == wanted
        for value in (
            obstruction.get("parent_pool"),
            obstruction.get("grand_pool"),
            obstruction.get("conditioned_pool"),
            obstruction.get("source_pool"),
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    parent_path = Path(certificate["parent_pool"])
    child_path = Path(certificate["conditioned_pool"])
    obstruction_path = Path(certificate["conditioned_obstruction"])
    parent = json.loads(parent_path.read_text())
    child = json.loads(child_path.read_text())
    obstruction = json.loads(obstruction_path.read_text())
    period = int(certificate["period"])

    pure = []
    residual = []
    for row in parent["choices"]:
        (pure if geometry(row, period)[0] == 1 else residual).append(row)

    child_rows = {int(row["p"]): row for row in child["choices"]}
    exact_residual_set = sorted(child_rows) == sorted(
        int(row["p"]) for row in residual
    )
    structural_match = exact_residual_set
    if structural_match:
        for row in residual:
            expected = geometry(row, period)
            actual_row = child_rows[int(row["p"])]
            h = int(actual_row["h"])
            actual = (
                h,
                int(actual_row["a"]) % h,
                int(actual_row["b"]) % h,
            )
            if expected != actual:
                structural_match = False
                break

    capacities = []
    pure_moduli_divide_period = True
    for row in pure:
        h = int(row["h"])
        if period % h:
            pure_moduli_divide_period = False
            capacity = period * period
        else:
            coefficient_gcd = math.gcd(
                int(row["a"]) % h,
                int(row["b"]) % h,
                h,
            )
            capacity = period * period * coefficient_gcd // h
        capacities.append((int(row["p"]), h, capacity))

    capacity_sum = sum(item[2] for item in capacities)
    certified_capacities = sorted(
        (
            int(item["p"]),
            int(item["h"]),
            int(item["cell_capacity_upper_bound"]),
        )
        for item in certificate["pure_period_rows"]
    )
    child_fixed = sorted(
        int(prime) for prime in child.get("fixed_inactive_primes", ())
    )
    pure_primes = sorted(int(row["p"]) for row in pure)

    checks = {
        "child_period_matches": (
            int(child.get("cell_condition", {}).get("period", 0)) == period
        ),
        "conditioned_obstruction_proved": bool(
            obstruction.get("proved_no_cover")
        ),
        "obstruction_references_child": references_pool(
            obstruction, child_path
        ),
        "child_has_no_full_cell_rows": not child.get("full_cell_primes"),
        "exact_residual_prime_set": exact_residual_set,
        "conditioned_geometry_matches": structural_match,
        "pure_moduli_divide_period": pure_moduli_divide_period,
        "pure_capacities_match": sorted(capacities) == certified_capacities,
        "pure_count_matches": (
            len(pure) == int(certificate["pure_period_row_count"])
        ),
        "residual_count_matches": (
            len(residual) == int(certificate["residual_row_count"])
        ),
        "child_fixed_inactive_is_exact_pure_set": child_fixed == pure_primes,
        "period_cell_count_matches": (
            period * period == int(certificate["period_cell_count"])
        ),
        "capacity_sum_matches": (
            capacity_sum
            == int(certificate["pure_union_capacity_upper_bound"])
        ),
        "union_capacity_forces_missed_cell": capacity_sum < period * period,
        "missed_cell_bound_matches": (
            period * period - capacity_sum
            == int(certificate["forced_missed_cells_lower_bound"])
        ),
        "certificate_claims_no_cover": bool(
            certificate.get("proved_no_cover")
        ),
    }
    passed = all(checks.values())
    result = {
        "certificate": str(args.certificate),
        "checks": checks,
        "recomputed": {
            "parent_rows": len(parent["choices"]),
            "pure_rows": len(pure),
            "residual_rows": len(residual),
            "period_cells": period * period,
            "pure_union_capacity_upper_bound": capacity_sum,
            "forced_missed_cells_lower_bound": period * period - capacity_sum,
        },
        "verified": passed,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS" if passed else "FAIL", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
