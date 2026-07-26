#!/usr/bin/env python3
"""Replay the 18-case translation symmetry split for density obstructions."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


ANCHORS = (19, 37, 73, 97, 109, 15121)
EXPECTED_MODULI = {
    19: 3,
    37: 3,
    73: 3,
    97: 4,
    109: 9,
    15121: 9,
}


def q_part(value: int, prime: int) -> int:
    part = 1
    while value % prime == 0:
        value //= prime
        part *= prime
    return part


def shifted(row: dict, target: int, dk: int, dl: int) -> int:
    return (
        target + int(row["a"]) * dk + int(row["b"]) * dl
    ) % int(row["h"])


def result_path(
    pattern: str,
    orbit: str,
    r109: int,
    r15121: int,
) -> Path:
    return Path(
        pattern.format(
            orbit=orbit,
            r109=r109,
            r15121=r15121,
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--result-pattern", required=True)
    parser.add_argument("--density-prime", type=int, default=2)
    parser.add_argument("--expected-engine", default="ortools-cp-sat")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    choices = payload["choices"]
    rows = {int(row["p"]): row for row in choices}
    if not set(ANCHORS) <= rows.keys():
        raise AssertionError("one or more symmetry anchors are missing")
    for prime, modulus in EXPECTED_MODULI.items():
        if int(rows[prime]["h"]) != modulus:
            raise AssertionError(
                f"p={prime} has h={rows[prime]['h']}, expected {modulus}"
            )

    residual_modulus = 1
    coarse_period = 1
    for row in choices:
        h = int(row["h"])
        residual = q_part(h, args.density_prime)
        residual_modulus = max(residual_modulus, residual)
        coarse_period = math.lcm(coarse_period, h // residual)

    normalize_3 = {}
    for c19 in range(3):
        for c37 in range(3):
            solutions = [
                (dk, dl)
                for dk in range(3)
                for dl in range(3)
                if shifted(rows[19], c19, dk, dl) == 0
                and shifted(rows[37], c37, dk, dl) == 0
            ]
            if len(solutions) != 1:
                raise AssertionError(
                    "order-3 anchors do not normalize uniquely"
                )
            normalize_3[c19, c37] = solutions[0]

    normalize_4 = {}
    for c97 in range(4):
        solutions = [
            (dk, dl)
            for dk in range(12)
            for dl in range(12)
            if dk % 3 == 0
            and dl % 3 == 0
            and shifted(rows[97], c97, dk, dl) == 0
        ]
        if not solutions:
            raise AssertionError("p=97 cannot be normalized")
        normalize_4[c97] = solutions[0]

    order9_checks = 0
    for c109 in range(9):
        for c15121 in range(9):
            solutions = [
                (12 * u, 3 * v)
                for u in range(3)
                for v in range(3)
                if shifted(rows[109], c109, 12 * u, 3 * v)
                == c109 % 3
                and shifted(rows[15121], c15121, 12 * u, 3 * v)
                == c15121 % 3
            ]
            if len(solutions) != 1:
                raise AssertionError(
                    "the two order-9 anchors do not normalize independently"
                )
            order9_checks += 1

    canonical = set()
    assignments = 0
    for c19 in range(3):
        for c37 in range(3):
            shift3 = normalize_3[c19, c37]
            for c73 in range(3):
                normalized73 = shifted(rows[73], c73, *shift3)
                for c97 in range(4):
                    shifted97 = shifted(rows[97], c97, *shift3)
                    shift4 = normalize_4[shifted97]
                    for c109 in range(9):
                        normalized109 = shifted(
                            rows[109], c109, *shift3
                        )
                        normalized109 = shifted(
                            rows[109], normalized109, *shift4
                        )
                        for c15121 in range(9):
                            normalized15121 = shifted(
                                rows[15121], c15121, *shift3
                            )
                            normalized15121 = shifted(
                                rows[15121], normalized15121, *shift4
                            )
                            orbit = (
                                "concurrent"
                                if normalized73 == 0
                                else "triangle"
                            )
                            if normalized73 == 2:
                                normalized109 = (-normalized109) % 9
                                normalized15121 = (
                                    -normalized15121
                                ) % 9
                            canonical.add(
                                (
                                    orbit,
                                    normalized109 % 3,
                                    normalized15121 % 3,
                                )
                            )
                            assignments += 1

    expected = {
        (orbit, r109, r15121)
        for orbit in ("concurrent", "triangle")
        for r109 in range(3)
        for r15121 in range(3)
    }
    if canonical != expected:
        raise AssertionError("the canonical case split is not exhaustive")

    files = []
    total_solve_seconds = 0.0
    for orbit, r109, r15121 in sorted(expected):
        path = result_path(
            args.result_pattern,
            orbit,
            r109,
            r15121,
        )
        result = json.loads(path.read_text())
        expected_fixed = {
            "19": 0,
            "37": 0,
            "73": 0 if orbit == "concurrent" else 1,
            "97": 0,
            "109": r109,
            "15121": r15121,
        }
        observed_fixed = {
            str(prime): int(target)
            for prime, target in result["fixed_targets"].items()
        }
        if observed_fixed != expected_fixed:
            raise AssertionError(f"fixed-target mismatch in {path}")
        checks = {
            "pool": str(args.pool),
            "prime": args.density_prime,
            "rows": len(choices),
            "residual_modulus": residual_modulus,
            "coarse_period": coarse_period,
            "cells": coarse_period * coarse_period,
            "engine": args.expected_engine,
            "status": "INFEASIBLE",
            "finite_density_core_unsat": True,
        }
        for key, expected_value in checks.items():
            if result.get(key) != expected_value:
                raise AssertionError(
                    f"{key} mismatch in {path}: "
                    f"{result.get(key)!r} != {expected_value!r}"
                )
        if int(result["constraints"]) <= 0:
            raise AssertionError(f"no density constraints in {path}")
        total_solve_seconds += float(result["solve_seconds"])
        files.append(str(path))

    period = payload.get("period_filter")
    audit = {
        "pool": str(args.pool),
        "period": int(period) if period is not None else None,
        "density_prime": args.density_prime,
        "rows": len(choices),
        "residual_modulus": residual_modulus,
        "coarse_period": coarse_period,
        "anchor_primes": list(ANCHORS),
        "anchor_assignments_checked": assignments,
        "first_anchor_normalizations": len(normalize_3),
        "p97_normalizations": len(normalize_4),
        "order9_normalizations": order9_checks,
        "canonical_cases": [
            {"orbit": orbit, "r109": r109, "r15121": r15121}
            for orbit, r109, r15121 in sorted(expected)
        ],
        "result_files": files,
        "engine": args.expected_engine,
        "total_solve_seconds": total_solve_seconds,
        "all_phase_assignments_fail_necessary_density_condition": True,
        "scope": (
            f"all phase assignments of the {len(choices)} rows in "
            "the supplied period-filtered component core"
        ),
    }
    args.output.write_text(json.dumps(audit, indent=2) + "\n")
    print(
        f"PASS assignments={assignments} cases={len(expected)} "
        f"engine={args.expected_engine}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
