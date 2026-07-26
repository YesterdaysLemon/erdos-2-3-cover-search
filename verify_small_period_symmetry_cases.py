#!/usr/bin/env python3
"""Verify the translation/negation reduction for small 2,3-period covers."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument("--concurrent-result", type=Path, required=True)
    parser.add_argument("--triangle-result", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.period % 12:
        raise SystemExit("period must be divisible by 12")

    rows = {
        int(row["p"]): row
        for row in json.loads(args.pool.read_text())["choices"]
        if args.period % int(row["h"]) == 0
    }
    required_primes = {19, 37, 73, 97}
    if not required_primes <= set(rows):
        raise AssertionError("normalization rows are missing")
    for prime in (19, 37, 73):
        if int(rows[prime]["h"]) != 3:
            raise AssertionError(f"prime {prime} is not an order-3 row")
    if int(rows[97]["h"]) != 4:
        raise AssertionError("prime 97 is not an order-4 row")

    # Exhaustively verify that every pair of phases for the independent
    # p=19 and p=37 rows can be translated to zero modulo 3.
    normalization_count = 0
    third_targets = set()
    for target_19 in range(3):
        for target_37 in range(3):
            solutions = []
            for shift_k in range(3):
                for shift_l in range(3):
                    transformed_19 = (
                        target_19
                        + int(rows[19]["a"]) * shift_k
                        + int(rows[19]["b"]) * shift_l
                    ) % 3
                    transformed_37 = (
                        target_37
                        + int(rows[37]["a"]) * shift_k
                        + int(rows[37]["b"]) * shift_l
                    ) % 3
                    if transformed_19 == transformed_37 == 0:
                        solutions.append((shift_k, shift_l))
            if len(solutions) != 1:
                raise AssertionError(
                    "first two order-3 anchors do not normalize uniquely"
                )
            normalization_count += 1
            shift_k, shift_l = solutions[0]
            for target_73 in range(3):
                third_targets.add(
                    (
                        target_73
                        + int(rows[73]["a"]) * shift_k
                        + int(rows[73]["b"]) * shift_l
                    )
                    % 3
                )
    if third_targets != {0, 1, 2}:
        raise AssertionError("third-anchor targets are not exhaustive")

    # With all mod-3 anchors fixed, a further translation that is zero mod 3
    # can set the p=97 phase to zero modulo 4.  Check every old phase.
    p97_normalizations = {}
    for target_97 in range(4):
        solutions = []
        for shift_k in range(12):
            for shift_l in range(12):
                if shift_k % 3 or shift_l % 3:
                    continue
                transformed = (
                    target_97
                    + int(rows[97]["a"]) * shift_k
                    + int(rows[97]["b"]) * shift_l
                ) % 4
                if transformed == 0:
                    solutions.append((shift_k, shift_l))
        if not solutions:
            raise AssertionError("p=97 cannot be normalized")
        p97_normalizations[str(target_97)] = solutions[0]

    # Negation preserves zero targets and interchanges third-anchor targets
    # 1 and 2.  Therefore target 0 and target 1 are complete orbit reps.
    if (-0) % 3 != 0 or (-1) % 3 != 2 or (-2) % 3 != 1:
        raise AssertionError("negation orbit computation failed")

    concurrent = json.loads(args.concurrent_result.read_text())
    triangle = json.loads(args.triangle_result.read_text())
    expected_cases = (
        (concurrent, {"19": 0, "37": 0, "73": 0, "97": 0}),
        (triangle, {"19": 0, "37": 0, "73": 1, "97": 0}),
    )
    for result, fixed in expected_cases:
        if bool(result["sat"]):
            raise AssertionError("a declared UNSAT orbit is SAT")
        if int(result["grid_period"]) != args.period:
            raise AssertionError("grid period mismatch")
        if int(result["point_count"]) != args.period * args.period:
            raise AssertionError("grid point count mismatch")
        if {
            str(prime): int(target)
            for prime, target in result["fixed_targets"].items()
        } != fixed:
            raise AssertionError("fixed-target orbit mismatch")

    output = {
        "pool": str(args.pool),
        "period": args.period,
        "selected_rows": len(rows),
        "first_anchor_phase_pairs_checked": normalization_count,
        "third_anchor_orbits": [[0], [1, 2]],
        "p97_phase_normalizations": p97_normalizations,
        "representative_fixed_targets": [
            expected_cases[0][1],
            expected_cases[1][1],
        ],
        "concurrent_unsat": True,
        "triangle_unsat": True,
        "all_phase_assignments_excluded": True,
    }
    if args.output:
        args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"PASS period={args.period} rows={len(rows)} "
        "orbits=2 both_unsat=True"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
