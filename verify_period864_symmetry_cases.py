#!/usr/bin/env python3
"""Verify the exhaustive symmetry split used for period-864 covers."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ANCHOR_PRIMES = (19, 37, 73, 97, 109, 15121)


def transformed_target(
    row: dict[str, object],
    target: int,
    shift_k: int,
    shift_l: int,
) -> int:
    modulus = int(row["h"])
    return (
        target
        + int(row["a"]) * shift_k
        + int(row["b"]) * shift_l
    ) % modulus


def result_path(pattern: str, orbit: str, r109: int, r15121: int) -> Path:
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
    parser.add_argument("--period", type=int, default=864)
    parser.add_argument(
        "--result-pattern",
        required=True,
        help=(
            "path pattern with {orbit}, {r109}, and {r15121} placeholders"
        ),
    )
    parser.add_argument("--expected-solver")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.period != 864:
        raise SystemExit("this audit is specific to period 864")

    payload = json.loads(args.pool.read_text())
    rows = {
        int(row["p"]): row
        for row in payload["choices"]
        if args.period % int(row["h"]) == 0
    }
    if not set(ANCHOR_PRIMES) <= set(rows):
        raise AssertionError("one or more symmetry anchors are missing")
    expected_moduli = {
        19: 3,
        37: 3,
        73: 3,
        97: 4,
        109: 9,
        15121: 9,
    }
    for prime, modulus in expected_moduli.items():
        if int(rows[prime]["h"]) != modulus:
            raise AssertionError(
                f"prime {prime} has modulus {rows[prime]['h']}, "
                f"expected {modulus}"
            )

    # First normalize the independent p=19 and p=37 order-3 anchors.
    first_anchor_normalizations: dict[tuple[int, int], tuple[int, int]] = {}
    for target_19 in range(3):
        for target_37 in range(3):
            solutions = []
            for shift_k in range(3):
                for shift_l in range(3):
                    if (
                        transformed_target(
                            rows[19], target_19, shift_k, shift_l
                        )
                        == 0
                        and transformed_target(
                            rows[37], target_37, shift_k, shift_l
                        )
                        == 0
                    ):
                        solutions.append((shift_k, shift_l))
            if len(solutions) != 1:
                raise AssertionError(
                    "p=19 and p=37 do not give a unique mod-3 translation"
                )
            first_anchor_normalizations[(target_19, target_37)] = solutions[0]

    # A translation divisible by 3 preserves the three order-3 anchors.
    # Verify that p=97 can then be set to zero.
    p97_normalizations: dict[int, tuple[int, int]] = {}
    for target_97 in range(4):
        solutions = []
        for shift_k in range(12):
            for shift_l in range(12):
                if shift_k % 3 or shift_l % 3:
                    continue
                if (
                    transformed_target(
                        rows[97], target_97, shift_k, shift_l
                    )
                    == 0
                ):
                    solutions.append((shift_k, shift_l))
            if solutions:
                break
        if not solutions:
            raise AssertionError("p=97 cannot be normalized after mod-3 anchors")
        p97_normalizations[target_97] = solutions[0]

    # A translation k=12u, l=3v preserves the order-3 anchors and p=97.
    # Exhaustively check that its action is transitive on the three lifts of
    # each fixed residue modulo 3 for both order-9 anchors simultaneously.
    order9_normalizations = 0
    for target_109 in range(9):
        for target_15121 in range(9):
            residue_109 = target_109 % 3
            residue_15121 = target_15121 % 3
            solutions = []
            for u in range(3):
                for v in range(3):
                    shift_k = 12 * u
                    shift_l = 3 * v
                    if (
                        transformed_target(
                            rows[109], target_109, shift_k, shift_l
                        )
                        == residue_109
                        and transformed_target(
                            rows[15121],
                            target_15121,
                            shift_k,
                            shift_l,
                        )
                        == residue_15121
                    ):
                        solutions.append((shift_k, shift_l))
            if len(solutions) != 1:
                raise AssertionError(
                    "order-9 anchors are not independently normalized"
                )
            order9_normalizations += 1

    # Exhaust every possible phase assignment on all six anchors and replay
    # the normalization sequence.  Negation maps p=73 target 2 to target 1;
    # target 0 stays in the concurrent orbit.
    canonical_cases: set[tuple[str, int, int]] = set()
    assignments_checked = 0
    for target_19 in range(3):
        for target_37 in range(3):
            shift_3 = first_anchor_normalizations[(target_19, target_37)]
            for target_73 in range(3):
                normalized_73 = transformed_target(
                    rows[73], target_73, *shift_3
                )
                for target_97 in range(4):
                    shifted_97 = transformed_target(
                        rows[97], target_97, *shift_3
                    )
                    shift_4 = p97_normalizations[shifted_97]
                    for target_109 in range(9):
                        shifted_109 = transformed_target(
                            rows[109], target_109, *shift_3
                        )
                        shifted_109 = transformed_target(
                            rows[109], shifted_109, *shift_4
                        )
                        for target_15121 in range(9):
                            shifted_15121 = transformed_target(
                                rows[15121], target_15121, *shift_3
                            )
                            shifted_15121 = transformed_target(
                                rows[15121], shifted_15121, *shift_4
                            )
                            orbit = (
                                "concurrent"
                                if normalized_73 == 0
                                else "triangle"
                            )
                            if normalized_73 == 2:
                                shifted_109 = (-shifted_109) % 9
                                shifted_15121 = (-shifted_15121) % 9
                            canonical_cases.add(
                                (
                                    orbit,
                                    shifted_109 % 3,
                                    shifted_15121 % 3,
                                )
                            )
                            assignments_checked += 1

    expected_cases = {
        (orbit, r109, r15121)
        for orbit in ("concurrent", "triangle")
        for r109 in range(3)
        for r15121 in range(3)
    }
    if canonical_cases != expected_cases:
        raise AssertionError("canonical case split is not exhaustive")

    result_files: list[str] = []
    solvers: set[str] = set()
    for orbit, r109, r15121 in sorted(expected_cases):
        path = result_path(args.result_pattern, orbit, r109, r15121)
        result = json.loads(path.read_text())
        fixed = {
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
        if observed_fixed != fixed:
            raise AssertionError(f"fixed-target mismatch in {path}")
        if bool(result["sat"]):
            raise AssertionError(f"canonical case is SAT: {path}")
        if int(result["grid_period"]) != args.period:
            raise AssertionError(f"grid-period mismatch in {path}")
        if int(result["point_count"]) != args.period * args.period:
            raise AssertionError(f"point-count mismatch in {path}")
        solver = str(result["solver"])
        if args.expected_solver and solver != args.expected_solver:
            raise AssertionError(
                f"solver mismatch in {path}: {solver} != "
                f"{args.expected_solver}"
            )
        solvers.add(solver)
        result_files.append(str(path))

    output = {
        "pool": str(args.pool),
        "period": args.period,
        "selected_rows": len(rows),
        "anchor_primes": list(ANCHOR_PRIMES),
        "anchor_assignments_checked": assignments_checked,
        "first_anchor_normalizations": len(first_anchor_normalizations),
        "p97_normalizations": len(p97_normalizations),
        "order9_normalizations": order9_normalizations,
        "canonical_cases": [
            {"orbit": orbit, "r109": r109, "r15121": r15121}
            for orbit, r109, r15121 in sorted(expected_cases)
        ],
        "result_files": result_files,
        "solvers": sorted(solvers),
        "all_phase_assignments_excluded": True,
    }
    if args.output:
        args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"PASS period={args.period} assignments={assignments_checked} "
        f"canonical_cases={len(expected_cases)} solvers={sorted(solvers)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
