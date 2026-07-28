#!/usr/bin/env python3
"""Build a rigorously screened one-axis affine-cover pool.

If ``ord_p(2)`` divides ``B``, then the fibre supplied by ``p`` is periodic
in the k-coordinate modulo B.  After fixing one admissible k residue class,
the fibre restricts to a one-dimensional residue class in l.  This script
uses a small MILP only to discover an allocation of those active k classes
whose one-dimensional reciprocal capacity is at least one in every column.
It then replays the capacities with exact rational arithmetic and writes the
corresponding target congruences for an exact whole-lattice CEGIS search.

The l-layered construction is the transposed analogue.

An optional unimodular basis change allows the same construction along any
primitive lattice direction.  If

``(k, l) = x * direction + y * transverse``,

the script replaces the stored coefficients by the exact coefficients of
``x`` and ``y`` before applying the k-layered construction.  Determinant
``+1`` or ``-1`` makes this a bijection of the complete integer lattice.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def transform_basis_row(
    raw: dict,
    direction: tuple[int, int],
    transverse: tuple[int, int],
) -> dict:
    """Express one affine row in a unimodular coordinate basis."""
    du, dv = direction
    tu, tv = transverse
    determinant = du * tv - dv * tu
    if abs(determinant) != 1:
        raise ValueError("coordinate basis must have determinant +1 or -1")
    h = int(raw["h"])
    source_a = int(raw["a"]) % h
    source_b = int(raw["b"]) % h
    if math.gcd(source_a, source_b, h) != 1:
        raise ValueError(f"source row is not surjective for p={raw['p']}")
    a = (source_a * du + source_b * dv) % h
    b = (source_a * tu + source_b * tv) % h
    if math.gcd(a, b, h) != 1:
        raise AssertionError("unimodular basis change lost surjectivity")
    row = dict(raw)
    row.update(
        {
            "a": a,
            "b": b,
            "ord2": h // math.gcd(a, h),
            "ord3": h // math.gcd(b, h),
            "source_a": source_a,
            "source_b": source_b,
            "source_ord2": int(raw["ord2"]),
            "source_ord3": int(raw["ord3"]),
        }
    )
    return row


def transform_basis_rows(
    rows: list[dict],
    direction: tuple[int, int],
    transverse: tuple[int, int],
) -> list[dict]:
    """Transform a complete row pool without changing its primes or phases."""
    return [
        transform_basis_row(row, direction, transverse) for row in rows
    ]


def layer_data(row: dict, axis: str) -> tuple[int, int, int, int]:
    """Return period, active modulus, active coefficient, and line weight."""
    h = int(row["h"])
    a = int(row["a"]) % h
    b = int(row["b"]) % h
    if axis == "k":
        coordinate_period = int(row["ord2"])
        active_modulus = math.gcd(b, h)
        active_coefficient = a
        residual_period = int(row["ord3"])
    else:
        coordinate_period = int(row["ord3"])
        active_modulus = math.gcd(a, h)
        active_coefficient = b
        residual_period = int(row["ord2"])
    if h != active_modulus * residual_period:
        raise ValueError(f"stored order mismatch for p={row['p']}")
    if math.gcd(active_coefficient, active_modulus) != 1:
        raise ValueError(f"active class is not surjective for p={row['p']}")
    return (
        coordinate_period,
        active_modulus,
        active_coefficient,
        residual_period,
    )


def select_rows(payload: dict, axis: str, period: int) -> list[dict]:
    selected = []
    seen_primes = set()
    for raw in payload["choices"]:
        prime = int(raw["p"])
        if prime in seen_primes:
            raise ValueError(f"duplicate source prime {prime}")
        seen_primes.add(prime)
        coordinate_period, _modulus, _coefficient, _residual = layer_data(
            raw, axis
        )
        if period % coordinate_period == 0:
            selected.append(dict(raw))
    if not selected:
        raise ValueError("axis-layered pool is empty")
    return selected


def prime_factors(value: int) -> tuple[int, ...]:
    """Return the distinct prime factors of a positive integer."""
    if value < 1:
        raise ValueError("modulus must be positive")
    factors = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1 if divisor == 2 else 2
    if value > 1:
        factors.append(value)
    return tuple(factors)


def prune_prime_deficit_rows(
    rows: list[dict],
    axis: str,
) -> tuple[list[dict], list[dict]]:
    """Remove rows redundant by the prime-deficit covering lemma."""
    current = list(rows)
    rounds = []
    while True:
        counts = Counter(
            prime
            for row in current
            for prime in prime_factors(layer_data(row, axis)[3])
        )
        deficient = sorted(
            prime for prime, count in counts.items() if count < prime
        )
        if not deficient:
            return current, rounds
        prime = deficient[0]
        removed = [
            row
            for row in current
            if layer_data(row, axis)[3] % prime == 0
        ]
        rounds.append(
            {
                "prime": prime,
                "row_count_before": len(current),
                "divisible_row_count": len(removed),
                "divisible_count_below_prime": len(removed) < prime,
                "removed_primes": [int(row["p"]) for row in removed],
                "row_count_after": len(current) - len(removed),
            }
        )
        removed_primes = {int(row["p"]) for row in removed}
        current = [
            row
            for row in current
            if int(row["p"]) not in removed_primes
        ]


def discover_placement(
    rows: list[dict],
    axis: str,
    period: int,
    time_limit: float,
    pair_cut_rounds: int,
    maximize_minimum: bool,
) -> tuple[list[int], dict]:
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import lil_matrix, vstack

    data = [layer_data(row, axis) for row in rows]
    offsets = []
    choice_variable_count = 0
    for _coordinate_period, modulus, _coefficient, _residual in data:
        offsets.append(choice_variable_count)
        choice_variable_count += modulus
    minimum_variable = (
        choice_variable_count if maximize_minimum else None
    )
    variable_count = choice_variable_count + int(maximize_minimum)

    constraint_count = len(rows) + period
    matrix = lil_matrix(
        (constraint_count, variable_count),
        dtype=np.float64,
    )
    lower = np.empty(constraint_count, dtype=np.float64)
    upper = np.empty(constraint_count, dtype=np.float64)
    for index, (
        _coordinate_period,
        modulus,
        _coefficient,
        _residual,
    ) in enumerate(data):
        matrix[index, offsets[index] : offsets[index] + modulus] = 1.0
        lower[index] = 1.0
        upper[index] = 1.0
    for coordinate in range(period):
        constraint = len(rows) + coordinate
        for index, (
            _coordinate_period,
            modulus,
            _coefficient,
            residual_period,
        ) in enumerate(data):
            matrix[constraint, offsets[index] + coordinate % modulus] = (
                1.0 / residual_period
            )
        if minimum_variable is not None:
            matrix[constraint, minimum_variable] = -1.0
            lower[constraint] = 0.0
        else:
            lower[constraint] = 1.0
        upper[constraint] = np.inf

    constraint_matrix = matrix.tocsr()
    objective = np.zeros(variable_count, dtype=np.float64)
    integrality = np.ones(variable_count, dtype=np.int8)
    bound_lower = np.zeros(variable_count, dtype=np.float64)
    bound_upper = np.ones(variable_count, dtype=np.float64)
    if minimum_variable is not None:
        objective[minimum_variable] = -1.0
        integrality[minimum_variable] = 0
        bound_lower[minimum_variable] = 1.0
        bound_upper[minimum_variable] = sum(
            1.0 / residual_period
            for (
                _coordinate_period,
                _modulus,
                _coefficient,
                residual_period,
            ) in data
        )
    started = time.monotonic()
    cut_count = 0
    added_cut_keys = set()
    solve_records = []
    placement = None
    remaining_violations = []
    for cut_round in range(pair_cut_rounds + 1):
        remaining_time = max(
            0.001,
            time_limit - (time.monotonic() - started),
        )
        result = milp(
            objective,
            integrality=integrality,
            bounds=Bounds(bound_lower, bound_upper),
            constraints=LinearConstraint(
                constraint_matrix,
                lower,
                upper,
            ),
            options={"presolve": True, "time_limit": remaining_time},
        )
        solve_records.append(
            {
                "cut_round": cut_round,
                "status": int(result.status),
                "success": bool(result.success),
                "message": str(result.message),
                "optimized_minimum": (
                    float(result.x[minimum_variable])
                    if result.x is not None
                    and minimum_variable is not None
                    else None
                ),
                "seconds": time.monotonic() - started,
            }
        )
        if result.x is None:
            raise RuntimeError(
                "capacity placement was not found: "
                f"status={result.status} message={result.message}"
            )
        placement = [
            int(
                np.argmax(
                    result.x[offsets[index] : offsets[index] + modulus]
                )
            )
            for index, (
                _coordinate_period,
                modulus,
                _coefficient,
                _residual,
            ) in enumerate(data)
        ]
        remaining_violations = unavoidable_pair_violations(
            rows,
            axis,
            period,
            placement,
        )
        print(
            f"capacity_round={cut_round} "
            f"constraints={constraint_matrix.shape[0]} "
            f"optimized_minimum="
            f"{solve_records[-1]['optimized_minimum']} "
            f"pair_violations={len(remaining_violations)} "
            f"elapsed_s={time.monotonic() - started:.3f}",
            flush=True,
        )
        if not remaining_violations or cut_round == pair_cut_rounds:
            break

        new_violations = []
        for violation in remaining_violations:
            key = (
                violation["coordinate"],
                *violation["row_indices"],
            )
            if key in added_cut_keys:
                continue
            added_cut_keys.add(key)
            new_violations.append(violation)
        if not new_violations:
            raise RuntimeError(
                "the MILP candidate violates only previously added exact "
                "pair cuts; solver tolerances are too loose"
            )
        cut_matrix = lil_matrix(
            (len(new_violations), variable_count),
            dtype=np.float64,
        )
        cut_lower = np.empty(
            len(new_violations),
            dtype=np.float64,
        )
        cut_upper = np.full(
            len(new_violations),
            np.inf,
            dtype=np.float64,
        )
        for cut_index, violation in enumerate(new_violations):
            coordinate = violation["coordinate"]
            overlap = Fraction(
                violation["overlap_numerator"],
                violation["overlap_denominator"],
            )
            # Two residue classes with coprime moduli intersect on exactly
            # 1/(n*m) of Z, independent of their targets.  Since
            # sum(densities) - union_density is the total overcoverage,
            # a covering column must retain capacity at least one after
            # subtracting that unavoidable pair intersection.  The big-M
            # form below becomes the ordinary capacity cut unless both
            # active-class variables are selected.
            for row_index, (
                _coordinate_period,
                modulus,
                _coefficient,
                residual_period,
            ) in enumerate(data):
                variable = offsets[row_index] + coordinate % modulus
                cut_matrix[cut_index, variable] += 1.0 / residual_period
            for row_index in violation["row_indices"]:
                modulus = data[row_index][1]
                variable = offsets[row_index] + coordinate % modulus
                cut_matrix[cut_index, variable] -= float(overlap)
            cut_lower[cut_index] = 1.0 - float(overlap)
        constraint_matrix = vstack(
            (constraint_matrix, cut_matrix.tocsr()),
            format="csr",
        )
        lower = np.concatenate((lower, cut_lower))
        upper = np.concatenate((upper, cut_upper))
        cut_count += len(new_violations)

    assert placement is not None
    metadata = {
        "engine": "scipy-highs-binary-capacity-with-lazy-pair-cuts",
        "variables": variable_count,
        "choice_variables": choice_variable_count,
        "maximize_minimum": maximize_minimum,
        "base_constraints": constraint_count,
        "pair_cut_round_limit": pair_cut_rounds,
        "pair_cuts_added": cut_count,
        "remaining_pair_violations": len(remaining_violations),
        "pair_safe": not remaining_violations,
        "solves": solve_records,
        "seconds": time.monotonic() - started,
    }
    return placement, metadata


def unavoidable_pair_violations(
    rows: list[dict],
    axis: str,
    period: int,
    placement: list[int],
) -> list[dict]:
    """Find exact pair overlaps that exceed a column's capacity slack."""
    data = [layer_data(row, axis) for row in rows]
    violations = []
    for coordinate in range(period):
        active = [
            index
            for index, (
                _coordinate_period,
                modulus,
                _coefficient,
                _residual_period,
            ) in enumerate(data)
            if coordinate % modulus == placement[index]
        ]
        capacity = sum(
            (Fraction(1, data[index][3]) for index in active),
            Fraction(),
        )
        for left_position, left in enumerate(active):
            left_period = data[left][3]
            for right in active[left_position + 1 :]:
                right_period = data[right][3]
                if math.gcd(left_period, right_period) != 1:
                    continue
                overlap = Fraction(1, left_period * right_period)
                if capacity - overlap >= 1:
                    continue
                violations.append(
                    {
                        "coordinate": coordinate,
                        "row_indices": [left, right],
                        "primes": [
                            int(rows[left]["p"]),
                            int(rows[right]["p"]),
                        ],
                        "moduli": [left_period, right_period],
                        "capacity_numerator": capacity.numerator,
                        "capacity_denominator": capacity.denominator,
                        "overlap_numerator": overlap.numerator,
                        "overlap_denominator": overlap.denominator,
                    }
                )
    return violations


def exact_capacity_replay(
    rows: list[dict],
    axis: str,
    period: int,
    placement: list[int],
) -> dict:
    if len(placement) != len(rows):
        raise ValueError("placement length mismatch")
    capacities = []
    for coordinate in range(period):
        capacity = Fraction()
        for row, active_class in zip(rows, placement):
            (
                _coordinate_period,
                modulus,
                _coefficient,
                residual_period,
            ) = layer_data(row, axis)
            if not 0 <= active_class < modulus:
                raise ValueError("active class is outside its modulus")
            if coordinate % modulus == active_class:
                capacity += Fraction(1, residual_period)
        capacities.append(capacity)
    minimum = min(capacities)
    if minimum < 1:
        raise RuntimeError(
            "floating-point discovery failed exact capacity replay: "
            f"minimum={minimum}"
        )
    return {
        "engine": "python-fraction-complete-column-replay",
        "column_count": period,
        "minimum_numerator": minimum.numerator,
        "minimum_denominator": minimum.denominator,
        "minimum_decimal": float(minimum),
        "minimum_multiplicity": sum(value == minimum for value in capacities),
    }


def materialize_rows(
    rows: list[dict],
    axis: str,
    placement: list[int],
) -> list[dict]:
    materialized = []
    for raw, active_class in zip(rows, placement):
        (
            _coordinate_period,
            modulus,
            coefficient,
            _residual_period,
        ) = layer_data(raw, axis)
        row = dict(raw)
        row["target_modulus"] = modulus
        row["target_residue"] = coefficient * active_class % modulus
        row["layer_active_class"] = active_class
        materialized.append(row)
    return materialized


def materialize_column_relaxation(
    rows: list[dict],
    axis: str,
    coordinate: int,
) -> list[dict]:
    """Relax one layer column to independent 1D residue-class phases."""
    relaxed = []
    for raw in rows:
        (
            _coordinate_period,
            modulus,
            _coefficient,
            residual_period,
        ) = layer_data(raw, axis)
        active_class = int(raw["layer_active_class"])
        if coordinate % modulus != active_class:
            continue
        row = {
            "p": int(raw["p"]),
            "h": residual_period,
            "a": 0 if axis == "k" else 1,
            "b": 1 if axis == "k" else 0,
            "ord2": 1 if axis == "k" else residual_period,
            "ord3": residual_period if axis == "k" else 1,
            "c": 0,
            "target_residue": 0,
            "target_modulus": 1,
            "source_h": int(raw["h"]),
        }
        relaxed.append(row)
    if not relaxed:
        raise ValueError("selected column has no active rows")
    return relaxed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--axis", choices=("k", "l"))
    parser.add_argument(
        "--direction",
        nargs=2,
        type=int,
        metavar=("DK", "DL"),
        help=(
            "first vector of a unimodular basis; enables a sheared "
            "coordinate system"
        ),
    )
    parser.add_argument(
        "--transverse",
        nargs=2,
        type=int,
        metavar=("TK", "TL"),
        help="second vector of the same unimodular basis",
    )
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--pair-cut-rounds", type=int, default=0)
    parser.add_argument("--maximize-minimum", action="store_true")
    parser.add_argument(
        "--prime-deficit-prune",
        action="store_true",
        help=(
            "iteratively remove residual moduli that contain a prime q "
            "appearing in fewer than q current rows"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--column",
        type=int,
        help="also emit the independent one-dimensional relaxation here",
    )
    parser.add_argument("--column-output", type=Path)
    args = parser.parse_args()
    if (
        args.period < 1
        or args.time_limit <= 0
        or args.pair_cut_rounds < 0
    ):
        raise SystemExit(
            "period and time limit must be positive; cut rounds nonnegative"
        )
    if (args.column is None) != (args.column_output is None):
        raise SystemExit("--column and --column-output must be used together")
    if args.direction is None:
        if args.transverse is not None:
            raise SystemExit("--transverse requires --direction")
        if args.axis is None:
            raise SystemExit("--axis is required without --direction")
        axis = args.axis
        coordinate_basis = None
    else:
        if args.axis is not None:
            raise SystemExit("--axis and --direction are mutually exclusive")
        if args.transverse is None:
            raise SystemExit("--direction requires --transverse")
        direction = tuple(args.direction)
        transverse = tuple(args.transverse)
        determinant = (
            direction[0] * transverse[1]
            - direction[1] * transverse[0]
        )
        if abs(determinant) != 1:
            raise SystemExit(
                "direction and transverse must have determinant +1 or -1"
            )
        axis = "k"
        coordinate_basis = {
            "direction": list(direction),
            "transverse": list(transverse),
            "determinant": determinant,
            "mapping": (
                "(original_k, original_l) = "
                "x * direction + y * transverse"
            ),
        }

    payload = json.loads(args.pool.read_text())
    source_rows = payload["choices"]
    if coordinate_basis is not None:
        source_rows = transform_basis_rows(
            source_rows,
            tuple(coordinate_basis["direction"]),
            tuple(coordinate_basis["transverse"]),
        )
    selected_rows = select_rows(
        {"choices": source_rows},
        axis,
        args.period,
    )
    pruning_rounds = []
    rows = selected_rows
    if args.prime_deficit_prune:
        rows, pruning_rounds = prune_prime_deficit_rows(rows, axis)
        if not rows:
            raise RuntimeError("prime-deficit pruning removed every row")
    capacity_pattern_period = math.lcm(
        *(layer_data(row, axis)[1] for row in rows)
    )
    placement, discovery = discover_placement(
        rows,
        axis,
        capacity_pattern_period,
        args.time_limit,
        args.pair_cut_rounds,
        args.maximize_minimum,
    )
    replay = exact_capacity_replay(
        rows,
        axis,
        args.period,
        placement,
    )
    materialized = materialize_rows(rows, axis, placement)
    pattern_repeat_count = args.period // capacity_pattern_period
    full_pair_violation_count = (
        discovery["remaining_pair_violations"] * pattern_repeat_count
    )
    result = {
        "schema": "axis_layered_pool_v1",
        "source": str(args.pool),
        "source_sha256": sha256(args.pool),
        "layer_axis": axis,
        "coordinate_basis": coordinate_basis,
        "layer_period": args.period,
        "capacity_pattern_period": capacity_pattern_period,
        "scope": (
            "rows whose selected-coordinate multiplicative order divides "
            "the declared layer period"
        ),
        "source_selection_row_count": len(selected_rows),
        "prime_deficit_pruning": {
            "enabled": args.prime_deficit_prune,
            "theorem": (
                "if fewer than q residue-class moduli are divisible by "
                "prime q, those classes are redundant in any cover of "
                "the integers"
            ),
            "rounds": pruning_rounds,
            "removed_row_count": len(selected_rows) - len(rows),
            "surviving_row_count": len(rows),
        },
        "row_count": len(materialized),
        "raw_reciprocal_density": sum(
            1.0 / int(row["h"]) for row in materialized
        ),
        "capacity_discovery": discovery,
        "pair_overlap_screen": {
            "complete": True,
            "pair_safe": full_pair_violation_count == 0,
            "remaining_violations": full_pair_violation_count,
            "pattern_violation_count": discovery[
                "remaining_pair_violations"
            ],
            "pattern_repeat_count": pattern_repeat_count,
            "claim": (
                "for each selected column, no two active coprime residual "
                "classes have forced intersection larger than the exact "
                "capacity slack"
            ),
        },
        "exact_capacity_replay": replay,
        "choices": materialized,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    if args.column is not None:
        coordinate = args.column % args.period
        column_rows = materialize_column_relaxation(
            materialized,
            axis,
            coordinate,
        )
        column_capacity = sum(
            (Fraction(1, int(row["h"])) for row in column_rows),
            Fraction(),
        )
        column_payload = {
            "schema": "axis_layered_column_relaxation_v1",
            "source": str(args.output),
            "layer_axis": axis,
            "coordinate_basis": coordinate_basis,
            "layer_period": args.period,
            "coordinate": coordinate,
            "scope": (
                "necessary one-dimensional relaxation with independent "
                "residual phases; UNSAT rules out this layered placement, "
                "while SAT does not produce a whole-lattice cover"
            ),
            "row_count": len(column_rows),
            "capacity_numerator": column_capacity.numerator,
            "capacity_denominator": column_capacity.denominator,
            "choices": column_rows,
        }
        args.column_output.write_text(
            json.dumps(column_payload, indent=2) + "\n"
        )
    print(
        f"axis={axis} period={args.period} rows={len(rows)} "
        f"variables={discovery['variables']} "
        f"minimum={replay['minimum_numerator']}/"
        f"{replay['minimum_denominator']} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
