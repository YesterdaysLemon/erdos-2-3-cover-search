#!/usr/bin/env python3
"""Independently replay an axis-layered affine-pool artifact.

This verifier does not rerun the discovery MILP.  It authenticates the source
pool, reconstructs an optional unimodular basis change, checks that the
artifact contains exactly every source row whose coordinate order divides the
declared period, and replays all column capacities and unavoidable coprime
pair overlaps with exact rational arithmetic.

Passing this verifier proves only those preprocessing statements.  It does
not choose the remaining row phases or prove a cover of the exponent lattice.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def transformed_signature(
    raw: dict,
    basis: dict | None,
) -> dict:
    h = int(raw["h"])
    source_a = int(raw["a"]) % h
    source_b = int(raw["b"]) % h
    if math.gcd(source_a, source_b, h) != 1:
        raise ValueError(f"source row is not surjective for p={raw['p']}")
    source_ord2 = int(raw["ord2"])
    source_ord3 = int(raw["ord3"])
    if source_ord2 != h // math.gcd(source_a, h):
        raise ValueError(f"source ord2 is inconsistent for p={raw['p']}")
    if source_ord3 != h // math.gcd(source_b, h):
        raise ValueError(f"source ord3 is inconsistent for p={raw['p']}")
    if basis is None:
        return {
            "h": h,
            "p": int(raw["p"]),
            "a": source_a,
            "b": source_b,
            "ord2": source_ord2,
            "ord3": source_ord3,
            "c": int(raw["c"]),
        }
    direction = tuple(map(int, basis["direction"]))
    transverse = tuple(map(int, basis["transverse"]))
    determinant = (
        direction[0] * transverse[1]
        - direction[1] * transverse[0]
    )
    if abs(determinant) != 1:
        raise ValueError("coordinate basis is not unimodular")
    if determinant != int(basis["determinant"]):
        raise ValueError("stored coordinate determinant is incorrect")
    a = (source_a * direction[0] + source_b * direction[1]) % h
    b = (
        source_a * transverse[0] + source_b * transverse[1]
    ) % h
    if math.gcd(a, b, h) != 1:
        raise AssertionError("basis change lost row surjectivity")
    return {
        "h": h,
        "p": int(raw["p"]),
        "a": a,
        "b": b,
        "ord2": h // math.gcd(a, h),
        "ord3": h // math.gcd(b, h),
        "c": int(raw["c"]),
        "source_a": source_a,
        "source_b": source_b,
        "source_ord2": source_ord2,
        "source_ord3": source_ord3,
    }


def layer_data(row: dict, axis: str) -> tuple[int, int, int, int]:
    h = int(row["h"])
    a = int(row["a"]) % h
    b = int(row["b"]) % h
    if axis == "k":
        coordinate_order = h // math.gcd(a, h)
        active_modulus = math.gcd(b, h)
        active_coefficient = a
    elif axis == "l":
        coordinate_order = h // math.gcd(b, h)
        active_modulus = math.gcd(a, h)
        active_coefficient = b
    else:
        raise ValueError("layer axis must be k or l")
    residual_period = h // active_modulus
    if math.gcd(active_coefficient, active_modulus) != 1:
        raise ValueError(f"active coefficient is not invertible for p={row['p']}")
    return (
        coordinate_order,
        active_modulus,
        active_coefficient,
        residual_period,
    )


def prime_factors(value: int) -> tuple[int, ...]:
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


def replay_prime_deficit_pruning(
    rows: list[dict],
    axis: str,
) -> tuple[list[dict], list[dict]]:
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


def verify_payload(source: dict, artifact: dict) -> dict:
    if artifact.get("schema") != "axis_layered_pool_v1":
        raise ValueError("unsupported layered-pool schema")
    axis = str(artifact["layer_axis"])
    basis = artifact.get("coordinate_basis")
    if basis is not None and axis != "k":
        raise ValueError("a transformed artifact must use its first coordinate")
    period = int(artifact["layer_period"])
    if period < 1:
        raise ValueError("layer period must be positive")

    expected = []
    seen_source_primes = set()
    for raw in source["choices"]:
        prime = int(raw["p"])
        if prime in seen_source_primes:
            raise ValueError(f"duplicate source prime {prime}")
        seen_source_primes.add(prime)
        row = transformed_signature(raw, basis)
        coordinate_order, _modulus, _coefficient, _residual = layer_data(
            row, axis
        )
        if period % coordinate_order == 0:
            expected.append(row)

    pruning = artifact.get("prime_deficit_pruning")
    pruning_metadata_matches = True
    if pruning is not None and bool(pruning.get("enabled")):
        source_selection_count = int(
            artifact.get("source_selection_row_count", -1)
        )
        selected_count_matches = source_selection_count == len(expected)
        expected, replayed_rounds = replay_prime_deficit_pruning(
            expected,
            axis,
        )
        pruning_metadata_matches = (
            selected_count_matches
            and pruning.get("rounds") == replayed_rounds
            and int(pruning.get("removed_row_count", -1))
            == source_selection_count - len(expected)
            and int(pruning.get("surviving_row_count", -1))
            == len(expected)
        )
    elif pruning is not None:
        pruning_metadata_matches = (
            pruning.get("rounds") == []
            and int(pruning.get("removed_row_count", -1)) == 0
            and int(pruning.get("surviving_row_count", -1))
            == len(expected)
        )

    rows = artifact["choices"]
    if len(rows) != len(expected):
        raise ValueError("artifact row count does not match source selection")
    if len(rows) != int(artifact["row_count"]):
        raise ValueError("stored artifact row count is incorrect")
    seen_artifact_primes = set()
    placements = []
    for actual, source_row in zip(rows, expected, strict=True):
        prime = int(actual["p"])
        if prime in seen_artifact_primes:
            raise ValueError(f"duplicate artifact prime {prime}")
        seen_artifact_primes.add(prime)
        for key, value in source_row.items():
            if int(actual[key]) != int(value):
                raise ValueError(f"transformed field {key} differs for p={prime}")
        (
            _coordinate_order,
            active_modulus,
            active_coefficient,
            _residual_period,
        ) = layer_data(actual, axis)
        active_class = int(actual["layer_active_class"])
        if not 0 <= active_class < active_modulus:
            raise ValueError(f"active class is invalid for p={prime}")
        if int(actual["target_modulus"]) != active_modulus:
            raise ValueError(f"target modulus is invalid for p={prime}")
        target_residue = active_coefficient * active_class % active_modulus
        if int(actual["target_residue"]) != target_residue:
            raise ValueError(f"target residue is invalid for p={prime}")
        placements.append(active_class)

    capacity_pattern_period = math.lcm(
        *(layer_data(row, axis)[1] for row in rows)
    )
    if period % capacity_pattern_period:
        raise AssertionError("capacity pattern does not divide layer period")
    repeat_count = period // capacity_pattern_period
    capacities = []
    pattern_pair_violations = 0
    for coordinate in range(capacity_pattern_period):
        active = []
        capacity = Fraction()
        for index, (row, active_class) in enumerate(
            zip(rows, placements, strict=True)
        ):
            (
                _coordinate_order,
                active_modulus,
                _active_coefficient,
                residual_period,
            ) = layer_data(row, axis)
            if coordinate % active_modulus != active_class:
                continue
            active.append((index, residual_period))
            capacity += Fraction(1, residual_period)
        capacities.append(capacity)
        for left_position, (_left_index, left_period) in enumerate(active):
            for _right_index, right_period in active[left_position + 1 :]:
                if math.gcd(left_period, right_period) != 1:
                    continue
                if capacity - Fraction(1, left_period * right_period) < 1:
                    pattern_pair_violations += 1

    minimum = min(capacities)
    pair_violations = pattern_pair_violations * repeat_count
    minimum_multiplicity = (
        sum(value == minimum for value in capacities) * repeat_count
    )
    pattern_metadata_matches = int(
        artifact.get(
            "capacity_pattern_period",
            capacity_pattern_period,
        )
    ) == capacity_pattern_period
    replay = artifact["exact_capacity_replay"]
    capacity_metadata_matches = (
        int(replay["column_count"]) == period
        and int(replay["minimum_numerator"]) == minimum.numerator
        and int(replay["minimum_denominator"]) == minimum.denominator
        and int(replay["minimum_multiplicity"])
        == minimum_multiplicity
    )
    pair_screen = artifact["pair_overlap_screen"]
    pair_metadata_matches = (
        bool(pair_screen["complete"])
        and bool(pair_screen["pair_safe"]) == (pair_violations == 0)
        and int(pair_screen["remaining_violations"]) == pair_violations
    )
    checks = {
        "exact_source_row_selection": True,
        "prime_deficit_pruning": pruning_metadata_matches,
        "target_restrictions": True,
        "capacity_pattern_period": pattern_metadata_matches,
        "all_columns_have_capacity_at_least_one": minimum >= 1,
        "capacity_metadata": capacity_metadata_matches,
        "complete_coprime_pair_screen": pair_metadata_matches,
    }
    return {
        "verified": all(checks.values()),
        "checks": checks,
        "row_count": len(rows),
        "layer_period": period,
        "minimum_capacity_numerator": minimum.numerator,
        "minimum_capacity_denominator": minimum.denominator,
        "minimum_capacity_decimal": float(minimum),
        "minimum_capacity_multiplicity": minimum_multiplicity,
        "unavoidable_pair_violations": pair_violations,
        "scope": (
            "independent source-selection, unimodular-transform, exact "
            "column-capacity, and unavoidable-coprime-pair replay; this "
            "does not choose residual phases or prove a lattice cover"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    artifact = json.loads(args.artifact.read_text())
    source_path = args.artifact.parent / str(artifact["source"])
    source_hash_matches = sha256(source_path) == artifact["source_sha256"]
    result = verify_payload(
        json.loads(source_path.read_text()),
        artifact,
    )
    result["artifact"] = str(args.artifact)
    result["artifact_sha256"] = sha256(args.artifact)
    result["source"] = str(source_path)
    result["source_sha256_matches"] = source_hash_matches
    result["checks"]["source_sha256"] = source_hash_matches
    result["verified"] = all(result["checks"].values())
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"verified={result['verified']} rows={result['row_count']} "
        f"period={result['layer_period']} "
        f"minimum={result['minimum_capacity_numerator']}/"
        f"{result['minimum_capacity_denominator']} "
        f"pair_violations={result['unavoidable_pair_violations']}",
        flush=True,
    )
    return 0 if result["verified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
