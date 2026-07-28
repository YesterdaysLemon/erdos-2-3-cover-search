#!/usr/bin/env python3
"""Deterministically phase an arithmetic finite-group line inventory.

For unrestricted line targets, independently uniform phases leave

    |G| * product_i (1 - 1 / h_i)

cells uncovered in expectation. Conditional expectation makes this proof
constructive: process the actual-prime rows one at a time and choose a target
covering the most currently uncovered cells. If the initial exact expectation
is below one, the final integer number of misses must be zero.

A zero-miss output is still a candidate until a structurally independent
full-domain checker validates it and ``build_crt_m.py`` validates every
arithmetic prime signature.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

from build_axis_layered_pool import sha256
from inventory_reverse_group_lines import (
    canonical_projective_direction,
    transformed_coefficients,
    transformed_target_restriction,
)


def source_target(
    h: int,
    source_a: int,
    source_b: int,
    canonical_a: int,
    canonical_b: int,
    canonical_target: int,
) -> int:
    multiplier = next(
        value
        for value in range(1, h + 1)
        if math.gcd(value, h) == 1
        and (
            value * source_a % h,
            value * source_b % h,
        )
        == (canonical_a, canonical_b)
    )
    if h == 1:
        return 0
    return pow(multiplier, -1, h) * canonical_target % h


def construct(inventory: dict, pool: dict, pool_path: Path) -> dict:
    if inventory.get("schema") != "reverse_group_arithmetic_line_inventory_v3":
        raise ValueError("inventory schema mismatch")
    if sha256(pool_path) != inventory["source_sha256"]:
        raise ValueError("inventory source SHA-256 mismatch")
    if json.loads(pool_path.read_text()) != pool:
        raise ValueError("pool payload differs from authenticated file")
    width = int(inventory["group"]["width"])
    height = int(inventory["group"]["height"])
    direction = tuple(int(value) for value in inventory["basis"]["direction"])
    transverse = tuple(
        int(value) for value in inventory["basis"]["transverse"]
    )
    if (
        width < 1
        or height < 1
        or abs(
            direction[0] * transverse[1]
            - direction[1] * transverse[0]
        )
        != 1
    ):
        raise ValueError("invalid inventory group or basis")
    raw_by_prime = {}
    for row in pool["choices"]:
        prime = int(row["p"])
        if prime in raw_by_prime:
            raise ValueError(f"duplicate source prime {prime}")
        raw_by_prime[prime] = row

    slots = []
    slot_primes = set()
    for line_type in inventory["line_types"]:
        if int(line_type["target_modulus"]) != 1:
            continue
        h = int(line_type["h"])
        canonical_a = int(line_type["a"]) % h
        canonical_b = int(line_type["b"]) % h
        if canonical_projective_direction(
            h,
            canonical_a,
            canonical_b,
        ) != (canonical_a, canonical_b):
            raise ValueError("inventory direction is not canonical")
        signature_primes = sorted(
            int(signature["p"])
            for signature in line_type["source_signatures"]
        )
        if (
            signature_primes
            != sorted(int(prime) for prime in line_type["primes"])
            or len(signature_primes) != int(line_type["prime_capacity"])
        ):
            raise ValueError("inventory line-type prime capacity mismatch")
        for signature in line_type["source_signatures"]:
            prime = int(signature["p"])
            if prime in slot_primes:
                raise ValueError(f"inventory reuses prime {prime}")
            slot_primes.add(prime)
            raw = raw_by_prime.get(prime)
            if raw is None:
                raise ValueError(f"inventory prime is absent: {prime}")
            source_h, source_a, source_b = transformed_coefficients(
                raw,
                direction,
                transverse,
            )
            if (
                source_h != h
                or source_a != int(signature["a"]) % h
                or source_b != int(signature["b"]) % h
                or (source_a * width) % h
                or (source_b * height) % h
            ):
                raise ValueError(f"inventory signature mismatch p={prime}")
            raw_canonical = canonical_projective_direction(
                h,
                source_a,
                source_b,
            )
            if raw_canonical != (canonical_a, canonical_b):
                raise ValueError(
                    f"inventory canonical direction mismatch p={prime}"
                )
            target_modulus, target_residue = (
                transformed_target_restriction(
                    raw,
                    h,
                    source_a,
                    source_b,
                    raw_canonical,
                )
            )
            if (
                target_modulus != int(line_type["target_modulus"])
                or target_residue
                != int(line_type["target_residue"]) % target_modulus
            ):
                raise ValueError(
                    f"inventory target restriction mismatch p={prime}"
                )
            if target_modulus != 1:
                raise ValueError(
                    f"restricted row entered unrestricted slots p={prime}"
                )
            slots.append(
                {
                    "p": prime,
                    "h": h,
                    "canonical_a": canonical_a,
                    "canonical_b": canonical_b,
                    "source_a": source_a,
                    "source_b": source_b,
                    "source_target_modulus": int(
                        raw.get("target_modulus", 1)
                    ),
                    "source_target_residue": int(
                        raw.get("target_residue", 0)
                    ),
                    "raw": raw,
                }
            )
    slots.sort(key=lambda item: (item["h"], item["p"]))
    expected = Fraction(width * height)
    for slot in slots:
        expected *= Fraction(slot["h"] - 1, slot["h"])
    if expected != Fraction(
        int(inventory["expected_uncovered_cell_count_numerator"]),
        int(inventory["expected_uncovered_cell_count_denominator"]),
    ):
        raise ValueError("inventory expected-uncovered value mismatch")

    uncovered = [
        (x, y)
        for x in range(width)
        for y in range(height)
    ]
    choices = []
    for slot in slots:
        h = slot["h"]
        counts = [0] * h
        for x, y in uncovered:
            target = (
                slot["canonical_a"] * x
                + slot["canonical_b"] * y
            ) % h
            counts[target] += 1
        best_count = max(counts)
        canonical_target = counts.index(best_count)
        uncovered = [
            (x, y)
            for x, y in uncovered
            if (
                slot["canonical_a"] * x
                + slot["canonical_b"] * y
            )
            % h
            != canonical_target
        ]
        target = source_target(
            h,
            slot["source_a"],
            slot["source_b"],
            slot["canonical_a"],
            slot["canonical_b"],
            canonical_target,
        )
        if (
            target - slot["source_target_residue"]
        ) % slot["source_target_modulus"]:
            raise AssertionError(
                f"constructed target violates source restriction "
                f"p={slot['p']}"
            )
        materialized = dict(slot["raw"])
        materialized["c"] = target
        choices.append(materialized)
    if expected < 1 and uncovered:
        raise AssertionError(
            "conditional expectation below one left an uncovered cell"
        )

    replay_misses = []
    for x in range(width):
        for y in range(height):
            if not any(
                (
                    transformed_coefficients(
                        row,
                        direction,
                        transverse,
                    )[1]
                    * x
                    + transformed_coefficients(
                        row,
                        direction,
                        transverse,
                    )[2]
                    * y
                    - int(row["c"])
                )
                % int(row["h"])
                == 0
                for row in choices
            ):
                replay_misses.append([x, y])
    if replay_misses != [list(point) for point in uncovered]:
        raise AssertionError("direct finite-group replay differs")
    return {
        "schema": "reverse_group_affine_cover_candidate_v1",
        "inventory_source": inventory["source"],
        "inventory_source_sha256": inventory["source_sha256"],
        "basis": inventory["basis"],
        "group": inventory["group"],
        "row_count": len(choices),
        "initial_expected_uncovered_numerator": expected.numerator,
        "initial_expected_uncovered_denominator": expected.denominator,
        "initial_expected_uncovered_decimal": float(expected),
        "constructor": "deterministic-conditional-expectation-greedy",
        "choices": choices,
        "finite_group_miss_count": len(replay_misses),
        "finite_group_misses": replay_misses,
        "claim": {
            "finite_group_cover_candidate": not replay_misses,
            "independent_full_domain_verification_passed": False,
            "integer_m_found": False,
        },
        "scope": (
            "deterministic phase assignment for unrestricted actual-prime "
            "rows descending to the declared quotient; a zero-miss result "
            "still requires independent arithmetic and full-domain replay"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inventory", type=Path)
    parser.add_argument("pool", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    inventory = json.loads(args.inventory.read_text())
    pool = json.loads(args.pool.read_text())
    result = construct(inventory, pool, args.pool)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"rows={result['row_count']} "
        f"expected={result['initial_expected_uncovered_numerator']}/"
        f"{result['initial_expected_uncovered_denominator']} "
        f"misses={result['finite_group_miss_count']} "
        f"candidate={result['claim']['finite_group_cover_candidate']} "
        f"output={args.output}",
        flush=True,
    )
    return 0 if not result["finite_group_misses"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
