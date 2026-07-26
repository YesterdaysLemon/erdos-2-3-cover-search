#!/usr/bin/env python3
"""Independent verifier for a component-density peeling certificate."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import exact_uncovered


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.pool.read_text())
    certificate = json.loads(args.certificate.read_text())
    threshold = int(certificate["min_parent_common"])
    rows = [
        row
        for row in payload["choices"]
        if int(row["h"]) > 1
        and int(row.get("parent_common", 1)) >= threshold
    ]
    if len(rows) != int(certificate["initial_row_count"]):
        raise AssertionError("initial row count mismatch")

    for index, step in enumerate(certificate["steps"], start=1):
        prime = int(step["prime"])
        incident = [
            row for row in rows if int(row["h"]) % prime == 0
        ]
        if not incident:
            raise AssertionError(f"step {index} has no incident rows")
        density = Fraction(0)
        directions = Counter()
        homogeneous_coefficients = [1]
        for row in incident:
            h = int(row["h"])
            factors = exact_uncovered.factor(h)
            if prime not in factors:
                raise AssertionError("incident factor disappeared")
            if int(row["a"]) % prime == int(row["b"]) % prime == 0:
                raise AssertionError(
                    f"row p={row['p']} is nonprimitive at {prime}"
                )
            density += Fraction(1, prime ** factors[prime])
            a = int(row["a"]) % prime
            b = int(row["b"]) % prime
            scale = pow(a, -1, prime) if a else pow(b, -1, prime)
            directions[
                (a * scale % prime, b * scale % prime)
            ] += 1
            next_coefficients = [0] * (
                len(homogeneous_coefficients) + 1
            )
            for exponent, coefficient in enumerate(
                homogeneous_coefficients
            ):
                next_coefficients[exponent] = (
                    next_coefficients[exponent] + coefficient * b
                ) % prime
                next_coefficients[exponent + 1] = (
                    next_coefficients[exponent + 1] + coefficient * a
                ) % prime
            homogeneous_coefficients = next_coefficients
        recorded = Fraction(
            int(step["density_numerator"]),
            int(step["density_denominator"]),
        )
        if density != recorded:
            raise AssertionError(f"step {index} density mismatch")
        maximum_direction_capacity = max(directions.values())
        reason = step.get("reason", "density")
        if reason == "density":
            if density >= 1:
                raise AssertionError(
                    f"step {index} is not density-peelable"
                )
        elif reason == "prime_projective_blocking_bound":
            if prime <= 2:
                raise AssertionError("blocking bound requires an odd prime")
            threshold_bound = (3 * prime + 1) // 2
            if maximum_direction_capacity >= prime:
                raise AssertionError(
                    f"step {index} permits a full parallel class"
                )
            if len(incident) >= threshold_bound:
                raise AssertionError(
                    f"step {index} does not beat the blocking bound"
                )
        elif reason == "top_degree_polynomial":
            witness = step.get("polynomial_witness")
            if not witness:
                raise AssertionError("missing polynomial witness")
            exponent_x = int(witness["exponent_x"])
            exponent_y = int(witness["exponent_y"])
            coefficient = int(witness["coefficient_mod_prime"])
            degree = len(incident)
            if exponent_x + exponent_y != degree:
                raise AssertionError("polynomial witness degree mismatch")
            if exponent_x >= prime or exponent_y >= prime:
                raise AssertionError(
                    "polynomial witness is not in the unreduced middle range"
                )
            actual = homogeneous_coefficients[exponent_x] % prime
            if actual != coefficient % prime or actual == 0:
                raise AssertionError(
                    f"step {index} polynomial coefficient mismatch"
                )
        else:
            raise AssertionError(f"unknown peel reason {reason}")
        if "maximum_direction_capacity" in step and int(
            step["maximum_direction_capacity"]
        ) != maximum_direction_capacity:
            raise AssertionError(
                f"step {index} direction capacity mismatch"
            )
        expected_threshold = (3 * prime + 1) // 2 if prime > 2 else None
        if (
            "prime_blocking_line_threshold" in step
            and step["prime_blocking_line_threshold"]
            != expected_threshold
        ):
            raise AssertionError(
                f"step {index} blocking threshold mismatch"
            )
        if int(step["row_count_before"]) != len(rows):
            raise AssertionError(f"step {index} row count mismatch")
        if int(step["removed_row_count"]) != len(incident):
            raise AssertionError(f"step {index} removed count mismatch")
        if sorted(int(row["p"]) for row in incident) != [
            int(value) for value in step["removed_fibre_primes"]
        ]:
            raise AssertionError(f"step {index} row identity mismatch")
        rows = [
            row for row in rows if int(row["h"]) % prime != 0
        ]

    recorded_remaining = [
        int(value) for value in certificate["remaining_fibre_primes"]
    ]
    if sorted(int(row["p"]) for row in rows) != recorded_remaining:
        raise AssertionError("remaining row identity mismatch")
    complete = not rows
    if complete != bool(certificate["complete_peeling"]):
        raise AssertionError("completion flag mismatch")
    if len(rows) != int(certificate["remaining_row_count"]):
        raise AssertionError("remaining row count mismatch")
    print(
        f"verified steps={len(certificate['steps'])} "
        f"remaining={len(rows)} complete={complete}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
