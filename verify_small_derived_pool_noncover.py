#!/usr/bin/env python3
"""Independently enumerate a small derived-pool noncover certificate."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import time
from pathlib import Path


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for prime in small:
        if value % prime == 0:
            return value == prime
    exponent = value - 1
    power = 0
    while exponent % 2 == 0:
        power += 1
        exponent //= 2
    for base in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if base % value == 0:
            continue
        witness = pow(base, exponent, value)
        if witness in (1, value - 1):
            continue
        for _ in range(power - 1):
            witness = witness * witness % value
            if witness == value - 1:
                break
        else:
            return False
    return True


def factor(value: int) -> dict[int, int]:
    result: dict[int, int] = {}
    divisor = 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            result[divisor] = result.get(divisor, 0) + 1
            value //= divisor
        divisor += 1
    if value > 1:
        result[value] = result.get(value, 0) + 1
    return result


def multiplicative_order(base: int, prime: int) -> int:
    if not is_prime(prime) or prime <= base or base % prime == 0:
        raise ValueError("multiplicative order requires a valid prime modulus")
    order = prime - 1
    for divisor in factor(order):
        while order % divisor == 0 and pow(base, order // divisor, prime) == 1:
            order //= divisor
    return order


def crt(congruences: list[tuple[int, int]]) -> int:
    value = 0
    modulus = 1
    for residue, next_modulus in congruences:
        common = math.gcd(modulus, next_modulus)
        if (residue - value) % common:
            raise ValueError("incompatible CRT congruences")
        left = modulus // common
        right = next_modulus // common
        step = 0 if right == 1 else (
            (residue - value) // common * pow(left, -1, right)
        ) % right
        value = (value + modulus * step) % (modulus * right)
        modulus *= right
    return value


def merge_congruences(
    first_residue: int,
    first_modulus: int,
    second_residue: int,
    second_modulus: int,
) -> tuple[int, int]:
    common = math.gcd(first_modulus, second_modulus)
    if (second_residue - first_residue) % common:
        raise ValueError("conditioned row has incompatible congruences")
    left = first_modulus // common
    right = second_modulus // common
    step = 0 if right == 1 else (
        (second_residue - first_residue)
        // common
        * pow(left, -1, right)
    ) % right
    modulus = first_modulus * right
    return (first_residue + first_modulus * step) % modulus, modulus


def power_target_congruence(
    h: int,
    prime: int,
    power: int,
) -> tuple[int, int]:
    """Independently derive targets making -g^c an n-th power mod p."""
    divisor = math.gcd(power, prime - 1)
    subgroup_step = (prime - 1) // h
    common = math.gcd(subgroup_step, divisor)
    half = (prime - 1) // 2
    if half % common:
        raise ValueError("row has no power-compatible phase")
    modulus = divisor // common
    if modulus == 1:
        return 0, 1
    residue = (
        half // common
        * pow(subgroup_step // common, -1, modulus)
    ) % modulus
    return residue, modulus


def validate_conditioned_rows(scope: dict, rows: list[dict]) -> None:
    """Reconstruct every derived row without consulting the source pool."""
    power = int(scope["power"])
    condition = scope["cell_condition"]
    period = int(condition["period"])
    k0 = int(condition["k_residue"]) % period
    l0 = int(condition["l_residue"]) % period
    shift_x = int(condition["coordinate_shift"]["x"])
    shift_y = int(condition["coordinate_shift"]["y"])

    canonical_primes = [
        prime
        for prime in factor(power)
        if prime % 2 and math.gcd(prime, period) == 1
    ]
    for prime in factor(power):
        if (
            prime % 2
            and math.gcd(prime, period) != 1
            and k0 % prime == 0
            and l0 % prime == 0
        ):
            raise ValueError("declared cell is wholly algebraically covered")
    if canonical_primes != [
        int(prime) for prime in scope["algebraic_primes"]
    ]:
        raise ValueError("canonical algebraic-prime list mismatch")
    expected_shift_x = crt(
        [
            ((-k0 * pow(period, -1, prime)) % prime, prime)
            for prime in canonical_primes
        ]
    )
    expected_shift_y = crt(
        [
            ((-l0 * pow(period, -1, prime)) % prime, prime)
            for prime in canonical_primes
        ]
    )
    if (shift_x, shift_y) != (expected_shift_x, expected_shift_y):
        raise ValueError("canonical coordinate shift mismatch")

    for row in rows:
        prime = int(row["p"])
        original_h = int(row["original_h"])
        original_a = int(row["original_a"]) % original_h
        original_b = int(row["original_b"]) % original_h
        order2 = multiplicative_order(2, prime)
        order3 = multiplicative_order(3, prime)
        if math.lcm(order2, order3) != original_h:
            raise ValueError(f"original subgroup order mismatch for p={prime}")
        if original_h // math.gcd(original_a, original_h) != order2:
            raise ValueError(f"original coefficient a mismatch for p={prime}")
        if original_h // math.gcd(original_b, original_h) != order3:
            raise ValueError(f"original coefficient b mismatch for p={prime}")
        if pow(2, original_b, prime) != pow(3, original_a, prime):
            raise ValueError(f"original affine-log relation fails for p={prime}")

        power_residue, power_modulus = power_target_congruence(
            original_h,
            prime,
            power,
        )
        base = (original_a * k0 + original_b * l0) % original_h
        common = math.gcd(
            period * original_a,
            period * original_b,
            original_h,
        )
        combined, combined_modulus = merge_congruences(
            power_residue,
            power_modulus,
            base,
            common,
        )
        h = original_h // common
        a = (period * original_a // common) % h
        b = (period * original_b // common) % h
        target_modulus = combined_modulus // common
        target_shift = (a * shift_x + b * shift_y) % h
        target_residue = (
            (combined - base) // common - target_shift
        ) % target_modulus
        expected = {
            "h": h,
            "a": a,
            "b": b,
            "ord2": h // math.gcd(a, h),
            "ord3": h // math.gcd(b, h),
            "c": target_residue,
            "target_residue": target_residue,
            "target_modulus": target_modulus,
            "coordinate_shift": target_shift,
        }
        for field, value in expected.items():
            if int(row[field]) != value:
                raise ValueError(
                    f"conditioned field {field} mismatch for p={prime}"
                )


def allowed_targets(row: dict) -> tuple[int, ...]:
    h = int(row["h"])
    modulus = int(row.get("target_modulus", 1))
    if modulus < 1 or h % modulus:
        raise ValueError("target modulus must be positive and divide h")
    residue = int(row.get("target_residue", 0)) % modulus
    return tuple(range(residue, h, modulus))


def enumerate_phase_cover(
    rows: list[dict],
    points: list[tuple[int, int]],
) -> tuple[bool, int, int, list[int] | None]:
    full = (1 << len(points)) - 1
    row_options = []
    assignment_count = 1
    for row in rows:
        h = int(row["h"])
        a = int(row["a"])
        b = int(row["b"])
        options = []
        for target in allowed_targets(row):
            mask = 0
            for index, (k, l) in enumerate(points):
                if (a * k + b * l) % h == target:
                    mask |= 1 << index
            options.append((target, mask))
        assignment_count *= len(options)
        row_options.append(options)

    enumerated = 0
    for selection in itertools.product(*row_options):
        enumerated += 1
        covered = 0
        for _target, mask in selection:
            covered |= mask
            if covered == full:
                return (
                    True,
                    enumerated,
                    assignment_count,
                    [target for target, _mask in selection],
                )
    return False, enumerated, assignment_count, None


def verify(payload: dict) -> dict:
    if payload.get("schema") != "small_derived_pool_noncover_v1":
        raise ValueError("unsupported certificate schema")
    scope = payload["scope"]
    rows = payload["rows"]
    points = [(int(k), int(l)) for k, l in payload["points"]]
    if len(rows) != int(scope["row_count"]):
        raise ValueError("row count mismatch")
    if len(points) != int(scope["point_count"]):
        raise ValueError("point count mismatch")
    if len(points) != len(set(points)):
        raise ValueError("duplicate certificate point")
    if len({int(row["p"]) for row in rows}) != len(rows):
        raise ValueError("duplicate source prime")
    if any(int(row["h"]) > int(scope["max_h"]) for row in rows):
        raise ValueError("certificate row exceeds the declared modulus bound")
    excluded_primes = {
        int(prime) for prime in scope["cell_condition"]["excluded_primes"]
    }
    if any(int(row["p"]) in excluded_primes for row in rows):
        raise ValueError("certificate contains an excluded source prime")
    if not all(is_prime(int(row["p"])) for row in rows):
        raise ValueError("certificate contains a composite source prime")
    if any(
        math.gcd(
            int(row["a"]),
            int(row["b"]),
            int(row["h"]),
        )
        != 1
        for row in rows
    ):
        raise ValueError("certificate contains a non-surjective row")
    if math.lcm(*(int(row["h"]) for row in rows)) != int(scope["period"]):
        raise ValueError("derived torus period mismatch")
    if any(
        not 0 <= coordinate < int(scope["period"])
        for point in points
        for coordinate in point
    ):
        raise ValueError("certificate point lies outside the declared torus")
    validate_conditioned_rows(scope, rows)
    algebraic_primes = tuple(
        int(value) for value in scope["algebraic_primes"]
    )
    if any(
        k % prime == 0 and l % prime == 0
        for k, l in points
        for prime in algebraic_primes
    ):
        raise ValueError("certificate contains an algebraically covered point")

    started = time.monotonic()
    cover_exists, enumerated, assignment_count, covering_phases = (
        enumerate_phase_cover(rows, points)
    )
    expected = int(scope["phase_assignment_count"])
    if assignment_count != expected:
        raise ValueError("phase-assignment count mismatch")
    verified = not cover_exists and enumerated == assignment_count
    return {
        "verified": verified,
        "row_count": len(rows),
        "point_count": len(points),
        "algebraic_primes": list(algebraic_primes),
        "phase_assignments_declared": expected,
        "phase_assignments_enumerated": enumerated,
        "cover_exists": cover_exists,
        "covering_phases": covering_phases,
        "verification_seconds": time.monotonic() - started,
        "engine": "independent-python-exhaustive-bitset-enumeration",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    raw = args.certificate.read_bytes()
    payload = json.loads(raw)
    result = verify(payload)
    result["certificate"] = str(args.certificate)
    result["certificate_sha256"] = hashlib.sha256(raw).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)
    return 0 if result["verified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
