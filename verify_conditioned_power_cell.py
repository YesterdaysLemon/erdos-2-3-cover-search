#!/usr/bin/env python3
"""Independently replay power compatibility and rectangular cell conditioning."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


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


def crt(congruences: list[tuple[int, int]]) -> int:
    value = 0
    modulus = 1
    for residue, new_modulus in congruences:
        gcd = math.gcd(modulus, new_modulus)
        if (residue - value) % gcd:
            raise AssertionError("incompatible CRT congruences")
        left = modulus // gcd
        right = new_modulus // gcd
        step = (
            (residue - value) // gcd * pow(left, -1, right)
        ) % right
        value += modulus * step
        modulus *= right
        value %= modulus
    return value


def merge_congruences(
    residue_a: int,
    modulus_a: int,
    residue_b: int,
    modulus_b: int,
) -> tuple[int, int] | None:
    gcd = math.gcd(modulus_a, modulus_b)
    if (residue_b - residue_a) % gcd:
        return None
    left = modulus_a // gcd
    right = modulus_b // gcd
    if right == 1:
        step = 0
    else:
        step = (
            (residue_b - residue_a) // gcd
            * pow(left, -1, right)
        ) % right
    modulus = modulus_a * right
    return (residue_a + modulus_a * step) % modulus, modulus


def power_target_congruence(h: int, p: int, power: int) -> tuple[int, int]:
    divisor = math.gcd(power, p - 1)
    subgroup_step = (p - 1) // h
    common = math.gcd(subgroup_step, divisor)
    half = (p - 1) // 2
    if half % common:
        raise RuntimeError("row has no power-compatible target")
    modulus = divisor // common
    if modulus == 1:
        return 0, 1
    residue = (
        (half // common)
        * pow(subgroup_step // common, -1, modulus)
    ) % modulus
    return residue, modulus


def row_key(row: dict[str, object]) -> tuple[int, ...]:
    return tuple(
        int(row[field])
        for field in (
            "h",
            "p",
            "a",
            "b",
            "ord2",
            "ord3",
            "c",
            "target_residue",
            "target_modulus",
            "original_h",
            "original_a",
            "original_b",
            "coordinate_shift",
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("conditioned", type=Path)
    parser.add_argument("--max-component", type=int, default=0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    source = json.loads(args.source.read_text())
    conditioned = json.loads(args.conditioned.read_text())
    power = int(conditioned["power"])
    condition = conditioned["cell_condition"]
    period = int(condition["period"])
    k0 = int(condition["k_residue"]) % period
    l0 = int(condition["l_residue"]) % period
    excluded = {int(prime) for prime in condition["excluded_primes"]}
    recorded_shift_x = int(condition["coordinate_shift"]["x"])
    recorded_shift_y = int(condition["coordinate_shift"]["y"])

    algebraic_primes = [
        prime for prime in factor(power) if prime % 2
    ]
    canonical_primes = []
    x_congruences = []
    y_congruences = []
    for prime in algebraic_primes:
        common = math.gcd(period, prime)
        if common != 1:
            if k0 % prime == 0 and l0 % prime == 0:
                raise AssertionError(
                    "selected cell is wholly algebraically covered"
                )
            continue
        inverse = pow(period, -1, prime)
        canonical_primes.append(prime)
        x_congruences.append(((-k0 * inverse) % prime, prime))
        y_congruences.append(((-l0 * inverse) % prime, prime))
    shift_x = crt(x_congruences) if x_congruences else 0
    shift_y = crt(y_congruences) if y_congruences else 0
    if (shift_x, shift_y) != (recorded_shift_x, recorded_shift_y):
        raise AssertionError("canonical coordinate shift mismatch")
    if canonical_primes != [
        int(prime) for prime in conditioned["algebraic_primes"]
    ]:
        raise AssertionError("canonical algebraic-prime list mismatch")

    expected = []
    incompatible = 0
    full_cell = 0
    filtered = 0
    for raw in source["choices"]:
        p = int(raw["p"])
        if p in excluded:
            continue
        h = int(raw["h"])
        if args.max_component and max(
            (
                prime**exponent
                for prime, exponent in factor(h).items()
            ),
            default=1,
        ) > args.max_component:
            filtered += 1
            continue
        a = int(raw["a"]) % h
        b = int(raw["b"]) % h
        power_residue, power_modulus = power_target_congruence(
            h, p, power
        )
        base = (a * k0 + b * l0) % h
        common = math.gcd(period * a, period * b, h)
        merged = merge_congruences(
            power_residue,
            power_modulus,
            base,
            common,
        )
        if merged is None:
            incompatible += 1
            continue
        combined, combined_modulus = merged
        new_h = h // common
        new_a = (period * a // common) % new_h
        new_b = (period * b // common) % new_h
        target_modulus = combined_modulus // common
        target_shift = (new_a * shift_x + new_b * shift_y) % new_h
        target_residue = (
            (combined - base) // common - target_shift
        ) % target_modulus
        if new_h == 1:
            full_cell += 1
        expected.append(
            {
                "h": new_h,
                "p": p,
                "a": new_a,
                "b": new_b,
                "ord2": new_h // math.gcd(new_a, new_h),
                "ord3": new_h // math.gcd(new_b, new_h),
                "c": target_residue,
                "target_residue": target_residue,
                "target_modulus": target_modulus,
                "original_h": h,
                "original_a": a,
                "original_b": b,
                "coordinate_shift": target_shift,
            }
        )

    expected.sort(key=lambda row: (int(row["h"]), int(row["p"])))
    actual = conditioned["choices"]
    if [row_key(row) for row in expected] != [
        row_key(row) for row in actual
    ]:
        raise AssertionError("conditioned rows differ from independent replay")
    if incompatible != int(conditioned["incompatible_rows"]):
        raise AssertionError("incompatible-row count mismatch")
    if full_cell != int(conditioned["full_cell_rows"]):
        raise AssertionError("full-cell row count mismatch")

    result = {
        "source": str(args.source),
        "conditioned": str(args.conditioned),
        "power": power,
        "period": period,
        "source_rows": len(source["choices"]),
        "conditioned_rows": len(actual),
        "filtered_by_component_bound": filtered,
        "incompatible_rows": incompatible,
        "full_cell_rows": full_cell,
        "canonical_algebraic_primes": canonical_primes,
        "verified": True,
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"PASS source_rows={len(source['choices'])} "
        f"conditioned_rows={len(actual)} filtered={filtered} "
        f"incompatible={incompatible} full_cell={full_cell}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
