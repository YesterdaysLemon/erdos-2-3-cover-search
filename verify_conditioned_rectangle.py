#!/usr/bin/env python3
"""Independently replay a rectangular derived-pool conditioning."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def merge_congruences(
    left_residue: int,
    left_modulus: int,
    right_residue: int,
    right_modulus: int,
) -> tuple[int, int] | None:
    common = math.gcd(left_modulus, right_modulus)
    if (right_residue - left_residue) % common:
        return None
    left_reduced = left_modulus // common
    right_reduced = right_modulus // common
    multiplier = (
        0
        if right_reduced == 1
        else (
            (right_residue - left_residue)
            // common
            * pow(left_reduced, -1, right_reduced)
        )
        % right_reduced
    )
    modulus = math.lcm(left_modulus, right_modulus)
    return (
        left_residue + left_modulus * multiplier
    ) % modulus, modulus


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("result", type=Path)
    args = parser.parse_args()

    source = json.loads(args.source.read_text())
    result = json.loads(args.result.read_text())
    condition = result["rectangle_condition"]
    k_period = int(condition["k_period"])
    l_period = int(condition["l_period"])
    k0 = int(condition["k_residue"])
    l0 = int(condition["l_residue"])
    shift_x = int(condition["coordinate_shift"]["x"])
    shift_y = int(condition["coordinate_shift"]["y"])
    excluded = {int(prime) for prime in result["excluded_primes"]}
    children = {int(row["p"]): row for row in result["choices"]}
    if len(children) != len(result["choices"]):
        raise AssertionError("result contains duplicate primes")

    for prime in result.get("algebraic_primes", ()):
        prime = int(prime)
        if (k0 + k_period * shift_x) % prime:
            raise AssertionError(
                f"x translation does not normalize algebraic prime {prime}"
            )
        if (l0 + l_period * shift_y) % prime:
            raise AssertionError(
                f"y translation does not normalize algebraic prime {prime}"
            )

    incompatible = 0
    full_cell = []
    expected_primes = set()
    for parent in source["choices"]:
        prime = int(parent["p"])
        if prime in excluded:
            if prime in children:
                raise AssertionError("excluded prime appears in result")
            continue
        h = int(parent["h"])
        a = int(parent["a"]) % h
        b = int(parent["b"]) % h
        residue = int(parent["target_residue"]) % h
        modulus = int(parent["target_modulus"])
        base = (a * k0 + b * l0) % h
        common = math.gcd(k_period * a, l_period * b, h)
        merged = merge_congruences(
            residue,
            modulus,
            base,
            common,
        )
        if merged is None:
            incompatible += 1
            if prime in children:
                raise AssertionError(
                    f"incompatible prime {prime} appears in result"
                )
            continue
        combined, combined_modulus = merged
        new_h = h // common
        if new_h == 1:
            full_cell.append(prime)
            if prime in children:
                raise AssertionError(
                    f"full-cell prime {prime} appears as an ordinary row"
                )
            continue

        expected_primes.add(prime)
        if prime not in children:
            raise AssertionError(f"compatible prime {prime} is missing")
        child = children[prime]
        new_a = (k_period * a // common) % new_h
        new_b = (l_period * b // common) % new_h
        target_shift = (new_a * shift_x + new_b * shift_y) % new_h
        new_modulus = combined_modulus // common
        new_residue = (
            (combined - base) // common - target_shift
        ) % new_modulus
        expected = {
            "h": new_h,
            "p": prime,
            "a": new_a,
            "b": new_b,
            "ord2": new_h // math.gcd(new_a, new_h),
            "ord3": new_h // math.gcd(new_b, new_h),
            "c": new_residue,
            "target_residue": new_residue,
            "target_modulus": new_modulus,
            "parent_h": h,
            "parent_a": a,
            "parent_b": b,
            "parent_base": base,
            "parent_common": common,
            "coordinate_shift": target_shift,
        }
        if child != expected:
            raise AssertionError(f"row mismatch for prime {prime}")

        # Check the divided affine equation at independent deterministic
        # coordinates and at several legal child targets.
        target_count = new_h // new_modulus
        target_indices = {
            0,
            target_count // 2,
            max(0, target_count - 1),
        }
        for target_index in target_indices:
            child_target = new_residue + new_modulus * target_index
            parent_target = (
                base
                + common * (child_target + target_shift)
            ) % h
            if parent_target % modulus != residue % modulus:
                raise AssertionError(
                    f"lifted target violates parent restriction for {prime}"
                )
            for x, y in ((0, 0), (1, 2), (7, 11), (new_h - 1, 3)):
                parent_k = k0 + k_period * (x + shift_x)
                parent_l = l0 + l_period * (y + shift_y)
                parent_hit = (
                    a * parent_k + b * parent_l - parent_target
                ) % h == 0
                child_hit = (
                    new_a * x + new_b * y - child_target
                ) % new_h == 0
                if parent_hit != child_hit:
                    raise AssertionError(
                        f"line equivalence fails for prime {prime}"
                    )

    if set(children) != expected_primes:
        raise AssertionError("result has unexpected child primes")
    if incompatible != int(result["incompatible_rows"]):
        raise AssertionError("incompatible-row count mismatch")
    if sorted(full_cell) != sorted(
        int(prime) for prime in result["full_cell_primes"]
    ):
        raise AssertionError("full-cell prime list mismatch")

    print(
        f"PASS source_rows={len(source['choices'])} "
        f"child_rows={len(children)} incompatible={incompatible} "
        f"full_cell={len(full_cell)} excluded={len(excluded)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
