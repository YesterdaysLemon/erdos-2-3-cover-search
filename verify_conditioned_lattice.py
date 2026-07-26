#!/usr/bin/env python3
"""Independently replay an affine-lattice derived-pool conditioning."""

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
    left = left_modulus // common
    right = right_modulus // common
    multiplier = (
        0
        if right == 1
        else (
            (right_residue - left_residue)
            // common
            * pow(left, -1, right)
        )
        % right
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
    condition = result["lattice_condition"]
    matrix = condition["matrix"]
    matrix_a = int(matrix["a"])
    matrix_b = int(matrix["b"])
    matrix_c = int(matrix["c"])
    matrix_d = int(matrix["d"])
    k0 = int(condition["offset"]["k"])
    l0 = int(condition["offset"]["l"])
    shift_x = int(condition["coordinate_shift"]["x"])
    shift_y = int(condition["coordinate_shift"]["y"])
    determinant = matrix_a * matrix_d - matrix_b * matrix_c
    if determinant != int(condition["determinant"]):
        raise AssertionError("determinant mismatch")
    excluded = {int(prime) for prime in result["excluded_primes"]}
    children = {int(row["p"]): row for row in result["choices"]}
    if len(children) != len(result["choices"]):
        raise AssertionError("duplicate child prime")

    for prime in result.get("algebraic_primes", ()):
        prime = int(prime)
        parent_k = (
            k0 + matrix_a * shift_x + matrix_b * shift_y
        ) % prime
        parent_l = (
            l0 + matrix_c * shift_x + matrix_d * shift_y
        ) % prime
        if parent_k or parent_l:
            raise AssertionError(
                f"algebraic translation fails modulo {prime}"
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
        coeff_x = a * matrix_a + b * matrix_c
        coeff_y = a * matrix_b + b * matrix_d
        common = math.gcd(coeff_x, coeff_y, h)
        merged = merge_congruences(
            residue, modulus, base, common
        )
        if merged is None:
            incompatible += 1
            if prime in children:
                raise AssertionError("incompatible prime appears in result")
            continue
        combined, combined_modulus = merged
        new_h = h // common
        if new_h == 1:
            full_cell.append(prime)
            if prime in children:
                raise AssertionError("full-cell prime appears as child row")
            continue
        expected_primes.add(prime)
        if prime not in children:
            raise AssertionError(f"compatible prime {prime} is missing")
        new_a = (coeff_x // common) % new_h
        new_b = (coeff_y // common) % new_h
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
            "parent_coeff_x": coeff_x,
            "parent_coeff_y": coeff_y,
            "parent_common": common,
            "coordinate_shift": target_shift,
        }
        if children[prime] != expected:
            raise AssertionError(f"row mismatch for prime {prime}")
        targets = {
            new_residue,
            (new_residue + new_modulus * (new_h // new_modulus // 2))
            % new_h,
            (new_residue - new_modulus) % new_h,
        }
        for child_target in targets:
            if child_target % new_modulus != new_residue:
                raise AssertionError("invalid deterministic child target")
            parent_target = (
                base + common * (child_target + target_shift)
            ) % h
            if parent_target % modulus != residue % modulus:
                raise AssertionError("lift violates parent target restriction")
            for x, y in ((0, 0), (1, 2), (7, 11), (new_h - 1, 3)):
                parent_k = (
                    k0
                    + matrix_a * (x + shift_x)
                    + matrix_b * (y + shift_y)
                )
                parent_l = (
                    l0
                    + matrix_c * (x + shift_x)
                    + matrix_d * (y + shift_y)
                )
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
        raise AssertionError("unexpected child prime")
    if incompatible != int(result["incompatible_rows"]):
        raise AssertionError("incompatible count mismatch")
    if sorted(full_cell) != sorted(
        int(prime) for prime in result["full_cell_primes"]
    ):
        raise AssertionError("full-cell list mismatch")
    print(
        f"PASS source_rows={len(source['choices'])} "
        f"child_rows={len(children)} incompatible={incompatible} "
        f"full_cell={len(full_cell)} excluded={len(excluded)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
