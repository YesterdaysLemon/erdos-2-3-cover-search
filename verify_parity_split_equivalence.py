#!/usr/bin/env python3
"""Certify the p=41 parity split and canonical branch equivalence."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


CANONICAL_FIELDS = (
    "h",
    "p",
    "a",
    "b",
    "ord2",
    "ord3",
    "target_residue",
    "target_modulus",
)


def verify_branch(parent, branch, parity: int, split_phase: int) -> None:
    condition = branch["rectangle_condition"]
    if (
        int(condition["k_period"]) != 1
        or int(condition["l_period"]) != 2
        or int(condition["k_residue"]) != 0
        or int(condition["l_residue"]) != parity
    ):
        raise AssertionError("unexpected rectangle condition")
    if [int(value) for value in branch["excluded_primes"]] != [41]:
        raise AssertionError("branch does not exclude exactly p=41")
    parent_by_prime = {
        int(row["p"]): row for row in parent["choices"]
    }
    split = parent_by_prime[41]
    if (
        int(split["h"]) != 2
        or int(split["a"]) % 2 != 0
        or int(split["b"]) % 2 != 1
        or int(split["target_modulus"]) != 1
    ):
        raise AssertionError("p=41 is not the unique parity splitter")
    # The split phase must cover the parity complementary to the residual
    # branch, and must miss every point of the residual branch.
    if split_phase % 2 != 1 - parity:
        raise AssertionError("split phase does not cover complementary parity")
    for l in range(4):
        hit = (l - split_phase) % 2 == 0
        if hit != (l % 2 == 1 - parity):
            raise AssertionError("p=41 parity coverage is incorrect")

    shift_x = int(condition["coordinate_shift"]["x"])
    shift_y = int(condition["coordinate_shift"]["y"])
    children = branch["choices"]
    if set(parent_by_prime) - {41} != {
        int(row["p"]) for row in children
    }:
        raise AssertionError("branch prime set mismatch")
    for child in children:
        prime = int(child["p"])
        row = parent_by_prime[prime]
        h = int(row["h"])
        a = int(row["a"]) % h
        b = int(row["b"]) % h
        base = b * parity % h
        common = math.gcd(a, 2 * b, h)
        new_h = h // common
        new_a = (a // common) % new_h
        new_b = (2 * b // common) % new_h
        shift = (new_a * shift_x + new_b * shift_y) % new_h
        if int(child["target_modulus"]) != 1:
            raise AssertionError("branch target is not unrestricted")
        for child_phase in {0, new_h // 2, new_h - 1}:
            parent_phase = (
                base + common * (child_phase + shift)
            ) % h
            if (
                parent_phase % int(row["target_modulus"])
                != int(row["target_residue"])
            ):
                raise AssertionError(
                    f"arbitrary child phase does not lift for prime {prime}"
                )
            for x, y in ((0, 0), (1, 2), (7, 11), (new_h - 1, 3)):
                parent_k = x + shift_x
                parent_l = parity + 2 * (y + shift_y)
                parent_hit = (
                    a * parent_k + b * parent_l - parent_phase
                ) % h == 0
                child_hit = (
                    new_a * x + new_b * y - child_phase
                ) % new_h == 0
                if parent_hit != child_hit:
                    raise AssertionError(
                        f"symbolic lift sample fails for prime {prime}"
                    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("parent", type=Path)
    parser.add_argument("even_branch", type=Path)
    parser.add_argument("odd_branch", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    parent = json.loads(args.parent.read_text())
    even = json.loads(args.even_branch.read_text())
    odd = json.loads(args.odd_branch.read_text())
    verify_branch(parent, even, parity=0, split_phase=1)
    verify_branch(parent, odd, parity=1, split_phase=0)
    even_rows = [
        tuple(int(row[field]) for field in CANONICAL_FIELDS)
        for row in even["choices"]
    ]
    odd_rows = [
        tuple(int(row[field]) for field in CANONICAL_FIELDS)
        for row in odd["choices"]
    ]
    if even_rows != odd_rows:
        raise AssertionError("canonical parity branches are not identical")
    result = {
        "parent": str(args.parent),
        "even_branch": str(args.even_branch),
        "odd_branch": str(args.odd_branch),
        "parent_rows": len(parent["choices"]),
        "canonical_branch_rows": len(even_rows),
        "split_prime": 41,
        "even_residual_split_phase": 1,
        "odd_residual_split_phase": 0,
        "canonical_rows_identical": True,
        "all_child_targets_unrestricted": True,
        "arbitrary_child_phases_lift": True,
        "cover_existence_equivalent": True,
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"PASS parent={len(parent['choices'])} branch={len(even_rows)} "
        "canonical_identical=True arbitrary_lift=True"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
