#!/usr/bin/env python3
"""Independently replay an affine-lattice phase lift."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("parent", type=Path)
    parser.add_argument("lattice", type=Path)
    parser.add_argument("child_phases", type=Path)
    parser.add_argument("parent_phases", type=Path)
    args = parser.parse_args()

    parent = json.loads(args.parent.read_text())
    lattice = json.loads(args.lattice.read_text())
    children = {int(row["p"]): row for row in lattice["choices"]}
    parents = {int(row["p"]): row for row in parent["choices"]}
    child_phases = {
        int(prime): int(value)
        for prime, value in json.loads(args.child_phases.read_text()).items()
    }
    parent_phases = {
        int(prime): int(value)
        for prime, value in json.loads(args.parent_phases.read_text()).items()
    }
    if set(parent_phases) != set(parents):
        raise AssertionError("parent phase prime set mismatch")
    if set(child_phases) != set(children):
        raise AssertionError("child phase prime set mismatch")

    condition = lattice["lattice_condition"]
    matrix = condition["matrix"]
    matrix_a = int(matrix["a"])
    matrix_b = int(matrix["b"])
    matrix_c = int(matrix["c"])
    matrix_d = int(matrix["d"])
    k0 = int(condition["offset"]["k"])
    l0 = int(condition["offset"]["l"])
    shift_x = int(condition["coordinate_shift"]["x"])
    shift_y = int(condition["coordinate_shift"]["y"])

    for prime, child in children.items():
        row = parents[prime]
        h = int(row["h"])
        a = int(row["a"]) % h
        b = int(row["b"]) % h
        base = (a * k0 + b * l0) % h
        coeff_x = a * matrix_a + b * matrix_c
        coeff_y = a * matrix_b + b * matrix_d
        common = math.gcd(coeff_x, coeff_y, h)
        new_h = h // common
        new_a = (coeff_x // common) % new_h
        new_b = (coeff_y // common) % new_h
        shift = (new_a * shift_x + new_b * shift_y) % new_h
        child_phase = child_phases[prime] % new_h
        expected_parent = (
            base + common * (child_phase + shift)
        ) % h
        actual_parent = parent_phases[prime] % h
        if actual_parent != expected_parent:
            raise AssertionError(f"lifted phase mismatch for {prime}")
        if (
            actual_parent % int(row["target_modulus"])
            != int(row["target_residue"])
        ):
            raise AssertionError(f"illegal parent phase for {prime}")
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
                a * parent_k + b * parent_l - actual_parent
            ) % h == 0
            child_hit = (
                new_a * x + new_b * y - child_phase
            ) % new_h == 0
            if parent_hit != child_hit:
                raise AssertionError(
                    f"line lift equivalence fails for {prime}"
                )

    extra = set(parents) - set(children)
    for prime in extra:
        row = parents[prime]
        phase = parent_phases[prime] % int(row["h"])
        if (
            phase % int(row["target_modulus"])
            != int(row["target_residue"])
        ):
            raise AssertionError(f"illegal extra phase for {prime}")
    print(
        f"PASS lifted={len(children)} extra={len(extra)} "
        f"parent={len(parents)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
