#!/usr/bin/env python3
"""Independently verify phases and points projected to a lattice child."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("lattice", type=Path)
    parser.add_argument("--parent-phase-file", type=Path, required=True)
    parser.add_argument("--child-phase-file", type=Path, required=True)
    parser.add_argument(
        "--parent-points-file", type=Path, action="append", default=[]
    )
    parser.add_argument("--child-points-file", type=Path)
    args = parser.parse_args()
    if bool(args.parent_points_file) != bool(args.child_points_file):
        raise SystemExit("parent and child point arguments must be paired")

    source = json.loads(args.source.read_text())
    lattice = json.loads(args.lattice.read_text())
    parents = {int(row["p"]): row for row in source["choices"]}
    children = {int(row["p"]): row for row in lattice["choices"]}
    parent_phases = {
        int(prime): int(value)
        for prime, value in json.loads(
            args.parent_phase_file.read_text()
        ).items()
    }
    child_phases = {
        int(prime): int(value)
        for prime, value in json.loads(
            args.child_phase_file.read_text()
        ).items()
    }
    if set(child_phases) != set(children):
        raise AssertionError("child phase prime set mismatch")
    condition = lattice["lattice_condition"]
    matrix = condition["matrix"]
    matrix_a = int(matrix["a"])
    matrix_b = int(matrix["b"])
    matrix_c = int(matrix["c"])
    matrix_d = int(matrix["d"])
    determinant = int(condition["determinant"])
    k0 = int(condition["offset"]["k"])
    l0 = int(condition["offset"]["l"])
    shift_x = int(condition["coordinate_shift"]["x"])
    shift_y = int(condition["coordinate_shift"]["y"])

    active = 0
    inactive = 0
    for prime, child in children.items():
        parent = parents[prime]
        h = int(parent["h"])
        a = int(parent["a"]) % h
        b = int(parent["b"]) % h
        base = (a * k0 + b * l0) % h
        coeff_x = a * matrix_a + b * matrix_c
        coeff_y = a * matrix_b + b * matrix_d
        common = math.gcd(coeff_x, coeff_y, h)
        new_h = h // common
        new_a = (coeff_x // common) % new_h
        new_b = (coeff_y // common) % new_h
        shift = (new_a * shift_x + new_b * shift_y) % new_h
        phase = parent_phases.get(prime)
        if phase is not None:
            phase %= h
            if (
                phase % int(parent["target_modulus"])
                != int(parent["target_residue"])
            ):
                raise AssertionError(f"illegal parent phase for {prime}")
        if phase is None or (phase - base) % common:
            expected = int(child["target_residue"])
            inactive += 1
        else:
            expected = ((phase - base) // common - shift) % new_h
            active += 1
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
                    a * parent_k + b * parent_l - phase
                ) % h == 0
                child_hit = (
                    new_a * x + new_b * y - expected
                ) % new_h == 0
                if parent_hit != child_hit:
                    raise AssertionError(
                        f"phase projection fails for {prime}"
                    )
        if child_phases[prime] % new_h != expected:
            raise AssertionError(f"child phase mismatch for {prime}")

    point_count = 0
    if args.parent_points_file:
        expected_points = []
        seen = set()
        for path in args.parent_points_file:
            for raw_k, raw_l in json.loads(path.read_text()):
                k = int(raw_k)
                l = int(raw_l)
                relative_k = k - k0
                relative_l = l - l0
                numerator_x = (
                    matrix_d * relative_k
                    - matrix_b * relative_l
                )
                numerator_y = (
                    -matrix_c * relative_k
                    + matrix_a * relative_l
                )
                if (
                    numerator_x % determinant
                    or numerator_y % determinant
                ):
                    continue
                point = (
                    numerator_x // determinant - shift_x,
                    numerator_y // determinant - shift_y,
                )
                if point not in seen:
                    seen.add(point)
                    expected_points.append([point[0], point[1]])
        actual_points = json.loads(args.child_points_file.read_text())
        if actual_points != expected_points:
            raise AssertionError("projected point list mismatch")
        point_count = len(actual_points)
    print(
        f"PASS phases={len(children)} active={active} "
        f"inactive_defaulted={inactive} points={point_count}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
