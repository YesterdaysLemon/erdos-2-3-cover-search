#!/usr/bin/env python3
"""Independently verify projected phases and lesson points."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("rectangle", type=Path)
    parser.add_argument("--parent-phase-file", type=Path, required=True)
    parser.add_argument("--child-phase-file", type=Path, required=True)
    parser.add_argument(
        "--parent-points-file",
        type=Path,
        action="append",
        default=[],
    )
    parser.add_argument("--child-points-file", type=Path)
    args = parser.parse_args()
    if bool(args.parent_points_file) != bool(args.child_points_file):
        raise SystemExit(
            "parent point files and --child-points-file must be used together"
        )

    source = json.loads(args.source.read_text())
    rectangle = json.loads(args.rectangle.read_text())
    parents = {int(row["p"]): row for row in source["choices"]}
    children = {int(row["p"]): row for row in rectangle["choices"]}
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
        raise AssertionError("projected phase prime set does not match rectangle")

    condition = rectangle["rectangle_condition"]
    k_period = int(condition["k_period"])
    l_period = int(condition["l_period"])
    k0 = int(condition["k_residue"])
    l0 = int(condition["l_residue"])
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
        common = math.gcd(k_period * a, l_period * b, h)
        new_h = h // common
        new_a = (k_period * a // common) % new_h
        new_b = (l_period * b // common) % new_h
        target_shift = (new_a * shift_x + new_b * shift_y) % new_h
        parent_phase = parent_phases.get(prime)
        if parent_phase is not None:
            parent_phase %= h
            parent_residue = int(parent["target_residue"])
            parent_modulus = int(parent["target_modulus"])
            if parent_phase % parent_modulus != parent_residue:
                raise AssertionError(
                    f"illegal parent phase for prime {prime}"
                )
        if parent_phase is None or (parent_phase - base) % common:
            expected = int(child["target_residue"])
            inactive += 1
        else:
            expected = (
                (parent_phase - base) // common - target_shift
            ) % new_h
            active += 1
            for x, y in ((0, 0), (1, 2), (7, 11), (new_h - 1, 3)):
                parent_k = k0 + k_period * (x + shift_x)
                parent_l = l0 + l_period * (y + shift_y)
                parent_hit = (
                    a * parent_k + b * parent_l - parent_phase
                ) % h == 0
                child_hit = (
                    new_a * x + new_b * y - expected
                ) % new_h == 0
                if parent_hit != child_hit:
                    raise AssertionError(
                        f"phase line projection fails for prime {prime}"
                    )
        actual = child_phases[prime] % new_h
        if actual != expected:
            raise AssertionError(
                f"projected phase mismatch for prime {prime}"
            )
        if (
            actual % int(child["target_modulus"])
            != int(child["target_residue"])
        ):
            raise AssertionError(
                f"illegal child phase for prime {prime}"
            )

    point_count = 0
    if args.parent_points_file:
        expected_points = []
        seen = set()
        for path in args.parent_points_file:
            for raw_k, raw_l in json.loads(path.read_text()):
                k = int(raw_k)
                l = int(raw_l)
                if (k - k0) % k_period or (l - l0) % l_period:
                    continue
                point = (
                    (k - k0) // k_period - shift_x,
                    (l - l0) // l_period - shift_y,
                )
                if point not in seen:
                    seen.add(point)
                    expected_points.append([point[0], point[1]])
        actual_points = json.loads(args.child_points_file.read_text())
        if actual_points != expected_points:
            raise AssertionError("projected lesson list mismatch")
        point_count = len(actual_points)
        for x, y in (
            actual_points[:3]
            + actual_points[len(actual_points) // 2 : len(actual_points) // 2 + 3]
            + actual_points[-3:]
        ):
            parent_k = k0 + k_period * (int(x) + shift_x)
            parent_l = l0 + l_period * (int(y) + shift_y)
            if (
                (parent_k - k0) % k_period
                or (parent_l - l0) % l_period
            ):
                raise AssertionError("lesson does not lift into rectangle")

    print(
        f"PASS phases={len(child_phases)} active={active} "
        f"inactive_defaulted={inactive} points={point_count}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
