#!/usr/bin/env python3
"""Project phases and lesson points into an affine-lattice child pool."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("lattice", type=Path)
    parser.add_argument("--phase-file", type=Path)
    parser.add_argument("--phase-output", type=Path)
    parser.add_argument(
        "--points-file", type=Path, action="append", default=[]
    )
    parser.add_argument("--points-output", type=Path)
    args = parser.parse_args()
    if bool(args.phase_file) != bool(args.phase_output):
        raise SystemExit("--phase-file and --phase-output must be paired")
    if bool(args.points_file) != bool(args.points_output):
        raise SystemExit("--points-file and --points-output must be paired")

    source = json.loads(args.source.read_text())
    lattice = json.loads(args.lattice.read_text())
    parents = {int(row["p"]): row for row in source["choices"]}
    children = {int(row["p"]): row for row in lattice["choices"]}
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

    if args.phase_file:
        parent_phases = {
            int(prime): int(value)
            for prime, value in json.loads(
                args.phase_file.read_text()
            ).items()
        }
        phases = {}
        active = 0
        inactive = 0
        for prime, child in children.items():
            parent = parents[prime]
            h = int(parent["h"])
            phase = parent_phases.get(prime)
            if phase is None:
                child_phase = int(child["target_residue"])
                inactive += 1
            else:
                phase %= h
                if (
                    phase % int(parent["target_modulus"])
                    != int(parent["target_residue"])
                ):
                    raise RuntimeError(f"illegal parent phase for {prime}")
                base = int(child["parent_base"])
                common = int(child["parent_common"])
                if (phase - base) % common:
                    child_phase = int(child["target_residue"])
                    inactive += 1
                else:
                    child_phase = (
                        (phase - base) // common
                        - int(child["coordinate_shift"])
                    ) % int(child["h"])
                    active += 1
            if (
                child_phase % int(child["target_modulus"])
                != int(child["target_residue"])
            ):
                raise AssertionError(f"illegal child phase for {prime}")
            phases[str(prime)] = child_phase
        args.phase_output.write_text(json.dumps(phases) + "\n")
        print(
            f"phases={len(phases)} active={active} "
            f"inactive_defaulted={inactive} output={args.phase_output}",
            flush=True,
        )

    if args.points_file:
        points = []
        seen = set()
        input_count = 0
        matching = 0
        for path in args.points_file:
            for raw_k, raw_l in json.loads(path.read_text()):
                input_count += 1
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
                matching += 1
                x = numerator_x // determinant - shift_x
                y = numerator_y // determinant - shift_y
                if (
                    k
                    != k0
                    + matrix_a * (x + shift_x)
                    + matrix_b * (y + shift_y)
                    or l
                    != l0
                    + matrix_c * (x + shift_x)
                    + matrix_d * (y + shift_y)
                ):
                    raise AssertionError("inverse lattice map failed")
                point = (x, y)
                if point not in seen:
                    seen.add(point)
                    points.append([x, y])
        args.points_output.write_text(json.dumps(points) + "\n")
        print(
            f"input_points={input_count} matching={matching} "
            f"unique={len(points)} output={args.points_output}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
