#!/usr/bin/env python3
"""Project phases and lesson points into a conditioned rectangle.

The rectangle constructor writes enough parent arithmetic into each child row
to make this projection exact.  A parent phase that is inactive on the chosen
rectangle is replaced by the child's first legal target; such rows remain
mutable during subsequent synthesis.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("rectangle", type=Path)
    parser.add_argument("--phase-file", type=Path)
    parser.add_argument("--phase-output", type=Path)
    parser.add_argument(
        "--points-file",
        type=Path,
        action="append",
        default=[],
        help="parent lesson file; may be supplied more than once",
    )
    parser.add_argument("--points-output", type=Path)
    args = parser.parse_args()
    if bool(args.phase_file) != bool(args.phase_output):
        raise SystemExit("--phase-file and --phase-output must be used together")
    if bool(args.points_file) != bool(args.points_output):
        raise SystemExit("--points-file and --points-output must be used together")

    source = json.loads(args.source.read_text())
    rectangle = json.loads(args.rectangle.read_text())
    source_by_prime = {
        int(row["p"]): row for row in source["choices"]
    }
    child_by_prime = {
        int(row["p"]): row for row in rectangle["choices"]
    }
    condition = rectangle["rectangle_condition"]
    k_period = int(condition["k_period"])
    l_period = int(condition["l_period"])
    k0 = int(condition["k_residue"])
    l0 = int(condition["l_residue"])
    shift_x = int(condition["coordinate_shift"]["x"])
    shift_y = int(condition["coordinate_shift"]["y"])

    if args.phase_file:
        parent_phases = {
            int(prime): int(value)
            for prime, value in json.loads(
                args.phase_file.read_text()
            ).items()
        }
        child_phases = {}
        active = 0
        inactive = 0
        for prime, child in child_by_prime.items():
            parent = source_by_prime[prime]
            h = int(parent["h"])
            parent_phase = parent_phases.get(prime)
            if parent_phase is None:
                child_phase = int(child["target_residue"])
                inactive += 1
            else:
                parent_phase %= h
                parent_residue = int(parent["target_residue"])
                parent_modulus = int(parent["target_modulus"])
                if parent_phase % parent_modulus != parent_residue:
                    raise RuntimeError(
                        f"parent phase violates target restriction for {prime}"
                    )
                base = int(child["parent_base"])
                common = int(child["parent_common"])
                target_shift = int(child["coordinate_shift"])
                if (parent_phase - base) % common:
                    child_phase = int(child["target_residue"])
                    inactive += 1
                else:
                    child_phase = (
                        (parent_phase - base) // common - target_shift
                    ) % int(child["h"])
                    active += 1
            if (
                child_phase % int(child["target_modulus"])
                != int(child["target_residue"])
            ):
                raise AssertionError(
                    f"projected child phase is illegal for {prime}"
                )
            child_phases[str(prime)] = child_phase
        args.phase_output.write_text(
            json.dumps(child_phases) + "\n"
        )
        print(
            f"phases={len(child_phases)} active={active} "
            f"inactive_defaulted={inactive} output={args.phase_output}",
            flush=True,
        )

    if args.points_file:
        points = []
        seen = set()
        input_count = 0
        matching_count = 0
        for points_file in args.points_file:
            for raw_k, raw_l in json.loads(points_file.read_text()):
                input_count += 1
                k = int(raw_k)
                l = int(raw_l)
                if (k - k0) % k_period or (l - l0) % l_period:
                    continue
                matching_count += 1
                point = (
                    (k - k0) // k_period - shift_x,
                    (l - l0) // l_period - shift_y,
                )
                if point not in seen:
                    seen.add(point)
                    points.append(point)
        args.points_output.write_text(
            json.dumps(points) + "\n"
        )
        print(
            f"input_points={input_count} matching={matching_count} "
            f"unique={len(points)} output={args.points_output}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
