#!/usr/bin/env python3
"""Lift a phase assignment from rectangle coordinates to its parent pool."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("parent", type=Path)
    parser.add_argument("rectangle", type=Path)
    parser.add_argument("child_phases", type=Path)
    parser.add_argument("--extra-phase", action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    parent = json.loads(args.parent.read_text())
    rectangle = json.loads(args.rectangle.read_text())
    child_phases = {
        int(prime): int(value)
        for prime, value in json.loads(args.child_phases.read_text()).items()
    }
    children = {
        int(row["p"]): row for row in rectangle["choices"]
    }
    parents = {int(row["p"]): row for row in parent["choices"]}
    if set(child_phases) != set(children):
        raise RuntimeError("child phase prime set does not match rectangle")

    phases = {}
    for prime, child in children.items():
        child_h = int(child["h"])
        child_phase = child_phases[prime] % child_h
        if (
            child_phase % int(child["target_modulus"])
            != int(child["target_residue"])
        ):
            raise RuntimeError(f"illegal child phase for prime {prime}")
        parent_row = parents[prime]
        parent_h = int(parent_row["h"])
        parent_phase = (
            int(child["parent_base"])
            + int(child["parent_common"])
            * (child_phase + int(child["coordinate_shift"]))
        ) % parent_h
        if (
            parent_phase % int(parent_row["target_modulus"])
            != int(parent_row["target_residue"])
        ):
            raise AssertionError(
                f"lifted phase violates parent target for prime {prime}"
            )
        phases[prime] = parent_phase

    for item in args.extra_phase:
        prime_text, separator, phase_text = item.partition(":")
        if not separator:
            raise SystemExit("--extra-phase must have form prime:phase")
        prime = int(prime_text)
        phase = int(phase_text)
        if prime in phases:
            raise RuntimeError(f"duplicate lifted/extra phase for {prime}")
        if prime not in parents:
            raise RuntimeError(f"extra phase prime {prime} is not in parent")
        row = parents[prime]
        phase %= int(row["h"])
        if (
            phase % int(row["target_modulus"])
            != int(row["target_residue"])
        ):
            raise RuntimeError(f"illegal extra phase for prime {prime}")
        phases[prime] = phase

    if set(phases) != set(parents):
        missing = sorted(set(parents) - set(phases))
        raise RuntimeError(
            f"lift does not assign every parent row; missing {missing[:10]}"
        )
    args.output.write_text(
        json.dumps({str(prime): phases[prime] for prime in phases}) + "\n"
    )
    print(
        f"lifted={len(children)} extra={len(args.extra_phase)} "
        f"parent={len(phases)} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
