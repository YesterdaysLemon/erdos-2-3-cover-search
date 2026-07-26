#!/usr/bin/env python3
"""Project original perfect-power phases into a conditioned cell."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("conditioned_pool", type=Path)
    parser.add_argument("original_phases", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.conditioned_pool.read_text())
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(args.original_phases.read_text()).items()
    }
    cell = payload["cell_condition"]
    k0 = int(cell["k_residue"])
    l0 = int(cell["l_residue"])
    projected = {}
    retained = 0
    empty = 0
    absent = 0
    for row in payload["choices"]:
        prime = int(row["p"])
        h = int(row["h"])
        original_h = int(row["original_h"])
        common = original_h // h
        original_a = int(row["original_a"])
        original_b = int(row["original_b"])
        base = (original_a * k0 + original_b * l0) % original_h
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        shift = int(row["coordinate_shift"])
        if prime not in phases:
            target = residue
            absent += 1
        else:
            original_target = phases[prime] % original_h
            delta = (original_target - base) % original_h
            if delta % common:
                target = residue
                empty += 1
            else:
                target = (delta // common - shift) % h
                if target % modulus != residue:
                    raise AssertionError(
                        f"projected target violates restriction for prime {prime}"
                    )
                retained += 1
        projected[str(prime)] = target

    args.output.write_text(json.dumps(projected) + "\n")
    print(
        f"rows={len(projected)} retained={retained} empty={empty} "
        f"absent={absent} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
