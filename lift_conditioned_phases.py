#!/usr/bin/env python3
"""Lift phases from conditioned cell coordinates back to original exponents."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from power_anchor_capacity_lp import power_target_congruence


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("conditioned_pool", type=Path)
    parser.add_argument("conditioned_phases", type=Path)
    parser.add_argument("baseline_original_phases", type=Path)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path)
    args = parser.parse_args()

    payload = json.loads(args.conditioned_pool.read_text())
    local_phases = {
        int(prime): int(target)
        for prime, target in json.loads(
            args.conditioned_phases.read_text()
        ).items()
    }
    output_phases = {
        int(prime): int(target)
        for prime, target in json.loads(
            args.baseline_original_phases.read_text()
        ).items()
    }
    cell = payload["cell_condition"]
    k0 = int(cell["k_residue"])
    l0 = int(cell["l_residue"])
    lifted = []
    for row in payload["choices"]:
        prime = int(row["p"])
        if prime not in local_phases:
            raise RuntimeError(f"missing conditioned phase for p={prime}")
        h = int(row["h"])
        target = local_phases[prime] % h
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        if target % modulus != residue:
            raise RuntimeError(
                f"conditioned phase violates target restriction for p={prime}"
            )

        original_h = int(row["original_h"])
        common = original_h // h
        original_a = int(row["original_a"])
        original_b = int(row["original_b"])
        base = (original_a * k0 + original_b * l0) % original_h
        shift = int(row["coordinate_shift"])
        original_target = (
            base + common * ((target + shift) % h)
        ) % original_h
        power_residue, power_modulus = power_target_congruence(
            original_h,
            prime,
            args.power,
        )
        if original_target % power_modulus != power_residue:
            raise AssertionError(
                f"lifted phase is not power-compatible for p={prime}"
            )
        # Reproject exactly, independently of any baseline value.
        delta = (original_target - base) % original_h
        if delta % common:
            raise AssertionError("lifted row misses the conditioned cell")
        reprojected = (delta // common - shift) % h
        if reprojected != target:
            raise AssertionError("conditioned lift/reprojection mismatch")
        output_phases[prime] = original_target
        lifted.append(
            {
                "p": prime,
                "conditioned_target": target,
                "original_target": original_target,
            }
        )

    args.output.write_text(
        json.dumps(
            {str(prime): target for prime, target in output_phases.items()}
        )
        + "\n"
    )
    if args.audit_output:
        args.audit_output.write_text(
            json.dumps(
                {
                    "conditioned_pool": str(args.conditioned_pool),
                    "conditioned_phases": str(args.conditioned_phases),
                    "baseline_original_phases": str(
                        args.baseline_original_phases
                    ),
                    "power": args.power,
                    "lifted_count": len(lifted),
                    "lifted": lifted,
                    "output": str(args.output),
                },
                indent=2,
            )
            + "\n"
        )
    print(
        f"lifted={len(lifted)} total_phases={len(output_phases)} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
