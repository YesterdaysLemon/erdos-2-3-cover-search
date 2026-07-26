#!/usr/bin/env python3
"""Lift a certified child-cell phase core to its parent affine-fibre pool.

For a conditioned row, condition_derived_cell.py records

    child_c = (parent_c - parent_base) / parent_common
              - coordinate_shift        (mod child_h).

Consequently the parent phase is recovered exactly as

    parent_c = parent_base
               + parent_common * (child_c + coordinate_shift)
               (mod parent_h).

Only rows named by the supplied core are changed.  This makes the resulting
phase assignment suitable for monotone certificate locking: phases outside
the new proof core remain byte-for-byte identical to the base assignment.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def unique_by_prime(rows: list[dict], label: str) -> dict[int, dict]:
    indexed: dict[int, dict] = {}
    for row in rows:
        prime = int(row["p"])
        if prime in indexed:
            raise RuntimeError(f"duplicate prime {prime} in {label}")
        indexed[prime] = row
    return indexed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("parent_pool", type=Path)
    parser.add_argument("child_pool", type=Path)
    parser.add_argument("child_cover", type=Path)
    parser.add_argument("core", type=Path)
    parser.add_argument("base_phases", type=Path)
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--rows-output", type=Path, required=True)
    args = parser.parse_args()

    parent_payload = json.loads(args.parent_pool.read_text())
    child_payload = json.loads(args.child_pool.read_text())
    cover_payload = json.loads(args.child_cover.read_text())
    core_payload = json.loads(args.core.read_text())
    base_phases = {
        int(prime): int(phase)
        for prime, phase in json.loads(args.base_phases.read_text()).items()
    }

    parent_rows = unique_by_prime(parent_payload["choices"], "parent pool")
    child_rows = unique_by_prime(child_payload["choices"], "child pool")
    cover_rows = cover_payload["choices"]
    core_indices = [int(index) for index in core_payload["row_indices"]]
    if len(core_indices) != len(set(core_indices)):
        raise RuntimeError("core row indices are not unique")

    lifted_rows = []
    lifted_phases: dict[int, int] = {}
    for index in core_indices:
        try:
            cover_row = cover_rows[index]
        except IndexError as exc:
            raise RuntimeError(f"core index {index} is outside child cover") from exc
        prime = int(cover_row["p"])
        child = child_rows.get(prime)
        parent = parent_rows.get(prime)
        if child is None or parent is None:
            raise RuntimeError(f"prime {prime} is missing from a source pool")

        for field in ("h", "a", "b"):
            if int(cover_row[field]) != int(child[field]):
                raise RuntimeError(
                    f"child cover/pool mismatch for prime {prime}, field {field}"
                )
        for child_field, parent_field in (
            ("parent_h", "h"),
            ("parent_a", "a"),
            ("parent_b", "b"),
        ):
            if int(child[child_field]) != int(parent[parent_field]):
                raise RuntimeError(
                    f"child/parent mismatch for prime {prime}, "
                    f"field {child_field}"
                )

        child_h = int(child["h"])
        child_phase = int(cover_row["c"]) % child_h
        parent_h = int(child["parent_h"])
        parent_common = int(child["parent_common"])
        if parent_h != parent_common * child_h:
            raise RuntimeError(f"invalid modulus quotient for prime {prime}")
        parent_phase = (
            int(child["parent_base"])
            + parent_common
            * (child_phase + int(child["coordinate_shift"]))
        ) % parent_h
        if prime not in base_phases:
            raise RuntimeError(f"prime {prime} is absent from base phase file")

        lifted_phases[prime] = parent_phase
        lifted = dict(parent)
        lifted["c"] = parent_phase
        lifted_rows.append(lifted)

    output_phases = dict(base_phases)
    output_phases.update(lifted_phases)
    args.phase_output.write_text(
        json.dumps(
            {str(prime): phase for prime, phase in output_phases.items()}
        )
        + "\n"
    )
    rows_payload = {
        "parent_pool": str(args.parent_pool),
        "child_pool": str(args.child_pool),
        "child_cover": str(args.child_cover),
        "child_core": str(args.core),
        "base_phases": str(args.base_phases),
        "cell_condition": child_payload["cell_condition"],
        "algebraic_primes": parent_payload.get("algebraic_primes", []),
        "sophie_germain": bool(parent_payload.get("sophie_germain", False)),
        "rows": lifted_rows,
        "lifted_phases": {
            str(prime): phase for prime, phase in lifted_phases.items()
        },
    }
    args.rows_output.write_text(json.dumps(rows_payload, indent=2) + "\n")
    print(
        f"core_rows={len(lifted_rows)} changed_primes={len(lifted_phases)} "
        f"phase_output={args.phase_output} rows_output={args.rows_output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
