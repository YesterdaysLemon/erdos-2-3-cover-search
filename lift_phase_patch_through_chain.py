#!/usr/bin/env python3
"""Lift a small phase patch through a chain of derived conditioned pools."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def rows_by_prime(payload: dict) -> dict[int, dict]:
    return {int(row["p"]): row for row in payload["choices"]}


def read_locked(path: Path | None) -> set[int]:
    if path is None:
        return set()
    payload = json.loads(path.read_text())
    values = payload.get("primes", ()) if isinstance(payload, dict) else payload
    return {int(value) for value in values}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("leaf_pool", type=Path)
    parser.add_argument("leaf_phases", type=Path)
    parser.add_argument("root_pool", type=Path)
    parser.add_argument("base_root_phases", type=Path)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--primes")
    selection.add_argument(
        "--all-unlocked-leaf-primes",
        action="store_true",
        help=(
            "lift every prime in the leaf phase file except primes declared "
            "by --locked-primes; avoids command-line length limits for a "
            "whole-pool patch"
        ),
    )
    parser.add_argument("--locked-primes", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path, required=True)
    args = parser.parse_args()

    locked = read_locked(args.locked_primes)
    leaf_targets = {
        int(prime): int(target)
        for prime, target in json.loads(args.leaf_phases.read_text()).items()
    }
    selected = (
        sorted(set(leaf_targets) - locked)
        if args.all_unlocked_leaf_primes
        else [int(value) for value in args.primes.split(",") if value]
    )
    if len(selected) != len(set(selected)):
        raise RuntimeError("selected primes are not unique")
    overlap = set(selected) & locked
    if overlap:
        raise RuntimeError(f"phase patch changes locked primes: {sorted(overlap)}")

    missing = set(selected) - set(leaf_targets)
    if missing:
        raise RuntimeError(f"leaf targets missing primes: {sorted(missing)}")
    phases = {prime: leaf_targets[prime] for prime in selected}

    current_path = args.leaf_pool
    steps = []
    seen = set()
    while current_path.name != args.root_pool.name:
        resolved = current_path.resolve()
        if resolved in seen:
            raise RuntimeError("cycle in conditioned-pool source chain")
        seen.add(resolved)
        child_payload = json.loads(current_path.read_text())
        source = child_payload.get("source")
        if not source:
            raise RuntimeError(
                f"root pool not reached before source ended at {current_path}"
            )
        parent_path = Path(source)
        if not parent_path.exists():
            parent_path = current_path.parent / source
        parent_payload = json.loads(parent_path.read_text())
        child_rows = rows_by_prime(child_payload)
        parent_rows = rows_by_prime(parent_payload)

        lifted = {}
        records = []
        for prime, child_target in phases.items():
            child = child_rows.get(prime)
            parent = parent_rows.get(prime)
            if child is None or parent is None:
                raise RuntimeError(
                    f"prime {prime} missing at {current_path.name}"
                )
            for child_field, parent_field in (
                ("parent_h", "h"),
                ("parent_a", "a"),
                ("parent_b", "b"),
            ):
                if int(child[child_field]) != int(parent[parent_field]):
                    raise RuntimeError(
                        f"row mapping mismatch p={prime} at "
                        f"{current_path.name}"
                    )
            child_h = int(child["h"])
            parent_h = int(child["parent_h"])
            common = int(child["parent_common"])
            if parent_h != child_h * common:
                raise RuntimeError(
                    f"invalid modulus quotient p={prime} at "
                    f"{current_path.name}"
                )
            child_target %= child_h
            parent_target = (
                int(child["parent_base"])
                + common
                * (child_target + int(child["coordinate_shift"]))
            ) % parent_h
            lifted[prime] = parent_target
            records.append(
                {
                    "p": prime,
                    "child_target": child_target,
                    "parent_target": parent_target,
                }
            )
        steps.append(
            {
                "child_pool": str(current_path),
                "parent_pool": str(parent_path),
                "records": records,
            }
        )
        phases = lifted
        current_path = parent_path

    if current_path.name != args.root_pool.name:
        raise RuntimeError("root pool mismatch")
    root_payload = json.loads(args.root_pool.read_text())
    root_rows = rows_by_prime(root_payload)
    base = {
        int(prime): int(target)
        for prime, target in json.loads(args.base_root_phases.read_text()).items()
    }
    if set(selected) - set(base):
        raise RuntimeError("base root phase file misses selected primes")
    if set(selected) - set(root_rows):
        raise RuntimeError("root pool misses selected primes")
    before = {prime: base[prime] for prime in selected}
    for prime, target in phases.items():
        h = int(root_rows[prime]["h"])
        base[prime] = target % h
    changed = [
        prime for prime in selected if before[prime] != base[prime]
    ]
    if set(changed) & locked:
        raise AssertionError("a locked phase changed")

    args.output.write_text(
        json.dumps(
            {str(prime): target for prime, target in base.items()}
        )
        + "\n"
    )
    audit = {
        "leaf_pool": str(args.leaf_pool),
        "leaf_phases": str(args.leaf_phases),
        "root_pool": str(args.root_pool),
        "base_root_phases": str(args.base_root_phases),
        "locked_primes": (
            str(args.locked_primes) if args.locked_primes else None
        ),
        "selected_primes": selected,
        "changed_primes": changed,
        "root_targets_before": {
            str(prime): before[prime] for prime in selected
        },
        "root_targets_after": {
            str(prime): base[prime] for prime in selected
        },
        "steps_leaf_to_root": steps,
        "output": str(args.output),
    }
    args.audit_output.write_text(json.dumps(audit, indent=2) + "\n")
    print(
        f"levels={len(steps)} selected={len(selected)} "
        f"changed={len(changed)} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
