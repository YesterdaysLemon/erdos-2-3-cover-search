#!/usr/bin/env python3
"""Independently reproject and verify a lifted phase patch."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def rows(payload: dict) -> dict[int, dict]:
    return {int(row["p"]): row for row in payload["choices"]}


def locked_values(path_value: str | None) -> set[int]:
    if not path_value:
        return set()
    payload = json.loads(Path(path_value).read_text())
    values = payload.get("primes", ()) if isinstance(payload, dict) else payload
    return {int(value) for value in values}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("audit", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    audit = json.loads(args.audit.read_text())
    base = {
        int(prime): int(target)
        for prime, target in json.loads(
            Path(audit["base_root_phases"]).read_text()
        ).items()
    }
    patched = {
        int(prime): int(target)
        for prime, target in json.loads(
            Path(audit["output"]).read_text()
        ).items()
    }
    selected = [int(prime) for prime in audit["selected_primes"]]
    locked = locked_values(audit.get("locked_primes"))

    unchanged_outside_patch = all(
        patched.get(prime) == target
        for prime, target in base.items()
        if prime not in selected
    ) and set(patched) == set(base)
    selected_disjoint_locked = not (set(selected) & locked)
    changed = [
        prime for prime in selected if patched[prime] != base[prime]
    ]
    changed_matches = changed == [
        int(prime) for prime in audit["changed_primes"]
    ]

    phases = {prime: patched[prime] for prime in selected}
    chain_valid = True
    for step in reversed(audit["steps_leaf_to_root"]):
        child_payload = json.loads(Path(step["child_pool"]).read_text())
        child_rows = rows(child_payload)
        projected = {}
        for prime, parent_target in phases.items():
            child = child_rows.get(prime)
            if child is None:
                chain_valid = False
                break
            common = int(child["parent_common"])
            delta = (
                parent_target - int(child["parent_base"])
            ) % int(child["parent_h"])
            if delta % common:
                chain_valid = False
                break
            child_target = (
                delta // common - int(child["coordinate_shift"])
            ) % int(child["h"])
            projected[prime] = child_target
        if not chain_valid:
            break
        phases = projected

    leaf_targets = {
        int(prime): int(target)
        for prime, target in json.loads(
            Path(audit["leaf_phases"]).read_text()
        ).items()
    }
    leaf_targets_match = chain_valid and all(
        phases[prime] == leaf_targets[prime] for prime in selected
    )
    checks = {
        "output_prime_set_matches_base": set(patched) == set(base),
        "unchanged_outside_selected_patch": unchanged_outside_patch,
        "selected_primes_are_not_locked": selected_disjoint_locked,
        "changed_prime_list_matches": changed_matches,
        "conditioning_chain_reprojects": chain_valid,
        "leaf_targets_match": leaf_targets_match,
    }
    passed = all(checks.values())
    result = {
        "audit": str(args.audit),
        "checks": checks,
        "recomputed": {
            "levels": len(audit["steps_leaf_to_root"]),
            "selected_primes": len(selected),
            "changed_primes": changed,
        },
        "verified": passed,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS" if passed else "FAIL", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
