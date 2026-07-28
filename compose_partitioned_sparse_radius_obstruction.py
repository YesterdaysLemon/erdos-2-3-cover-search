#!/usr/bin/env python3
"""Compose sparse-radius leaf certificates into an exhaustive phase tree."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def legal_targets(row: dict) -> list[int]:
    h = int(row["h"])
    modulus = int(row.get("target_modulus", 1))
    residue = int(row.get("target_residue", 0)) % modulus
    if h < 1 or modulus < 1 or h % modulus:
        raise RuntimeError(f"malformed row p={row['p']}")
    return list(range(residue, h, modulus))


def is_legal_target(row: dict, target: int) -> bool:
    h = int(row["h"])
    modulus = int(row.get("target_modulus", 1))
    residue = int(row.get("target_residue", 0)) % modulus
    return (
        h >= 1
        and modulus >= 1
        and h % modulus == 0
        and 0 <= target < h
        and target % modulus == residue
    )


def effective_phases(rows: list[dict], phase_path: Path) -> dict[int, int]:
    raw = {
        int(prime): int(target)
        for prime, target in json.loads(phase_path.read_text()).items()
    }
    result = {}
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        target = int(raw.get(prime, residue)) % h
        if not is_legal_target(row, target):
            raise RuntimeError(f"illegal phase p={prime} target={target}")
        result[prime] = target
    return result


def compose_manifest(
    pool_path: Path,
    base_phase_path: Path,
    tree_spec_path: Path,
    fixed_primes: set[int],
    max_changes: int,
) -> dict:
    if max_changes < 0:
        raise ValueError("negative change radius")
    pool_sha256 = sha256_file(pool_path)
    pool = json.loads(pool_path.read_text())
    rows = pool["choices"]
    row_by_prime = {int(row["p"]): row for row in rows}
    if len(row_by_prime) != len(rows):
        raise RuntimeError("pool repeats a prime")
    if fixed_primes - set(row_by_prime):
        raise RuntimeError("a root fixed prime is absent from the pool")
    global_phases = effective_phases(rows, base_phase_path)
    spec = json.loads(tree_spec_path.read_text())
    leaf_count = 0

    def visit(
        node: dict,
        path_targets: dict[int, int],
        changed_count: int,
    ) -> dict:
        nonlocal leaf_count
        if "certificate" in node:
            certificate_path = Path(node["certificate"])
            certificate = json.loads(certificate_path.read_text())
            expected_fixed = fixed_primes | set(path_targets)
            remaining = max_changes - changed_count
            if remaining < 0:
                raise RuntimeError("an over-budget path has a leaf")
            if certificate.get("pool_sha256") != pool_sha256:
                raise RuntimeError("leaf pool hash differs from root")
            if set(map(int, certificate["fixed_primes"])) != expected_fixed:
                raise RuntimeError("leaf fixed primes differ from tree path")
            if int(certificate["max_changes"]) != remaining:
                raise RuntimeError("leaf radius differs from remaining budget")
            leaf_phase_path = Path(certificate["base_phase"])
            if (
                sha256_file(leaf_phase_path)
                != certificate["base_phase_sha256"]
            ):
                raise RuntimeError("leaf base-phase hash mismatch")
            leaf_phases = effective_phases(rows, leaf_phase_path)
            expected_phases = dict(global_phases)
            expected_phases.update(path_targets)
            if leaf_phases != expected_phases:
                raise RuntimeError("leaf phase map differs from tree path")
            if (
                certificate.get("complete") is not True
                or certificate["discovery"].get("status") != "UNSAT"
                or certificate["discovery"].get("complete_negative")
                is not True
            ):
                raise RuntimeError("leaf is not a complete obstruction")
            leaf_count += 1
            return {
                "type": "leaf",
                "complete": True,
                "fixed_primes": sorted(expected_fixed),
                "max_changes": remaining,
                "points": certificate["points"],
                "discovery": certificate["discovery"],
                "source_certificate_sha256": sha256_file(certificate_path),
            }

        prime = int(node["prime"])
        if prime in fixed_primes or prime in path_targets:
            raise RuntimeError(f"tree repeats or pre-fixes p={prime}")
        row = row_by_prime.get(prime)
        if row is None:
            raise RuntimeError(f"tree prime p={prime} is absent")
        targets = legal_targets(row)
        branches = node["branches"]
        if set(branches) != {str(target) for target in targets}:
            raise RuntimeError(f"tree does not partition every phase of p={prime}")
        composed_branches = {}
        for target in targets:
            child = branches[str(target)]
            next_changed = changed_count + int(
                target != global_phases[prime]
            )
            if next_changed > max_changes:
                if child != {"over_budget": True}:
                    raise RuntimeError("over-budget branch is not marked")
                composed_branches[str(target)] = {
                    "type": "over_budget"
                }
                continue
            next_targets = dict(path_targets)
            next_targets[prime] = target
            composed_branches[str(target)] = visit(
                child,
                next_targets,
                next_changed,
            )
        return {
            "type": "partition",
            "prime": prime,
            "h": int(row["h"]),
            "legal_targets": targets,
            "base_target": global_phases[prime],
            "branches": composed_branches,
        }

    tree = visit(spec, {}, 0)
    return {
        "problem": (
            "finite affine phase repair obstruction partitioned by "
            "selected row phases"
        ),
        "complete": True,
        "scope": (
            "the embedded finite row pool, global base phase, root fixed "
            "rows, and declared Hamming radius; every legal target of each "
            "partition row is represented"
        ),
        "pool": str(pool_path),
        "pool_sha256": pool_sha256,
        "base_phase": str(base_phase_path),
        "base_phase_sha256": sha256_file(base_phase_path),
        "fixed_primes": sorted(fixed_primes),
        "max_changes": max_changes,
        "tree_spec_sha256": sha256_file(tree_spec_path),
        "leaf_count": leaf_count,
        "tree": tree,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--base-phase", type=Path, required=True)
    parser.add_argument("--tree-spec", type=Path, required=True)
    parser.add_argument("--fixed-primes", required=True)
    parser.add_argument("--max-changes", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = compose_manifest(
        args.pool,
        args.base_phase,
        args.tree_spec,
        {
            int(value)
            for value in args.fixed_primes.split(",")
            if value
        },
        args.max_changes,
    )
    args.output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"COMPOSED_PARTITIONED_RADIUS_OBSTRUCTION "
        f"radius={manifest['max_changes']} "
        f"leaves={manifest['leaf_count']} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
