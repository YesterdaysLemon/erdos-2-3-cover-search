#!/usr/bin/env python3
"""Independent verifier for a phase-partitioned sparse-radius obstruction."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from verify_sparse_radius_obstruction import verify_certificate


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


def verify_manifest(certificate_path: Path) -> dict:
    started = time.monotonic()
    certificate = json.loads(certificate_path.read_text())
    pool_path = Path(certificate["pool"])
    base_phase_path = Path(certificate["base_phase"])
    if sha256_file(pool_path) != certificate["pool_sha256"]:
        raise RuntimeError("pool hash mismatch")
    if (
        sha256_file(base_phase_path)
        != certificate["base_phase_sha256"]
    ):
        raise RuntimeError("base-phase hash mismatch")

    pool = json.loads(pool_path.read_text())
    rows = pool["choices"]
    row_by_prime = {int(row["p"]): row for row in rows}
    if len(row_by_prime) != len(rows):
        raise RuntimeError("pool repeats a prime")
    raw_phases = {
        int(prime): int(target)
        for prime, target in json.loads(
            base_phase_path.read_text()
        ).items()
    }
    global_phases = {}
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        target = int(raw_phases.get(prime, residue)) % h
        if not is_legal_target(row, target):
            raise RuntimeError(f"illegal root phase p={prime}")
        global_phases[prime] = target

    root_fixed = set(map(int, certificate["fixed_primes"]))
    if root_fixed - set(row_by_prime):
        raise RuntimeError("a root fixed prime is absent from the pool")
    max_changes = int(certificate["max_changes"])
    if max_changes < 0:
        raise RuntimeError("negative root radius")
    leaf_results = []
    partition_count = 0
    over_budget_count = 0

    def visit(
        node: dict,
        path_targets: dict[int, int],
        changed_count: int,
    ) -> None:
        nonlocal partition_count, over_budget_count
        node_type = node.get("type")
        if node_type == "leaf":
            expected_fixed = root_fixed | set(path_targets)
            remaining = max_changes - changed_count
            if (
                node.get("complete") is not True
                or set(map(int, node["fixed_primes"]))
                != expected_fixed
                or int(node["max_changes"]) != remaining
            ):
                raise RuntimeError("leaf metadata differs from tree path")
            phases = dict(global_phases)
            phases.update(path_targets)
            embedded = {
                "complete": node["complete"],
                "points": node["points"],
                "fixed_primes": node["fixed_primes"],
                "max_changes": node["max_changes"],
                "discovery": node["discovery"],
            }
            result = verify_certificate(
                None,
                certificate_data=embedded,
                pool_data=pool,
                phase_data=phases,
            )
            if not result["verified"]:
                raise RuntimeError("embedded leaf failed scalar replay")
            leaf_results.append(result)
            return
        if node_type == "over_budget":
            if changed_count <= max_changes:
                raise RuntimeError("reachable branch marked over budget")
            over_budget_count += 1
            return
        if node_type != "partition":
            raise RuntimeError("unknown tree node type")

        partition_count += 1
        prime = int(node["prime"])
        if prime in root_fixed or prime in path_targets:
            raise RuntimeError(f"tree repeats or pre-fixes p={prime}")
        row = row_by_prime.get(prime)
        if row is None:
            raise RuntimeError(f"tree prime p={prime} is absent")
        targets = legal_targets(row)
        if (
            int(node["h"]) != int(row["h"])
            or list(map(int, node["legal_targets"])) != targets
            or int(node["base_target"]) != global_phases[prime]
        ):
            raise RuntimeError("partition metadata differs from pool")
        branches = node["branches"]
        if set(branches) != {str(target) for target in targets}:
            raise RuntimeError("partition omits a legal target")
        for target in targets:
            next_targets = dict(path_targets)
            next_targets[prime] = target
            next_changed = changed_count + int(
                target != global_phases[prime]
            )
            visit(
                branches[str(target)],
                next_targets,
                next_changed,
            )

    visit(certificate["tree"], {}, 0)
    count_match = len(leaf_results) == int(certificate["leaf_count"])
    verified = (
        certificate.get("complete") is True
        and count_match
        and all(result["verified"] for result in leaf_results)
    )
    return {
        "certificate": str(certificate_path),
        "verified": verified,
        "repair_exists": any(
            result["repair_exists"] for result in leaf_results
        ),
        "row_count": len(rows),
        "max_changes": max_changes,
        "partition_count": partition_count,
        "leaf_count": len(leaf_results),
        "declared_leaf_count": int(certificate["leaf_count"]),
        "over_budget_branch_count": over_budget_count,
        "leaf_count_match": count_match,
        "total_leaf_point_count": sum(
            result["point_count"] for result in leaf_results
        ),
        "total_full_skeleton_count": sum(
            result["full_skeleton_count"]
            for result in leaf_results
        ),
        "elapsed_seconds": time.monotonic() - started,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify_manifest(args.certificate)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"verified={result['verified']} "
        f"repair_exists={result['repair_exists']} "
        f"partitions={result['partition_count']} "
        f"leaves={result['leaf_count']} "
        f"leaf_points={result['total_leaf_point_count']} "
        f"elapsed_s={result['elapsed_seconds']:.3f}",
        flush=True,
    )
    return 0 if result["verified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
