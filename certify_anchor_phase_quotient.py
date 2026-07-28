#!/usr/bin/env python3
"""Build a finite exact obstruction for a frozen-remainder phase quotient.

The certificate exhausts every legal target assignment for a selected set of
anchor rows.  All other row phases remain fixed.  For each anchor branch it
records one supplied exponent pair that is outside every affine fibre and
outside every algebraically covered origin sublattice.

This is a local statement about the embedded finite row family and frozen
nonanchor phases.  It is not a no-cover theorem for arbitrary phase maps.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path


CELL_KEYS = ("cells", "exact_misses", "misses", "new_cells")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_cells(path: Path) -> list[tuple[int, int]]:
    payload = json.loads(path.read_text())
    if isinstance(payload, dict):
        key = next((name for name in CELL_KEYS if name in payload), None)
        if key is None:
            raise ValueError(f"{path} has no supported point-list field")
        payload = payload[key]
    cells = []
    for raw_cell in payload:
        if not isinstance(raw_cell, list) or len(raw_cell) != 2:
            raise ValueError(f"{path} contains a malformed point")
        k, l = map(int, raw_cell)
        if k < 0 or l < 0:
            raise ValueError(f"{path} contains a negative exponent")
        cells.append((k, l))
    return cells


def legal_targets(row: dict) -> tuple[int, ...]:
    h = int(row["h"])
    modulus = int(row.get("target_modulus", 1))
    if h < 1 or modulus < 1 or h % modulus:
        raise ValueError("invalid row target restriction")
    residue = int(row.get("target_residue", 0)) % modulus
    return tuple(target for target in range(h) if target % modulus == residue)


def row_covers(
    row: dict,
    target: int,
    point: tuple[int, int],
) -> bool:
    k, l = point
    h = int(row["h"])
    return (
        int(row["a"]) * k + int(row["b"]) * l - target
    ) % h == 0


def build_certificate(
    pool_path: Path,
    phase_path: Path,
    point_paths: list[Path],
    anchor_primes: tuple[int, ...],
) -> dict:
    pool = json.loads(pool_path.read_text())
    rows = pool["choices"]
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(phase_path.read_text()).items()
    }
    if len(set(anchor_primes)) != len(anchor_primes):
        raise ValueError("anchor primes must be distinct")

    row_by_prime = {}
    for row in rows:
        prime = int(row["p"])
        if prime in row_by_prime:
            raise ValueError(f"pool repeats row prime {prime}")
        row_by_prime[prime] = row
    if set(row_by_prime) - phases.keys():
        raise ValueError("phase file omits at least one pool row")
    if set(anchor_primes) - row_by_prime.keys():
        raise ValueError("an anchor prime is absent from the pool")

    normalized_phases = {}
    compact_rows = []
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        target = phases[prime] % h
        if target not in legal_targets(row):
            raise ValueError(f"base phase for p={prime} is forbidden")
        normalized_phases[prime] = target
        compact_rows.append(
            [
                prime,
                h,
                int(row["a"]) % h,
                int(row["b"]) % h,
                target,
                int(row.get("target_modulus", 1)),
                int(row.get("target_residue", 0)),
            ]
        )

    points = []
    seen_points = set()
    for point_path in point_paths:
        for point in load_cells(point_path):
            if point not in seen_points:
                seen_points.add(point)
                points.append(point)

    algebraic_primes = tuple(
        int(prime) for prime in pool.get("algebraic_primes", ())
    )
    anchor_set = set(anchor_primes)
    nonanchor_rows = [
        row for row in rows if int(row["p"]) not in anchor_set
    ]
    eligible = []
    algebraically_covered = 0
    nonanchor_covered = 0
    for point in points:
        k, l = point
        if any(
            k % prime == 0 and l % prime == 0
            for prime in algebraic_primes
        ):
            algebraically_covered += 1
            continue
        if any(
            row_covers(
                row,
                normalized_phases[int(row["p"])],
                point,
            )
            for row in nonanchor_rows
        ):
            nonanchor_covered += 1
            continue
        signature = tuple(
            (
                int(row_by_prime[prime]["a"]) * k
                + int(row_by_prime[prime]["b"]) * l
            )
            % int(row_by_prime[prime]["h"])
            for prime in anchor_primes
        )
        eligible.append((point, signature))

    target_spaces = [
        legal_targets(row_by_prime[prime])
        for prime in anchor_primes
    ]
    branches = []
    open_branches = []
    for targets in itertools.product(*target_spaces):
        witnesses = [
            point
            for point, signature in eligible
            if all(
                target != covered_target
                for target, covered_target in zip(targets, signature)
            )
        ]
        if not witnesses:
            open_branches.append(list(targets))
            continue
        branches.append(
            {
                "targets": list(targets),
                "witness": list(witnesses[0]),
                "available_witnesses": len(witnesses),
            }
        )
    if open_branches:
        preview = open_branches[:10]
        raise RuntimeError(
            f"{len(open_branches)} anchor branches lack a witness: {preview}"
        )

    source_points = [
        {
            "path": path.name,
            "sha256": sha256_file(path),
        }
        for path in point_paths
    ]
    anchor_metadata = []
    for prime in anchor_primes:
        row = row_by_prime[prime]
        anchor_metadata.append(
            {
                "p": prime,
                "h": int(row["h"]),
                "a": int(row["a"]) % int(row["h"]),
                "b": int(row["b"]) % int(row["h"]),
                "target_modulus": int(row.get("target_modulus", 1)),
                "target_residue": (
                    int(row.get("target_residue", 0))
                    % int(row.get("target_modulus", 1))
                ),
                "legal_targets": list(legal_targets(row)),
            }
        )

    return {
        "schema_version": 1,
        "claim": (
            "Every legal target assignment for the listed anchor rows has "
            "an explicit uncovered exponent pair when every nonanchor row "
            "remains at its embedded base phase."
        ),
        "scope_warning": (
            "This frozen-remainder finite quotient is not a proof that the "
            "full row pool, a larger pool, or the Erdos problem has no cover."
        ),
        "source": {
            "pool_path": pool_path.name,
            "pool_sha256": sha256_file(pool_path),
            "phase_path": phase_path.name,
            "phase_sha256": sha256_file(phase_path),
            "point_files": source_points,
            "distinct_input_points": len(points),
            "eligible_points": len(eligible),
            "algebraically_covered_input_points": algebraically_covered,
            "nonanchor_covered_input_points": nonanchor_covered,
        },
        "row_columns": [
            "p",
            "h",
            "a",
            "b",
            "base_phase",
            "target_modulus",
            "target_residue",
        ],
        "rows": compact_rows,
        "algebraic_primes": list(algebraic_primes),
        "anchors": anchor_metadata,
        "branches": branches,
        "summary": {
            "rows": len(rows),
            "nonanchor_rows": len(nonanchor_rows),
            "anchor_rows": len(anchor_primes),
            "legal_anchor_branches": len(branches),
            "open_anchor_branches": 0,
            "minimum_available_witnesses": min(
                branch["available_witnesses"] for branch in branches
            ),
            "maximum_available_witnesses": max(
                branch["available_witnesses"] for branch in branches
            ),
            "engine": "python-exact-frozen-remainder-incidence-enumeration",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument("--points-file", type=Path, action="append", required=True)
    parser.add_argument("--anchors", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    anchors = tuple(
        int(value) for value in args.anchors.split(",") if value
    )
    if not anchors:
        parser.error("--anchors must list at least one row prime")
    certificate = build_certificate(
        args.pool,
        args.phase_file,
        args.points_file,
        anchors,
    )
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    summary = certificate["summary"]
    print(
        f"rows={summary['rows']} anchors={summary['anchor_rows']} "
        f"branches={summary['legal_anchor_branches']} "
        f"open={summary['open_anchor_branches']} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
