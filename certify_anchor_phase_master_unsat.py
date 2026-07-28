#!/usr/bin/env python3
"""Build a self-contained finite UNSAT certificate for an anchor quotient."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import anchor_phase_master
from certify_anchor_phase_quotient import load_cells


ROW_COLUMNS = [
    "p",
    "h",
    "a",
    "b",
    "base_phase",
    "target_modulus",
    "target_residue",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_certificate(
    pool_path: Path,
    phase_path: Path,
    points_path: Path,
    anchor_primes: tuple[int, ...],
) -> dict:
    payload = json.loads(pool_path.read_text())
    rows = payload["choices"]
    base_phases = {
        int(prime): int(target)
        for prime, target in json.loads(phase_path.read_text()).items()
    }
    points = list(dict.fromkeys(load_cells(points_path)))
    if len(points) != len(load_cells(points_path)):
        raise ValueError("core point file contains duplicates")
    algebraic_primes = tuple(
        int(prime) for prime in payload.get("algebraic_primes", ())
    )
    solved_phases, master = anchor_phase_master.solve_anchor_master(
        rows,
        base_phases,
        anchor_primes,
        points,
        algebraic_primes,
    )
    if solved_phases is not None or master["sat"]:
        raise RuntimeError("supplied anchor quotient is SAT")
    if master["eligible_points"] != len(points):
        raise RuntimeError(
            "every certificate point must require an anchor-row cover"
        )

    compact_rows, normalized_phases = anchor_phase_master.normalize_rows(
        rows,
        base_phases,
    )
    compact_with_phases = [
        [
            prime,
            h,
            a,
            b,
            normalized_phases[prime],
            modulus,
            residue,
        ]
        for prime, h, a, b, modulus, residue in compact_rows
    ]
    row_by_prime = {row[0]: row for row in compact_rows}
    anchor_metadata = []
    target_space_sizes = []
    for prime in anchor_primes:
        _prime, h, a, b, modulus, residue = row_by_prime[prime]
        legal = [
            target for target in range(h) if target % modulus == residue
        ]
        target_space_sizes.append(len(legal))
        anchor_metadata.append(
            {
                "p": prime,
                "h": h,
                "a": a,
                "b": b,
                "target_modulus": modulus,
                "target_residue": residue,
                "legal_target_count": len(legal),
            }
        )
    stable_master = {
        key: master[key]
        for key in (
            "sat",
            "rows",
            "anchor_rows",
            "frozen_rows",
            "input_points",
            "eligible_points",
            "algebraically_discharged_points",
            "frozen_discharged_points",
            "empty_clause_point",
            "target_options",
            "variables",
            "clauses",
            "engine",
        )
    }
    return {
        "schema_version": 1,
        "claim": (
            "No legal joint target assignment for the listed anchor rows "
            "covers all recorded exponent pairs when every nonanchor row "
            "remains at its embedded base phase."
        ),
        "scope_warning": (
            "This is an exact finite obstruction for one frozen remainder. "
            "It is not a no-cover theorem for arbitrary phases, the full "
            "candidate universe, or the Erdos problem."
        ),
        "source": {
            "pool_path": pool_path.name,
            "pool_sha256": sha256_file(pool_path),
            "phase_path": phase_path.name,
            "phase_sha256": sha256_file(phase_path),
            "points_path": points_path.name,
            "points_sha256": sha256_file(points_path),
        },
        "row_columns": ROW_COLUMNS,
        "rows": compact_with_phases,
        "algebraic_primes": list(algebraic_primes),
        "anchors": anchor_metadata,
        "points": [list(point) for point in points],
        "master": stable_master,
        "summary": {
            "rows": len(rows),
            "anchor_rows": len(anchor_primes),
            "frozen_rows": len(rows) - len(anchor_primes),
            "core_points": len(points),
            "legal_target_space_sizes": target_space_sizes,
            "joint_legal_target_assignments": math.prod(target_space_sizes),
            "pysat_unsat": True,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--base-phases", type=Path, required=True)
    parser.add_argument("--points", type=Path, required=True)
    parser.add_argument("--anchors", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    anchors = tuple(
        int(value) for value in args.anchors.split(",") if value
    )
    certificate = build_certificate(
        args.pool,
        args.base_phases,
        args.points,
        anchors,
    )
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    summary = certificate["summary"]
    print(
        f"UNSAT anchors={summary['anchor_rows']} "
        f"points={summary['core_points']} "
        f"assignments={summary['joint_legal_target_assignments']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
