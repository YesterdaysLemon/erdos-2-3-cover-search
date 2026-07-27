#!/usr/bin/env python3
"""Package a complete finite sparse-radius obstruction certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import exact_greedy
from certify_anchor_phase_quotient import load_cells
from finite_sample_mask_repair import search_mask_repair


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path, nargs="+")
    parser.add_argument("--base-phase", type=Path, required=True)
    parser.add_argument("--fixed-primes", required=True)
    parser.add_argument("--max-changes", type=int, required=True)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    candidates = exact_greedy.load_candidates(args.pool, True)
    points = []
    seen_points = set()
    for path in args.points:
        for point in load_cells(path):
            if point not in seen_points:
                seen_points.add(point)
                points.append(point)
    initial = {
        int(prime): int(target)
        for prime, target in json.loads(
            args.base_phase.read_text()
        ).items()
    }
    fixed_primes = {
        int(value) for value in args.fixed_primes.split(",") if value
    }
    result = search_mask_repair(
        rows,
        candidates,
        points,
        initial,
        fixed_primes,
        args.max_changes,
        args.solver,
        0,
        0.0,
    )
    if result["status"] != "UNSAT" or not result["complete_negative"]:
        raise RuntimeError(
            "the supplied finite instance is not a complete obstruction"
        )
    certificate = {
        "problem": (
            "finite affine phase repair obstruction within a Hamming ball"
        ),
        "complete": True,
        "scope": (
            "the embedded finite row pool, base phase, fixed rows, supplied "
            "points, and declared maximum number of changed nonfixed rows"
        ),
        "pool": str(args.pool),
        "pool_sha256": sha256_file(args.pool),
        "base_phase": str(args.base_phase),
        "base_phase_sha256": sha256_file(args.base_phase),
        "point_sources": [
            {
                "path": str(path),
                "sha256": sha256_file(path),
            }
            for path in args.points
        ],
        "points": [[k, l] for k, l in points],
        "fixed_primes": sorted(fixed_primes),
        "max_changes": args.max_changes,
        "discovery": {
            key: value
            for key, value in result.items()
            if key not in {"repaired"}
        },
    }
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    print(
        f"CERTIFIED_SPARSE_RADIUS_OBSTRUCTION "
        f"points={len(points)} initial_misses={result['initial_misses']} "
        f"radius={args.max_changes} masks={result['gain_mask_count']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
