#!/usr/bin/env python3
"""Package a compact finite-point obstruction for a small derived row pool."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--points", type=Path, required=True)
    parser.add_argument("--primary-result", type=Path, required=True)
    parser.add_argument("--max-h", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    rows = sorted(
        (
            row
            for row in source["choices"]
            if int(row["h"]) <= args.max_h
        ),
        key=lambda row: (int(row["h"]), int(row["p"])),
    )
    if not rows:
        raise RuntimeError("small derived pool is empty")
    points = [
        (int(k), int(l))
        for k, l in json.loads(args.points.read_text())
    ]
    if len(points) != len(set(points)):
        raise RuntimeError("point certificate contains duplicates")
    primary = json.loads(args.primary_result.read_text())
    if not primary.get("master_unsat"):
        raise RuntimeError("primary result does not report master UNSAT")
    if int(primary["master_points"]) != len(points):
        raise RuntimeError("primary point count does not match certificate")

    algebraic_primes = tuple(
        int(value) for value in source.get("algebraic_primes", ())
    )
    for point in points:
        if any(
            point[0] % prime == 0 and point[1] % prime == 0
            for prime in algebraic_primes
        ):
            raise RuntimeError(
                f"certificate point {point} is algebraically covered"
            )

    phase_assignment_count = 1
    for row in rows:
        h = int(row["h"])
        modulus = int(row.get("target_modulus", 1))
        if modulus < 1 or h % modulus:
            raise RuntimeError("invalid target restriction")
        phase_assignment_count *= h // modulus
    period = math.lcm(*(int(row["h"]) for row in rows))
    certificate = {
        "schema": "small_derived_pool_noncover_v1",
        "claim": (
            "No allowed phase assignment on the embedded derived rows covers "
            "all embedded points; every point also avoids the declared "
            "algebraic sublattices."
        ),
        "scope": {
            "source_pool": str(args.pool),
            "source_pool_sha256": file_sha256(args.pool),
            "selection": f"h <= {args.max_h}",
            "max_h": args.max_h,
            "period": period,
            "power": int(source["power"]),
            "cell_condition": source["cell_condition"],
            "row_count": len(rows),
            "point_count": len(points),
            "phase_assignment_count": phase_assignment_count,
            "algebraic_primes": list(algebraic_primes),
        },
        "primary_discovery": {
            "engine": "exact_cegis.py with PySAT master",
            "master_solver": primary.get("master_solver"),
            "master_unsat": True,
            "master_points": int(primary["master_points"]),
        },
        "rows": rows,
        "points": [[k, l] for k, l in points],
    }
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    print(
        f"rows={len(rows)} points={len(points)} "
        f"phase_assignments={phase_assignment_count} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
