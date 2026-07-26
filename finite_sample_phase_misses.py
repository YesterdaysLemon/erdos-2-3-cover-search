#!/usr/bin/env python3
"""Check a derived phase map against an explicit finite point set."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

from local_phase_cegis import build_targets


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument("--points-key", default="")
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="retain at most this many misses; zero retains every miss",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.limit < 0:
        raise SystemExit("--limit must be nonnegative")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    candidates = [
        (
            int(row["h"]),
            int(row["p"]),
            int(row["a"]),
            int(row["b"]),
            int(row.get("ord2", row["h"])),
            int(row.get("ord3", row["h"])),
        )
        for row in rows
    ]
    points_payload = json.loads(args.points.read_text())
    if args.points_key:
        points_payload = points_payload[args.points_key]
    points = [(int(k), int(l)) for k, l in points_payload]
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(args.phase_file.read_text()).items()
    }
    assignment = []
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        target = phases.get(prime, residue) % h
        if target % modulus != residue:
            raise RuntimeError(f"invalid phase for prime {prime}")
        assignment.append(target)
    assignment_array = np.asarray(assignment, dtype=np.uint32)

    started = time.monotonic()
    targets, build_seconds = build_targets(points, candidates, np)
    coverage = np.count_nonzero(
        targets == assignment_array, axis=1
    ).astype(np.int32)
    miss_indices = np.flatnonzero(coverage == 0)
    retained_indices = (
        miss_indices[: args.limit] if args.limit else miss_indices
    )
    histogram_values, histogram_counts = np.unique(
        np.minimum(coverage, 10), return_counts=True
    )
    result = {
        "pool": str(args.pool),
        "points": str(args.points),
        "phase_file": str(args.phase_file),
        "point_count": len(points),
        "miss_count": int(len(miss_indices)),
        "retained_miss_count": int(len(retained_indices)),
        "coverage_histogram_cap10": {
            str(int(value)): int(count)
            for value, count in zip(histogram_values, histogram_counts)
        },
        "misses": [
            [points[int(index)][0], points[int(index)][1]]
            for index in retained_indices
        ],
        "build_seconds": build_seconds,
        "elapsed_seconds": time.monotonic() - started,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"points={len(points)} misses={len(miss_indices)} "
        f"retained={len(retained_indices)} build_s={build_seconds:.3f} "
        f"output={args.output}",
        flush=True,
    )
    return 1 if len(miss_indices) else 0


if __name__ == "__main__":
    raise SystemExit(main())
