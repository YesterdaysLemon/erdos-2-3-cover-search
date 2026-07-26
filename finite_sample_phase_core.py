#!/usr/bin/env python3
"""Extract the low-coverage critical core of a finite affine-cover sample."""

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
    parser.add_argument("--max-coverage", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path)
    args = parser.parse_args()
    if args.max_coverage < 0:
        raise SystemExit("--max-coverage must be nonnegative")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    points = [
        (int(k), int(l)) for k, l in json.loads(args.points.read_text())
    ]
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(args.phase_file.read_text()).items()
    }
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
    assignment = np.empty(len(rows), dtype=np.uint32)
    for row_index, row in enumerate(rows):
        prime = int(row["p"])
        h = int(row["h"])
        residue = int(row.get("target_residue", 0))
        modulus = int(row.get("target_modulus", 1))
        target = phases.get(prime, residue) % h
        if target % modulus != residue:
            raise RuntimeError(f"invalid phase for p={prime}")
        assignment[row_index] = target

    started = time.monotonic()
    targets, build_seconds = build_targets(points, candidates, np)
    coverage = np.count_nonzero(targets == assignment, axis=1)
    selected_indices = np.flatnonzero(coverage <= args.max_coverage)
    selected = [
        [points[int(index)][0], points[int(index)][1]]
        for index in selected_indices
    ]
    args.output.write_text(json.dumps(selected) + "\n")
    values, counts = np.unique(np.minimum(coverage, 10), return_counts=True)
    report = {
        "pool": str(args.pool),
        "points": str(args.points),
        "phase_file": str(args.phase_file),
        "point_count": len(points),
        "max_coverage": args.max_coverage,
        "selected_count": len(selected),
        "miss_count": int(np.count_nonzero(coverage == 0)),
        "coverage_histogram_cap10": {
            str(int(value)): int(count)
            for value, count in zip(values, counts)
        },
        "build_seconds": build_seconds,
        "elapsed_seconds": time.monotonic() - started,
    }
    if args.report_output:
        args.report_output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"points={len(points)} selected={len(selected)} "
        f"misses={report['miss_count']} build_s={build_seconds:.3f} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
