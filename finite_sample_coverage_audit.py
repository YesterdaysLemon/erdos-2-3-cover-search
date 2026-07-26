#!/usr/bin/env python3
"""Independently recount affine-fibre coverage on a finite point file."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument("--required-coverage", type=int, default=1)
    parser.add_argument(
        "--prefix",
        type=int,
        default=0,
        help="audit only this many leading points; zero audits all",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.required_coverage < 1:
        raise SystemExit("--required-coverage must be positive")
    if args.prefix < 0:
        raise SystemExit("--prefix must be nonnegative")

    rows = json.loads(args.pool.read_text())["choices"]
    points = [
        (int(k), int(l)) for k, l in json.loads(args.points.read_text())
    ]
    if args.prefix:
        if args.prefix > len(points):
            raise SystemExit("--prefix exceeds point count")
        points = points[: args.prefix]
    phases = {
        int(prime): int(value)
        for prime, value in json.loads(args.phase_file.read_text()).items()
    }
    missing = [
        int(row["p"]) for row in rows if int(row["p"]) not in phases
    ]
    if missing:
        raise RuntimeError(f"phase file is missing primes {missing[:10]}")
    histogram = Counter()
    below = []
    deficit = 0
    for k, l in points:
        coverage = 0
        for row in rows:
            h = int(row["h"])
            prime = int(row["p"])
            phase = phases[prime] % h
            if (
                int(row["a"]) * (k % h)
                + int(row["b"]) * (l % h)
                - phase
            ) % h == 0:
                coverage += 1
        histogram[coverage] += 1
        if coverage < args.required_coverage:
            deficit += args.required_coverage - coverage
            below.append([k, l, coverage])
    result = {
        "pool": str(args.pool),
        "points": str(args.points),
        "phase_file": str(args.phase_file),
        "point_count": len(points),
        "row_count": len(rows),
        "required_coverage": args.required_coverage,
        "minimum_coverage": min(histogram) if histogram else None,
        "coverage_histogram": {
            str(coverage): histogram[coverage]
            for coverage in sorted(histogram)
        },
        "below_required": len(below),
        "total_deficit": deficit,
        "below_points": below,
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"PASS points={len(points)} rows={len(rows)} "
        f"minimum={result['minimum_coverage']} "
        f"below={len(below)} deficit={deficit}",
        flush=True,
    )
    return 1 if below else 0


if __name__ == "__main__":
    raise SystemExit(main())
