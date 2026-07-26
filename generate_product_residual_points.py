#!/usr/bin/env python3
"""Expand the holes of one CRT component block across another component."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def merge_coprime(first: int, first_modulus: int, second: int, second_modulus: int) -> int:
    if math.gcd(first_modulus, second_modulus) != 1:
        raise ValueError("component moduli must be coprime")
    step = (
        (second - first)
        * pow(first_modulus, -1, second_modulus)
    ) % second_modulus
    return (first + first_modulus * step) % (
        first_modulus * second_modulus
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument("--base-period", type=int, required=True)
    parser.add_argument("--extension-period", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    args = parser.parse_args()

    if math.gcd(args.base_period, args.extension_period) != 1:
        raise SystemExit("base and extension periods must be coprime")
    payload = json.loads(args.pool.read_text())
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(
            args.phase_file.read_text()
        ).items()
    }
    rows = [
        row
        for row in payload["choices"]
        if args.base_period % int(row["h"]) == 0
    ]
    if not rows:
        raise RuntimeError("no rows divide the base period")

    base_holes = []
    for k in range(args.base_period):
        for l in range(args.base_period):
            if all(
                (
                    int(row["a"]) * k
                    + int(row["b"]) * l
                    - phases[int(row["p"])]
                )
                % int(row["h"])
                != 0
                for row in rows
            ):
                base_holes.append((k, l))

    k_lifts = {
        (k, residue): merge_coprime(
            k,
            args.base_period,
            residue,
            args.extension_period,
        )
        for k in {point[0] for point in base_holes}
        for residue in range(args.extension_period)
    }
    l_lifts = {
        (l, residue): merge_coprime(
            l,
            args.base_period,
            residue,
            args.extension_period,
        )
        for l in {point[1] for point in base_holes}
        for residue in range(args.extension_period)
    }
    points = [
        [
            k_lifts[k, extension_k],
            l_lifts[l, extension_l],
        ]
        for k, l in base_holes
        for extension_k in range(args.extension_period)
        for extension_l in range(args.extension_period)
    ]
    args.output.write_text(json.dumps(points) + "\n")
    summary = {
        "pool": str(args.pool),
        "phase_file": str(args.phase_file),
        "base_period": args.base_period,
        "extension_period": args.extension_period,
        "base_row_count": len(rows),
        "base_hole_count": len(base_holes),
        "base_holes": [list(point) for point in base_holes],
        "expanded_point_count": len(points),
        "combined_period": (
            args.base_period * args.extension_period
        ),
    }
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"base_rows={len(rows)} base_holes={len(base_holes)} "
        f"expanded={len(points)} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
