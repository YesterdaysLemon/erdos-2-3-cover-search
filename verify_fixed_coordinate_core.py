#!/usr/bin/env python3
"""Dual-verify a lifted core on one or more fixed coordinate cells."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import exact_uncovered
import exact_uncovered_z3_bv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("core", type=Path)
    parser.add_argument(
        "--cells",
        required=True,
        help="comma-separated modulus:k:l coordinate restrictions",
    )
    parser.add_argument("--max-component", type=int, default=16384)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.core.read_text())
    rows = payload["rows"]
    algebraic_primes = tuple(
        int(prime) for prime in payload.get("algebraic_primes", ())
    )
    sophie_germain = bool(payload.get("sophie_germain", False))
    cells = []
    for raw in args.cells.split(","):
        modulus, k_residue, l_residue = map(int, raw.split(":"))
        cells.append((modulus, k_residue % modulus, l_residue % modulus))

    checks = []
    for cell in cells:
        fixed = (cell,)
        for name, engine in (
            ("pysat", exact_uncovered),
            ("z3bv", exact_uncovered_z3_bv),
        ):
            misses, meta = engine.find_uncovered(
                rows,
                max_component=args.max_component,
                limit=1,
                algebraic_primes=algebraic_primes,
                sophie_germain=sophie_germain,
                fixed_coordinate_residues=fixed,
            )
            if misses or meta["sat"]:
                raise RuntimeError(f"{name} failed fixed cell {cell}")
            checks.append(
                {
                    "cell": list(cell),
                    "checker": name,
                    "sat": bool(meta["sat"]),
                }
            )
            print(f"{name} cell={cell} UNSAT", flush=True)
    result = {
        "core": str(args.core),
        "checks": checks,
        "all_passed": True,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"PASS checks={len(checks)} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
