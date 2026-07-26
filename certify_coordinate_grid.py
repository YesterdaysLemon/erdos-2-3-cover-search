#!/usr/bin/env python3
"""Classify and independently certify every cell of one coordinate grid.

The primary checker extracts a small row core for each UNSAT cell.  Every
closed core is then replayed with the independent Z3 bit-vector checker.
Open cells are reported as search targets and are never treated as evidence.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import exact_uncovered
import exact_uncovered_z3_bv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("phase_file", type=Path)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument("--max-component", type=int, default=16384)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--primary-output", type=Path, required=True)
    parser.add_argument("--z3-output", type=Path, required=True)
    args = parser.parse_args()
    if args.period < 1:
        raise SystemExit("--period must be positive")

    source = json.loads(args.pool.read_text())
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(args.phase_file.read_text()).items()
    }
    rows = []
    for raw in source["choices"]:
        row = {
            key: int(raw[key])
            for key in ("h", "p", "a", "b", "ord2", "ord3")
        }
        row["c"] = phases[row["p"]] % row["h"]
        rows.append(row)
    algebraic_primes = tuple(
        int(prime) for prime in source.get("algebraic_primes", ())
    )
    sophie_germain = bool(source.get("sophie_germain", False))
    cells = tuple(
        ((args.period, k, l),)
        for k in range(args.period)
        for l in range(args.period)
    )

    primary_started = time.monotonic()
    _witnesses, meta = exact_uncovered.find_uncovered(
        rows,
        max_component=args.max_component,
        limit=1,
        algebraic_primes=algebraic_primes,
        sophie_germain=sophie_germain,
        solver_name=args.solver,
        core_coordinate_cells=cells,
    )
    records = meta["core_coordinate_cells"]
    open_cells = [
        [int(record["cell"][0][1]), int(record["cell"][0][2])]
        for record in records
        if not record["closed"]
    ]
    closed_records = [record for record in records if record["closed"]]
    union_indices = sorted(
        {
            int(index)
            for record in closed_records
            for index in record["row_indices"]
        }
    )
    primary_payload = {
        "pool": str(args.pool),
        "phase_file": str(args.phase_file),
        "period": args.period,
        "open_cells": open_cells,
        "closed_cell_count": len(closed_records),
        "checker": meta,
        "union_row_indices": union_indices,
        "union_primes": [int(rows[index]["p"]) for index in union_indices],
        "elapsed_seconds": time.monotonic() - primary_started,
    }
    args.primary_output.write_text(
        json.dumps(primary_payload, indent=2) + "\n"
    )
    print(
        f"PRIMARY period={args.period} open={len(open_cells)} "
        f"closed={len(closed_records)} union_rows={len(union_indices)}",
        flush=True,
    )

    z3_started = time.monotonic()
    checks = []
    for number, record in enumerate(closed_records, 1):
        core_rows = [rows[int(index)] for index in record["row_indices"]]
        modulus, k_residue, l_residue = (
            int(value) for value in record["cell"][0]
        )
        misses, z3_meta = exact_uncovered_z3_bv.find_uncovered(
            core_rows,
            max_component=args.max_component,
            limit=1,
            algebraic_primes=algebraic_primes,
            sophie_germain=sophie_germain,
            fixed_coordinate_residues=(
                (modulus, k_residue, l_residue),
            ),
        )
        if misses or z3_meta["sat"]:
            raise RuntimeError(
                "independent replay failed for "
                f"{(k_residue, l_residue)}"
            )
        checks.append(
            {
                "cell": [k_residue, l_residue],
                "row_indices": [
                    int(index) for index in record["row_indices"]
                ],
            }
        )
        if number % 100 == 0 or number == len(closed_records):
            print(
                f"Z3 {number}/{len(closed_records)}",
                flush=True,
            )
    z3_payload = {
        "source": str(args.primary_output),
        "checks": checks,
        "elapsed_seconds": time.monotonic() - z3_started,
    }
    args.z3_output.write_text(json.dumps(z3_payload, indent=2) + "\n")
    print(
        f"PASS period={args.period} closed={len(checks)} "
        f"open={len(open_cells)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
