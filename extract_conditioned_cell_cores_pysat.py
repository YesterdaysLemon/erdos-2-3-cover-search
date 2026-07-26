#!/usr/bin/env python3
"""Extract primary-SAT row cores for closed (mod 16, mod 27) cells."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import exact_uncovered


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_pool", type=Path)
    parser.add_argument("phase_file", type=Path)
    parser.add_argument("cells", type=Path)
    parser.add_argument("--max-component", type=int, default=16384)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.candidate_pool.read_text())
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
    raw_cells = json.loads(args.cells.read_text())
    core_cells = tuple(
        (
            (16, int(k16), int(l16)),
            (27, int(k27), int(l27)),
        )
        for k16, l16, k27, l27 in raw_cells
    )
    _witnesses, meta = exact_uncovered.find_uncovered(
        rows,
        max_component=args.max_component,
        limit=1,
        algebraic_primes=tuple(
            int(prime) for prime in source.get("algebraic_primes", ())
        ),
        sophie_germain=bool(source.get("sophie_germain", False)),
        solver_name=args.solver,
        core_coordinate_cells=core_cells,
    )
    records = meta["core_coordinate_cells"]
    unexpectedly_open = [
        record for record in records if not record["closed"]
    ]
    if unexpectedly_open:
        raise RuntimeError(
            f"{len(unexpectedly_open)} supplied cells are not closed"
        )
    union_indices = sorted(
        {
            index
            for record in records
            for index in record["row_indices"]
        }
    )
    payload = {
        "candidate_pool": str(args.candidate_pool),
        "phase_file": str(args.phase_file),
        "cells_file": str(args.cells),
        "checker": meta,
        "union_row_indices": union_indices,
        "union_primes": [rows[index]["p"] for index in union_indices],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"PASS cells={len(records)} union_rows={len(union_indices)} "
        f"core_sizes={[len(record['row_indices']) for record in records]} "
        f"output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
