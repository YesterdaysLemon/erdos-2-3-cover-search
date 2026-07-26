#!/usr/bin/env python3
"""Extract Z3 UNSAT row cores for closed low-coordinate cells.

Every guarded row asserts that an uncovered point lies outside that row's
selected affine fibre.  If a subset of guards is already UNSAT after fixing a
low CRT cell, those rows cover the whole cell (apart from the declared
algebraic identities).  Keeping their phases fixed therefore preserves that
closed cell during later CEGIS rounds.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import exact_uncovered


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_pool", type=Path)
    parser.add_argument("phase_file", type=Path)
    parser.add_argument("cells", type=Path)
    parser.add_argument("--max-component", type=int, default=16384)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import z3  # type: ignore

    source = json.loads(args.candidate_pool.read_text())
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(args.phase_file.read_text()).items()
    }
    rows = []
    factors_by_row = []
    maximal: dict[int, int] = {}
    for raw in source["choices"]:
        row = {
            key: int(raw[key])
            for key in ("h", "p", "a", "b")
        }
        row["c"] = phases[row["p"]] % row["h"]
        factors = exact_uncovered.factor(row["h"])
        rows.append(row)
        factors_by_row.append(factors)
        for prime, exponent in factors.items():
            maximal[prime] = max(maximal.get(prime, 0), exponent)
    algebraic_primes = tuple(
        int(prime) for prime in source.get("algebraic_primes", ())
    )
    for prime in algebraic_primes:
        maximal[prime] = max(maximal.get(prime, 0), 1)
    components = {
        prime: prime**exponent for prime, exponent in maximal.items()
    }
    if max(components.values(), default=1) > args.max_component:
        raise RuntimeError(
            f"component guard exceeded: {max(components.values())}"
        )
    if components.get(2, 1) < 16 or components.get(3, 1) < 27:
        raise RuntimeError("the requested 16-by-27 cell coordinates are absent")

    solver = z3.Solver()
    kval = {prime: z3.Int(f"k_{prime}") for prime in components}
    lval = {prime: z3.Int(f"l_{prime}") for prime in components}
    for prime, modulus in components.items():
        solver.add(kval[prime] >= 0, kval[prime] < modulus)
        solver.add(lval[prime] >= 0, lval[prime] < modulus)

    activations = []
    activation_index = {}
    for index, (row, factors) in enumerate(zip(rows, factors_by_row)):
        activation = z3.Bool(f"use_row_{index}")
        activations.append(activation)
        activation_index[activation.decl().name()] = index
        failures = []
        for prime, exponent in factors.items():
            modulus = prime**exponent
            failures.append(
                (
                    row["a"] * kval[prime]
                    + row["b"] * lval[prime]
                    - row["c"]
                )
                % modulus
                != 0
            )
        solver.add(z3.Implies(activation, z3.Or(*failures)))
    for prime in algebraic_primes:
        solver.add(
            z3.Or(kval[prime] % prime != 0, lval[prime] % prime != 0)
        )
    if source.get("sophie_germain", False):
        solver.add(z3.Or(kval[2] % 4 != 2, lval[2] % 4 != 0))

    raw_cells = json.loads(args.cells.read_text())
    records = []
    union_indices: set[int] = set()
    started = time.monotonic()
    for cell_no, raw_cell in enumerate(raw_cells, 1):
        k16, l16, k27, l27 = map(int, raw_cell)
        solver.push()
        solver.add(kval[2] % 16 == k16, lval[2] % 16 == l16)
        solver.add(kval[3] % 27 == k27, lval[3] % 27 == l27)
        status = solver.check(*activations)
        if status != z3.unsat:
            raise RuntimeError(
                f"cell {raw_cell} is not closed under the supplied phases: "
                f"{status}"
            )
        core_indices = sorted(
            activation_index[value.decl().name()]
            for value in solver.unsat_core()
        )
        solver.pop()
        if not core_indices:
            raise RuntimeError(f"empty row core for cell {raw_cell}")
        union_indices.update(core_indices)
        records.append(
            {
                "cell": [k16, l16, k27, l27],
                "row_indices": core_indices,
                "primes": [rows[index]["p"] for index in core_indices],
            }
        )
        print(
            f"cell={cell_no}/{len(raw_cells)} core={len(core_indices)} "
            f"union={len(union_indices)}",
            flush=True,
        )

    union_indices_sorted = sorted(union_indices)
    payload = {
        "candidate_pool": str(args.candidate_pool),
        "phase_file": str(args.phase_file),
        "algebraic_primes": list(algebraic_primes),
        "components": components,
        "cell_coordinate_moduli": [16, 27],
        "cells": records,
        "union_row_indices": union_indices_sorted,
        "union_primes": [rows[index]["p"] for index in union_indices_sorted],
        "elapsed_seconds": time.monotonic() - started,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"PASS cells={len(records)} union_rows={len(union_indices_sorted)} "
        f"output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
