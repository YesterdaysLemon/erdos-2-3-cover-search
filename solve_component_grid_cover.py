#!/usr/bin/env python3
"""Solve an exact affine-line cover of F_q^2 for rows with modulus q."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path


def parse_targets(text: str) -> dict[int, int]:
    result = {}
    for item in text.split(","):
        if not item:
            continue
        row_prime, target = item.split(":", 1)
        result[int(row_prime)] = int(target)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--fixed-targets", default="")
    parser.add_argument("--initial-phases", type=Path)
    parser.add_argument(
        "--required-cells-file",
        type=Path,
        help=(
            "optional JSON list of [k,l] cells (or object with open_cells); "
            "only these cells must be covered"
        ),
    )
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType
    from pysat.formula import IDPool
    from pysat.solvers import Solver

    payload = json.loads(args.pool.read_text())
    rows = [
        row
        for row in payload["choices"]
        if int(row["h"]) == args.prime
    ]
    if not rows:
        raise RuntimeError("pool has no rows of the requested modulus")
    for row in rows:
        if (
            int(row["a"]) % args.prime == 0
            and int(row["b"]) % args.prime == 0
        ):
            raise RuntimeError("degenerate component row")
    fixed = parse_targets(args.fixed_targets)
    row_primes = {int(row["p"]) for row in rows}
    if fixed.keys() - row_primes:
        raise RuntimeError("a fixed target is not a component row")
    initial = {}
    if args.initial_phases:
        initial = {
            int(row_prime): int(target)
            for row_prime, target in json.loads(
                args.initial_phases.read_text()
            ).items()
        }
    if args.required_cells_file:
        required_payload = json.loads(args.required_cells_file.read_text())
        raw_cells = (
            required_payload.get("open_cells", ())
            if isinstance(required_payload, dict)
            else required_payload
        )
        required_cells = sorted(
            {
                (
                    int(cell[0]) % args.prime,
                    int(cell[1]) % args.prime,
                )
                for cell in raw_cells
            }
        )
    else:
        required_cells = [
            (k, l)
            for k in range(args.prime)
            for l in range(args.prime)
        ]

    mutable = [
        row for row in rows if int(row["p"]) not in fixed
    ]
    fixed_rows = [
        row for row in rows if int(row["p"]) in fixed
    ]
    vpool = IDPool()
    variable = {}
    solver = Solver(name=args.solver)
    clause_count = 0
    for row_index, row in enumerate(mutable):
        literals = []
        for target in range(args.prime):
            literal = vpool.id((row_index, target))
            variable[(row_index, target)] = literal
            literals.append(literal)
        encoding = CardEnc.equals(
            literals,
            bound=1,
            vpool=vpool,
            encoding=EncType.seqcounter,
        )
        solver.append_formula(encoding.clauses)
        clause_count += len(encoding.clauses)

    uncovered_by_fixed = []
    for k, l in required_cells:
        if any(
            (
                int(row["a"]) * k
                + int(row["b"]) * l
                - fixed[int(row["p"])]
            )
            % args.prime
            == 0
            for row in fixed_rows
        ):
            continue
        clause = [
            variable[
                (
                    row_index,
                    (
                        int(row["a"]) * k
                        + int(row["b"]) * l
                    )
                    % args.prime,
                )
            ]
            for row_index, row in enumerate(mutable)
        ]
        solver.add_clause(clause)
        clause_count += 1
        uncovered_by_fixed.append((k, l))

    hints = []
    for row_index, row in enumerate(mutable):
        row_prime = int(row["p"])
        if row_prime in initial:
            hints.append(
                variable[
                    (row_index, initial[row_prime] % args.prime)
                ]
            )
    try:
        solver.set_phases(hints)
    except NotImplementedError:
        hints = []
    print(
        f"prime={args.prime} rows={len(rows)} fixed={len(fixed_rows)} "
        f"mutable={len(mutable)} cells={len(uncovered_by_fixed)} "
        f"variables={vpool.top} clauses={clause_count}",
        flush=True,
    )
    started = time.monotonic()
    sat = solver.solve()
    elapsed = time.monotonic() - started
    phases = dict(fixed)
    if sat:
        positive = {
            literal for literal in solver.get_model() if literal > 0
        }
        for row_index, row in enumerate(mutable):
            phases[int(row["p"])] = next(
                target
                for target in range(args.prime)
                if variable[(row_index, target)] in positive
            )
    solver.delete()

    misses = []
    if sat:
        for k, l in required_cells:
            if not any(
                (
                    int(row["a"]) * k
                    + int(row["b"]) * l
                    - phases[int(row["p"])]
                )
                % args.prime
                == 0
                for row in rows
            ):
                misses.append((k, l))
        if misses:
            raise AssertionError("decoded SAT assignment leaves a cell")
        args.phase_output.write_text(
            json.dumps(
                {
                    str(row_prime): target
                    for row_prime, target in sorted(phases.items())
                }
            )
            + "\n"
        )

    result = {
        "pool": str(args.pool),
        "prime": args.prime,
        "sat": sat,
        "fixed_rows": len(fixed_rows),
        "mutable_rows": len(mutable),
        "cells_not_covered_by_fixed": len(uncovered_by_fixed),
        "required_cells": len(required_cells),
        "required_cells_file": (
            str(args.required_cells_file)
            if args.required_cells_file
            else None
        ),
        "variables": vpool.top,
        "clauses": clause_count,
        "solver": args.solver,
        "solve_seconds": elapsed,
        "phase_file": str(args.phase_output) if sat else None,
        "verified_misses": len(misses) if sat else None,
    }
    args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"{'SAT' if sat else 'UNSAT'} elapsed_s={elapsed:.3f}",
        flush=True,
    )
    return 0 if sat else 2


if __name__ == "__main__":
    raise SystemExit(main())
