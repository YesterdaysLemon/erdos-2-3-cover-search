#!/usr/bin/env python3
"""Synthesize a cover of a finite period plane using divisor-modulus rows.

This is intended for exact repair of an already conditioned branch.  Rows
whose moduli divide P define full affine fibres on (Z/PZ)^2.  Some targets may
be fixed while a named set of rows is retargeted.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def parse_targets(text: str) -> dict[int, int]:
    result = {}
    for item in text.split(","):
        if not item:
            continue
        prime, target = item.split(":", 1)
        result[int(prime)] = int(target)
    return result


def parse_primes(text: str) -> list[int]:
    return [int(item) for item in text.split(",") if item]


def covers(row: dict, target: int, x: int, y: int) -> bool:
    h = int(row["h"])
    return (
        int(row["a"]) * x + int(row["b"]) * y - target
    ) % h == 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument("--fixed-targets", required=True)
    parser.add_argument("--mutable-primes", required=True)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.period <= 0:
        raise RuntimeError("period must be positive")
    fixed = parse_targets(args.fixed_targets)
    mutable_primes = parse_primes(args.mutable_primes)
    if len(mutable_primes) != len(set(mutable_primes)):
        raise RuntimeError("mutable primes are not unique")
    if set(fixed) & set(mutable_primes):
        raise RuntimeError("fixed and mutable sets overlap")

    payload = json.loads(args.pool.read_text())
    rows_by_prime = {
        int(row["p"]): row
        for row in payload["choices"]
        if args.period % int(row["h"]) == 0
    }
    selected_primes = list(fixed) + mutable_primes
    missing = set(selected_primes) - set(rows_by_prime)
    if missing:
        raise RuntimeError(f"selected rows are absent: {sorted(missing)}")
    for prime, target in fixed.items():
        row = rows_by_prime[prime]
        h = int(row["h"])
        if target % int(row["target_modulus"]) != int(
            row["target_residue"]
        ):
            raise RuntimeError(f"fixed target violates restriction p={prime}")
        fixed[prime] = target % h

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType
    from pysat.formula import IDPool
    from pysat.solvers import Solver

    vpool = IDPool()
    variables = {}
    solver = Solver(name=args.solver)
    clause_count = 0
    for index, prime in enumerate(mutable_primes):
        row = rows_by_prime[prime]
        h = int(row["h"])
        literals = []
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        for target in range(residue, h, modulus):
            literal = vpool.id((index, target))
            variables[index, target] = literal
            literals.append(literal)
        encoding = CardEnc.equals(
            literals,
            bound=1,
            vpool=vpool,
            encoding=EncType.seqcounter,
        )
        solver.append_formula(encoding.clauses)
        clause_count += len(encoding.clauses)

    cells_covered_by_fixed = 0
    for x in range(args.period):
        for y in range(args.period):
            if any(
                covers(rows_by_prime[prime], target, x, y)
                for prime, target in fixed.items()
            ):
                cells_covered_by_fixed += 1
                continue
            clause = []
            for index, prime in enumerate(mutable_primes):
                row = rows_by_prime[prime]
                h = int(row["h"])
                target = (
                    int(row["a"]) * x + int(row["b"]) * y
                ) % h
                literal = variables.get((index, target))
                if literal is not None:
                    clause.append(literal)
            solver.add_clause(clause)
            clause_count += 1

    sat = solver.solve()
    if not sat:
        solver.delete()
        raise RuntimeError("selected mutable rows cannot complete the cover")
    positive = {literal for literal in solver.get_model() if literal > 0}
    mutable_targets = {}
    for index, prime in enumerate(mutable_primes):
        row = rows_by_prime[prime]
        h = int(row["h"])
        mutable_targets[prime] = next(
            target
            for target in range(h)
            if variables.get((index, target)) in positive
        )
    solver.delete()

    phases = dict(fixed)
    phases.update(mutable_targets)
    active_rows = [rows_by_prime[prime] for prime in selected_primes]
    misses = [
        [x, y]
        for x in range(args.period)
        for y in range(args.period)
        if not any(
            covers(row, phases[int(row["p"])], x, y)
            for row in active_rows
        )
    ]
    if misses:
        raise AssertionError("decoded assignment leaves period cells")

    core = list(active_rows)
    removed = []
    changed = True
    while changed:
        changed = False
        for row in list(core):
            prime = int(row["p"])
            others = [item for item in core if int(item["p"]) != prime]
            if all(
                any(
                    covers(
                        other,
                        phases[int(other["p"])],
                        x,
                        y,
                    )
                    for other in others
                )
                for x in range(args.period)
                for y in range(args.period)
            ):
                core.remove(row)
                removed.append(prime)
                changed = True
                break

    cover_rows = []
    for row in core:
        item = dict(row)
        item["c"] = phases[int(row["p"])]
        cover_rows.append(item)
    result = {
        "pool": str(args.pool),
        "period": args.period,
        "fixed_targets": {
            str(prime): target for prime, target in sorted(fixed.items())
        },
        "mutable_primes": mutable_primes,
        "mutable_targets": {
            str(prime): target
            for prime, target in sorted(mutable_targets.items())
        },
        "all_selected_primes": selected_primes,
        "cells": args.period * args.period,
        "cells_covered_by_fixed": cells_covered_by_fixed,
        "sat_solver": args.solver,
        "variables": vpool.top,
        "clauses": clause_count,
        "verified_misses": 0,
        "removed_redundant_primes": removed,
        "choices": cover_rows,
        "core_row_count": len(cover_rows),
        "proved_cover": True,
        "scope": (
            "the complete affine period plane in the coordinates of the "
            "conditioned source pool"
        ),
    }
    args.phase_output.write_text(
        json.dumps(
            {str(prime): target for prime, target in sorted(phases.items())}
        )
        + "\n"
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"period={args.period} fixed={len(fixed)} "
        f"mutable={len(mutable_primes)} core={len(core)} "
        f"cells={args.period * args.period} PASS",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
