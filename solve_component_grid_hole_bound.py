#!/usr/bin/env python3
"""Decide whether prescribed affine-line directions can leave few grid holes."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--max-holes", type=int, required=True)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument(
        "--normalize-two",
        action="store_true",
        help=(
            "fix two nonparallel row targets to zero, WLOG under coordinate "
            "translation"
        ),
    )
    parser.add_argument(
        "--order-parallel-groups",
        action="store_true",
        help=(
            "strictly order normalized offsets of interchangeable parallel "
            "rows; duplicate lines can always be moved to an unused offset"
        ),
    )
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    args = parser.parse_args()

    if args.max_holes < 0:
        raise SystemExit("max-holes must be nonnegative")
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType
    from pysat.formula import IDPool
    from pysat.solvers import Solver

    rows = [
        row
        for row in json.loads(args.pool.read_text())["choices"]
        if int(row["h"]) == args.prime
    ]
    if not rows:
        raise RuntimeError("pool has no rows of the requested modulus")

    vpool = IDPool()
    variable = {}
    solver = Solver(name=args.solver)
    clause_count = 0
    for row_index, _row in enumerate(rows):
        literals = []
        for target in range(args.prime):
            literal = vpool.id(("phase", row_index, target))
            variable[row_index, target] = literal
            literals.append(literal)
        encoding = CardEnc.equals(
            literals,
            bound=1,
            vpool=vpool,
            encoding=EncType.seqcounter,
        )
        solver.append_formula(encoding.clauses)
        clause_count += len(encoding.clauses)

    normalized_rows = []
    if args.normalize_two:
        first = 0
        second = next(
            index
            for index in range(1, len(rows))
            if (
                int(rows[first]["a"]) * int(rows[index]["b"])
                - int(rows[first]["b"]) * int(rows[index]["a"])
            )
            % args.prime
            != 0
        )
        normalized_rows = [first, second]
        for index in normalized_rows:
            solver.add_clause([variable[index, 0]])
            clause_count += 1

    ordered_groups = []
    if args.order_parallel_groups:
        direction_groups = {}
        scales = {}
        for index, row in enumerate(rows):
            a = int(row["a"]) % args.prime
            b = int(row["b"]) % args.prime
            if a:
                scale = pow(a, -1, args.prime)
                direction = (1, b * scale % args.prime)
            elif b:
                scale = pow(b, -1, args.prime)
                direction = (0, 1)
            else:
                raise RuntimeError("degenerate affine form")
            scales[index] = scale
            direction_groups.setdefault(direction, []).append(index)
        for direction, indices in sorted(direction_groups.items()):
            if len(indices) < 2:
                continue
            ordered_groups.append(
                {
                    "direction": list(direction),
                    "row_indices": indices,
                }
            )
            for left, right in zip(indices, indices[1:]):
                left_inverse_scale = pow(
                    scales[left],
                    -1,
                    args.prime,
                )
                right_inverse_scale = pow(
                    scales[right],
                    -1,
                    args.prime,
                )
                for left_offset in range(args.prime):
                    left_target = (
                        left_offset * left_inverse_scale
                    ) % args.prime
                    for right_offset in range(
                        left_offset + 1
                    ):
                        right_target = (
                            right_offset * right_inverse_scale
                        ) % args.prime
                        solver.add_clause(
                            [
                                -variable[left, left_target],
                                -variable[right, right_target],
                            ]
                        )
                        clause_count += 1

    hole_variables = []
    cells = []
    for k in range(args.prime):
        for l in range(args.prime):
            hole = vpool.id(("hole", k, l))
            hole_variables.append(hole)
            cells.append((k, l))
            cover_literals = [
                variable[
                    row_index,
                    (
                        int(row["a"]) * k
                        + int(row["b"]) * l
                    )
                    % args.prime,
                ]
                for row_index, row in enumerate(rows)
            ]
            solver.add_clause([*cover_literals, hole])
            clause_count += 1

    bound = CardEnc.atmost(
        hole_variables,
        bound=args.max_holes,
        vpool=vpool,
        encoding=EncType.seqcounter,
    )
    solver.append_formula(bound.clauses)
    clause_count += len(bound.clauses)
    print(
        f"prime={args.prime} rows={len(rows)} cells={len(cells)} "
        f"max_holes={args.max_holes} normalized={normalized_rows} "
        f"variables={vpool.top} clauses={clause_count}",
        flush=True,
    )

    started = time.monotonic()
    sat = solver.solve()
    elapsed = time.monotonic() - started
    phases = {}
    misses = []
    if sat:
        positive = {
            literal for literal in solver.get_model() if literal > 0
        }
        for row_index, row in enumerate(rows):
            phases[int(row["p"])] = next(
                target
                for target in range(args.prime)
                if variable[row_index, target] in positive
            )
        misses = [
            [k, l]
            for k, l in cells
            if all(
                (
                    int(row["a"]) * k
                    + int(row["b"]) * l
                    - phases[int(row["p"])]
                )
                % args.prime
                != 0
                for row in rows
            )
        ]
        if len(misses) > args.max_holes:
            raise AssertionError("decoded model violates the hole bound")
        args.phase_output.write_text(
            json.dumps(
                {
                    str(prime): target
                    for prime, target in phases.items()
                }
            )
            + "\n"
        )
    solver.delete()

    result = {
        "pool": str(args.pool),
        "prime": args.prime,
        "row_count": len(rows),
        "max_holes": args.max_holes,
        "sat": sat,
        "solver": args.solver,
        "normalized_row_indices": normalized_rows,
        "normalized_row_primes": [
            int(rows[index]["p"]) for index in normalized_rows
        ],
        "ordered_parallel_groups": ordered_groups,
        "parallel_ordering_argument": (
            "Rows with proportional affine forms are interchangeable after "
            "normalizing their offsets. A duplicate selected line can be "
            "moved to an unused parallel offset without losing coverage, "
            "so a strictly increasing order is WLOG."
            if ordered_groups
            else None
        ),
        "normalization_argument": (
            "Two nonparallel affine forms give an invertible translation "
            "map, so any phase assignment can be translated to make both "
            "selected targets zero without changing the number of holes."
            if normalized_rows
            else None
        ),
        "variables": vpool.top,
        "clauses": clause_count,
        "solve_seconds": elapsed,
        "decoded_holes": misses if sat else None,
        "phase_file": str(args.phase_output) if sat else None,
    }
    args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"{'SAT' if sat else 'UNSAT'} holes="
        f"{len(misses) if sat else '-'} elapsed_s={elapsed:.3f}",
        flush=True,
    )
    return 0 if sat else 2


if __name__ == "__main__":
    raise SystemExit(main())
