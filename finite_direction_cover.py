#!/usr/bin/env python3
"""Exact symmetry-reduced cover solver for a prime affine-plane component.

Rows with the same projective normal direction are interchangeable when their
target modulus is one.  The ordinary row-by-row encoding has a factorial
symmetry among those rows.  This solver instead chooses a subset of parallel
lines for each direction, bounded by the number of available prime fibres.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument(
        "--fill-capacities",
        action="store_true",
        help=(
            "select exactly the available number of distinct lines in every "
            "direction; this is without loss for the monotone cover problem"
        ),
    )
    parser.add_argument(
        "--anchor-first-two",
        action="store_true",
        help=(
            "with filled capacities, use translation symmetry to require "
            "target zero in the first two projective directions"
        ),
    )
    parser.add_argument("--phase-output", type=Path)
    parser.add_argument("--result-output", type=Path)
    args = parser.parse_args()
    q = args.prime
    if q < 2:
        raise SystemExit("--prime must be at least two")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = [row for row in payload["choices"] if int(row["h"]) == q]
    groups: dict[tuple[int, int], list[tuple[dict, int]]] = defaultdict(list)
    for row in rows:
        if int(row["target_modulus"]) != 1:
            raise RuntimeError(
                "symmetry reduction currently requires target_modulus=1"
            )
        a = int(row["a"]) % q
        b = int(row["b"]) % q
        if a:
            scale = pow(a, -1, q)
        elif b:
            scale = pow(b, -1, q)
        else:
            raise RuntimeError("zero affine normal")
        direction = (a * scale % q, b * scale % q)
        groups[direction].append((row, scale))
    directions = sorted(groups)
    if args.anchor_first_two and (
        not args.fill_capacities or len(directions) < 2
    ):
        raise SystemExit(
            "--anchor-first-two requires filled capacities and two directions"
        )

    vpool = IDPool()
    solver = Solver(name=args.solver)
    line_var = {}
    clause_count = 0
    started = time.monotonic()
    for direction in directions:
        members = groups[direction]
        variables = []
        for target in range(q):
            variable = vpool.id(("line", direction, target))
            line_var[(direction, target)] = variable
            variables.append(variable)
        if len(members) < q:
            if args.fill_capacities:
                encoding = CardEnc.equals(
                    variables,
                    bound=len(members),
                    vpool=vpool,
                    encoding=EncType.seqcounter,
                )
            else:
                encoding = CardEnc.atmost(
                    variables,
                    bound=len(members),
                    vpool=vpool,
                    encoding=EncType.seqcounter,
                )
            solver.append_formula(encoding.clauses)
            clause_count += len(encoding.clauses)
    if args.anchor_first_two:
        for direction in directions[:2]:
            solver.add_clause([line_var[(direction, 0)]])
            clause_count += 1

    for k in range(q):
        for l in range(q):
            clause = [
                line_var[
                    (
                        direction,
                        (direction[0] * k + direction[1] * l) % q,
                    )
                ]
                for direction in directions
            ]
            solver.add_clause(clause)
            clause_count += 1
    build_seconds = time.monotonic() - started
    print(
        f"prime={q} rows={len(rows)} directions={len(groups)} "
        f"variables={vpool.top} clauses={clause_count} "
        f"build_s={build_seconds:.3f}",
        flush=True,
    )

    solve_started = time.monotonic()
    sat = solver.solve()
    solve_seconds = time.monotonic() - solve_started
    selected = defaultdict(list)
    phases = {}
    if sat:
        model = {literal for literal in solver.get_model() if literal > 0}
        for (direction, target), variable in line_var.items():
            if variable in model:
                selected[direction].append(target)
        for direction in directions:
            members = groups[direction]
            targets = selected[direction]
            if len(targets) > len(members):
                raise AssertionError("direction capacity exceeded")
            if not targets:
                targets = [0]
            for index, (row, scale) in enumerate(members):
                canonical_target = targets[min(index, len(targets) - 1)]
                target = canonical_target * pow(scale, -1, q) % q
                phases[str(int(row["p"]))] = target
        for k in range(q):
            for l in range(q):
                if not any(
                    (
                        int(row["a"]) * k
                        + int(row["b"]) * l
                        - phases[str(int(row["p"]))]
                    )
                    % q
                    == 0
                    for row in rows
                ):
                    raise AssertionError("decoded model misses a grid point")
        if args.phase_output:
            args.phase_output.write_text(json.dumps(phases) + "\n")

    result = {
        "pool": str(args.pool),
        "prime": q,
        "row_count": len(rows),
        "direction_count": len(groups),
        "direction_capacities": {
            f"{a},{b}": len(members)
            for (a, b), members in groups.items()
        },
        "variable_count": vpool.top,
        "clause_count": clause_count,
        "solver": args.solver,
        "fill_capacities": args.fill_capacities,
        "anchor_first_two": args.anchor_first_two,
        "sat": sat,
        "selected_line_counts": {
            f"{a},{b}": len(selected[(a, b)])
            for a, b in groups
        }
        if sat
        else {},
        "build_seconds": build_seconds,
        "solve_seconds": solve_seconds,
    }
    if args.result_output:
        args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"result={'SAT' if sat else 'UNSAT'} solve_s={solve_seconds:.3f}",
        flush=True,
    )
    solver.delete()
    return 0 if sat else 2


if __name__ == "__main__":
    raise SystemExit(main())
