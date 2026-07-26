#!/usr/bin/env python3
"""Exact compact SAT master for a finite set of affine-cover points.

For each row, only target values that occur on the declared points can help.
The encoding creates one variable for each such valid (row, target) pair,
allows at most one target per row, and requires every point to select at
least one incident pair.  A row selecting no variable represents any target
that covers none of the finite points.  Consequently SAT is equivalent to a
phase assignment covering the entire finite point set, while UNSAT is a
rigorous finite obstruction.
"""

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
    parser.add_argument("points", type=Path)
    parser.add_argument(
        "--points-key",
        default="",
        help="read the point list from this top-level key of a JSON object",
    )
    parser.add_argument(
        "--grid-period",
        type=int,
        default=0,
        help="replace the points file contents by the full square residue grid",
    )
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument(
        "--native-atmost",
        action="store_true",
        help=(
            "send at-most-one phase constraints directly to a solver with "
            "native cardinality support (for example MiniCard)"
        ),
    )
    parser.add_argument(
        "--period",
        type=int,
        default=0,
        help=(
            "if nonzero, retain only rows whose modulus divides this period; "
            "coverage of one representative then certifies its whole period cell"
        ),
    )
    parser.add_argument(
        "--min-parent-common",
        type=int,
        default=0,
        help=(
            "retain conditioned rows whose recorded parent_common is at "
            "least this value"
        ),
    )
    parser.add_argument(
        "--fixed-targets",
        default="",
        help="comma-separated prime:target phases fixed by an external symmetry proof",
    )
    parser.add_argument(
        "--initial-phases",
        type=Path,
        help="optional satisfying-subset phases used only as solver hints",
    )
    parser.add_argument("--phase-output", type=Path)
    parser.add_argument("--result-output", type=Path)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    payload = json.loads(args.pool.read_text())
    if args.period < 0:
        raise SystemExit("--period must be nonnegative")
    rows = [
        row
        for row in payload["choices"]
        if (not args.period or args.period % int(row["h"]) == 0)
        and int(row.get("parent_common", 1)) >= args.min_parent_common
    ]
    fixed_targets = {}
    for item in args.fixed_targets.split(","):
        if not item:
            continue
        prime, target = item.split(":", 1)
        fixed_targets[int(prime)] = int(target)
    row_by_prime = {int(row["p"]): row for row in rows}
    if set(fixed_targets) - set(row_by_prime):
        raise RuntimeError("a fixed-target prime is absent from the pool")
    for prime, target in fixed_targets.items():
        row = row_by_prime[prime]
        h = int(row["h"])
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        if not 0 <= target < h or target % modulus != residue:
            raise RuntimeError(f"invalid fixed target {prime}:{target}")
    if args.grid_period < 0:
        raise SystemExit("--grid-period must be nonnegative")
    if args.grid_period:
        points = [
            (k, l)
            for k in range(args.grid_period)
            for l in range(args.grid_period)
        ]
    else:
        points_payload = json.loads(args.points.read_text())
        if args.points_key:
            if not isinstance(points_payload, dict):
                raise RuntimeError("--points-key requires a JSON object")
            points_payload = points_payload[args.points_key]
        points = [(int(k), int(l)) for k, l in points_payload]
    if len(set(points)) != len(points):
        raise RuntimeError("points file contains duplicates")

    pool = IDPool()
    solver = Solver(name=args.solver)
    point_clauses = [[] for _ in points]
    variable_target: dict[int, tuple[int, int]] = {}
    clause_count = 0
    native_atmost_count = 0
    option_count = 0
    started = time.monotonic()

    for row_index, row in enumerate(rows):
        h = int(row["h"])
        p = int(row["p"])
        a = int(row["a"])
        b = int(row["b"])
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        target_variable: dict[int, int] = {}
        for point_index, (k, l) in enumerate(points):
            target = (a * k + b * l) % h
            if target % modulus != residue:
                continue
            if p in fixed_targets and target != fixed_targets[p]:
                continue
            variable = target_variable.get(target)
            if variable is None:
                variable = pool.id(("target", row_index, target))
                target_variable[target] = variable
                variable_target[variable] = (row_index, target)
            point_clauses[point_index].append(variable)
        variables = list(target_variable.values())
        option_count += len(variables)
        if p in fixed_targets and variables:
            if len(variables) != 1:
                raise AssertionError(
                    f"fixed target for p={p} created multiple variables"
                )
            solver.add_clause([variables[0]])
            clause_count += 1
        elif len(variables) > 1:
            if args.native_atmost:
                solver.add_atmost(variables, 1)
                native_atmost_count += 1
            else:
                encoding = CardEnc.atmost(
                    variables,
                    bound=1,
                    vpool=pool,
                    encoding=EncType.seqcounter,
                )
                solver.append_formula(encoding.clauses)
                clause_count += len(encoding.clauses)
            valid_target_count = h // modulus
            if len(variables) == valid_target_count:
                # Every allowable target occurs on the finite point set, so
                # the usual "unselected means an unseen target" relaxation
                # is unnecessary.  A real phase assignment must select one.
                solver.add_clause(variables)
                clause_count += 1

    for point_index, clause in enumerate(point_clauses):
        if not clause:
            raise RuntimeError(f"point {points[point_index]} has no valid row")
        solver.add_clause(clause)
        clause_count += 1
    build_seconds = time.monotonic() - started
    phase_hint_count = 0
    if args.initial_phases:
        initial_phases = {
            int(prime): int(target)
            for prime, target in json.loads(
                args.initial_phases.read_text()
            ).items()
        }
        hints = []
        for variable, (row_index, target) in variable_target.items():
            prime = int(rows[row_index]["p"])
            if initial_phases.get(prime) == target:
                hints.append(variable)
        try:
            solver.set_phases(hints)
            phase_hint_count = len(hints)
        except NotImplementedError:
            phase_hint_count = 0
    print(
        f"rows={len(rows)} points={len(points)} options={option_count} "
        f"variables={pool.top} clauses={clause_count} "
        f"native_atmost={native_atmost_count} "
        f"phase_hints={phase_hint_count} build_s={build_seconds:.3f}",
        flush=True,
    )

    solve_started = time.monotonic()
    sat = solver.solve()
    solve_seconds = time.monotonic() - solve_started
    phases = {}
    if sat:
        model = {literal for literal in solver.get_model() if literal > 0}
        selected_by_row = {}
        for variable in model & variable_target.keys():
            row_index, target = variable_target[variable]
            if row_index in selected_by_row:
                raise AssertionError("model selects two targets for one row")
            selected_by_row[row_index] = target
        for row_index, row in enumerate(rows):
            prime = int(row["p"])
            target = fixed_targets.get(
                prime,
                selected_by_row.get(row_index, int(row["target_residue"])),
            )
            phases[str(int(row["p"]))] = target
        for k, l in points:
            if not any(
                (
                    int(row["a"]) * k
                    + int(row["b"]) * l
                    - phases[str(int(row["p"]))]
                )
                % int(row["h"])
                == 0
                for row in rows
            ):
                raise AssertionError("SAT model does not cover a declared point")
        if args.phase_output:
            args.phase_output.write_text(json.dumps(phases) + "\n")

    result = {
        "pool": str(args.pool),
        "points": str(args.points),
        "points_key": args.points_key or None,
        "grid_period": args.grid_period or None,
        "row_count": len(rows),
        "period": args.period,
        "point_count": len(points),
        "option_count": option_count,
        "variable_count": pool.top,
        "clause_count": clause_count,
        "native_atmost_count": native_atmost_count,
        "solver": args.solver,
        "fixed_targets": {
            str(prime): target for prime, target in fixed_targets.items()
        },
        "initial_phases": (
            str(args.initial_phases) if args.initial_phases else None
        ),
        "phase_hint_count": phase_hint_count,
        "sat": sat,
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
