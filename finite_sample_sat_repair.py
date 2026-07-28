#!/usr/bin/env python3
"""Exact SAT search for a bounded-change finite affine-cover repair.

The current phase is represented explicitly for every movable row, while each
legal noncurrent target observed on the supplied point corpus is an alternate
literal.  Exactly one phase is selected per row and a global cardinality
constraint bounds the number of alternates.  Point clauses are guarded by
assumptions so an UNSAT run can report a finite lesson core.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import exact_greedy
from certify_anchor_phase_quotient import load_cells
from local_phase_cegis import build_targets


def solve_repair(
    rows: list[dict],
    candidates: list[tuple],
    points: list[tuple[int, int]],
    initial: dict[int, int],
    fixed_primes: set[int],
    max_changes: int,
    solver_name: str,
    time_limit: float,
    minimize_changes: bool = False,
):
    if time_limit > 0:
        raise ValueError(
            "in-process SAT timeouts are unsupported; run this solver "
            "through sparse_anchor_quotient_sweep.py for a hard subprocess "
            "deadline"
        )
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool, WCNF  # type: ignore
    from pysat.examples.rc2 import RC2  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    assignment = []
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        residue = int(row.get("target_residue", 0))
        modulus = int(row.get("target_modulus", 1))
        target = int(initial.get(prime, residue)) % h
        if target % modulus != residue % modulus:
            raise RuntimeError(f"invalid initial phase for p={prime}")
        assignment.append(target)
    targets, build_seconds = build_targets(points, candidates, np)
    full_cover = np.count_nonzero(
        targets == np.asarray(assignment, dtype=targets.dtype),
        axis=1,
    )

    vpool = IDPool()
    hard_clauses = []
    option_variables: list[dict[int, int] | None] = []
    alternate_variables = []
    preferred_variables = []
    hard_clause_count = 0
    started = time.monotonic()
    for row_index, row in enumerate(rows):
        prime = int(row["p"])
        current = int(assignment[row_index])
        if prime in fixed_primes:
            option_variables.append(None)
            continue
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        observed = {
            int(target)
            for target in targets[:, row_index]
            if int(target) % modulus == residue
        }
        alternatives = sorted(observed - {current})
        if not alternatives:
            option_variables.append(None)
            continue
        mapping = {
            target: vpool.id(("phase", row_index, target))
            for target in (current, *alternatives)
        }
        option_variables.append(mapping)
        variables = list(mapping.values())
        encoding = CardEnc.equals(
            variables,
            bound=1,
            vpool=vpool,
            encoding=EncType.seqcounter,
        )
        hard_clauses.extend(encoding.clauses)
        hard_clause_count += len(encoding.clauses)
        preferred_variables.append(mapping[current])
        alternate_variables.extend(
            mapping[target] for target in alternatives
        )
    if not minimize_changes and len(alternate_variables) > max_changes:
        encoding = CardEnc.atmost(
            alternate_variables,
            bound=max_changes,
            vpool=vpool,
            encoding=EncType.seqcounter,
        )
        hard_clauses.extend(encoding.clauses)
        hard_clause_count += len(encoding.clauses)

    point_activations = []
    activation_to_index = {}
    constant_covered_points = 0
    for point_index, point_targets in enumerate(targets):
        clause = []
        constant_covered = False
        for row_index, target in enumerate(point_targets):
            target = int(target)
            mapping = option_variables[row_index]
            if mapping is None:
                if target == int(assignment[row_index]):
                    constant_covered = True
                    break
                continue
            variable = mapping.get(target)
            if variable is not None:
                clause.append(variable)
        if constant_covered:
            constant_covered_points += 1
            continue
        if minimize_changes:
            hard_clauses.append(clause)
        else:
            activation = vpool.id(("point", point_index))
            point_activations.append(activation)
            activation_to_index[activation] = point_index
            hard_clauses.append([-activation, *clause])
        hard_clause_count += 1

    encode_seconds = time.monotonic() - started
    solve_started = time.monotonic()
    solver = None
    optimizer = None
    minimum_changes = None
    if minimize_changes:
        formula = WCNF()
        formula.extend(hard_clauses)
        for variable in preferred_variables:
            formula.append([variable], weight=1)
        optimizer = RC2(formula, solver=solver_name)
        model_literals = optimizer.compute()
        sat = model_literals is not None
        minimum_changes = optimizer.cost if sat else None
    else:
        solver = Solver(name=solver_name, bootstrap_with=hard_clauses)
        sat = solver.solve(assumptions=point_activations)
        model_literals = solver.get_model() if sat else None
    solve_seconds = time.monotonic() - solve_started

    repaired = None
    core_indices = []
    changed = None
    full_misses = None
    if sat is True:
        model = {literal for literal in model_literals if literal > 0}
        repaired = list(assignment)
        for row_index, mapping in enumerate(option_variables):
            if mapping is None:
                continue
            chosen = [
                target
                for target, variable in mapping.items()
                if variable in model
            ]
            if len(chosen) != 1:
                raise AssertionError("SAT model does not select one phase")
            repaired[row_index] = chosen[0]
        repaired_array = np.asarray(repaired, dtype=targets.dtype)
        repaired_cover = np.count_nonzero(
            targets == repaired_array,
            axis=1,
        )
        miss_indices = np.flatnonzero(repaired_cover == 0)
        full_misses = int(len(miss_indices))
        if full_misses:
            raise AssertionError("SAT repair misses a supplied point")
        changed = sum(
            old != new for old, new in zip(assignment, repaired)
        )
        if not minimize_changes and changed > max_changes:
            raise AssertionError("SAT repair exceeds the change bound")
    elif sat is False:
        if solver is not None:
            core = set(solver.get_core() or ())
            core_indices = sorted(
                activation_to_index[literal]
                for literal in core
                if literal in activation_to_index
            )
    if solver is not None:
        solver.delete()
    if optimizer is not None:
        optimizer.delete()
    return {
        "sat": sat,
        "assignment": assignment,
        "repaired": repaired,
        "changed_phases": changed,
        "minimum_changes": minimum_changes,
        "within_change_bound": (
            changed is not None and changed <= max_changes
        ),
        "full_misses": full_misses,
        "core_indices": core_indices,
        "initial_misses": int(np.count_nonzero(full_cover == 0)),
        "constant_covered_points": constant_covered_points,
        "variable_count": vpool.top,
        "alternate_variable_count": len(alternate_variables),
        "clause_count": hard_clause_count,
        "build_seconds": build_seconds,
        "encode_seconds": encode_seconds,
        "solve_seconds": solve_seconds,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path, nargs="+")
    parser.add_argument("--initial-phases", type=Path, required=True)
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    parser.add_argument("--fixed-primes", default="")
    parser.add_argument("--max-changes", type=int, default=4)
    parser.add_argument("--solver", default="glucose42")
    parser.add_argument(
        "--minimize-changes",
        action="store_true",
        help="compute the exact minimum Hamming distance with RC2 MaxSAT",
    )
    parser.add_argument(
        "--time-limit",
        type=float,
        default=0.0,
        help=(
            "must be zero for this in-process solver; the quotient sweep "
            "enforces hard subprocess deadlines"
        ),
    )
    args = parser.parse_args()
    if args.max_changes < 0:
        raise SystemExit("--max-changes must be nonnegative")

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    candidates = exact_greedy.load_candidates(args.pool, True)
    points = []
    seen_points = set()
    for path in args.points:
        for point in load_cells(path):
            if point not in seen_points:
                seen_points.add(point)
                points.append(point)
    initial = {
        int(prime): int(target)
        for prime, target in json.loads(
            args.initial_phases.read_text()
        ).items()
    }
    fixed_primes = {
        int(value) for value in args.fixed_primes.split(",") if value
    }
    row_primes = {int(row["p"]) for row in rows}
    missing_fixed = fixed_primes - row_primes
    if missing_fixed:
        raise RuntimeError(f"fixed primes absent from pool: {missing_fixed}")

    result = solve_repair(
        rows,
        candidates,
        points,
        initial,
        fixed_primes,
        args.max_changes,
        args.solver,
        args.time_limit,
        args.minimize_changes,
    )
    if result["sat"] is True:
        status = (
            "INTEGER_MODEL"
            if result["within_change_bound"]
            else "OVER_BUDGET"
        )
        phase_map = {
            str(int(row["p"])): int(target)
            for row, target in zip(rows, result["repaired"])
        }
        args.phase_output.write_text(json.dumps(phase_map) + "\n")
        return_code = 0 if status == "INTEGER_MODEL" else 1
    elif result["sat"] is False:
        status = "UNSAT"
        return_code = 2
    else:
        status = "TIME_LIMIT"
        return_code = 1
    output = {
        "pool": str(args.pool),
        "points": [str(path) for path in args.points],
        "initial_phases": str(args.initial_phases),
        "result": status,
        "solver": args.solver,
        "fixed_primes": sorted(fixed_primes),
        "max_changes": args.max_changes,
        "minimize_changes": args.minimize_changes,
        "point_count": len(points),
        "core_point_count": len(result["core_indices"]),
        "core_points": [
            list(points[index]) for index in result["core_indices"]
        ],
        **{
            key: value
            for key, value in result.items()
            if key not in {"sat", "assignment", "repaired", "core_indices"}
        },
    }
    args.result_output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"rows={len(rows)} points={len(points)} "
        f"initial_misses={result['initial_misses']} "
        f"variables={result['variable_count']} "
        f"clauses={result['clause_count']} "
        f"result={status} changed={result['changed_phases']} "
        f"core={len(result['core_indices'])} "
        f"solve_s={result['solve_seconds']:.3f}",
        flush=True,
    )
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
