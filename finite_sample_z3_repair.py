#!/usr/bin/env python3
"""Exact bounded-change repair using sparse Z3 pseudo-Boolean constraints.

Each Boolean variable means that one row moves from its base phase to one
legal target observed on the finite corpus.  Keeping the base phase is
implicit.  Per-row and global at-most constraints enforce a genuine Hamming
ball, while signed point-incidence changes are encoded as pseudo-Boolean
coverage inequalities.
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
    time_limit: float,
):
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    import z3  # type: ignore

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
    assignment_array = np.asarray(assignment, dtype=targets.dtype)
    base_cover = np.count_nonzero(
        targets == assignment_array,
        axis=1,
    )

    solver = z3.Solver()
    if time_limit > 0:
        solver.set(timeout=max(1, round(1000 * time_limit)))
    move_variables: list[dict[int, object] | None] = []
    all_moves = []
    build_started = time.monotonic()
    for row_index, row in enumerate(rows):
        prime = int(row["p"])
        current = int(assignment[row_index])
        if prime in fixed_primes:
            move_variables.append(None)
            continue
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        alternatives = sorted(
            {
                int(target)
                for target in targets[:, row_index]
                if (
                    int(target) != current
                    and int(target) % modulus == residue
                )
            }
        )
        if not alternatives:
            move_variables.append(None)
            continue
        mapping = {
            target: z3.Bool(f"move_{row_index}_{option_index}")
            for option_index, target in enumerate(alternatives)
        }
        variables = list(mapping.values())
        move_variables.append(mapping)
        all_moves.extend(variables)
        if len(variables) > 1:
            solver.add(z3.AtMost(*variables, 1))
    if len(all_moves) > max_changes:
        solver.add(z3.AtMost(*all_moves, max_changes))

    point_activations = []
    activation_to_index = {}
    for point_index, point_targets in enumerate(targets):
        terms = []
        loss_count = 0
        for row_index, observed_target in enumerate(point_targets):
            mapping = move_variables[row_index]
            if mapping is None:
                continue
            base_covers = (
                int(observed_target) == int(assignment[row_index])
            )
            for target, variable in mapping.items():
                alternate_covers = int(observed_target) == target
                if alternate_covers and not base_covers:
                    terms.append((variable, 1))
                elif base_covers and not alternate_covers:
                    terms.append((z3.Not(variable), 1))
                    loss_count += 1
        lower = 1 - int(base_cover[point_index]) + loss_count
        if lower <= 0:
            continue
        activation = z3.Bool(f"point_{point_index}")
        point_activations.append(activation)
        activation_to_index[activation] = point_index
        constraint = z3.PbGe(terms, lower) if terms else z3.BoolVal(False)
        solver.add(z3.Implies(activation, constraint))
    encode_seconds = time.monotonic() - build_started

    solve_started = time.monotonic()
    check = solver.check(*point_activations)
    solve_seconds = time.monotonic() - solve_started
    repaired = None
    changed = None
    full_misses = None
    core_indices = []
    if check == z3.sat:
        model = solver.model()
        repaired = list(assignment)
        for row_index, mapping in enumerate(move_variables):
            if mapping is None:
                continue
            chosen = [
                target
                for target, variable in mapping.items()
                if z3.is_true(model.eval(variable, model_completion=True))
            ]
            if len(chosen) > 1:
                raise AssertionError("Z3 selected two targets for one row")
            if chosen:
                repaired[row_index] = chosen[0]
        repaired_array = np.asarray(repaired, dtype=targets.dtype)
        repaired_cover = np.count_nonzero(
            targets == repaired_array,
            axis=1,
        )
        full_misses = int(np.count_nonzero(repaired_cover == 0))
        if full_misses:
            raise AssertionError("Z3 repair misses a supplied point")
        changed = sum(
            old != new for old, new in zip(assignment, repaired)
        )
        if changed > max_changes:
            raise AssertionError("Z3 repair exceeds the change bound")
    elif check == z3.unsat:
        core = set(solver.unsat_core())
        core_indices = sorted(
            activation_to_index[activation]
            for activation in core
            if activation in activation_to_index
        )
    return {
        "check": str(check),
        "assignment": assignment,
        "repaired": repaired,
        "changed_phases": changed,
        "full_misses": full_misses,
        "core_indices": core_indices,
        "initial_misses": int(np.count_nonzero(base_cover == 0)),
        "move_variable_count": len(all_moves),
        "point_constraint_count": len(point_activations),
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
    parser.add_argument("--time-limit", type=float, default=240.0)
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
        args.time_limit,
    )
    if result["check"] == "sat":
        status = "INTEGER_MODEL"
        phase_map = {
            str(int(row["p"])): int(target)
            for row, target in zip(rows, result["repaired"])
        }
        args.phase_output.write_text(json.dumps(phase_map) + "\n")
        return_code = 0
    elif result["check"] == "unsat":
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
        "engine": "z3-pseudo-boolean-bounded-repair",
        "fixed_primes": sorted(fixed_primes),
        "max_changes": args.max_changes,
        "point_count": len(points),
        "core_point_count": len(result["core_indices"]),
        "core_points": [
            list(points[index]) for index in result["core_indices"]
        ],
        **{
            key: value
            for key, value in result.items()
            if key not in {
                "check",
                "assignment",
                "repaired",
                "core_indices",
            }
        },
    }
    args.result_output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"rows={len(rows)} points={len(points)} "
        f"initial_misses={result['initial_misses']} "
        f"moves={result['move_variable_count']} "
        f"constraints={result['point_constraint_count']} "
        f"result={status} changed={result['changed_phases']} "
        f"core={len(result['core_indices'])} "
        f"solve_s={result['solve_seconds']:.3f}",
        flush=True,
    )
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
