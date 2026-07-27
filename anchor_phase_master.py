#!/usr/bin/env python3
"""Exact SAT master for a selected phase-coordinate quotient.

Selected anchor rows may choose any legal target.  Every other row remains
at a supplied base phase.  Supplied exact exponent pairs that are already
covered by a frozen row or an algebraic identity are discharged; each
remaining point becomes one short clause over the anchor target variables.

SAT yields a coordinated anchor assignment for an exact checker to attack.
UNSAT is a finite obstruction for this anchor set and frozen remainder only.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

from certify_anchor_phase_quotient import load_cells


def normalize_rows(
    rows: list[dict],
    phases: dict[int, int],
) -> tuple[list[tuple[int, int, int, int, int, int]], dict[int, int]]:
    compact = []
    normalized_phases = {}
    seen_primes = set()
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        modulus = int(row.get("target_modulus", 1))
        if prime in seen_primes:
            raise ValueError(f"row prime {prime} is repeated")
        if h < 1 or modulus < 1 or h % modulus:
            raise ValueError(f"row p={prime} has an invalid restriction")
        a = int(row["a"]) % h
        b = int(row["b"]) % h
        residue = int(row.get("target_residue", 0)) % modulus
        if prime not in phases:
            raise ValueError(f"base phase omits row p={prime}")
        phase = int(phases[prime]) % h
        if phase % modulus != residue:
            raise ValueError(f"base phase for p={prime} is forbidden")
        seen_primes.add(prime)
        normalized_phases[prime] = phase
        compact.append((prime, h, a, b, modulus, residue))
    return compact, normalized_phases


def solve_anchor_master(
    rows: list[dict],
    base_phases: dict[int, int],
    anchor_primes: tuple[int, ...],
    points: list[tuple[int, int]],
    algebraic_primes: tuple[int, ...] = (),
    solver_name: str = "cadical195",
) -> tuple[dict[int, int] | None, dict]:
    if not anchor_primes:
        raise ValueError("at least one anchor row is required")
    if len(set(anchor_primes)) != len(anchor_primes):
        raise ValueError("anchor primes must be distinct")
    compact_rows, normalized_phases = normalize_rows(rows, base_phases)
    row_by_prime = {row[0]: row for row in compact_rows}
    if set(anchor_primes) - row_by_prime.keys():
        raise ValueError("an anchor prime is absent from the row pool")
    if len(set(points)) != len(points):
        raise ValueError("master points must be distinct")
    if any(k < 0 or l < 0 for k, l in points):
        raise ValueError("master points must be nonnegative")

    anchor_set = set(anchor_primes)
    frozen_rows = [
        row for row in compact_rows if row[0] not in anchor_set
    ]
    eligible_points = []
    algebraically_discharged = 0
    frozen_discharged = 0
    filter_started = time.monotonic()
    for k, l in points:
        if any(
            k % prime == 0 and l % prime == 0
            for prime in algebraic_primes
        ):
            algebraically_discharged += 1
            continue
        if any(
            (a * k + b * l - normalized_phases[prime]) % h == 0
            for prime, h, a, b, _modulus, _residue in frozen_rows
        ):
            frozen_discharged += 1
            continue
        eligible_points.append((k, l))
    filter_seconds = time.monotonic() - filter_started

    dependency_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dependency_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    build_started = time.monotonic()
    variable_pool = IDPool()
    solver = Solver(name=solver_name)
    target_variable = {}
    option_count = 0
    for prime in anchor_primes:
        _prime, h, _a, _b, modulus, residue = row_by_prime[prime]
        variables = []
        for target in range(h):
            if target % modulus != residue:
                continue
            variable = variable_pool.id(("target", prime, target))
            target_variable[prime, target] = variable
            variables.append(variable)
        option_count += len(variables)
        encoding = CardEnc.equals(
            variables,
            bound=1,
            vpool=variable_pool,
            encoding=EncType.seqcounter,
        )
        solver.append_formula(encoding.clauses)

    empty_point = None
    for k, l in eligible_points:
        clause = []
        for prime in anchor_primes:
            _prime, h, a, b, modulus, residue = row_by_prime[prime]
            target = (a * k + b * l) % h
            if target % modulus == residue:
                clause.append(target_variable[prime, target])
        if not clause:
            empty_point = (k, l)
            solver.add_clause([])
        else:
            solver.add_clause(clause)
    build_seconds = time.monotonic() - build_started

    solve_started = time.monotonic()
    sat = solver.solve()
    solve_seconds = time.monotonic() - solve_started
    phases = None
    selected_targets = None
    if sat:
        model = {literal for literal in solver.get_model() if literal > 0}
        selected_targets = {}
        for prime in anchor_primes:
            targets = [
                target
                for (row_prime, target), variable in target_variable.items()
                if row_prime == prime and variable in model
            ]
            if len(targets) != 1:
                raise AssertionError("master did not select one anchor target")
            selected_targets[prime] = targets[0]
        phases = dict(normalized_phases)
        phases.update(selected_targets)
        for k, l in points:
            algebraic_cover = any(
                k % prime == 0 and l % prime == 0
                for prime in algebraic_primes
            )
            row_cover = any(
                (a * k + b * l - phases[prime]) % h == 0
                for prime, h, a, b, _modulus, _residue in compact_rows
            )
            if not algebraic_cover and not row_cover:
                raise AssertionError("SAT phase misses a supplied point")

    metadata = {
        "sat": sat,
        "rows": len(compact_rows),
        "anchor_rows": len(anchor_primes),
        "frozen_rows": len(frozen_rows),
        "input_points": len(points),
        "eligible_points": len(eligible_points),
        "algebraically_discharged_points": algebraically_discharged,
        "frozen_discharged_points": frozen_discharged,
        "empty_clause_point": list(empty_point) if empty_point else None,
        "target_options": option_count,
        "variables": variable_pool.top,
        "clauses": solver.nof_clauses(),
        "filter_seconds": filter_seconds,
        "build_seconds": build_seconds,
        "solve_seconds": solve_seconds,
        "selected_targets": (
            {str(prime): target for prime, target in selected_targets.items()}
            if selected_targets
            else None
        ),
        "engine": "pysat-exact-anchor-phase-master",
    }
    solver.delete()
    return phases, metadata


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--base-phases", type=Path, required=True)
    parser.add_argument("--anchors", required=True)
    parser.add_argument("--points-file", type=Path, action="append", required=True)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(args.base_phases.read_text()).items()
    }
    anchors = tuple(
        int(value) for value in args.anchors.split(",") if value
    )
    points = []
    seen = set()
    for path in args.points_file:
        for point in load_cells(path):
            if point not in seen:
                seen.add(point)
                points.append(point)
    solved_phases, metadata = solve_anchor_master(
        payload["choices"],
        phases,
        anchors,
        points,
        tuple(int(q) for q in payload.get("algebraic_primes", ())),
        args.solver,
    )
    if solved_phases is not None:
        args.phase_output.write_text(
            json.dumps(
                {
                    str(prime): target
                    for prime, target in solved_phases.items()
                }
            )
            + "\n"
        )
    result = {
        "pool": str(args.pool),
        "base_phases": str(args.base_phases),
        "anchors": list(anchors),
        "points_files": [str(path) for path in args.points_file],
        "complete": False,
        "scope": "finite points with every nonanchor phase frozen",
        "master": metadata,
        "phase_file": str(args.phase_output) if solved_phases else None,
    }
    args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    label = "ANCHOR_MASTER_SAT" if solved_phases else "ANCHOR_MASTER_UNSAT"
    print(
        f"{label} anchors={len(anchors)} "
        f"points={metadata['eligible_points']} "
        f"options={metadata['target_options']} "
        f"solve_s={metadata['solve_seconds']} "
        f"output={args.result_output}",
        flush=True,
    )
    return 0 if solved_phases else 2


if __name__ == "__main__":
    raise SystemExit(main())
