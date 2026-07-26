#!/usr/bin/env python3
"""Exact symmetry-reduced cover over several independent prime components."""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from collections import defaultdict
from math import prod
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--primes", required=True)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--phase-output", type=Path)
    parser.add_argument("--result-output", type=Path)
    args = parser.parse_args()
    primes = tuple(int(value) for value in args.primes.split(",") if value)
    if not primes or len(set(primes)) != len(primes):
        raise SystemExit("--primes must be a nonempty unique list")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows_by_prime = {
        q: [row for row in payload["choices"] if int(row["h"]) == q]
        for q in primes
    }
    groups = defaultdict(list)
    for q, rows in rows_by_prime.items():
        for row in rows:
            if int(row["target_modulus"]) != 1:
                raise RuntimeError(
                    "symmetry reduction requires target_modulus=1"
                )
            a = int(row["a"]) % q
            b = int(row["b"]) % q
            scale = pow(a, -1, q) if a else pow(b, -1, q)
            direction = (a * scale % q, b * scale % q)
            groups[(q, direction)].append((row, scale))

    vpool = IDPool()
    solver = Solver(name=args.solver)
    line_var = {}
    clause_count = 0
    started = time.monotonic()
    for (q, direction), members in groups.items():
        variables = []
        for target in range(q):
            variable = vpool.id(("line", q, direction, target))
            line_var[(q, direction, target)] = variable
            variables.append(variable)
        if len(members) < q:
            encoding = CardEnc.atmost(
                variables,
                bound=len(members),
                vpool=vpool,
                encoding=EncType.seqcounter,
            )
            solver.append_formula(encoding.clauses)
            clause_count += len(encoding.clauses)

    local_points = {
        q: tuple((k, l) for k in range(q) for l in range(q))
        for q in primes
    }
    for point_tuple in itertools.product(
        *(local_points[q] for q in primes)
    ):
        clause = []
        for q, (k, l) in zip(primes, point_tuple):
            for component_q, direction in groups:
                if component_q != q:
                    continue
                target = (direction[0] * k + direction[1] * l) % q
                clause.append(line_var[(q, direction, target)])
        solver.add_clause(clause)
        clause_count += 1
    build_seconds = time.monotonic() - started
    grid_points = prod(q * q for q in primes)
    print(
        f"primes={primes} rows={sum(map(len, rows_by_prime.values()))} "
        f"directions={len(groups)} grid_points={grid_points} "
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
        for (q, direction, target), variable in line_var.items():
            if variable in model:
                selected[(q, direction)].append(target)
        for (q, direction), members in groups.items():
            targets = selected[(q, direction)]
            if len(targets) > len(members):
                raise AssertionError("direction capacity exceeded")
            if not targets:
                targets = [0]
            for index, (row, scale) in enumerate(members):
                canonical_target = targets[min(index, len(targets) - 1)]
                phases[str(int(row["p"]))] = (
                    canonical_target * pow(scale, -1, q) % q
                )
        for point_tuple in itertools.product(
            *(local_points[q] for q in primes)
        ):
            covered = False
            for q, (k, l) in zip(primes, point_tuple):
                for row in rows_by_prime[q]:
                    if (
                        int(row["a"]) * k
                        + int(row["b"]) * l
                        - phases[str(int(row["p"]))]
                    ) % q == 0:
                        covered = True
                        break
                if covered:
                    break
            if not covered:
                raise AssertionError("decoded model misses a grid point")
        if args.phase_output:
            args.phase_output.write_text(json.dumps(phases) + "\n")

    result = {
        "pool": str(args.pool),
        "primes": list(primes),
        "row_counts": {
            str(q): len(rows_by_prime[q]) for q in primes
        },
        "direction_count": len(groups),
        "direction_capacities": {
            f"{q}:{a},{b}": len(members)
            for (q, (a, b)), members in groups.items()
        },
        "grid_point_count": grid_points,
        "variable_count": vpool.top,
        "clause_count": clause_count,
        "solver": args.solver,
        "sat": sat,
        "selected_line_counts": {
            f"{q}:{a},{b}": len(selected[(q, (a, b))])
            for q, (a, b) in groups
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
