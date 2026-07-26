#!/usr/bin/env python3
"""Exact symmetry-reduced cover on a squarefree CRT coordinate product.

Every retained row modulus divides the requested squarefree period.  Rows with
the same componentwise projective normal directions are interchangeable, so
the solver chooses bounded subsets of affine CRT flats rather than assigning
labelled targets to indistinguishable prime fibres.
"""

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

import exact_uncovered


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--primes", required=True)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument(
        "--fill-capacities",
        action="store_true",
        help=(
            "select exactly every group capacity; monotonicity makes this "
            "equivalent to the ordinary at-most-capacity cover problem"
        ),
    )
    parser.add_argument("--phase-output", type=Path)
    parser.add_argument("--result-output", type=Path)
    args = parser.parse_args()
    primes = tuple(int(value) for value in args.primes.split(",") if value)
    if not primes or len(set(primes)) != len(primes):
        raise SystemExit("--primes must be a nonempty unique list")
    period = prod(primes)

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = [
        row
        for row in payload["choices"]
        if int(row["h"]) > 1 and period % int(row["h"]) == 0
    ]
    groups = defaultdict(list)
    for row in rows:
        h = int(row["h"])
        factors = exact_uncovered.factor(h)
        if any(exponent != 1 for exponent in factors.values()):
            raise RuntimeError(f"row modulus {h} is not squarefree")
        if int(row["target_modulus"]) != 1:
            raise RuntimeError(
                "symmetry reduction requires target_modulus=1"
            )
        component_data = []
        for q in sorted(factors):
            a = int(row["a"]) % q
            b = int(row["b"]) % q
            scale = pow(a, -1, q) if a else pow(b, -1, q)
            direction = (a * scale % q, b * scale % q)
            component_data.append((q, direction, scale))
        signature = tuple(
            (q, direction) for q, direction, _scale in component_data
        )
        groups[signature].append((row, tuple(component_data)))

    vpool = IDPool()
    solver = Solver(name=args.solver)
    flat_var = {}
    target_spaces = {}
    clause_count = 0
    started = time.monotonic()
    for signature, members in groups.items():
        target_space = tuple(
            itertools.product(*(range(q) for q, _direction in signature))
        )
        target_spaces[signature] = target_space
        variables = []
        for targets in target_space:
            variable = vpool.id(("flat", signature, targets))
            flat_var[(signature, targets)] = variable
            variables.append(variable)
        if len(members) < len(variables):
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

    local_points = {
        q: tuple((k, l) for k in range(q) for l in range(q))
        for q in primes
    }
    prime_index = {q: index for index, q in enumerate(primes)}
    for point_tuple in itertools.product(
        *(local_points[q] for q in primes)
    ):
        clause = []
        for signature in groups:
            targets = []
            for q, direction in signature:
                k, l = point_tuple[prime_index[q]]
                targets.append(
                    (direction[0] * k + direction[1] * l) % q
                )
            clause.append(flat_var[(signature, tuple(targets))])
        solver.add_clause(clause)
        clause_count += 1
    build_seconds = time.monotonic() - started
    grid_points = prod(q * q for q in primes)
    print(
        f"primes={primes} period={period} rows={len(rows)} "
        f"groups={len(groups)} grid_points={grid_points} "
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
        for key, variable in flat_var.items():
            if variable in model:
                signature, targets = key
                selected[signature].append(targets)
        for signature, members in groups.items():
            targets_list = selected[signature]
            if len(targets_list) > len(members):
                raise AssertionError("group capacity exceeded")
            if not targets_list:
                targets_list = [
                    tuple(0 for _q, _direction in signature)
                ]
            for index, (row, component_data) in enumerate(members):
                canonical_targets = targets_list[
                    min(index, len(targets_list) - 1)
                ]
                congruences = []
                for (
                    canonical_target,
                    (q, _direction, scale),
                ) in zip(canonical_targets, component_data):
                    congruences.append(
                        (
                            canonical_target * pow(scale, -1, q) % q,
                            q,
                        )
                    )
                phases[str(int(row["p"]))] = exact_uncovered.crt(
                    congruences
                )

        for point_tuple in itertools.product(
            *(local_points[q] for q in primes)
        ):
            covered = False
            for row in rows:
                h = int(row["h"])
                row_factors = exact_uncovered.factor(h)
                residues = {
                    q: point_tuple[prime_index[q]]
                    for q in row_factors
                }
                if all(
                    (
                        int(row["a"]) * residues[q][0]
                        + int(row["b"]) * residues[q][1]
                        - phases[str(int(row["p"]))]
                    )
                    % q
                    == 0
                    for q in row_factors
                ):
                    covered = True
                    break
            if not covered:
                raise AssertionError("decoded model misses a grid point")
        if args.phase_output:
            args.phase_output.write_text(json.dumps(phases) + "\n")

    result = {
        "pool": str(args.pool),
        "primes": list(primes),
        "period": period,
        "row_count": len(rows),
        "group_count": len(groups),
        "group_capacities": {
            repr(signature): len(members)
            for signature, members in groups.items()
        },
        "grid_point_count": grid_points,
        "variable_count": vpool.top,
        "clause_count": clause_count,
        "solver": args.solver,
        "fill_capacities": args.fill_capacities,
        "sat": sat,
        "selected_flat_counts": {
            repr(signature): len(selected[signature])
            for signature in groups
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
