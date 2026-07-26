#!/usr/bin/env python3
"""Exact all-cell residual-density master using OR-Tools CP-SAT.

For a selected prime q, each row h=q^e*s chooses one target modulo s.
On every (k,l) cell modulo lcm(s), the sum of the compatible q-primary
line densities must be at least one.  This is a necessary condition for an
affine cover.  The complete finite grid below is exact; SAT is only a phase
hint, while INFEASIBLE is a finite-family obstruction for the supplied rows
and fixed targets.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path


def q_part(value: int, prime: int) -> int:
    part = 1
    while value % prime == 0:
        value //= prime
        part *= prime
    return part


def crt_pair(
    residue_a: int,
    modulus_a: int,
    residue_b: int,
    modulus_b: int,
) -> int:
    if math.gcd(modulus_a, modulus_b) != 1:
        raise ValueError("CRT moduli must be coprime")
    step = (
        (residue_b - residue_a)
        * pow(modulus_a, -1, modulus_b)
    ) % modulus_b
    return (residue_a + modulus_a * step) % (modulus_a * modulus_b)


def parse_fixed(raw: str) -> dict[int, int]:
    fixed = {}
    for item in raw.split(","):
        if not item:
            continue
        prime, target = item.split(":", 1)
        fixed[int(prime)] = int(target)
    return fixed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--initial-phases", type=Path, required=True)
    parser.add_argument("--fixed-targets", default="")
    parser.add_argument("--max-time", type=float, default=0.0)
    parser.add_argument("--workers", type=int, default=0)
    parser.add_argument("--log-search", action="store_true")
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    args = parser.parse_args()

    ortools_path = Path(os.environ.get("TEMP", ".")) / "erdos203-ortools"
    sys.path.insert(0, str(ortools_path))
    from ortools.sat.python import cp_model  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    if any(int(row.get("target_modulus", 1)) != 1 for row in rows):
        raise RuntimeError("this master requires unrestricted row targets")
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(
            args.initial_phases.read_text()
        ).items()
    }
    row_primes = {int(row["p"]) for row in rows}
    if row_primes - phases.keys():
        raise RuntimeError("initial phase file omits a pool prime")
    fixed = parse_fixed(args.fixed_targets)
    if fixed.keys() - row_primes:
        raise RuntimeError("a fixed prime is absent from the pool")

    row_data = []
    residual_modulus = 1
    coarse_period = 1
    for row in rows:
        h = int(row["h"])
        residual = q_part(h, args.prime)
        other = h // residual
        residual_modulus = max(residual_modulus, residual)
        coarse_period = math.lcm(coarse_period, other)
        row_data.append((other, residual))

    model = cp_model.CpModel()
    variables: dict[tuple[int, int], object] = {}
    variable_count = 0
    for index, (row, (other, _residual)) in enumerate(
        zip(rows, row_data)
    ):
        row_prime = int(row["p"])
        if other == 1 or row_prime in fixed:
            continue
        options = []
        for target in range(other):
            variable = model.new_bool_var(f"r{index}_t{target}")
            variables[index, target] = variable
            options.append(variable)
            variable_count += 1
        model.add_exactly_one(options)
        model.add_hint(
            variables[index, phases[row_prime] % other],
            1,
        )

    constraint_count = 0
    nonzeros = 0
    impossible_cell = None
    build_started = time.monotonic()
    for k in range(coarse_period):
        for l in range(coarse_period):
            constant = 0
            terms = []
            coefficients = []
            for index, (row, (other, residual)) in enumerate(
                zip(rows, row_data)
            ):
                weight = residual_modulus // residual
                required = (
                    int(row["a"]) * k + int(row["b"]) * l
                ) % other
                row_prime = int(row["p"])
                if other == 1:
                    constant += weight
                elif row_prime in fixed:
                    if fixed[row_prime] % other == required:
                        constant += weight
                else:
                    terms.append(variables[index, required])
                    coefficients.append(weight)
            bound = residual_modulus - constant
            if bound <= 0:
                continue
            if sum(coefficients) < bound:
                impossible_cell = (k, l)
                break
            model.add(
                cp_model.LinearExpr.weighted_sum(terms, coefficients)
                >= bound
            )
            constraint_count += 1
            nonzeros += len(terms)
        if impossible_cell is not None:
            break
    build_seconds = time.monotonic() - build_started

    common = {
        "pool": str(args.pool),
        "prime": args.prime,
        "fixed_targets": {
            str(prime): target
            for prime, target in sorted(fixed.items())
        },
        "rows": len(rows),
        "residual_modulus": residual_modulus,
        "coarse_period": coarse_period,
        "cells": coarse_period * coarse_period,
        "variables": variable_count,
        "constraints": constraint_count,
        "nonzeros": nonzeros,
        "build_seconds": build_seconds,
        "engine": "ortools-cp-sat",
    }
    if impossible_cell is not None:
        result = {
            **common,
            "status": "INFEASIBLE",
            "finite_density_core_unsat": True,
            "reason": "one cell lacks enough total available weight",
            "impossible_cell": list(impossible_cell),
        }
        args.result_output.write_text(
            json.dumps(result, indent=2) + "\n"
        )
        print(
            f"INFEASIBLE immediate_cell={impossible_cell}",
            flush=True,
        )
        return 2

    solver = cp_model.CpSolver()
    if args.max_time > 0:
        solver.parameters.max_time_in_seconds = args.max_time
    if args.workers > 0:
        solver.parameters.num_workers = args.workers
    solver.parameters.log_search_progress = args.log_search
    solve_started = time.monotonic()
    status = solver.solve(model)
    solve_seconds = time.monotonic() - solve_started
    status_name = solver.status_name(status)
    stats = {
        **common,
        "status": status_name,
        "solve_seconds": solve_seconds,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time": solver.wall_time,
    }

    if status == cp_model.INFEASIBLE:
        result = {
            **stats,
            "finite_density_core_unsat": True,
        }
        args.result_output.write_text(
            json.dumps(result, indent=2) + "\n"
        )
        print("INFEASIBLE exact all-cell density master", flush=True)
        return 2
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        result = {
            **stats,
            "complete": False,
        }
        args.result_output.write_text(
            json.dumps(result, indent=2) + "\n"
        )
        print(f"INCOMPLETE status={status_name}", flush=True)
        return 1

    answer_phases = dict(phases)
    for index, (row, (other, residual)) in enumerate(
        zip(rows, row_data)
    ):
        row_prime = int(row["p"])
        if row_prime in fixed:
            answer_phases[row_prime] = fixed[row_prime] % int(row["h"])
        elif other > 1:
            target = next(
                candidate
                for candidate in range(other)
                if solver.boolean_value(variables[index, candidate])
            )
            answer_phases[row_prime] = crt_pair(
                target,
                other,
                phases[row_prime] % residual,
                residual,
            )

    minimum = residual_modulus
    minimum_cell = None
    for k in range(coarse_period):
        for l in range(coarse_period):
            density = sum(
                residual_modulus // residual
                for row, (other, residual) in zip(rows, row_data)
                if answer_phases[int(row["p"])] % other
                == (
                    int(row["a"]) * k + int(row["b"]) * l
                ) % other
            )
            if density < minimum:
                minimum = density
                minimum_cell = (k, l)
    if minimum < residual_modulus:
        raise AssertionError(
            f"CP-SAT answer fails cell {minimum_cell}: "
            f"{minimum} < {residual_modulus}"
        )

    args.phase_output.write_text(
        json.dumps(
            {
                str(prime): target
                for prime, target in answer_phases.items()
            }
        )
        + "\n"
    )
    result = {
        **stats,
        "complete": True,
        "necessary_density_condition": True,
        "phase_file": str(args.phase_output),
        "minimum_scaled_density": minimum,
        "minimum_cell": (
            list(minimum_cell) if minimum_cell is not None else None
        ),
    }
    args.result_output.write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(
        "FEASIBLE exact all-cell residual-density condition; "
        "this is not a cover certificate",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
