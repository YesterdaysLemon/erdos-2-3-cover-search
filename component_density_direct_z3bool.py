#!/usr/bin/env python3
"""Exact all-cell residual-density master using Z3 Boolean PB constraints.

For a selected prime q, each row h=q^e*s chooses one target modulo s.
On every (k,l) cell modulo lcm(s), the sum of compatible q-primary line
densities must be at least one.  This is a necessary condition for an affine
cover.  The formulation here is independent of the OR-Tools CP-SAT encoding:
one Z3 Boolean is created for every row-target pair, with native exactly-one
and weighted pseudo-Boolean constraints.
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
    parser.add_argument("--timeout-ms", type=int, default=0)
    parser.add_argument(
        "--solver-logic",
        choices=("generic", "QF_FD", "QF_BV"),
        default="generic",
    )
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import z3  # type: ignore

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

    solver = (
        z3.Solver()
        if args.solver_logic == "generic"
        else z3.SolverFor(args.solver_logic)
    )
    if args.timeout_ms > 0:
        solver.set(timeout=args.timeout_ms)
    variables: dict[tuple[int, int], object] = {}
    variable_count = 0
    exactly_one_constraints = 0
    for index, (row, (other, _residual)) in enumerate(
        zip(rows, row_data)
    ):
        row_prime = int(row["p"])
        if other == 1 or row_prime in fixed:
            continue
        options = []
        for target in range(other):
            variable = z3.Bool(f"r{index}_t{target}")
            variables[index, target] = variable
            options.append((variable, 1))
            variable_count += 1
        solver.add(z3.PbEq(options, 1))
        exactly_one_constraints += 1

    constraint_count = 0
    nonzeros = 0
    impossible_cell = None
    build_started = time.monotonic()
    for k in range(coarse_period):
        for l in range(coarse_period):
            constant = 0
            terms = []
            available = 0
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
                    terms.append((variables[index, required], weight))
                    available += weight
            bound = residual_modulus - constant
            if bound <= 0:
                continue
            if available < bound:
                impossible_cell = (k, l)
                break
            solver.add(z3.PbGe(terms, bound))
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
        "exactly_one_constraints": exactly_one_constraints,
        "density_constraints": constraint_count,
        "nonzeros": nonzeros,
        "build_seconds": build_seconds,
        "engine": (
            "z3-bool-pb"
            if args.solver_logic == "generic"
            else f"z3-{args.solver_logic.lower()}-bool-pb"
        ),
        "solver_logic": args.solver_logic,
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

    solve_started = time.monotonic()
    status = solver.check()
    solve_seconds = time.monotonic() - solve_started
    status_name = str(status)
    stats = {
        **common,
        "status": status_name,
        "solve_seconds": solve_seconds,
        "statistics": str(solver.statistics()),
    }

    if status == z3.unsat:
        result = {
            **stats,
            "status": "INFEASIBLE",
            "finite_density_core_unsat": True,
        }
        args.result_output.write_text(
            json.dumps(result, indent=2) + "\n"
        )
        print("INFEASIBLE exact all-cell density master", flush=True)
        return 2
    if status != z3.sat:
        result = {
            **stats,
            "status": "UNKNOWN",
            "complete": False,
            "reason": solver.reason_unknown(),
        }
        args.result_output.write_text(
            json.dumps(result, indent=2) + "\n"
        )
        print(
            f"INCOMPLETE status={status_name} "
            f"reason={solver.reason_unknown()}",
            flush=True,
        )
        return 1

    model = solver.model()
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
                if z3.is_true(
                    model.eval(
                        variables[index, candidate],
                        model_completion=True,
                    )
                )
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
            scaled_density = 0
            for row, (other, residual) in zip(rows, row_data):
                if (
                    int(row["a"]) * k
                    + int(row["b"]) * l
                    - answer_phases[int(row["p"])]
                ) % other == 0:
                    scaled_density += residual_modulus // residual
            if scaled_density < minimum:
                minimum = scaled_density
                minimum_cell = (k, l)
    if minimum < residual_modulus:
        raise AssertionError("extracted Z3 phase violates a density cell")

    args.phase_output.write_text(
        json.dumps(
            {
                str(prime): target
                for prime, target in sorted(answer_phases.items())
            },
            indent=2,
        )
        + "\n"
    )
    result = {
        **stats,
        "status": "SAT",
        "finite_density_core_unsat": False,
        "minimum_scaled_density": minimum,
        "minimum_cell": (
            list(minimum_cell) if minimum_cell is not None else None
        ),
        "phase_output": str(args.phase_output),
    }
    args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"SAT minimum_scaled_density={minimum}/{residual_modulus}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
