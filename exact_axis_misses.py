#!/usr/bin/env python3
"""Exactly enumerate missed residues on k- or l-axis for a saved phase.

Restricting an affine fibre a*k+b*l=c (mod h) to one coordinate axis gives
either an empty set or one ordinary residue class.  This script reduces every
eligible fibre to that class, adds the perfect-power algebraic classes, and
uses a prime-power CRT SAT encoding to decide whether their union covers Z.
Every reported integer is therefore a genuine missed exponent pair.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
from pathlib import Path

import exact_uncovered
from power_anchor_capacity_lp import power_target_congruence


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--axis", choices=("k", "l"), required=True)
    parser.add_argument("--row-max-component", type=int, default=0)
    parser.add_argument("--max-component", type=int, default=1024)
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument(
        "--diversity-modulus",
        type=int,
        default=0,
        help="if positive, enumerate a quota in each x residue modulo this",
    )
    parser.add_argument("--diversity-quota", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.diversity_modulus < 0:
        raise SystemExit("--diversity-modulus must be nonnegative")
    if args.diversity_quota < 1:
        raise SystemExit("--diversity-quota must be positive")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    source = json.loads(args.pool.read_text())
    phases = json.loads(args.phase_file.read_text())
    classes: set[tuple[int, int]] = set()
    used_rows = 0
    for raw in source["choices"]:
        h = int(raw["h"])
        if args.row_max_component and max(
            (
                prime**exponent
                for prime, exponent in exact_uncovered.factor(h).items()
            ),
            default=1,
        ) > args.row_max_component:
            continue
        p = int(raw["p"])
        try:
            power_residue, power_modulus = power_target_congruence(
                h, p, args.power
            )
        except RuntimeError:
            continue
        c = int(phases.get(str(p), power_residue)) % h
        if c % power_modulus != power_residue:
            raise RuntimeError(f"saved target for prime {p} is incompatible")
        coefficient = int(raw["a" if args.axis == "k" else "b"]) % h
        common = math.gcd(coefficient, h)
        if c % common:
            continue
        modulus = h // common
        if modulus == 1:
            residue = 0
        else:
            residue = (
                (c // common)
                * pow(coefficient // common, -1, modulus)
            ) % modulus
        classes.add((residue, modulus))
        used_rows += 1

    algebraic_primes = tuple(
        prime
        for prime in exact_uncovered.factor(args.power)
        if prime % 2
    )
    classes.update((0, prime) for prime in algebraic_primes)
    if args.axis == "k" and args.power % 4 == 0:
        classes.add((2, 4))

    maximal: dict[int, int] = {}
    for _residue, modulus in classes:
        for prime, exponent in exact_uncovered.factor(modulus).items():
            maximal[prime] = max(maximal.get(prime, 0), exponent)
    if args.diversity_modulus:
        for prime, exponent in exact_uncovered.factor(
            args.diversity_modulus
        ).items():
            maximal[prime] = max(maximal.get(prime, 0), exponent)
    components = {
        prime: prime**exponent for prime, exponent in maximal.items()
    }
    if components and max(components.values()) > args.max_component:
        raise ValueError(
            f"largest prime-power component {max(components.values())} "
            f"exceeds guard {args.max_component}"
        )

    pool = IDPool()
    solver = Solver(name="cadical195")
    values: dict[tuple[int, int], int] = {}
    for prime, modulus in components.items():
        variables = [
            pool.id(f"x_{prime}_{residue}")
            for residue in range(modulus)
        ]
        for residue, variable in enumerate(variables):
            values[prime, residue] = variable
        solver.append_formula(
            CardEnc.equals(
                variables,
                1,
                vpool=pool,
                encoding=EncType.seqcounter,
            ).clauses
        )

    equality_cache: dict[tuple[int, int, int], int] = {}

    def equality(prime: int, modulus: int, residue: int) -> int:
        key = (prime, modulus, residue % modulus)
        if key in equality_cache:
            return equality_cache[key]
        maximal_modulus = components[prime]
        if modulus == maximal_modulus:
            variable = values[prime, residue % modulus]
        else:
            variable = pool.id(
                f"eq_{prime}_{modulus}_{residue % modulus}"
            )
            refinements = [
                values[prime, fine]
                for fine in range(maximal_modulus)
                if fine % modulus == residue % modulus
            ]
            solver.add_clause([-variable, *refinements])
            for refinement in refinements:
                solver.add_clause([-refinement, variable])
        equality_cache[key] = variable
        return variable

    for residue, modulus in classes:
        if modulus == 1:
            solver.add_clause([])
            continue
        local_equalities = []
        for prime, exponent in exact_uncovered.factor(modulus).items():
            prime_power = prime**exponent
            local_equalities.append(
                equality(prime, prime_power, residue % prime_power)
            )
        # A missed axis point must fail this entire residue class.
        solver.add_clause([-variable for variable in local_equalities])

    seed = sum(
        (residue + 1) * (modulus + 1)
        for residue, modulus in classes
    )
    rng = random.Random(seed)
    misses = []
    previous: dict[int, int] = {}
    fingerprint_counts: dict[int, int] = {}
    for _ in range(args.limit):
        preferred = []
        for prime, modulus in components.items():
            if prime in previous:
                preferred.append(-values[prime, previous[prime]])
            residue = rng.randrange(modulus)
            preferred.append(values[prime, residue])
            previous[prime] = residue
        solver.set_phases(preferred)
        if not solver.solve():
            break
        model = {literal for literal in solver.get_model() if literal > 0}
        residues = []
        block = []
        for prime, modulus in components.items():
            residue = next(
                value
                for value in range(modulus)
                if values[prime, value] in model
            )
            residues.append((residue, modulus))
            block.append(-values[prime, residue])
        coordinate = exact_uncovered.crt(residues)
        assert all(
            (coordinate - residue) % modulus
            for residue, modulus in classes
        )
        misses.append(
            (coordinate, 0)
            if args.axis == "k"
            else (0, coordinate)
        )
        if not args.diversity_modulus:
            solver.add_clause(block)
            continue
        fingerprint = coordinate % args.diversity_modulus
        fingerprint_counts[fingerprint] = (
            fingerprint_counts.get(fingerprint, 0) + 1
        )
        if fingerprint_counts[fingerprint] < args.diversity_quota:
            solver.add_clause(block)
            continue
        fingerprint_equalities = []
        for prime, exponent in exact_uncovered.factor(
            args.diversity_modulus
        ).items():
            prime_power = prime**exponent
            fingerprint_equalities.append(
                equality(
                    prime,
                    prime_power,
                    fingerprint % prime_power,
                )
            )
        solver.add_clause(
            [-variable for variable in fingerprint_equalities]
        )
    solver.delete()

    args.output.write_text(json.dumps(misses) + "\n")
    print(
        f"axis={args.axis} rows={used_rows} classes={len(classes)} "
        f"components={len(components)} misses={len(misses)} "
        f"fingerprints={len(fingerprint_counts)} "
        f"diversity_modulus={args.diversity_modulus} "
        f"diversity_quota={args.diversity_quota} "
        f"covered={not misses} output={args.output}",
        flush=True,
    )
    return 1 if misses else 0


if __name__ == "__main__":
    raise SystemExit(main())
