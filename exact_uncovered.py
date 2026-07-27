#!/usr/bin/env python3
"""Exact counterexample finder for an affine-line covering assignment.

The moduli are decomposed into their maximal prime-power components.  Residues
of k and l on those components are independent by CRT.  A SAT model therefore
corresponds to a genuine exponent pair missed by every selected line; UNSAT is
an exact proof that the lines cover Z^2.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
from pathlib import Path

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def factor(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def crt(residues: list[tuple[int, int]]) -> int:
    value = 0
    modulus = 1
    for residue, new_modulus in residues:
        step = ((residue - value) * pow(modulus, -1, new_modulus)) % new_modulus
        value += modulus * step
        modulus *= new_modulus
    return value % modulus


def target_allowed(row: dict, target: int) -> bool:
    """Whether a synthesized phase target obeys this row's restrictions."""
    h = int(row["h"])
    modulus = int(row.get("target_modulus", 1))
    if modulus < 1 or h % modulus:
        raise ValueError("target modulus must be positive and divide h")
    residue = int(row.get("target_residue", 0)) % modulus
    return target % modulus == residue


def find_uncovered(
    rows: list[dict],
    max_component: int = 256,
    limit: int = 1,
    algebraic_primes: tuple[int, ...] = (),
    solver_name: str = "cadical195",
    diversity_primes: tuple[int, ...] = (),
    diversity_coordinate_moduli: tuple[int, ...] = (),
    diversity_quota: int = 1,
    diversity_target_cap: int = 0,
    sophie_germain: bool = False,
    core_coordinate_cells: tuple[
        tuple[tuple[int, int, int], ...], ...
    ] = (),
    fixed_coordinate_residues: tuple[tuple[int, int, int], ...] = (),
) -> tuple[list[tuple[int, int]], dict]:
    if diversity_quota < 1:
        raise ValueError("diversity_quota must be positive")
    if diversity_target_cap < 0:
        raise ValueError("diversity_target_cap must be nonnegative")
    if diversity_target_cap and not diversity_primes:
        raise ValueError(
            "diversity_target_cap requires at least one diversity prime"
        )
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    maximal: dict[int, int] = {}
    row_factors = []
    for row in rows:
        if not target_allowed(row, int(row["c"])):
            raise ValueError("selected row target violates its restriction")
        factors = factor(int(row["h"]))
        row_factors.append(factors)
        for prime, exponent in factors.items():
            maximal[prime] = max(maximal.get(prime, 0), exponent)
    for prime in algebraic_primes:
        maximal[prime] = max(maximal.get(prime, 0), 1)
    if sophie_germain:
        maximal[2] = max(maximal.get(2, 0), 2)
    requested_coordinate_moduli = tuple(
        dict.fromkeys(
            (
                *diversity_coordinate_moduli,
                *(
                    modulus
                    for cell in core_coordinate_cells
                    for modulus, _kr, _lr in cell
                ),
                *(
                    modulus
                    for modulus, _kr, _lr in fixed_coordinate_residues
                ),
            )
        )
    )
    # A coordinate restriction is part of the ambient CRT domain even when
    # the supplied row core only uses a proper divisor of that component.
    # This occurs naturally when a row of modulus 2 certifies one cell of a
    # mod-32 grid.  Add the requested prime-power level before allocating the
    # component variables; the independent Z3 checker already does the same.
    for modulus in requested_coordinate_moduli:
        factors = factor(modulus)
        if len(factors) != 1:
            raise ValueError(
                f"coordinate modulus {modulus} is not a prime power"
            )
        prime, exponent = next(iter(factors.items()))
        maximal[prime] = max(maximal.get(prime, 0), exponent)
    components = {prime: prime**exponent for prime, exponent in maximal.items()}
    if components and max(components.values()) > max_component:
        raise ValueError(
            f"largest prime-power component {max(components.values())} "
            f"exceeds guard {max_component}"
        )

    pool = IDPool()
    solver = Solver(name=solver_name)
    kval = {}
    lval = {}
    for prime, modulus in components.items():
        ks = [pool.id(f"k_{prime}_{r}") for r in range(modulus)]
        ls = [pool.id(f"l_{prime}_{r}") for r in range(modulus)]
        for r, variable in enumerate(ks):
            kval[prime, r] = variable
        for r, variable in enumerate(ls):
            lval[prime, r] = variable
        solver.append_formula(
            CardEnc.equals(ks, 1, vpool=pool, encoding=EncType.seqcounter).clauses
        )
        solver.append_formula(
            CardEnc.equals(ls, 1, vpool=pool, encoding=EncType.seqcounter).clauses
        )

    # Rows may use a smaller p-power than the maximal component.  Build linked
    # one-hot residue selectors at each used level once, so individual line
    # constraints do not need quadratic tables over maximal residues.
    used_moduli: dict[int, set[int]] = {prime: set() for prime in components}
    for factors in row_factors:
        for prime, exponent in factors.items():
            used_moduli[prime].add(prime**exponent)
    for prime in algebraic_primes:
        used_moduli[prime].add(prime)
    if sophie_germain:
        used_moduli[2].add(4)
    coordinate_spec_by_modulus = {}
    for modulus in requested_coordinate_moduli:
        factors = factor(modulus)
        if len(factors) != 1:
            raise ValueError(
                f"coordinate modulus {modulus} is not a prime power"
            )
        prime, exponent = next(iter(factors.items()))
        if prime not in components or components[prime] % modulus:
            raise ValueError(
                f"coordinate modulus {modulus} is absent"
            )
        used_moduli[prime].add(modulus)
        coordinate_spec_by_modulus[modulus] = (prime, modulus)
    diversity_coordinate_specs = [
        coordinate_spec_by_modulus[modulus]
        for modulus in diversity_coordinate_moduli
    ]
    klevel = {}
    llevel = {}
    for prime, moduli in used_moduli.items():
        maximal_modulus = components[prime]
        for modulus in moduli:
            if modulus == maximal_modulus:
                for residue in range(modulus):
                    klevel[prime, modulus, residue] = kval[prime, residue]
                    llevel[prime, modulus, residue] = lval[prime, residue]
                continue
            coarse_k = [pool.id(f"k_{prime}_{modulus}_{r}") for r in range(modulus)]
            coarse_l = [pool.id(f"l_{prime}_{modulus}_{r}") for r in range(modulus)]
            for residue, variable in enumerate(coarse_k):
                klevel[prime, modulus, residue] = variable
                refinements = [
                    kval[prime, fine]
                    for fine in range(maximal_modulus)
                    if fine % modulus == residue
                ]
                solver.add_clause([-variable, *refinements])
            for residue, variable in enumerate(coarse_l):
                llevel[prime, modulus, residue] = variable
                refinements = [
                    lval[prime, fine]
                    for fine in range(maximal_modulus)
                    if fine % modulus == residue
                ]
                solver.add_clause([-variable, *refinements])
            for fine in range(maximal_modulus):
                solver.add_clause(
                    [-kval[prime, fine], coarse_k[fine % modulus]]
                )
                solver.add_clause(
                    [-lval[prime, fine], coarse_l[fine % modulus]]
                )

    # Reify each distinct local affine equality once.  The original encoding
    # introduced one match variable per satisfying (k,l) pair and repeated
    # the predicate independently for every global row.  Since both
    # coordinates are one-hot, a single equality variable needs only two
    # ternary clauses per possible k (or l):
    #
    #   equality & k_r -> l_f(r)
    #   k_r & l_f(r) -> equality.
    #
    # This is exact, halves the local clauses, removes all pair variables, and
    # lets rows sharing the same prime-power predicate reuse it.
    equality_cache: dict[tuple[int, int, int, int, int], int] = {}
    equality_pair_count = 0

    def component_equality(
        prime: int, modulus: int, a: int, b: int, c: int
    ) -> int:
        nonlocal equality_pair_count
        key = (
            prime,
            modulus,
            a % modulus,
            b % modulus,
            c % modulus,
        )
        cached = equality_cache.get(key)
        if cached is not None:
            return cached
        equality = pool.id(f"component_eq_{len(equality_cache)}")
        equality_cache[key] = equality
        if math.gcd(b, modulus) == 1:
            inverse = pow(b, -1, modulus)
            pairs = (
                (kr, ((c - a * kr) * inverse) % modulus)
                for kr in range(modulus)
            )
            parameter_is_k = True
        elif math.gcd(a, modulus) == 1:
            inverse = pow(a, -1, modulus)
            pairs = (
                (((c - b * lr) * inverse) % modulus, lr)
                for lr in range(modulus)
            )
            parameter_is_k = False
        else:
            raise AssertionError("line map is not surjective on a component")
        for kr, lr in pairs:
            kvar = klevel[prime, modulus, kr]
            lvar = llevel[prime, modulus, lr]
            if parameter_is_k:
                solver.add_clause([-equality, -kvar, lvar])
            else:
                solver.add_clause([-equality, -lvar, kvar])
            solver.add_clause([-kvar, -lvar, equality])
            equality_pair_count += 1
        return equality

    row_activations = []
    activation_row_index = {}
    for index, (row, factors) in enumerate(zip(rows, row_factors)):
        failures = []
        a, b, c = int(row["a"]), int(row["b"]), int(row["c"])
        h = int(row["h"])
        for prime, exponent in factors.items():
            modulus = prime**exponent
            equality = component_equality(prime, modulus, a, b, c)
            failures.append(-equality)
        # To be outside this affine line, at least one local prime-power
        # component of its congruence must fail.
        if core_coordinate_cells:
            activation = pool.id(f"activate_row_{index}")
            row_activations.append(activation)
            activation_row_index[activation] = index
            solver.add_clause([-activation, *failures])
        else:
            solver.add_clause(failures)

    # If m=M^D and an odd prime q|D divides both exponents, then
    # 2^k 3^l m + 1 is X^q+1 and is algebraically composite.  A genuine
    # uncovered witness must therefore avoid (k,l)=(0,0) modulo every such q.
    for prime in algebraic_primes:
        solver.add_clause(
            [-klevel[prime, prime, 0], -llevel[prime, prime, 0]]
        )
    # If 4|D and m=M^D, then k=2 (mod 4), l=0 (mod 4) gives
    # 2^k 3^l m + 1 = 4X^4+1, composite by Sophie Germain's identity.
    if sophie_germain:
        solver.add_clause([-klevel[2, 4, 2], -llevel[2, 4, 0]])

    # A guarded row core is a reusable preservation certificate.  Under the
    # fixed coordinate cell, UNSAT means the rows whose activation literals
    # occur in the core already cover every admissible refinement.  Their
    # phases can therefore be frozen without freezing the rest of the pool.
    if core_coordinate_cells:
        core_records = []
        for cell in core_coordinate_cells:
            cell_assumptions = []
            normalized_cell = []
            for modulus, kr, lr in cell:
                prime, _ = coordinate_spec_by_modulus[modulus]
                kr %= modulus
                lr %= modulus
                cell_assumptions.extend(
                    (
                        klevel[prime, modulus, kr],
                        llevel[prime, modulus, lr],
                    )
                )
                normalized_cell.append((modulus, kr, lr))
            assumptions = [*row_activations, *cell_assumptions]
            if solver.solve(assumptions=assumptions):
                core_records.append(
                    {
                        "cell": normalized_cell,
                        "closed": False,
                        "row_indices": [],
                    }
                )
                continue
            core = set(solver.get_core() or ())
            core_indices = sorted(
                activation_row_index[literal]
                for literal in core
                if literal in activation_row_index
            )
            core_records.append(
                {
                    "cell": normalized_cell,
                    "closed": True,
                    "row_indices": core_indices,
                }
            )
        meta = {
            "components": components,
            "period": math.prod(components.values()),
            "rows": len(rows),
            "algebraic_primes": list(algebraic_primes),
            "sophie_germain": sophie_germain,
            "solver": solver_name,
            "component_equalities": len(equality_cache),
            "component_equality_pairs": equality_pair_count,
            "core_coordinate_cells": core_records,
            "sat": any(not record["closed"] for record in core_records),
        }
        solver.delete()
        return [], meta

    fixed_coordinate_assumptions = []
    normalized_fixed_coordinates = []
    for modulus, kr, lr in fixed_coordinate_residues:
        prime, _ = coordinate_spec_by_modulus[modulus]
        kr %= modulus
        lr %= modulus
        fixed_coordinate_assumptions.extend(
            (
                klevel[prime, modulus, kr],
                llevel[prime, modulus, lr],
            )
        )
        normalized_fixed_coordinates.append((modulus, kr, lr))

    # CaDiCaL otherwise tends to return the same low-residue corner (often
    # k=0) on successive rebuilt checker instances.  A deterministic phase
    # derived from the current targets yields diverse, reproducible exact
    # counterexamples without changing satisfiability.
    phase_seed = sum(
        (index + 1) * (int(row["c"]) + 1) * int(row["p"])
        for index, row in enumerate(rows)
    )
    rng = random.Random(phase_seed)
    diversity_indices = []
    if diversity_primes:
        row_by_prime = {int(row["p"]): index for index, row in enumerate(rows)}
        missing = [prime for prime in diversity_primes if prime not in row_by_prime]
        if missing:
            raise ValueError(f"diversity primes missing from rows: {missing}")
        diversity_indices = [row_by_prime[prime] for prime in diversity_primes]

    witnesses = []
    previous_preference: dict[int, tuple[int, int]] = {}
    fingerprint_counts: dict[tuple, int] = {}
    target_counts: dict[tuple[int, int], int] = {}
    for _ in range(limit):
        # Request a different deterministic CRT corner on every model.  When
        # enumerating a large batch, leaving one phase fixed makes CaDiCaL
        # vary only irrelevant high-order residues and yields many witnesses
        # with the same low-order hole.  Updating two preferred one-hot values
        # per component is cheap and does not alter satisfiability.
        phases = []
        for prime, modulus in components.items():
            if prime in previous_preference:
                old_k, old_l = previous_preference[prime]
                phases.extend((-kval[prime, old_k], -lval[prime, old_l]))
            preferred_k = rng.randrange(modulus)
            preferred_l = rng.randrange(modulus)
            phases.extend(
                (kval[prime, preferred_k], lval[prime, preferred_l])
            )
            previous_preference[prime] = (preferred_k, preferred_l)
        solver.set_phases(phases)
        if not solver.solve(assumptions=fixed_coordinate_assumptions):
            break
        model = {literal for literal in solver.get_model() if literal > 0}
        kres = []
        lres = []
        block = []
        for prime, modulus in components.items():
            kr = next(r for r in range(modulus) if kval[prime, r] in model)
            lr = next(r for r in range(modulus) if lval[prime, r] in model)
            kres.append((kr, modulus))
            lres.append((lr, modulus))
            block.extend([-kval[prime, kr], -lval[prime, lr]])
        k = crt(kres)
        l = crt(lres)
        assert all((int(r["a"]) * k + int(r["b"]) * l - int(r["c"])) % int(r["h"]) for r in rows)
        assert all(k % prime or l % prime for prime in algebraic_primes)
        assert all(
            k % modulus == kr and l % modulus == lr
            for modulus, kr, lr in normalized_fixed_coordinates
        )
        witnesses.append((k, l))
        if not diversity_indices and not diversity_coordinate_specs:
            solver.add_clause(block)
            continue

        # For synthesis, another CRT point with the same targets on a small
        # jointly-surjective anchor family often gives the master an identical
        # or nearly identical lesson.  Block that whole target fingerprint,
        # rather than merely this enormous full CRT assignment.  The clauses
        # below exactly reify each affine equality at the current target.
        row_fingerprint = tuple(
            (
                row_index,
                (
                    int(rows[row_index]["a"]) * k
                    + int(rows[row_index]["b"]) * l
                )
                % int(rows[row_index]["h"]),
            )
            for row_index in diversity_indices
            if target_allowed(
                rows[row_index],
                (
                    int(rows[row_index]["a"]) * k
                    + int(rows[row_index]["b"]) * l
                )
                % int(rows[row_index]["h"]),
            )
        )
        coordinate_fingerprint = tuple(
            (prime, modulus, k % modulus, l % modulus)
            for prime, modulus in diversity_coordinate_specs
        )
        if not row_fingerprint and not coordinate_fingerprint:
            # Every selected row target at this witness is illegal for phase
            # synthesis, so none should constrain diversity.  Exclude only
            # this exact CRT assignment and continue the honest enumeration.
            solver.add_clause(block)
            continue
        fingerprint = (row_fingerprint, coordinate_fingerprint)
        fingerprint_counts[fingerprint] = (
            fingerprint_counts.get(fingerprint, 0) + 1
        )
        if diversity_target_cap:
            for row_index, target in row_fingerprint:
                target_key = (row_index, target)
                target_counts[target_key] = target_counts.get(target_key, 0) + 1
                if target_counts[target_key] != diversity_target_cap:
                    continue
                row = rows[row_index]
                factors = row_factors[row_index]
                a, b = int(row["a"]), int(row["b"])
                target_equalities = [
                    component_equality(
                        prime,
                        prime**exponent,
                        a,
                        b,
                        target,
                    )
                    for prime, exponent in factors.items()
                ]
                # Once this row-target value has supplied the requested
                # number of witnesses, forbid that complete congruence on
                # later models.  Thus changing any one selected fibre phase
                # can absorb at most `diversity_target_cap` members of this
                # checker batch.
                solver.add_clause([-value for value in target_equalities])
        if fingerprint_counts[fingerprint] < diversity_quota:
            # Retain several genuinely different CRT representatives of each
            # coarse fingerprint before excluding the whole cell.  This gives
            # synthesis a controlled cross-section through every residual
            # anchor cell while still guaranteeing eventual enumeration.
            solver.add_clause(block)
            continue

        fingerprint_component_equalities = []
        for row_index, target in row_fingerprint:
            row = rows[row_index]
            factors = row_factors[row_index]
            a, b = int(row["a"]), int(row["b"])
            for prime, exponent in factors.items():
                modulus = prime**exponent
                fingerprint_component_equalities.append(
                    component_equality(
                        prime,
                        modulus,
                        a,
                        b,
                        target,
                    )
                )
        for prime, modulus, kr, lr in coordinate_fingerprint:
            fingerprint_component_equalities.extend(
                (
                    klevel[prime, modulus, kr],
                    llevel[prime, modulus, lr],
                )
            )
        # The whole fingerprint holds iff every local component equality of
        # every anchor holds.  One clause blocks that conjunction directly.
        solver.add_clause(
            [-value for value in fingerprint_component_equalities]
        )

    meta = {
        "components": components,
        "period": math.prod(components.values()),
        "rows": len(rows),
        "algebraic_primes": list(algebraic_primes),
        "sophie_germain": sophie_germain,
        "solver": solver_name,
        "diversity_primes": list(diversity_primes),
        "diversity_coordinate_moduli": list(
            diversity_coordinate_moduli
        ),
        "diversity_quota": diversity_quota,
        "diversity_target_cap": diversity_target_cap,
        "fixed_coordinate_residues": normalized_fixed_coordinates,
        "component_equalities": len(equality_cache),
        "component_equality_pairs": equality_pair_count,
        "sat": bool(witnesses),
    }
    solver.delete()
    return witnesses, meta


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--max-component", type=int, default=256)
    parser.add_argument("--limit", type=int, default=10)
    args = parser.parse_args()
    payload = json.loads(args.candidate.read_text())
    rows = payload.get("choices", payload.get("rows"))
    witnesses, meta = find_uncovered(rows, args.max_component, args.limit)
    print(json.dumps({"meta": meta, "uncovered": witnesses}, indent=2))
    return 1 if witnesses else 0


if __name__ == "__main__":
    raise SystemExit(main())
