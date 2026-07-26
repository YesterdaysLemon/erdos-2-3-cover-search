#!/usr/bin/env python3
"""Independent Z3 bit-vector checker for affine-line covers.

Every prime-power CRT coordinate is a bounded bit vector.  Row congruences are
evaluated in a wider vector that cannot overflow before the exact unsigned
remainder operation.  This avoids the difficult nonlinear integer-modulo
reasoning in exact_uncovered_z3.py while remaining structurally independent of
the primary PySAT one-hot encoding.
"""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path

import exact_uncovered


def find_uncovered(
    rows: list[dict],
    max_component: int = 256,
    limit: int = 1,
    algebraic_primes: tuple[int, ...] = (),
    sophie_germain: bool = False,
    fixed_coordinate_residues: tuple[tuple[int, int, int], ...] = (),
):
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import z3  # type: ignore

    maximal: dict[int, int] = {}
    row_factors = []
    for row in rows:
        factors = exact_uncovered.factor(int(row["h"]))
        row_factors.append(factors)
        for prime, exponent in factors.items():
            maximal[prime] = max(maximal.get(prime, 0), exponent)
    for prime in algebraic_primes:
        maximal[prime] = max(maximal.get(prime, 0), 1)
    if sophie_germain:
        maximal[2] = max(maximal.get(2, 0), 2)
    for modulus, _kr, _lr in fixed_coordinate_residues:
        factors = exact_uncovered.factor(modulus)
        if len(factors) != 1:
            raise ValueError(f"coordinate modulus {modulus} is not a prime power")
        prime, exponent = next(iter(factors.items()))
        maximal[prime] = max(maximal.get(prime, 0), exponent)
    components = {
        prime: prime**exponent for prime, exponent in maximal.items()
    }
    largest = max(components.values(), default=1)
    if largest > max_component:
        raise ValueError(
            f"largest prime-power component {largest} exceeds guard "
            f"{max_component}"
        )

    # For m <= largest, a*k+b*l+(m-c) is below 2*m^2+m.  This width therefore
    # makes every pre-remainder expression exact, with no bit-vector wrap.
    expression_width = max(8, 2 * largest.bit_length() + 2)
    solver = z3.SolverFor("QF_BV")
    kval = {}
    lval = {}
    widths = {}
    for prime, modulus in components.items():
        width = max(1, (modulus - 1).bit_length())
        widths[prime] = width
        kval[prime] = z3.BitVec(f"k_{prime}", width)
        lval[prime] = z3.BitVec(f"l_{prime}", width)
        if modulus != 1 << width:
            solver.add(z3.ULT(kval[prime], z3.BitVecVal(modulus, width)))
            solver.add(z3.ULT(lval[prime], z3.BitVecVal(modulus, width)))

    def widened(value, prime):
        return z3.ZeroExt(expression_width - widths[prime], value)

    def residue(value, prime, modulus):
        return z3.URem(
            widened(value, prime),
            z3.BitVecVal(modulus, expression_width),
        )

    for row, factors in zip(rows, row_factors):
        a, b, c = int(row["a"]), int(row["b"]), int(row["c"])
        failures = []
        for prime, exponent in factors.items():
            modulus = prime**exponent
            kr = residue(kval[prime], prime, modulus)
            lr = residue(lval[prime], prime, modulus)
            expression = (
                z3.BitVecVal(a % modulus, expression_width) * kr
                + z3.BitVecVal(b % modulus, expression_width) * lr
                + z3.BitVecVal((-c) % modulus, expression_width)
            )
            failures.append(
                z3.URem(
                    expression,
                    z3.BitVecVal(modulus, expression_width),
                )
                != 0
            )
        solver.add(z3.Or(*failures))
    for prime in algebraic_primes:
        solver.add(
            z3.Or(
                residue(kval[prime], prime, prime) != 0,
                residue(lval[prime], prime, prime) != 0,
            )
        )
    if sophie_germain:
        solver.add(
            z3.Or(
                residue(kval[2], 2, 4) != 2,
                residue(lval[2], 2, 4) != 0,
            )
        )
    normalized_fixed_coordinates = []
    for modulus, kr, lr in fixed_coordinate_residues:
        prime = next(iter(exact_uncovered.factor(modulus)))
        kr %= modulus
        lr %= modulus
        solver.add(residue(kval[prime], prime, modulus) == kr)
        solver.add(residue(lval[prime], prime, modulus) == lr)
        normalized_fixed_coordinates.append((modulus, kr, lr))

    witnesses = []
    for _ in range(limit):
        if solver.check() != z3.sat:
            break
        model = solver.model()
        kres = []
        lres = []
        blocks = []
        for prime, modulus in components.items():
            kr = model.eval(kval[prime], model_completion=True).as_long()
            lr = model.eval(lval[prime], model_completion=True).as_long()
            kres.append((kr, modulus))
            lres.append((lr, modulus))
            blocks.extend((kval[prime] != kr, lval[prime] != lr))
        k = exact_uncovered.crt(kres)
        l = exact_uncovered.crt(lres)
        assert all(
            (
                int(row["a"]) * k
                + int(row["b"]) * l
                - int(row["c"])
            )
            % int(row["h"])
            for row in rows
        )
        assert all(k % prime or l % prime for prime in algebraic_primes)
        assert not (sophie_germain and k % 4 == 2 and l % 4 == 0)
        assert all(
            k % modulus == kr and l % modulus == lr
            for modulus, kr, lr in normalized_fixed_coordinates
        )
        witnesses.append((k, l))
        solver.add(z3.Or(*blocks))
    return witnesses, {
        "components": components,
        "period": math.prod(components.values()),
        "rows": len(rows),
        "algebraic_primes": list(algebraic_primes),
        "sophie_germain": sophie_germain,
        "sat": bool(witnesses),
        "engine": "z3-bitvector",
        "expression_width": expression_width,
        "fixed_coordinate_residues": normalized_fixed_coordinates,
    }
