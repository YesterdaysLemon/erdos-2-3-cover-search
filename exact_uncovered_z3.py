#!/usr/bin/env python3
"""Z3 modular checker for an uncovered exponent pair."""

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
    components = {prime: prime**exponent for prime, exponent in maximal.items()}
    if components and max(components.values()) > max_component:
        raise ValueError(
            f"largest prime-power component {max(components.values())} exceeds guard {max_component}"
        )

    solver = z3.Solver()
    kval = {prime: z3.Int(f"k_{prime}") for prime in components}
    lval = {prime: z3.Int(f"l_{prime}") for prime in components}
    for prime, modulus in components.items():
        solver.add(kval[prime] >= 0, kval[prime] < modulus)
        solver.add(lval[prime] >= 0, lval[prime] < modulus)
    for row, factors in zip(rows, row_factors):
        a, b, c = int(row["a"]), int(row["b"]), int(row["c"])
        failures = []
        for prime, exponent in factors.items():
            modulus = prime**exponent
            failures.append((a * kval[prime] + b * lval[prime] - c) % modulus != 0)
        solver.add(z3.Or(*failures))
    for prime in algebraic_primes:
        solver.add(z3.Or(kval[prime] % prime != 0, lval[prime] % prime != 0))
    if sophie_germain:
        solver.add(z3.Or(kval[2] % 4 != 2, lval[2] % 4 != 0))

    witnesses = []
    for _ in range(limit):
        if solver.check() != z3.sat:
            break
        model = solver.model()
        kres = []
        lres = []
        block = []
        for prime, modulus in components.items():
            kr = model.eval(kval[prime]).as_long()
            lr = model.eval(lval[prime]).as_long()
            kres.append((kr, modulus))
            lres.append((lr, modulus))
            block.extend((kval[prime] != kr, lval[prime] != lr))
        k = exact_uncovered.crt(kres)
        l = exact_uncovered.crt(lres)
        assert all(
            (int(row["a"]) * k + int(row["b"]) * l - int(row["c"]))
            % int(row["h"])
            for row in rows
        )
        assert all(k % prime or l % prime for prime in algebraic_primes)
        assert not (sophie_germain and k % 4 == 2 and l % 4 == 0)
        witnesses.append((k, l))
        solver.add(z3.Or(*block))
    return witnesses, {
        "components": components,
        "period": math.prod(components.values()),
        "rows": len(rows),
        "algebraic_primes": list(algebraic_primes),
        "sophie_germain": sophie_germain,
        "sat": bool(witnesses),
        "engine": "z3",
    }
