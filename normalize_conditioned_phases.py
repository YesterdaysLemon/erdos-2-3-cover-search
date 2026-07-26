#!/usr/bin/env python3
"""Normalize several conditioned phases by an algebraic-preserving translation."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import exact_uncovered


def valuation(value: int, prime: int) -> int:
    exponent = 0
    while value % prime == 0:
        exponent += 1
        value //= prime
    return exponent


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("phases", type=Path)
    parser.add_argument("--zero-primes", required=True)
    parser.add_argument("--translation-multiple", type=int)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(args.phases.read_text()).items()
    }
    zero_primes = tuple(
        int(value) for value in args.zero_primes.split(",") if value
    )
    missing = set(zero_primes) - set(by_prime)
    if missing:
        raise RuntimeError(f"normalization primes absent from pool: {missing}")
    translation_multiple = (
        args.translation_multiple
        if args.translation_multiple is not None
        else math.prod(
            int(value) for value in payload.get("algebraic_primes", ())
        )
    )

    component_exponents = {}
    for prime in zero_primes:
        row = by_prime[prime]
        for base, exponent in exact_uncovered.factor(int(row["h"])).items():
            component_exponents[base] = max(
                component_exponents.get(base, 0),
                exponent,
            )

    local_solutions = []
    for base, max_exponent in sorted(component_exponents.items()):
        modulus = base**max_exponent
        equations = []
        for prime in zero_primes:
            row = by_prime[prime]
            exponent = valuation(int(row["h"]), base)
            if not exponent:
                continue
            row_modulus = base**exponent
            equations.append(
                (
                    int(row["a"]) * translation_multiple % row_modulus,
                    int(row["b"]) * translation_multiple % row_modulus,
                    -phases.get(prime, int(row["target_residue"]))
                    % row_modulus,
                    row_modulus,
                    prime,
                )
            )
        solution = None
        for u in range(modulus):
            for v in range(modulus):
                if all(
                    (a * u + b * v - target) % equation_modulus == 0
                    for a, b, target, equation_modulus, _prime in equations
                ):
                    solution = (u, v)
                    break
            if solution is not None:
                break
        if solution is None:
            raise RuntimeError(
                f"normalization equations inconsistent modulo {modulus}"
            )
        local_solutions.append((solution[0], solution[1], modulus, equations))

    u = exact_uncovered.crt(
        [(item[0], item[2]) for item in local_solutions]
    )
    v = exact_uncovered.crt(
        [(item[1], item[2]) for item in local_solutions]
    )
    shift_k = translation_multiple * u
    shift_l = translation_multiple * v

    normalized = {}
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        target = phases.get(prime, int(row["target_residue"])) % h
        new_target = (
            target
            + int(row["a"]) * shift_k
            + int(row["b"]) * shift_l
        ) % h
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        if new_target % modulus != residue:
            raise AssertionError(
                f"translation broke target restriction for p={prime}"
            )
        normalized[str(prime)] = new_target
    for prime in zero_primes:
        if normalized[str(prime)] != 0:
            raise AssertionError(f"failed to normalize p={prime}")

    args.output.write_text(json.dumps(normalized) + "\n")
    audit = {
        "pool": str(args.pool),
        "input_phases": str(args.phases),
        "zero_primes": list(zero_primes),
        "translation_multiple": translation_multiple,
        "u": u,
        "v": v,
        "shift_k": shift_k,
        "shift_l": shift_l,
        "component_solutions": [
            {
                "modulus": modulus,
                "u": local_u,
                "v": local_v,
                "equation_primes": [
                    equation[-1] for equation in equations
                ],
            }
            for local_u, local_v, modulus, equations in local_solutions
        ],
        "output": str(args.output),
    }
    if args.audit_output:
        args.audit_output.write_text(json.dumps(audit, indent=2) + "\n")
    print(
        f"rows={len(rows)} zeroed={len(zero_primes)} "
        f"translation_multiple={translation_multiple} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
