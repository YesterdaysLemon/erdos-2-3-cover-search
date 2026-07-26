#!/usr/bin/env python3
"""Build and independently verify the CRT integer from an exact line cover."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import cegis_cover
import exact_uncovered
import exact_uncovered_z3
import search_cover

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def merge(value: int, modulus: int, residue: int, prime: int) -> tuple[int, int]:
    step = ((residue - value) * pow(modulus, -1, prime)) % prime
    return value + modulus * step, modulus * prime


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cover", type=Path)
    parser.add_argument("--max-component", type=int, default=100000)
    parser.add_argument(
        "--materialize-max-digits",
        type=int,
        default=100000,
        help="store a decimal m only up to this estimated digit count",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.cover.read_text())
    rows = payload["choices"]
    power = int(payload.get("power", 1))
    if power < 1:
        raise RuntimeError("cover power must be a positive integer")
    algebraic_primes = tuple(
        prime for prime in exact_uncovered.factor(power) if prime % 2
    )
    sophie_germain = power % 4 == 0
    missed, meta = exact_uncovered.find_uncovered(
        rows,
        max_component=args.max_component,
        limit=1,
        algebraic_primes=algebraic_primes,
        sophie_germain=sophie_germain,
    )
    if missed:
        raise RuntimeError(f"input is not an exact cover; miss={missed[0]}")
    z3_missed, z3_meta = exact_uncovered_z3.find_uncovered(
        rows,
        max_component=args.max_component,
        limit=1,
        algebraic_primes=algebraic_primes,
        sophie_germain=sophie_germain,
    )
    if z3_missed:
        raise RuntimeError(
            f"independent Z3 checker rejects the cover; miss={z3_missed[0]}"
        )

    value, modulus = 1, 6
    congruences = []
    seen = set()
    for row in rows:
        p = int(row["p"])
        if p in seen or p in (2, 3):
            raise RuntimeError(f"invalid or repeated prime {p}")
        seen.add(p)
        ord2 = search_cover.multiplicative_order(2, p)
        ord3 = search_cover.multiplicative_order(3, p)
        h, a, b, ord2, ord3 = cegis_cover.signature(p, ord2, ord3)
        expected = (h, a, b, ord2, ord3)
        actual = tuple(int(row[key]) for key in ("h", "a", "b", "ord2", "ord3"))
        if actual != expected:
            raise RuntimeError(f"signature mismatch for prime {p}")
        root = cegis_cover.primitive_root(p)
        generator = pow(root, (p - 1) // h, p)
        c = int(row["c"]) % h
        m_residue = -pow(pow(generator, c, p), -1, p) % p
        if power == 1:
            base_residue = m_residue
        else:
            divisor = math.gcd(power, p - 1)
            exponent = ((p - 1) // 2 - ((p - 1) // h) * c) % (p - 1)
            if exponent % divisor:
                raise RuntimeError(f"target for prime {p} is not a {power}-th power")
            reduced_modulus = (p - 1) // divisor
            base_exponent = (
                (exponent // divisor)
                * pow(power // divisor, -1, reduced_modulus)
            ) % reduced_modulus
            base_residue = pow(root, base_exponent, p)
            if pow(base_residue, power, p) != m_residue:
                raise AssertionError(p)
        value, modulus = merge(value, modulus, base_residue, p)
        congruences.append(
            {
                "p": p,
                "base_residue": base_residue,
                "m_residue": m_residue,
                "h": h,
                "a": a,
                "b": b,
                "c": c,
            }
        )

    if value <= max(seen):
        value += modulus
    if math.gcd(value, 6) != 1 or value <= max(seen):
        raise AssertionError("bad CRT representative")
    for item in congruences:
        if value % item["p"] != item["base_residue"]:
            raise AssertionError(item["p"])
        if pow(value, power, item["p"]) != item["m_residue"]:
            raise AssertionError(item["p"])

    estimated_digits = (
        1 if value == 1 else math.floor(power * math.log10(value)) + 1
    )
    materialized = estimated_digits <= args.materialize_max_digits
    m_value = value**power if materialized else None
    result = {
        "m": str(m_value) if materialized else f"({value})^{power}",
        "m_materialized": materialized,
        "m_estimated_digits": estimated_digits,
        "base_M": str(value),
        "power": power,
        "crt_modulus": str(modulus),
        "rows": congruences,
        "checker": meta,
        "independent_checker": z3_meta,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"wrote={args.output} primes={len(rows)} power={power} "
        f"m_digits={estimated_digits} materialized={materialized} "
        f"exact_checkers=UNSAT/UNSAT"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
