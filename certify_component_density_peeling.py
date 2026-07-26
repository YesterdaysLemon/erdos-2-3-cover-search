#!/usr/bin/env python3
"""Certify a phase-independent CRT component-density obstruction.

For a prime p and a row whose modulus contains p^e, a primitive affine
congruence occupies exactly a 1/p^e fraction of the p-primary coordinate
plane.  If the incident-row density is below one, some p-coordinate avoids
every incident row, so all those rows may be deleted.  Iterating this peeling
rule to the empty row set proves that no assignment of row phases covers Z^2.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from math import gcd
from pathlib import Path

import exact_uncovered


def component_modulus(h: int, prime: int) -> int:
    exponent = exact_uncovered.factor(h)[prime]
    return prime**exponent


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--min-parent-common", type=int, default=0)
    parser.add_argument(
        "--use-prime-blocking-bound",
        action="store_true",
        help=(
            "also peel prime components using Blokhuis' lower bound "
            "3(p+1)/2 for nontrivial blocking sets in PG(2,p)"
        ),
    )
    parser.add_argument(
        "--use-top-degree-polynomial-obstruction",
        action="store_true",
        help=(
            "also peel prescribed-direction families whose top homogeneous "
            "line polynomial has a forbidden middle coefficient modulo p"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.pool.read_text())
    rows = [
        row
        for row in payload["choices"]
        if int(row["h"]) > 1
        and int(row.get("parent_common", 1)) >= args.min_parent_common
    ]
    initial_count = len(rows)
    steps = []

    while rows:
        primes = sorted(
            {
                prime
                for row in rows
                for prime in exact_uncovered.factor(int(row["h"]))
            }
        )
        candidates = []
        densities = []
        for prime in primes:
            incident = [
                row for row in rows if int(row["h"]) % prime == 0
            ]
            for row in incident:
                if (
                    int(row["a"]) % prime == 0
                    and int(row["b"]) % prime == 0
                ):
                    raise RuntimeError(
                        f"row p={row['p']} is nonprimitive at {prime}"
                    )
            density = sum(
                (
                    Fraction(
                        1,
                        component_modulus(int(row["h"]), prime),
                    )
                    for row in incident
                ),
                Fraction(0),
            )
            directions = Counter()
            homogeneous_coefficients = [1]
            for row in incident:
                a = int(row["a"]) % prime
                b = int(row["b"]) % prime
                scale = pow(a, -1, prime) if a else pow(b, -1, prime)
                directions[
                    (a * scale % prime, b * scale % prime)
                ] += 1
                next_coefficients = [0] * (
                    len(homogeneous_coefficients) + 1
                )
                for exponent, coefficient in enumerate(
                    homogeneous_coefficients
                ):
                    next_coefficients[exponent] = (
                        next_coefficients[exponent]
                        + coefficient * b
                    ) % prime
                    next_coefficients[exponent + 1] = (
                        next_coefficients[exponent + 1]
                        + coefficient * a
                    ) % prime
                homogeneous_coefficients = next_coefficients
            maximum_direction_capacity = max(directions.values())
            blocking_threshold = (
                (3 * prime + 1) // 2 if prime > 2 else None
            )
            densities.append((density, prime, incident))
            reason = None
            if density < 1:
                reason = "density"
            elif (
                args.use_prime_blocking_bound
                and prime > 2
                and maximum_direction_capacity < prime
                and len(incident) < blocking_threshold
            ):
                reason = "prime_projective_blocking_bound"
            polynomial_witness = None
            if args.use_top_degree_polynomial_obstruction:
                degree = len(incident)
                lower = max(0, degree - prime + 1)
                upper = min(degree, prime - 1)
                for exponent_x in range(lower, upper + 1):
                    coefficient = homogeneous_coefficients[exponent_x]
                    if coefficient:
                        polynomial_witness = (
                            exponent_x,
                            degree - exponent_x,
                            coefficient,
                        )
                        reason = "top_degree_polynomial"
                        break
            if reason:
                candidates.append(
                    (
                        len(incident),
                        prime,
                        density,
                        incident,
                        reason,
                        maximum_direction_capacity,
                        blocking_threshold,
                        polynomial_witness,
                    )
                )
        if not candidates:
            break
        if (
            args.use_prime_blocking_bound
            or args.use_top_degree_polynomial_obstruction
        ):
            (
                _incident_count,
                prime,
                density,
                incident,
                reason,
                maximum_direction_capacity,
                blocking_threshold,
                polynomial_witness,
            ) = min(candidates, key=lambda item: (item[0], item[1]))
        else:
            density, prime, incident = min(
                densities,
                key=lambda item: (item[0], item[1]),
            )
            if density >= 1:
                break
            reason = "density"
            directions = Counter()
            for row in incident:
                a = int(row["a"]) % prime
                b = int(row["b"]) % prime
                scale = pow(a, -1, prime) if a else pow(b, -1, prime)
                directions[
                    (a * scale % prime, b * scale % prime)
                ] += 1
            maximum_direction_capacity = max(directions.values())
            blocking_threshold = (
                (3 * prime + 1) // 2 if prime > 2 else None
            )
            polynomial_witness = None
        incident_primes = sorted(int(row["p"]) for row in incident)
        steps.append(
            {
                "prime": prime,
                "reason": reason,
                "density_numerator": density.numerator,
                "density_denominator": density.denominator,
                "density_decimal": float(density),
                "row_count_before": len(rows),
                "removed_row_count": len(incident),
                "maximum_direction_capacity": (
                    maximum_direction_capacity
                ),
                "prime_blocking_line_threshold": blocking_threshold,
                "polynomial_witness": (
                    {
                        "exponent_x": polynomial_witness[0],
                        "exponent_y": polynomial_witness[1],
                        "coefficient_mod_prime": polynomial_witness[2],
                    }
                    if polynomial_witness
                    else None
                ),
                "removed_fibre_primes": incident_primes,
            }
        )
        rows = [
            row for row in rows if int(row["h"]) % prime != 0
        ]

    certificate = {
        "pool": str(args.pool),
        "min_parent_common": args.min_parent_common,
        "use_prime_blocking_bound": args.use_prime_blocking_bound,
        "use_top_degree_polynomial_obstruction": (
            args.use_top_degree_polynomial_obstruction
        ),
        "initial_row_count": initial_count,
        "complete_peeling": not rows,
        "step_count": len(steps),
        "steps": steps,
        "remaining_row_count": len(rows),
        "remaining_fibre_primes": sorted(int(row["p"]) for row in rows),
        "theorem": (
            "At each step the union bound on the selected p-primary "
            "coordinate plane is strictly below one. A coordinate outside "
            "all incident fibres kills those rows. Pairwise CRT combines "
            "the successive component choices."
        ),
        "prime_blocking_reference": (
            "A. Blokhuis, On the size of a blocking set in PG(2,p), "
            "Combinatorica 14 (1994), 111-114, "
            "doi:10.1007/BF01305953"
            if args.use_prime_blocking_bound
            else None
        ),
        "polynomial_obstruction": (
            "If translated affine lines cover F_p^2, their product vanishes "
            "modulo (x^p-x,y^p-y). Its top homogeneous part is the product "
            "of the direction forms. Every coefficient x^j y^(N-j) with "
            "j<p and N-j<p must therefore vanish modulo p."
            if args.use_top_degree_polynomial_obstruction
            else None
        ),
    }
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    print(
        f"rows={initial_count} steps={len(steps)} "
        f"remaining={len(rows)} complete={not rows} output={args.output}",
        flush=True,
    )
    return 0 if not rows else 2


if __name__ == "__main__":
    raise SystemExit(main())
