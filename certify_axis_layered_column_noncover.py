#!/usr/bin/env python3
"""Certify a prime-deficit obstruction in one layered-pool column.

For a finite family of residue classes modulo ``n_i``, fix a prime ``q``.
If fewer than ``q`` of the moduli are divisible by ``q``, then those classes
are redundant in any putative cover of the integers.  Indeed, after fixing
all prime-to-``q`` coordinates, ``q`` suitable translates are indistinguish-
able to the other classes, while each ``q``-divisible class can hit at most
one translate.  If the reciprocal density of the remaining classes is below
one, no cover exists.

One failed column rules out the declared layered placement even after its
active rows are relaxed to choose their one-dimensional phases independently.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def prime_factors(value: int) -> tuple[int, ...]:
    if value < 1:
        raise ValueError("modulus must be positive")
    factors = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1 if divisor == 2 else 2
    if value > 1:
        factors.append(value)
    return tuple(factors)


def layer_row(row: dict, axis: str) -> dict:
    h = int(row["h"])
    a = int(row["a"]) % h
    b = int(row["b"]) % h
    if h < 1 or math.gcd(a, b, h) != 1:
        raise ValueError(f"invalid affine row for p={row.get('p')}")
    if axis == "k":
        active_modulus = math.gcd(b, h)
        active_coefficient = a
    elif axis == "l":
        active_modulus = math.gcd(a, h)
        active_coefficient = b
    else:
        raise ValueError("layer axis must be k or l")
    if math.gcd(active_coefficient, active_modulus) != 1:
        raise AssertionError("active coordinate is not surjective")
    active_class = int(row["layer_active_class"])
    if not 0 <= active_class < active_modulus:
        raise ValueError("layer active class is outside its modulus")
    expected_target = active_coefficient * active_class % active_modulus
    if (
        int(row["target_modulus"]) != active_modulus
        or int(row["target_residue"]) % active_modulus != expected_target
    ):
        raise ValueError("stored target restriction does not match placement")
    return {
        "p": int(row["p"]),
        "active_modulus": active_modulus,
        "active_class": active_class,
        "residual_modulus": h // active_modulus,
    }


def active_rows(payload: dict, coordinate: int) -> list[dict]:
    axis = payload["layer_axis"]
    rows = [layer_row(row, axis) for row in payload["choices"]]
    return [
        row
        for row in rows
        if coordinate % row["active_modulus"] == row["active_class"]
    ]


def prime_deficit_obstructions(rows: list[dict]) -> list[dict]:
    total = sum(
        (Fraction(1, row["residual_modulus"]) for row in rows),
        Fraction(),
    )
    primes = sorted(
        {
            prime
            for row in rows
            for prime in prime_factors(row["residual_modulus"])
        }
    )
    obstructions = []
    for prime in primes:
        divisible = [
            row for row in rows if row["residual_modulus"] % prime == 0
        ]
        if len(divisible) >= prime:
            continue
        surviving_density = total - sum(
            (
                Fraction(1, row["residual_modulus"])
                for row in divisible
            ),
            Fraction(),
        )
        if surviving_density >= 1:
            continue
        obstructions.append(
            {
                "prime": prime,
                "divisible_rows": divisible,
                "surviving_density": surviving_density,
                "deficit": 1 - surviving_density,
            }
        )
    return obstructions


def best_obstruction(payload: dict, coordinate: int | None) -> tuple[int, dict]:
    period = int(payload["layer_period"])
    if period < 1:
        raise ValueError("layer period must be positive")
    coordinates = range(period) if coordinate is None else (coordinate % period,)
    candidates = []
    for candidate_coordinate in coordinates:
        rows = active_rows(payload, candidate_coordinate)
        for obstruction in prime_deficit_obstructions(rows):
            candidates.append(
                (
                    obstruction["deficit"],
                    -candidate_coordinate,
                    -obstruction["prime"],
                    candidate_coordinate,
                    rows,
                    obstruction,
                )
            )
    if not candidates:
        raise RuntimeError("no prime-deficit column obstruction was found")
    *_rank, selected_coordinate, rows, obstruction = max(candidates)
    return selected_coordinate, {
        "rows": rows,
        **obstruction,
    }


def fraction_record(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def build_certificate(
    source_path: Path,
    payload: dict,
    coordinate: int | None,
) -> dict:
    selected_coordinate, obstruction = best_obstruction(payload, coordinate)
    rows = obstruction["rows"]
    divisible_rows = obstruction["divisible_rows"]
    total_density = sum(
        (Fraction(1, row["residual_modulus"]) for row in rows),
        Fraction(),
    )
    surviving_density = obstruction["surviving_density"]
    prime = obstruction["prime"]
    return {
        "schema": "axis_layered_column_prime_deficit_obstruction_v1",
        "source": str(source_path),
        "source_sha256": sha256(source_path),
        "layer_axis": payload["layer_axis"],
        "coordinate_basis": payload.get("coordinate_basis"),
        "layer_period": int(payload["layer_period"]),
        "column": selected_coordinate,
        "active_row_count": len(rows),
        "total_reciprocal_density": fraction_record(total_density),
        "obstruction_prime": prime,
        "divisible_row_count": len(divisible_rows),
        "divisible_rows": [
            {
                "p": row["p"],
                "residual_modulus": row["residual_modulus"],
            }
            for row in divisible_rows
        ],
        "surviving_row_count": len(rows) - len(divisible_rows),
        "surviving_reciprocal_density": fraction_record(surviving_density),
        "density_deficit": fraction_record(1 - surviving_density),
        "theorem": (
            "if fewer than q residue-class moduli are divisible by prime q, "
            "those classes are redundant in any cover of the integers"
        ),
        "proof": [
            f"only {len(divisible_rows)} < {prime} active residual moduli "
            f"are divisible by q={prime}",
            (
                "the q-divisible classes are therefore redundant by the "
                "q-translate argument"
            ),
            (
                "the exact reciprocal density of all surviving classes is "
                f"{surviving_density.numerator}/"
                f"{surviving_density.denominator} < 1"
            ),
        ],
        "proved_no_independent_column_cover": True,
        "proved_no_declared_layered_cover": True,
        "scope": (
            "the declared finite layered placement; the active residual "
            "phases were relaxed to be independent, so this does not rule "
            "out other placements, other finite pools, or the original "
            "infinite problem"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument(
        "--column",
        type=int,
        help="check only this column; otherwise scan the declared period",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text())
    if payload.get("schema") != "axis_layered_pool_v1":
        raise RuntimeError("input is not an axis-layered pool artifact")
    certificate = build_certificate(args.input, payload, args.column)
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    density = certificate["surviving_reciprocal_density"]
    print(
        f"column={certificate['column']} "
        f"active_rows={certificate['active_row_count']} "
        f"q={certificate['obstruction_prime']} "
        f"q_rows={certificate['divisible_row_count']} "
        f"surviving_density={density['numerator']}/{density['denominator']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
