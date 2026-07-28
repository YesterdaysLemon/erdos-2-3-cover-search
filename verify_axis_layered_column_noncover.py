#!/usr/bin/env python3
"""Independently replay an axis-layered prime-deficit obstruction."""

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


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return value == divisor
        divisor += 1 if divisor == 2 else 2
    return True


def exact_fraction(record: dict) -> Fraction:
    value = Fraction(
        int(record["numerator"]),
        int(record["denominator"]),
    )
    if float(record["decimal"]) != float(value):
        raise AssertionError("stored decimal does not replay")
    return value


def replay(source: dict, certificate: dict) -> dict:
    if source.get("schema") != "axis_layered_pool_v1":
        raise AssertionError("source schema mismatch")
    if (
        certificate.get("schema")
        != "axis_layered_column_prime_deficit_obstruction_v1"
    ):
        raise AssertionError("certificate schema mismatch")
    axis = source["layer_axis"]
    if axis not in {"k", "l"}:
        raise AssertionError("invalid layer axis")
    period = int(source["layer_period"])
    coordinate = int(certificate["column"])
    if period < 1 or not 0 <= coordinate < period:
        raise AssertionError("invalid layer period or column")
    if (
        certificate["layer_axis"] != axis
        or int(certificate["layer_period"]) != period
        or certificate.get("coordinate_basis")
        != source.get("coordinate_basis")
    ):
        raise AssertionError("certificate metadata does not match source")

    seen_primes = set()
    active = []
    for row in source["choices"]:
        row_prime = int(row["p"])
        if row_prime in seen_primes:
            raise AssertionError("source contains a duplicate row prime")
        seen_primes.add(row_prime)
        h = int(row["h"])
        if h < 1:
            raise AssertionError("row modulus must be positive")
        a = int(row["a"]) % h
        b = int(row["b"]) % h
        if math.gcd(a, b, h) != 1:
            raise AssertionError("row map is not surjective")
        if axis == "k":
            active_modulus = math.gcd(b, h)
            active_coefficient = a
        else:
            active_modulus = math.gcd(a, h)
            active_coefficient = b
        if math.gcd(active_coefficient, active_modulus) != 1:
            raise AssertionError("active-coordinate map is not surjective")
        active_class = int(row["layer_active_class"])
        if not 0 <= active_class < active_modulus:
            raise AssertionError("active class is outside its modulus")
        if int(row["target_modulus"]) != active_modulus:
            raise AssertionError("target modulus does not match placement")
        expected_residue = (
            active_coefficient * active_class % active_modulus
        )
        if int(row["target_residue"]) % active_modulus != expected_residue:
            raise AssertionError("target residue does not match placement")
        if coordinate % active_modulus == active_class:
            active.append(
                {
                    "p": row_prime,
                    "residual_modulus": h // active_modulus,
                }
            )

    obstruction_prime = int(certificate["obstruction_prime"])
    if not is_prime(obstruction_prime):
        raise AssertionError("obstruction factor is not prime")
    divisible = [
        row
        for row in active
        if row["residual_modulus"] % obstruction_prime == 0
    ]
    survivors = [
        row
        for row in active
        if row["residual_modulus"] % obstruction_prime
    ]
    total_density = sum(
        (Fraction(1, row["residual_modulus"]) for row in active),
        Fraction(),
    )
    surviving_density = sum(
        (Fraction(1, row["residual_modulus"]) for row in survivors),
        Fraction(),
    )
    divisible_summary = [
        {
            "p": row["p"],
            "residual_modulus": row["residual_modulus"],
        }
        for row in divisible
    ]
    if int(certificate["active_row_count"]) != len(active):
        raise AssertionError("active-row count mismatch")
    if int(certificate["divisible_row_count"]) != len(divisible):
        raise AssertionError("divisible-row count mismatch")
    if certificate["divisible_rows"] != divisible_summary:
        raise AssertionError("divisible-row list mismatch")
    if int(certificate["surviving_row_count"]) != len(survivors):
        raise AssertionError("surviving-row count mismatch")
    if exact_fraction(certificate["total_reciprocal_density"]) != total_density:
        raise AssertionError("total density mismatch")
    if (
        exact_fraction(certificate["surviving_reciprocal_density"])
        != surviving_density
    ):
        raise AssertionError("surviving density mismatch")
    if exact_fraction(certificate["density_deficit"]) != 1 - surviving_density:
        raise AssertionError("density deficit mismatch")
    if not len(divisible) < obstruction_prime:
        raise AssertionError("prime-deficit hypothesis does not hold")
    if not surviving_density < 1:
        raise AssertionError("surviving density is not below one")
    if not (
        certificate.get("proved_no_independent_column_cover") is True
        and certificate.get("proved_no_declared_layered_cover") is True
    ):
        raise AssertionError("certificate does not declare its proved claims")
    return {
        "verified": True,
        "source_row_count": len(source["choices"]),
        "column": coordinate,
        "active_row_count": len(active),
        "obstruction_prime": obstruction_prime,
        "divisible_row_count": len(divisible),
        "prime_deficit_holds": True,
        "surviving_density_numerator": surviving_density.numerator,
        "surviving_density_denominator": surviving_density.denominator,
        "surviving_density_below_one": True,
        "proved_no_independent_column_cover": True,
        "proved_no_declared_layered_cover": True,
        "engine": "independent-python-exact-rational-replay",
        "scope": certificate["scope"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    certificate = json.loads(args.certificate.read_text())
    source_path = Path(certificate["source"])
    if sha256(source_path) != certificate["source_sha256"]:
        raise AssertionError("source SHA-256 mismatch")
    source = json.loads(source_path.read_text())
    report = replay(source, certificate)
    report.update(
        {
            "certificate": str(args.certificate),
            "certificate_sha256": sha256(args.certificate),
            "source": str(source_path),
            "source_sha256": certificate["source_sha256"],
        }
    )
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"verified={report['verified']} column={report['column']} "
        f"q={report['obstruction_prime']} "
        f"q_rows={report['divisible_row_count']} "
        f"density={report['surviving_density_numerator']}/"
        f"{report['surviving_density_denominator']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
