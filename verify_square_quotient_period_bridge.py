#!/usr/bin/env python3
"""Verify that a square quotient is a sheared divisor-period family.

For a primitive affine row and a unimodular change of coordinates, descent
to (Z/n)^2 is equivalent to h dividing n.  This verifier authenticates that
bridge, transports an existing forced-pair obstruction through the shear,
and independently enumerates the small anchor target map.

The result concerns only the declared source pool and square quotient.  It
does not prove a global statement or construct an integer m.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

from build_axis_layered_pool import sha256


def transformed(
    row: dict,
    direction: tuple[int, int],
    transverse: tuple[int, int],
) -> dict:
    h = int(row["h"])
    source_a = int(row["a"]) % h
    source_b = int(row["b"]) % h
    if h <= 1 or math.gcd(source_a, source_b, h) != 1:
        raise AssertionError("source row is not primitive")
    a = (
        source_a * direction[0] + source_b * direction[1]
    ) % h
    b = (
        source_a * transverse[0] + source_b * transverse[1]
    ) % h
    if math.gcd(a, b, h) != 1:
        raise AssertionError("unimodular shear lost primitivity")
    return {
        "p": int(row["p"]),
        "h": h,
        "a": a,
        "b": b,
    }


def joint_index(left: dict, right: dict) -> int:
    h1 = int(left["h"])
    h2 = int(right["h"])
    a1 = int(left["a"]) % h1
    b1 = int(left["b"]) % h1
    a2 = int(right["a"]) % h2
    b2 = int(right["b"]) % h2
    values = (
        a1 * b2 - b1 * a2,
        h1 * a2,
        h1 * b2,
        h2 * a1,
        h2 * b1,
        h1 * h2,
    )
    result = 0
    for value in values:
        result = math.gcd(result, value)
    return result


def enumerate_joint_map(left: dict, right: dict) -> dict:
    h1 = int(left["h"])
    h2 = int(right["h"])
    period = math.lcm(h1, h2)
    counts = {}
    for x in range(period):
        for y in range(period):
            target = (
                (
                    int(left["a"]) * x + int(left["b"]) * y
                )
                % h1,
                (
                    int(right["a"]) * x + int(right["b"]) * y
                )
                % h2,
            )
            counts[target] = counts.get(target, 0) + 1
    target_count = len(counts)
    fibre_sizes = set(counts.values())
    surjective = target_count == h1 * h2
    uniform = len(fibre_sizes) == 1
    density = (
        Fraction(next(iter(fibre_sizes)), period * period)
        if surjective and uniform
        else Fraction()
    )
    return {
        "period": period,
        "cell_count": period * period,
        "target_count": target_count,
        "expected_target_count": h1 * h2,
        "uniform_fibres": uniform,
        "intersection_density": density,
    }


def read_fraction(payload: dict) -> Fraction:
    return Fraction(
        int(payload["numerator"]),
        int(payload["denominator"]),
    )


def verify(
    source_path: Path,
    period_pool_path: Path,
    certificate_path: Path,
    period: int,
    direction: tuple[int, int],
    transverse: tuple[int, int],
) -> dict:
    determinant = (
        direction[0] * transverse[1]
        - direction[1] * transverse[0]
    )
    if period < 1 or abs(determinant) != 1:
        raise AssertionError("invalid square quotient or basis")
    source = json.loads(source_path.read_text())
    period_pool = json.loads(period_pool_path.read_text())
    certificate = json.loads(certificate_path.read_text())
    source_rows = source["choices"]
    if len({int(row["p"]) for row in source_rows}) != len(source_rows):
        raise AssertionError("source pool repeats a prime")

    transformed_rows = [
        transformed(row, direction, transverse)
        for row in source_rows
    ]
    descending = [
        row
        for row in transformed_rows
        if (
            int(row["a"]) * period % int(row["h"]) == 0
            and int(row["b"]) * period % int(row["h"]) == 0
        )
    ]
    divisor_source_rows = [
        row for row in source_rows if period % int(row["h"]) == 0
    ]
    divisor_primes = sorted(int(row["p"]) for row in divisor_source_rows)
    descending_primes = sorted(int(row["p"]) for row in descending)
    if descending_primes != divisor_primes:
        raise AssertionError("square descent is not the divisor-period set")
    if period_pool["choices"] != divisor_source_rows:
        raise AssertionError("reference period pool differs from source")
    if int(period_pool["period_filter"]) != period:
        raise AssertionError("reference period metadata mismatch")
    if Path(period_pool["source"]).name != source_path.name:
        raise AssertionError("reference source filename mismatch")

    row_count = len(descending)
    density = sum(
        (Fraction(1, int(row["h"])) for row in descending),
        Fraction(),
    )
    if (
        Path(certificate["pool"]).name != period_pool_path.name
        or int(certificate["row_count"]) != row_count
        or read_fraction(certificate["total_reciprocal_density"])
        != density
    ):
        raise AssertionError("forced-pair certificate family mismatch")
    anchor_primes = [int(value) for value in certificate["anchor_primes"]]
    if len(anchor_primes) != 2 or anchor_primes[0] == anchor_primes[1]:
        raise AssertionError("forced-pair anchors are invalid")
    source_by_prime = {
        int(row["p"]): row for row in divisor_source_rows
    }
    recorded_anchors = certificate["anchor_rows"]
    if (
        len(recorded_anchors) != 2
        or [int(row["p"]) for row in recorded_anchors]
        != anchor_primes
        or any(
            int(recorded[key]) != int(source_by_prime[prime][key])
            for recorded, prime in zip(
                recorded_anchors,
                anchor_primes,
                strict=True,
            )
            for key in ("p", "h", "a", "b")
        )
    ):
        raise AssertionError("forced-pair anchor source rows mismatch")
    by_prime = {int(row["p"]): row for row in descending}
    if any(prime not in by_prime for prime in anchor_primes):
        raise AssertionError("forced-pair anchor is absent after shear")
    left = by_prime[anchor_primes[0]]
    right = by_prime[anchor_primes[1]]
    index = joint_index(left, right)
    enumeration = enumerate_joint_map(left, right)
    if (
        index != 1
        or int(certificate["joint_image_index"]) != 1
        or certificate["joint_target_map_surjective"] is not True
        or enumeration["target_count"]
        != enumeration["expected_target_count"]
        or not enumeration["uniform_fibres"]
    ):
        raise AssertionError("sheared anchor map is not surjective")
    overlap = enumeration["intersection_density"]
    upper_bound = density - overlap
    if (
        read_fraction(certificate["forced_pair_overlap_density"])
        != overlap
        or read_fraction(certificate["union_density_upper_bound"])
        != upper_bound
        or certificate["proved_no_cover"] is not True
        or upper_bound >= 1
    ):
        raise AssertionError("forced-pair noncover claim mismatch")
    return {
        "schema": "square_quotient_period_bridge_verification_v1",
        "verified": True,
        "source": str(source_path),
        "source_sha256": sha256(source_path),
        "period_pool": str(period_pool_path),
        "period_pool_sha256": sha256(period_pool_path),
        "forced_pair_certificate": str(certificate_path),
        "forced_pair_certificate_sha256": sha256(certificate_path),
        "period": period,
        "basis": {
            "direction": list(direction),
            "transverse": list(transverse),
            "determinant": determinant,
        },
        "row_count": row_count,
        "descending_primes_equal_divisor_period_primes": True,
        "raw_density_numerator": density.numerator,
        "raw_density_denominator": density.denominator,
        "anchor_primes": anchor_primes,
        "joint_image_index": index,
        "enumerated_anchor_period": enumeration["period"],
        "enumerated_anchor_cell_count": enumeration["cell_count"],
        "enumerated_target_count": enumeration["target_count"],
        "forced_overlap_numerator": overlap.numerator,
        "forced_overlap_denominator": overlap.denominator,
        "union_upper_bound_numerator": upper_bound.numerator,
        "union_upper_bound_denominator": upper_bound.denominator,
        "proved_no_declared_square_quotient_cover": True,
        "integer_m_found": False,
        "scope": (
            "all phase assignments of the declared source rows descending "
            "to the declared sheared square quotient only"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("period_pool", type=Path)
    parser.add_argument("forced_pair_certificate", type=Path)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument(
        "--direction",
        nargs=2,
        type=int,
        required=True,
        metavar=("K", "L"),
    )
    parser.add_argument(
        "--transverse",
        nargs=2,
        type=int,
        required=True,
        metavar=("K", "L"),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = verify(
        args.source,
        args.period_pool,
        args.forced_pair_certificate,
        args.period,
        tuple(args.direction),
        tuple(args.transverse),
    )
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"verified={report['verified']} rows={report['row_count']} "
        f"upper={report['union_upper_bound_numerator']}/"
        f"{report['union_upper_bound_denominator']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
