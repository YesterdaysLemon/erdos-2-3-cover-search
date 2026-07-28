#!/usr/bin/env python3
"""Independently replay a layered-column projection dual obstruction."""

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


def reconstruct_rows(
    source: dict,
    coordinate: int,
    projection: int,
) -> list[dict]:
    if source.get("schema") != "axis_layered_pool_v1":
        raise AssertionError("source schema mismatch")
    axis = source["layer_axis"]
    if axis not in {"k", "l"}:
        raise AssertionError("invalid layer axis")
    pattern_period = int(source["capacity_pattern_period"])
    if pattern_period < 1 or projection < 1:
        raise AssertionError("invalid period")
    coordinate %= pattern_period
    rows = []
    seen = set()
    for raw in source["choices"]:
        prime = int(raw["p"])
        if prime in seen:
            raise AssertionError("duplicate source prime")
        seen.add(prime)
        h = int(raw["h"])
        if h < 1:
            raise AssertionError("row modulus must be positive")
        a = int(raw["a"]) % h
        b = int(raw["b"]) % h
        if math.gcd(a, b, h) != 1:
            raise AssertionError("row map is not surjective")
        active_modulus = math.gcd(b if axis == "k" else a, h)
        coefficient = a if axis == "k" else b
        if math.gcd(coefficient, active_modulus) != 1:
            raise AssertionError("active-coordinate map is not surjective")
        active_class = int(raw["layer_active_class"])
        if (
            not 0 <= active_class < active_modulus
            or int(raw["target_modulus"]) != active_modulus
            or int(raw["target_residue"]) % active_modulus
            != coefficient * active_class % active_modulus
        ):
            raise AssertionError("source target restriction mismatch")
        if coordinate % active_modulus != active_class:
            continue
        residual_modulus = h // active_modulus
        projected_modulus = math.gcd(residual_modulus, projection)
        rows.append(
            {
                "p": prime,
                "residual_modulus": residual_modulus,
                "projected_modulus": projected_modulus,
                "conditional_denominator": (
                    residual_modulus // projected_modulus
                ),
            }
        )
    return rows


def replay(source: dict, certificate: dict) -> dict:
    if (
        certificate.get("schema")
        != "axis_layered_column_projection_dual_obstruction_v1"
    ):
        raise AssertionError("certificate schema mismatch")
    coordinate = int(certificate["column"])
    projection = int(certificate["projection_modulus"])
    rows = reconstruct_rows(source, coordinate, projection)
    if int(certificate["active_row_count"]) != len(rows):
        raise AssertionError("active row count mismatch")
    by_prime = {row["p"]: row for row in rows}
    base_primes = tuple(int(p) for p in certificate["base_anchor_primes"])
    if len(base_primes) != 2 or len(set(base_primes)) != 2:
        raise AssertionError("invalid base anchors")
    branch_prime = int(certificate["branch_anchor_prime"])
    anchor_primes = (*base_primes, branch_prime)
    if len(set(anchor_primes)) != 3:
        raise AssertionError("anchor primes are not distinct")
    if any(prime not in by_prime for prime in anchor_primes):
        raise AssertionError("anchor is not active")
    first, second = (by_prime[prime] for prime in base_primes)
    branch_anchor = by_prime[branch_prime]
    first_modulus = int(first["projected_modulus"])
    second_modulus = int(second["projected_modulus"])
    if (
        int(first["conditional_denominator"]) != 1
        or int(second["conditional_denominator"]) != 1
        or int(branch_anchor["conditional_denominator"]) != 1
        or math.gcd(first_modulus, second_modulus) != 1
        or math.lcm(first_modulus, second_modulus) != projection
        or int(branch_anchor["projected_modulus"]) != projection
        or list(certificate["base_anchor_moduli"])
        != [first_modulus, second_modulus]
        or int(certificate["branch_anchor_modulus"]) != projection
    ):
        raise AssertionError("anchor decomposition mismatch")
    tail = [row for row in rows if row["p"] not in anchor_primes]
    if int(certificate["tail_row_count"]) != len(tail):
        raise AssertionError("tail row count mismatch")

    branches = certificate["branches"]
    if [int(branch["branch_anchor_target"]) for branch in branches] != list(
        range(projection)
    ):
        raise AssertionError("branch targets are not exhaustive")
    minimum_gap = None
    for branch in branches:
        target = int(branch["branch_anchor_target"])
        covered = {
            cell
            for cell in range(projection)
            if cell % first_modulus == 0
            or cell % second_modulus == 0
            or cell == target
        }
        if branch["covered_cells"] != sorted(covered):
            raise AssertionError("covered-cell list mismatch")
        weights = [int(value) for value in branch["cell_weights"]]
        if (
            len(weights) != projection
            or any(value < 0 for value in weights)
            or any(weights[cell] for cell in covered)
        ):
            raise AssertionError("invalid branch cell weights")
        total_weight = sum(weights)
        if (
            total_weight <= 0
            or int(branch["total_uncovered_cell_weight"]) != total_weight
        ):
            raise AssertionError("invalid total cell weight")
        maximum = Fraction()
        for row in tail:
            modulus = int(row["projected_modulus"])
            denominator = int(row["conditional_denominator"])
            bucket_maximum = max(
                sum(
                    weights[cell]
                    for cell in range(residue, projection, modulus)
                )
                for residue in range(modulus)
            )
            maximum += Fraction(bucket_maximum, denominator)
        if maximum != Fraction(
            int(branch["maximum_tail_weight_numerator"]),
            int(branch["maximum_tail_weight_denominator"]),
        ):
            raise AssertionError("maximum tail weight mismatch")
        gap = Fraction(total_weight) - maximum
        if gap != Fraction(
            int(branch["strict_gap_numerator"]),
            int(branch["strict_gap_denominator"]),
        ):
            raise AssertionError("strict gap mismatch")
        if gap <= 0:
            raise AssertionError("branch dual does not have a strict gap")
        minimum_gap = gap if minimum_gap is None else min(minimum_gap, gap)
    if not (
        certificate.get("proved_no_independent_column_cover") is True
        and certificate.get("proved_no_declared_layered_cover") is True
    ):
        raise AssertionError("certificate does not declare proved claims")
    assert minimum_gap is not None
    return {
        "verified": True,
        "column": coordinate,
        "projection_modulus": projection,
        "active_row_count": len(rows),
        "tail_row_count": len(tail),
        "branch_count": len(branches),
        "minimum_strict_gap_numerator": minimum_gap.numerator,
        "minimum_strict_gap_denominator": minimum_gap.denominator,
        "proved_no_independent_column_cover": True,
        "proved_no_declared_layered_cover": True,
        "engine": "independent-python-exact-rational-weight-replay",
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
        f"projection={report['projection_modulus']} "
        f"branches={report['branch_count']} "
        f"minimum_gap={report['minimum_strict_gap_numerator']}/"
        f"{report['minimum_strict_gap_denominator']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
