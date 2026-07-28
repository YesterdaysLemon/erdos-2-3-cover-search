#!/usr/bin/env python3
"""Certify a projected weighted obstruction in one layered column.

Fix one layer column and project its residual coordinate modulo ``D``.  A
residue class modulo ``n`` is active in one class modulo ``g=gcd(n,D)`` and
has conditional density ``1/(n/g)`` there.  Hence every covered projection
cell must have total conditional density at least one.

Two coprime full-cell rows may be translated to target zero simultaneously.
A third full-cell row modulo ``D`` then has only ``D`` possible targets.  For
each target this script discovers a nonnegative integer weight on the still
uncovered cells for which the maximum possible contribution of every other
row sums to strictly less than the total cell weight.  Those exact weights,
not the floating-point discovery LP, form the certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from fractions import Fraction
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_primes(text: str) -> tuple[int, ...]:
    values = tuple(int(value) for value in text.split(",") if value)
    if len(values) != len(set(values)):
        raise ValueError("anchor primes must be distinct")
    return values


def active_projection_rows(
    payload: dict,
    coordinate: int,
    projection: int,
) -> list[dict]:
    if payload.get("schema") != "axis_layered_pool_v1":
        raise ValueError("input is not an axis-layered pool artifact")
    axis = payload["layer_axis"]
    if axis not in {"k", "l"}:
        raise ValueError("layer axis must be k or l")
    if projection < 1:
        raise ValueError("projection must be positive")
    period = int(payload["capacity_pattern_period"])
    coordinate %= period
    rows = []
    seen = set()
    for raw in payload["choices"]:
        prime = int(raw["p"])
        if prime in seen:
            raise ValueError(f"duplicate source prime {prime}")
        seen.add(prime)
        h = int(raw["h"])
        a = int(raw["a"]) % h
        b = int(raw["b"]) % h
        if h < 1 or math.gcd(a, b, h) != 1:
            raise ValueError(f"invalid affine row for p={prime}")
        active_modulus = math.gcd(b if axis == "k" else a, h)
        active_coefficient = a if axis == "k" else b
        if math.gcd(active_coefficient, active_modulus) != 1:
            raise AssertionError("active-coordinate map is not surjective")
        active_class = int(raw["layer_active_class"])
        if (
            not 0 <= active_class < active_modulus
            or int(raw["target_modulus"]) != active_modulus
            or int(raw["target_residue"]) % active_modulus
            != active_coefficient * active_class % active_modulus
        ):
            raise ValueError("stored target restriction does not match row")
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
    if not rows:
        raise ValueError("selected column has no active rows")
    return rows


def maximum_weighted_tail(
    rows: list[dict],
    weights: list[int],
) -> Fraction:
    period = len(weights)
    total = Fraction()
    for row in rows:
        modulus = int(row["projected_modulus"])
        denominator = int(row["conditional_denominator"])
        bucket_maximum = max(
            sum(weights[cell] for cell in range(residue, period, modulus))
            for residue in range(modulus)
        )
        total += Fraction(bucket_maximum, denominator)
    return total


def discover_branch_weights(
    tail: list[dict],
    covered: set[int],
    projection: int,
) -> list[int]:
    dependency_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dependency_path))
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import lil_matrix

    variable_count = projection + len(tail)
    constraint_count = sum(
        int(row["projected_modulus"]) for row in tail
    )
    matrix = lil_matrix((constraint_count, variable_count))
    constraint = 0
    for index, row in enumerate(tail):
        modulus = int(row["projected_modulus"])
        denominator = int(row["conditional_denominator"])
        for residue in range(modulus):
            for cell in range(residue, projection, modulus):
                matrix[constraint, cell] = 1.0 / denominator
            matrix[constraint, projection + index] = -1.0
            constraint += 1
    equality = np.zeros((1, variable_count))
    for cell in range(projection):
        if cell not in covered:
            equality[0, cell] = 1.0
    objective = np.zeros(variable_count)
    objective[projection:] = 1.0
    bounds = [
        (0.0, 0.0) if cell in covered else (0.0, None)
        for cell in range(projection)
    ] + [(0.0, None)] * len(tail)
    result = linprog(
        objective,
        A_ub=matrix.tocsr(),
        b_ub=np.zeros(constraint_count),
        A_eq=equality,
        b_eq=np.ones(1),
        bounds=bounds,
        method="highs",
    )
    if not result.success:
        raise RuntimeError(
            f"dual-weight discovery failed: {result.message}"
        )
    rational = [
        Fraction(float(value)).limit_denominator(1_000_000)
        for value in result.x[:projection]
    ]
    scale = math.lcm(*(value.denominator for value in rational))
    weights = [
        value.numerator * (scale // value.denominator)
        for value in rational
    ]
    divisor = math.gcd(*weights)
    if divisor:
        weights = [value // divisor for value in weights]
    total_weight = sum(weights)
    maximum = maximum_weighted_tail(tail, weights)
    if not maximum < total_weight:
        raise RuntimeError(
            "rationalized discovery weights do not give a strict gap"
        )
    return weights


def build_certificate(
    source_path: Path,
    payload: dict,
    coordinate: int,
    projection: int,
    base_anchor_primes: tuple[int, int],
    branch_anchor_prime: int,
) -> dict:
    rows = active_projection_rows(payload, coordinate, projection)
    by_prime = {row["p"]: row for row in rows}
    anchor_primes = (*base_anchor_primes, branch_anchor_prime)
    if any(prime not in by_prime for prime in anchor_primes):
        raise ValueError("an anchor prime is not active in the column")
    first, second = (by_prime[prime] for prime in base_anchor_primes)
    branch_anchor = by_prime[branch_anchor_prime]
    if (
        first["conditional_denominator"] != 1
        or second["conditional_denominator"] != 1
        or branch_anchor["conditional_denominator"] != 1
    ):
        raise ValueError("all anchors must cover complete projection cells")
    first_modulus = int(first["projected_modulus"])
    second_modulus = int(second["projected_modulus"])
    if (
        math.gcd(first_modulus, second_modulus) != 1
        or math.lcm(first_modulus, second_modulus) != projection
        or int(branch_anchor["projected_modulus"]) != projection
    ):
        raise ValueError(
            "base anchor moduli must be coprime with lcm D and branch "
            "anchor modulus D"
        )
    tail = [row for row in rows if row["p"] not in anchor_primes]
    branches = []
    for target in range(projection):
        covered = {
            cell
            for cell in range(projection)
            if cell % first_modulus == 0
            or cell % second_modulus == 0
            or cell == target
        }
        weights = discover_branch_weights(tail, covered, projection)
        total_weight = sum(weights)
        maximum = maximum_weighted_tail(tail, weights)
        branches.append(
            {
                "branch_anchor_target": target,
                "covered_cells": sorted(covered),
                "cell_weights": weights,
                "total_uncovered_cell_weight": total_weight,
                "maximum_tail_weight_numerator": maximum.numerator,
                "maximum_tail_weight_denominator": maximum.denominator,
                "strict_gap_numerator": (
                    Fraction(total_weight) - maximum
                ).numerator,
                "strict_gap_denominator": (
                    Fraction(total_weight) - maximum
                ).denominator,
            }
        )
    coordinate %= int(payload["capacity_pattern_period"])
    return {
        "schema": "axis_layered_column_projection_dual_obstruction_v1",
        "source": str(source_path),
        "source_sha256": sha256(source_path),
        "layer_axis": payload["layer_axis"],
        "coordinate_basis": payload.get("coordinate_basis"),
        "layer_period": int(payload["layer_period"]),
        "capacity_pattern_period": int(payload["capacity_pattern_period"]),
        "column": coordinate,
        "active_row_count": len(rows),
        "projection_modulus": projection,
        "base_anchor_primes": list(base_anchor_primes),
        "base_anchor_moduli": [first_modulus, second_modulus],
        "branch_anchor_prime": branch_anchor_prime,
        "branch_anchor_modulus": projection,
        "tail_row_count": len(tail),
        "branches": branches,
        "proof": [
            (
                "translate the residual coordinate to put both coprime "
                "base-anchor targets at zero"
            ),
            (
                "enumerate every target of the full-projection branch "
                "anchor"
            ),
            (
                "in each branch, exact nonnegative cell weights have total "
                "weight strictly larger than the sum of every tail row's "
                "maximum possible weighted contribution"
            ),
        ],
        "proved_no_independent_column_cover": True,
        "proved_no_declared_layered_cover": True,
        "scope": (
            "the declared finite layered placement; the selected column's "
            "residual phases were allowed to vary independently, so this "
            "does not rule out other placements, pools, or the original "
            "infinite problem"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--column", type=int, required=True)
    parser.add_argument("--projection", type=int, required=True)
    parser.add_argument("--base-anchor-primes", required=True)
    parser.add_argument("--branch-anchor-prime", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = parse_primes(args.base_anchor_primes)
    if len(base) != 2:
        raise SystemExit("--base-anchor-primes must contain exactly two primes")
    payload = json.loads(args.input.read_text())
    certificate = build_certificate(
        args.input,
        payload,
        args.column,
        args.projection,
        base,
        args.branch_anchor_prime,
    )
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    gaps = [
        Fraction(
            branch["strict_gap_numerator"],
            branch["strict_gap_denominator"],
        )
        for branch in certificate["branches"]
    ]
    print(
        f"column={certificate['column']} "
        f"projection={certificate['projection_modulus']} "
        f"branches={len(gaps)} minimum_gap={min(gaps)} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
