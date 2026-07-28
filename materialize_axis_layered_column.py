#!/usr/bin/env python3
"""Materialize one exact 1D column of an axis-layered affine family.

For a fixed layer coordinate, every active affine row restricts to one
residue class in the transverse coordinate.  Its phase is free modulo the
residual modulus, so replacing the restricted row by ``y = c (mod n)`` is a
bijection on its legal phases.  The resulting artifact is suitable for the
existing exact CEGIS and uncovered-point checkers as a necessary,
single-column covering-system problem.
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


def materialize_column(payload: dict, coordinate: int) -> dict:
    if payload.get("schema") != "axis_layered_pool_v1":
        raise ValueError("input is not an axis-layered pool artifact")
    axis = payload["layer_axis"]
    if axis not in {"k", "l"}:
        raise ValueError("layer axis must be k or l")
    period = int(payload["capacity_pattern_period"])
    if period < 1:
        raise ValueError("capacity pattern period must be positive")
    coordinate %= period

    choices = []
    seen_primes = set()
    for raw in payload["choices"]:
        prime = int(raw["p"])
        if prime in seen_primes:
            raise ValueError(f"duplicate source prime {prime}")
        seen_primes.add(prime)
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
        if not 0 <= active_class < active_modulus:
            raise ValueError("layer active class is outside its modulus")
        if (
            int(raw["target_modulus"]) != active_modulus
            or int(raw["target_residue"]) % active_modulus
            != active_coefficient * active_class % active_modulus
        ):
            raise ValueError("stored target restriction does not match row")
        if coordinate % active_modulus != active_class:
            continue
        residual_modulus = h // active_modulus
        choices.append(
            {
                "p": prime,
                "h": residual_modulus,
                "a": 0,
                "b": 1,
                "ord2": 1,
                "ord3": residual_modulus,
                "c": 0,
                "target_modulus": 1,
                "target_residue": 0,
                "source_h": h,
                "source_active_modulus": active_modulus,
                "source_active_class": active_class,
            }
        )
    if not choices:
        raise ValueError("selected column has no active rows")
    density = sum(
        (Fraction(1, int(row["h"])) for row in choices),
        Fraction(),
    )
    transverse_period = math.lcm(*(int(row["h"]) for row in choices))
    return {
        "schema": "axis_layered_column_relaxation_v1",
        "layer_axis": axis,
        "coordinate_basis": payload.get("coordinate_basis"),
        "layer_period": int(payload["layer_period"]),
        "capacity_pattern_period": period,
        "coordinate": coordinate,
        "row_count": len(choices),
        "transverse_period": transverse_period,
        "capacity_numerator": density.numerator,
        "capacity_denominator": density.denominator,
        "capacity_decimal": float(density),
        "scope": (
            "one exact layer column with independent residual phases; "
            "a no-cover result rules out the source placement, while a "
            "cover does not by itself couple phases across other columns"
        ),
        "choices": choices,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--column", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text())
    result = materialize_column(payload, args.column)
    result["source"] = str(args.input)
    result["source_sha256"] = sha256(args.input)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"column={result['coordinate']} rows={result['row_count']} "
        f"period={result['transverse_period']} "
        f"capacity={result['capacity_numerator']}/"
        f"{result['capacity_denominator']} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
