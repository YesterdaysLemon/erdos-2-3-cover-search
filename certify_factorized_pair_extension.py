#!/usr/bin/env python3
"""Extend an exact block by a CRT factor row and one lifted factor row.

Let a verified base block depend on a period P and have maximum union
density U.  Let F be a line modulo s, with gcd(P, s) = 1.  Let G be a line
modulo s*r, where gcd(r, P*s) = 1.  If the reductions of F and G modulo s
are jointly surjective, then on the independent s-by-r CRT component

    density(F union G) = 1/s + 1/(s*r) - 1/(s*s*r).

The full block maximum is therefore U + (1-U) times that density.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def recorded_anchor_rows(block: dict) -> list[dict]:
    if block.get("anchor_rows"):
        return list(block["anchor_rows"])
    if block.get("base_certificate"):
        base = json.loads(Path(block["base_certificate"]).read_text())
        return [*recorded_anchor_rows(base), block["extra_row"]]
    raise RuntimeError("base certificate does not expose its anchor rows")


def first_fraction(block: dict, *keys: str) -> Fraction:
    for key in keys:
        if block.get(key) is not None:
            return read_fraction(block[key])
    raise RuntimeError(f"base certificate omits all of {keys}")


def component_pair_union(shared: int, residual: int) -> Fraction:
    return (
        Fraction(1, shared)
        + Fraction(1, shared * residual)
        - Fraction(1, shared * shared * residual)
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("base_certificate", type=Path)
    parser.add_argument("--factor-prime", type=int, required=True)
    parser.add_argument("--lifted-prime", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    base = json.loads(args.base_certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    base_primes = tuple(int(p) for p in base["anchor_primes"])
    required = set(base_primes) | {
        args.factor_prime,
        args.lifted_prime,
    }
    missing = required - by_prime.keys()
    if missing:
        raise RuntimeError(f"anchor primes are absent: {sorted(missing)}")
    if (
        args.factor_prime in base_primes
        or args.lifted_prime in base_primes
        or args.factor_prime == args.lifted_prime
    ):
        raise RuntimeError("new anchor primes must be distinct and absent")

    base_anchor_rows = recorded_anchor_rows(base)
    for recorded in base_anchor_rows:
        row = by_prime[int(recorded["p"])]
        if any(
            int(row[key]) != int(recorded[key])
            for key in ("h", "a", "b")
        ):
            raise RuntimeError("base anchor row differs in the new pool")

    factor = by_prime[args.factor_prime]
    lifted = by_prime[args.lifted_prime]
    factor_h = int(factor["h"])
    lifted_h = int(lifted["h"])
    base_period = math.lcm(
        *(int(row["h"]) for row in base_anchor_rows)
    )
    if math.gcd(factor_h, base_period) != 1:
        raise RuntimeError("factor row modulus is not coprime to base period")
    if lifted_h % factor_h:
        raise RuntimeError("lifted row modulus is not a multiple of factor")
    residual = lifted_h // factor_h
    if residual <= 1 or math.gcd(residual, base_period * factor_h) != 1:
        raise RuntimeError("lifted residual is not a new coprime component")
    for row, label in ((factor, "factor"), (lifted, "lifted")):
        if int(row.get("target_modulus", 1)) != 1:
            raise RuntimeError(f"{label} row target is restricted")
        if math.gcd(int(row["a"]), int(row["b"]), int(row["h"])) != 1:
            raise RuntimeError(f"{label} row target map is not surjective")

    determinant = (
        int(factor["a"]) * int(lifted["b"])
        - int(factor["b"]) * int(lifted["a"])
    )
    if math.gcd(determinant, factor_h) != 1:
        raise RuntimeError(
            "factor and lifted rows are not jointly surjective "
            "on the shared component"
        )

    base_union = first_fraction(
        base,
        "maximum_block_union_density",
        "extended_maximum_union_density",
    )
    base_sum = first_fraction(
        base,
        "block_individual_density_sum",
        "extended_individual_density_sum",
    )
    component_union = component_pair_union(factor_h, residual)
    extended_union = base_union + (1 - base_union) * component_union
    extended_sum = (
        base_sum
        + Fraction(1, factor_h)
        + Fraction(1, lifted_h)
    )
    overlap_loss = extended_sum - extended_union
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    upper_bound = total_density - overlap_loss
    anchor_rows = [
        {
            key: int(row[key]) for key in ("p", "h", "a", "b")
        }
        for row in [*base_anchor_rows, factor, lifted]
    ]
    result = {
        "pool": str(args.pool),
        "row_count": len(rows),
        "base_certificate": str(args.base_certificate),
        "base_anchor_primes": list(base_primes),
        "base_period": base_period,
        "base_maximum_union_density": fraction_payload(base_union),
        "base_individual_density_sum": fraction_payload(base_sum),
        "factor_row": {
            key: int(factor[key]) for key in ("p", "h", "a", "b")
        },
        "lifted_row": {
            key: int(lifted[key]) for key in ("p", "h", "a", "b")
        },
        "shared_modulus": factor_h,
        "residual_modulus": residual,
        "shared_determinant": determinant,
        "shared_pair_jointly_surjective": True,
        "pair_component_union_density": fraction_payload(component_union),
        "anchor_primes": [
            *base_primes,
            args.factor_prime,
            args.lifted_prime,
        ],
        "anchor_rows": anchor_rows,
        "enumerated_period": base_period * factor_h * residual,
        "extended_maximum_union_density": fraction_payload(extended_union),
        "extended_individual_density_sum": fraction_payload(extended_sum),
        "maximum_block_union_density": fraction_payload(extended_union),
        "block_individual_density_sum": fraction_payload(extended_sum),
        "forced_overlap_loss": fraction_payload(overlap_loss),
        "total_pool_density": fraction_payload(total_density),
        "pool_union_density_upper_bound": fraction_payload(upper_bound),
        "proved_no_cover": upper_bound < 1,
        "argument": (
            "The base block is independent of the shared and residual CRT "
            "components. The factor line and lifted line have a jointly "
            "surjective two-target map on the shared component, so their "
            "intersection has density 1/shared^2 there. The lifted line "
            "also has density 1/residual on the new coprime component."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base_union={base_union} shared={factor_h} residual={residual} "
        f"component_union={component_union} extended_union={extended_union} "
        f"loss={overlap_loss} pool_upper={upper_bound} "
        f"proved_no_cover={upper_bound < 1}",
        flush=True,
    )
    return 0 if upper_bound < 1 else 1


if __name__ == "__main__":
    raise SystemExit(main())
