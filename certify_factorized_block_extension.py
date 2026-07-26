#!/usr/bin/env python3
"""Extend an exact anchor-block union bound by a coprime row."""

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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("base_certificate", type=Path)
    parser.add_argument("--extra-prime", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    base = json.loads(args.base_certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    base_primes = tuple(int(p) for p in base["anchor_primes"])
    missing = (set(base_primes) | {args.extra_prime}) - by_prime.keys()
    if missing:
        raise RuntimeError(f"anchor primes are absent: {sorted(missing)}")
    if args.extra_prime in base_primes:
        raise RuntimeError("extra prime is already in the base block")
    base_anchor_rows = recorded_anchor_rows(base)
    for recorded in base_anchor_rows:
        row = by_prime[int(recorded["p"])]
        if any(
            int(row[key]) != int(recorded[key])
            for key in ("h", "a", "b")
        ):
            raise RuntimeError("base anchor row differs in the new pool")
    extra = by_prime[args.extra_prime]
    extra_h = int(extra["h"])
    base_period = math.lcm(
        *(int(row["h"]) for row in base_anchor_rows)
    )
    if math.gcd(extra_h, base_period) != 1:
        raise RuntimeError("extra row modulus is not coprime to base period")
    if int(extra.get("target_modulus", 1)) != 1:
        raise RuntimeError("extra row target is restricted")
    if math.gcd(int(extra["a"]), int(extra["b"]), extra_h) != 1:
        raise RuntimeError("extra row target map is not surjective")

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
    extended_union = base_union + (1 - base_union) / extra_h
    extended_sum = base_sum + Fraction(1, extra_h)
    overlap_loss = extended_sum - extended_union
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    upper_bound = total_density - overlap_loss
    proved = upper_bound < 1
    anchor_rows = [
        {
            key: int(row[key]) for key in ("p", "h", "a", "b")
        }
        for row in [*base_anchor_rows, extra]
    ]
    enumerated_period = math.lcm(base_period, extra_h)
    result = {
        "pool": str(args.pool),
        "row_count": len(rows),
        "base_certificate": str(args.base_certificate),
        "base_anchor_primes": list(base_primes),
        "base_period": base_period,
        "base_maximum_union_density": fraction_payload(base_union),
        "base_individual_density_sum": fraction_payload(base_sum),
        "extra_row": {
            key: int(extra[key]) for key in ("p", "h", "a", "b")
        },
        "coprime_factorization": math.gcd(extra_h, base_period) == 1,
        "anchor_primes": [*base_primes, args.extra_prime],
        "anchor_rows": anchor_rows,
        "enumerated_period": enumerated_period,
        "extended_maximum_union_density": fraction_payload(extended_union),
        "extended_individual_density_sum": fraction_payload(extended_sum),
        "maximum_block_union_density": fraction_payload(extended_union),
        "block_individual_density_sum": fraction_payload(extended_sum),
        "forced_overlap_loss": fraction_payload(overlap_loss),
        "total_pool_density": fraction_payload(total_density),
        "pool_union_density_upper_bound": fraction_payload(upper_bound),
        "proved_no_cover": proved,
        "argument": (
            "The extra row depends on a CRT component coprime to the entire "
            "base block. If the base block covers density U, adding one "
            "surjective line of density 1/h on the independent component "
            "covers exactly U + (1-U)/h."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base_union={base_union} extra_h={extra_h} "
        f"extended_union={extended_union} loss={overlap_loss} "
        f"pool_upper={upper_bound} proved_no_cover={proved}",
        flush=True,
    )
    return 0 if proved else 1


if __name__ == "__main__":
    raise SystemExit(main())
