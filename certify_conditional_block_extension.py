#!/usr/bin/env python3
"""Extend a certified anchor block by one conditional-overlap row.

If a base block has forced overlap loss L and every placement of an extra
row meets the base union in density at least delta, then the enlarged block
has forced overlap loss at least L + delta.  The conditional intersection
certificate may use the full internal structure of the base block.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from certify_block_plus_overlap_forest import (
    read_fraction,
    recorded_anchor_rows,
)


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def recorded_row(row: dict) -> dict:
    return {key: int(row[key]) for key in ("p", "h", "a", "b")}


def same_path(left: str | Path, right: str | Path) -> bool:
    return Path(left).resolve() == Path(right).resolve()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("base_certificate", type=Path)
    parser.add_argument("conditional_certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    base = json.loads(args.base_certificate.read_text())
    conditional = json.loads(args.conditional_certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    if len(by_prime) != len(rows):
        raise RuntimeError("pool contains repeated fibre primes")

    base_anchor_rows = recorded_anchor_rows(base)
    base_primes = tuple(int(row["p"]) for row in base_anchor_rows)
    if base_primes != tuple(int(value) for value in base["anchor_primes"]):
        raise RuntimeError("base anchor order is inconsistent")
    for row in base_anchor_rows:
        prime = int(row["p"])
        if prime not in by_prime or recorded_row(by_prime[prime]) != (
            recorded_row(row)
        ):
            raise RuntimeError("base anchor differs from the pool")

    if not same_path(
        conditional["block_certificate"],
        args.base_certificate,
    ):
        raise RuntimeError("conditional certificate uses another base block")
    if tuple(int(value) for value in conditional["anchor_primes"]) != (
        base_primes
    ):
        raise RuntimeError(
            "conditional certificate does not use the complete base block"
        )
    if [
        recorded_row(row) for row in conditional["anchor_rows"]
    ] != [recorded_row(row) for row in base_anchor_rows]:
        raise RuntimeError("conditional anchor rows differ from the base")

    extra_prime = int(conditional["outside_prime"])
    if extra_prime in base_primes or extra_prime not in by_prime:
        raise RuntimeError("invalid extra prime")
    extra = by_prime[extra_prime]
    if recorded_row(extra) != recorded_row(conditional["outside_row"]):
        raise RuntimeError("conditional outside row differs from the pool")

    base_loss = read_fraction(base["forced_overlap_loss"])
    base_sum = sum(
        (Fraction(1, int(row["h"])) for row in base_anchor_rows),
        Fraction(0),
    )
    recorded_base_sum = base.get("block_individual_density_sum")
    if (
        recorded_base_sum is not None
        and read_fraction(recorded_base_sum) != base_sum
    ):
        raise RuntimeError("base individual-density sum is inconsistent")
    base_union_upper = base_sum - base_loss
    conditional_lower = read_fraction(
        conditional["forced_intersection_density"]
    )
    if not 0 <= conditional_lower <= Fraction(1, int(extra["h"])):
        raise RuntimeError("conditional intersection is out of range")

    extended_sum = base_sum + Fraction(1, int(extra["h"]))
    extended_loss = base_loss + conditional_lower
    extended_union_upper = extended_sum - extended_loss
    anchor_rows = [
        *[recorded_row(row) for row in base_anchor_rows],
        recorded_row(extra),
    ]
    total_density = sum(
        (Fraction(1, int(row["h"])) for row in rows),
        Fraction(0),
    )
    pool_upper = total_density - extended_loss
    result = {
        "schema": "conditional_block_extension_v1",
        "pool": str(args.pool),
        "row_count": len(rows),
        "base_certificate": str(args.base_certificate),
        "conditional_certificate": str(args.conditional_certificate),
        "base_anchor_primes": list(base_primes),
        "base_anchor_rows": [
            recorded_row(row) for row in base_anchor_rows
        ],
        "base_individual_density_sum": fraction_payload(base_sum),
        "base_forced_overlap_loss": fraction_payload(base_loss),
        "base_union_density_upper_bound": fraction_payload(base_union_upper),
        "extra_row": recorded_row(extra),
        "conditional_intersection_lower_bound": fraction_payload(
            conditional_lower
        ),
        "anchor_primes": [*base_primes, extra_prime],
        "anchor_rows": anchor_rows,
        "block_individual_density_sum": fraction_payload(extended_sum),
        "forced_overlap_loss": fraction_payload(extended_loss),
        "block_union_density_upper_bound": fraction_payload(
            extended_union_upper
        ),
        "total_pool_density": fraction_payload(total_density),
        "pool_union_density_upper_bound": fraction_payload(pool_upper),
        "proved_no_cover": pool_upper < 1,
        "argument": (
            "The base block contributes its certified forced overlap loss. "
            "The independently checkable conditional certificate proves "
            "that every target of the appended row intersects the complete "
            "base union by at least the recorded density. Adding that "
            "intersection to the base loss gives a valid upper bound for "
            "the enlarged block union."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base_anchors={len(base_primes)} extra={extra_prime} "
        f"conditional={conditional_lower} loss={extended_loss} "
        f"block_upper={extended_union_upper}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
