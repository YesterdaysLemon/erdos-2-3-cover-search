#!/usr/bin/env python3
"""Rank finite divisor-period subfamilies by exact reciprocal density."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import exact_uncovered


def primary_parts(value: int) -> dict[int, int]:
    return {
        prime: prime**exponent
        for prime, exponent in exact_uncovered.factor(value).items()
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--min-density", type=float, default=1.0)
    parser.add_argument("--limit", type=int, default=200)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.min_density < 0 or args.limit < 1:
        raise SystemExit("invalid arguments")

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    row_parts = [primary_parts(int(row["h"])) for row in rows]
    primes = sorted({prime for parts in row_parts for prime in parts})
    levels = {
        prime: [1]
        + sorted(
            {
                parts[prime]
                for parts in row_parts
                if prime in parts
            }
        )
        for prime in primes
    }

    families = []
    enumerated = 0
    for selected_levels in itertools.product(*(levels[p] for p in primes)):
        enumerated += 1
        if all(level == 1 for level in selected_levels):
            continue
        chosen = dict(zip(primes, selected_levels))
        selected_rows = [
            row
            for row, parts in zip(rows, row_parts)
            if all(parts.get(prime, 1) <= chosen[prime] for prime in primes)
        ]
        density = sum(
            (Fraction(1, int(row["h"])) for row in selected_rows),
            Fraction(0),
        )
        if float(density) < args.min_density:
            continue
        components = [level for level in selected_levels if level > 1]
        period = math.prod(components)
        families.append(
            {
                "components": components,
                "period": period,
                "cells": period * period,
                "rows": len(selected_rows),
                "distinct_moduli": len(
                    {int(row["h"]) for row in selected_rows}
                ),
                "density_numerator": density.numerator,
                "density_denominator": density.denominator,
                "density": float(density),
            }
        )

    by_period = sorted(
        families,
        key=lambda item: (
            item["period"],
            item["rows"],
            -item["density"],
        ),
    )
    by_rows = sorted(
        families,
        key=lambda item: (
            item["rows"],
            item["period"],
            -item["density"],
        ),
    )
    by_density = sorted(
        families,
        key=lambda item: (
            -item["density"],
            item["rows"],
            item["period"],
        ),
    )
    result = {
        "pool": str(args.pool),
        "prime_levels": {str(p): levels[p] for p in primes},
        "enumerated_families": enumerated,
        "eligible_families": len(families),
        "minimum_density": args.min_density,
        "smallest_period": by_period[: args.limit],
        "fewest_rows": by_rows[: args.limit],
        "highest_density": by_density[: args.limit],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"enumerated={enumerated} eligible={len(families)}",
        flush=True,
    )
    for family in by_period[:10]:
        print(
            f"period={family['period']} components={family['components']} "
            f"rows={family['rows']} density={family['density']:.12f}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
