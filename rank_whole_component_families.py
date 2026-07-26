#!/usr/bin/env python3
"""Rank exact whole-component subfamilies for symmetry-reduced CEGIS.

For each rational prime, choose at most one prime-power component level.
A row is retained only when each of its nontrivial primary parts is exactly
the chosen level.  Its modulus is then a product of whole pairwise-coprime
components, which is the input contract of ``finite_component_cegis.py``.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

import exact_uncovered


def primary_parts(value: int) -> dict[int, int]:
    return {
        prime: prime**exponent
        for prime, exponent in exact_uncovered.factor(value).items()
    }


def normalized_direction(row: dict, modulus: int) -> tuple[int, int]:
    a = int(row["a"]) % modulus
    b = int(row["b"]) % modulus
    if math.gcd(a, modulus) == 1:
        scale = pow(a, -1, modulus)
    elif math.gcd(b, modulus) == 1:
        scale = pow(b, -1, modulus)
    else:
        raise RuntimeError(
            f"row p={row['p']} has no unit coefficient modulo {modulus}"
        )
    return a * scale % modulus, b * scale % modulus


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--min-density", type=float, default=1.0)
    parser.add_argument("--max-flat-variables", type=int, default=5_000_000)
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.min_density < 0 or args.max_flat_variables < 1:
        raise SystemExit("invalid density or variable limit")
    if args.limit < 1:
        raise SystemExit("--limit must be positive")

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

    ranked = []
    enumerated = 0
    density_eligible = 0
    maximum_density = Fraction(0)
    maximum_density_summary = None
    for selected_levels in itertools.product(*(levels[p] for p in primes)):
        enumerated += 1
        if all(level == 1 for level in selected_levels):
            continue
        chosen = dict(zip(primes, selected_levels))
        selected_rows = []
        for row, parts in zip(rows, row_parts):
            if all(chosen[prime] == part for prime, part in parts.items()):
                selected_rows.append(row)
        density = sum(
            (Fraction(1, int(row["h"])) for row in selected_rows),
            Fraction(0),
        )
        if density > maximum_density:
            maximum_density = density
            maximum_density_summary = {
                "components": [
                    level for level in selected_levels if level > 1
                ],
                "period": math.prod(
                    level for level in selected_levels if level > 1
                ),
                "rows": len(selected_rows),
                "density_numerator": density.numerator,
                "density_denominator": density.denominator,
                "density": float(density),
            }
        if float(density) < args.min_density:
            continue
        density_eligible += 1
        components = tuple(level for level in selected_levels if level > 1)
        period = math.prod(components)
        groups = Counter()
        for row in selected_rows:
            h = int(row["h"])
            signature = tuple(
                (component, normalized_direction(row, component))
                for component in components
                if h % component == 0
            )
            if math.prod(component for component, _direction in signature) != h:
                raise AssertionError("row is not a product of whole components")
            groups[signature] += 1
        flat_variables = sum(
            math.prod(component for component, _direction in signature)
            for signature in groups
        )
        if flat_variables > args.max_flat_variables:
            continue
        saturated_groups = [
            {
                "signature": [
                    [component, list(direction)]
                    for component, direction in signature
                ],
                "members": member_count,
                "flats": math.prod(
                    component for component, _direction in signature
                ),
            }
            for signature, member_count in groups.items()
            if member_count
            >= math.prod(component for component, _direction in signature)
        ]
        ranked.append(
            {
                "components": list(components),
                "period": period,
                "rows": len(selected_rows),
                "groups": len(groups),
                "flat_variables": flat_variables,
                "density_numerator": density.numerator,
                "density_denominator": density.denominator,
                "density": float(density),
                "saturated_groups": saturated_groups,
            }
        )

    ranked.sort(
        key=lambda item: (
            bool(item["saturated_groups"]) is False,
            item["flat_variables"],
            item["period"],
            -item["density"],
        )
    )
    result = {
        "pool": str(args.pool),
        "prime_levels": {str(p): levels[p] for p in primes},
        "enumerated_families": enumerated,
        "density_eligible_families": density_eligible,
        "maximum_density_family": maximum_density_summary,
        "max_flat_variables": args.max_flat_variables,
        "reported_families": min(args.limit, len(ranked)),
        "families": ranked[: args.limit],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"enumerated={enumerated} density_eligible={density_eligible} "
        f"within_variable_guard={len(ranked)} reported={len(result['families'])}",
        flush=True,
    )
    if maximum_density_summary is not None:
        print(
            f"maximum_density={maximum_density_summary['density']:.12f} "
            f"components={maximum_density_summary['components']} "
            f"rows={maximum_density_summary['rows']}",
            flush=True,
        )
    for family in result["families"][:10]:
        print(
            f"components={family['components']} rows={family['rows']} "
            f"density={family['density']:.12f} "
            f"groups={family['groups']} flats={family['flat_variables']} "
            f"saturated={len(family['saturated_groups'])}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
