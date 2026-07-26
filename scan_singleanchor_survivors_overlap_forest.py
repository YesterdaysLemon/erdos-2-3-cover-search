#!/usr/bin/env python3
"""Apply exact maximum-weight overlap forests to remaining ranked families."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from certify_block_plus_overlap_forest import DisjointSet
from certify_forced_pair_overlap import joint_index


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("prior_scan", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    prior = json.loads(args.prior_scan.read_text())
    rows = source["choices"]
    moduli = [int(row["h"]) for row in rows]
    primes = [int(row["p"]) for row in rows]
    edges = sorted(
        (
            (
                moduli[first] * moduli[second],
                min(primes[first], primes[second]),
                max(primes[first], primes[second]),
                first,
                second,
            )
            for first in range(len(rows))
            for second in range(first)
            if joint_index(rows[first], rows[second]) == 1
        ),
        key=lambda edge: (edge[0], edge[1], edge[2]),
    )

    results = []
    proved = 0
    for family in prior["unresolved"]:
        period = int(family["period"])
        selected = [
            index
            for index, modulus in enumerate(moduli)
            if period % modulus == 0
        ]
        selected_set = set(selected)
        total_density = Fraction(
            sum(period // moduli[index] for index in selected),
            period,
        )
        dsu = DisjointSet(selected)
        overlap = Fraction(0)
        selected_edges = []
        no_cover = False
        for product, _first_prime, _second_prime, first, second in edges:
            if first not in selected_set or second not in selected_set:
                continue
            if not dsu.union(first, second):
                continue
            overlap += Fraction(1, product)
            selected_edges.append([primes[first], primes[second]])
            if total_density - overlap < 1:
                no_cover = True
                break
        upper = total_density - overlap
        proved += int(no_cover)
        results.append(
            {
                "period": period,
                "components": family["components"],
                "rows": len(selected),
                "total_density": fraction_payload(total_density),
                "selected_forest_edges": (
                    selected_edges if no_cover else []
                ),
                "forest_overlap_sum": fraction_payload(overlap),
                "union_upper_bound": fraction_payload(upper),
                "proved_no_cover": no_cover,
            }
        )

    unresolved = [
        result for result in results if not result["proved_no_cover"]
    ]
    unresolved.sort(
        key=lambda item: (
            item["union_upper_bound"]["decimal"],
            item["period"],
        )
    )
    output = {
        "pool": str(args.pool),
        "prior_scan": str(args.prior_scan),
        "families_checked": len(results),
        "proved_no_cover": proved,
        "unresolved_count": len(unresolved),
        "unresolved": unresolved,
        "results": results,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"families={len(results)} proved={proved} "
        f"unresolved={len(unresolved)}",
        flush=True,
    )
    for result in unresolved[:20]:
        print(
            f"period={result['period']} rows={result['rows']} "
            f"upper={result['union_upper_bound']['decimal']:.12f}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
