#!/usr/bin/env python3
"""Discover sixth-order conditional fibre-to-block overlap bounds.

For a fixed outside fibre, every compatible anchor has an exact pair
intersection.  For larger intersections, the target-map image index gives
an upper bound that is valid for every choice of fibre targets; when the
index is one, the intersection is target-independent and exact.

Applying the Bonferroni lower bound through order six to a small anchor
subset therefore gives a phase-uniform lower bound for the outside fibre's
intersection with the anchor union.  This is a discovery tool.  Any useful
winner still needs a dedicated certificate and independent replay before it
can enter a period proof.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

from certify_block_plus_overlap_forest import recorded_anchor_rows
from certify_forced_pair_overlap import joint_index
from scan_all_ranked_pairanchor_star import determinant_bareiss


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def target_image_index(rows: tuple[dict, ...]) -> int:
    """Index of the two-coordinate target map into cyclic row factors."""
    size = len(rows)
    moduli = tuple(int(row["h"]) for row in rows)
    columns = (
        *(
            tuple(moduli[index] if row == index else 0 for row in range(size))
            for index in range(size)
        ),
        tuple(int(row["a"]) for row in rows),
        tuple(int(row["b"]) for row in rows),
    )
    minors = (
        abs(
            determinant_bareiss(
                tuple(
                    tuple(columns[column][row] for column in selected)
                    for row in range(size)
                )
            )
        )
        for selected in itertools.combinations(range(size + 2), size)
    )
    return math.gcd(*minors)


def uniform_intersection_term(
    outside: dict,
    anchors: tuple[dict, ...],
    *,
    positive: bool,
) -> Fraction:
    """Return a target-uniform lower or upper intersection term."""
    rows = (outside, *anchors)
    index = target_image_index(rows)
    denominator = math.prod(int(row["h"]) for row in rows)
    if positive:
        return Fraction(1, denominator) if index == 1 else Fraction(0)
    return Fraction(index, denominator)


def sixterm_subset_lower(
    outside: dict,
    anchors: tuple[dict, ...],
) -> Fraction:
    """Bonferroni lower bound through the largest even order at most six."""
    if not anchors:
        return Fraction(0)
    cutoff = min(6, len(anchors))
    if cutoff % 2:
        cutoff -= 1
    if cutoff == 0:
        cutoff = 1
    value = Fraction(0)
    for order in range(1, cutoff + 1):
        positive = order % 2 == 1
        terms = (
            uniform_intersection_term(
                outside,
                selected,
                positive=positive,
            )
            for selected in itertools.combinations(anchors, order)
        )
        total = sum(terms, Fraction(0))
        value += total if positive else -total
    return value


def best_sixterm_lower(
    outside: dict,
    anchors: list[dict],
    anchor_limit: int,
) -> tuple[Fraction, tuple[int, ...], int]:
    """Maximize the 1-, 2-, 4-, and 6-anchor Bonferroni bounds."""
    compatible = sorted(
        (
            anchor
            for anchor in anchors
            if joint_index(outside, anchor) == 1
        ),
        key=lambda row: (int(row["h"]), int(row["p"])),
    )[:anchor_limit]
    terms: dict[tuple[int, ...], Fraction] = {}
    for size in range(1, min(6, len(compatible)) + 1):
        positive = size % 2 == 1
        for selected_indices in itertools.combinations(
            range(len(compatible)),
            size,
        ):
            selected = tuple(
                compatible[index] for index in selected_indices
            )
            terms[selected_indices] = uniform_intersection_term(
                outside,
                selected,
                positive=positive,
            )

    best = Fraction(0)
    best_primes: tuple[int, ...] = ()
    evaluated = 0
    for size in (1, 2, 4, 6):
        if size > len(compatible):
            continue
        for selected_indices in itertools.combinations(
            range(len(compatible)),
            size,
        ):
            evaluated += 1
            value = Fraction(0)
            for order in range(1, size + 1):
                subtotal = sum(
                    (
                        terms[subset]
                        for subset in itertools.combinations(
                            selected_indices,
                            order,
                        )
                    ),
                    Fraction(0),
                )
                value += subtotal if order % 2 else -subtotal
            if value > best:
                best = value
                best_primes = tuple(
                    int(compatible[index]["p"])
                    for index in selected_indices
                )
    return best, best_primes, evaluated


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("block_certificate", type=Path)
    parser.add_argument("period_certificate", type=Path)
    parser.add_argument("--outside-primes", required=True)
    parser.add_argument("--anchor-limit", type=int, default=14)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.anchor_limit < 1:
        raise SystemExit("--anchor-limit must be positive")

    source = json.loads(args.pool.read_text())
    block = json.loads(args.block_certificate.read_text())
    period = json.loads(args.period_certificate.read_text())
    rows_by_prime = {
        int(row["p"]): row for row in source["choices"]
    }
    anchors = recorded_anchor_rows(block)
    baselines = {
        int(record["outside_prime"]): read_fraction(
            record["used_intersection_lower_bound"]
        )
        for record in period["outside_row_lowers"]
    }
    outside_primes = tuple(
        int(item) for item in args.outside_primes.split(",") if item
    )
    missing = set(outside_primes) - rows_by_prime.keys()
    missing |= set(outside_primes) - baselines.keys()
    if missing:
        raise RuntimeError(f"outside primes are missing: {sorted(missing)}")

    results = []
    for prime in outside_primes:
        outside = rows_by_prime[prime]
        baseline = baselines[prime]
        value, selected, evaluated = best_sixterm_lower(
            outside,
            anchors,
            args.anchor_limit,
        )
        results.append(
            {
                "outside_prime": prime,
                "anchor_limit": args.anchor_limit,
                "subsets_evaluated": evaluated,
                "baseline": fraction_payload(baseline),
                "sixterm_lower_bound": fraction_payload(value),
                "improvement": fraction_payload(max(value - baseline, 0)),
                "selected_anchor_primes": list(selected),
            }
        )
        print(
            f"outside={prime} value={value} baseline={baseline} "
            f"improvement={max(value - baseline, 0)} "
            f"anchors={selected}",
            flush=True,
        )

    output = {
        "schema": "sixterm_conditional_bound_search_v1",
        "pool": str(args.pool),
        "block_certificate": str(args.block_certificate),
        "period_certificate": str(args.period_certificate),
        "outside_primes": list(outside_primes),
        "anchor_limit": args.anchor_limit,
        "results": results,
        "status": "discovery_only_requires_certificate_and_replay",
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
