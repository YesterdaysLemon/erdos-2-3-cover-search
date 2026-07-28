#!/usr/bin/env python3
"""Discover full-block pairwise bounds by meet-in-the-middle selection.

The pairwise Bonferroni lower bound for an outside fibre and an anchor subset
is a quadratic set function: exact outside/anchor intersections contribute
positive vertex weights and target-map image indices give phase-uniform
negative pair penalties.  The existing period scanner enumerates every
subset of a supplied anchor list, which becomes expensive beyond 20 anchors.

This discovery tool splits up to a few dozen compatible anchors in half and
uses chunked NumPy matrix products to optimize the same quadratic objective.
The selected subset is then evaluated again with exact ``Fraction``
arithmetic.  A proof need only record that subset and replay its exact
Bonferroni value; optimality of the floating-point discovery step is not a
proof dependency.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

try:
    import numpy as np
except ImportError:  # pragma: no cover - explicit runtime guard
    np = None

from certify_block_plus_overlap_forest import recorded_anchor_rows
from certify_forced_pair_overlap import joint_index
from scan_all_ranked_pairanchor_star import triple_image_index_minors


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def subset_incidence(size: int) -> object:
    masks = np.arange(1 << size, dtype=np.uint64)
    bits = np.arange(size, dtype=np.uint64)
    return ((masks[:, None] >> bits[None, :]) & 1).astype(np.float64)


def half_scores(
    incidence: object,
    weights: object,
    penalties: object,
) -> object:
    linear = incidence @ weights
    quadratic = 0.5 * np.sum(
        (incidence @ penalties) * incidence,
        axis=1,
    )
    return linear - quadratic


def mitm_quadratic_subset(
    weights: object,
    penalties: object,
    *,
    chunk_size: int = 256,
) -> tuple[tuple[int, ...], float]:
    """Maximize positive vertex weights minus pair penalties."""
    if np is None:
        raise RuntimeError("NumPy is required for MITM selection")
    count = len(weights)
    split = count // 2
    first_count = split
    second_count = count - split
    first_incidence = subset_incidence(first_count)
    second_incidence = subset_incidence(second_count)
    first_scores = half_scores(
        first_incidence,
        weights[:split],
        penalties[:split, :split],
    )
    second_scores = half_scores(
        second_incidence,
        weights[split:],
        penalties[split:, split:],
    )
    cross_penalties = penalties[:split, split:]

    best_score = 0.0
    best_first = 0
    best_second = 0
    for start in range(0, len(second_incidence), chunk_size):
        stop = min(start + chunk_size, len(second_incidence))
        second_chunk = second_incidence[start:stop]
        cross = first_incidence @ cross_penalties @ second_chunk.T
        scores = (
            first_scores[:, None]
            + second_scores[None, start:stop]
            - cross
        )
        flat = int(np.argmax(scores))
        first_index, chunk_index = np.unravel_index(flat, scores.shape)
        value = float(scores[first_index, chunk_index])
        if value > best_score:
            best_score = value
            best_first = int(first_index)
            best_second = start + int(chunk_index)

    selected = tuple(
        index
        for index in range(count)
        if (
            index < split
            and (best_first >> index) & 1
            or index >= split
            and (best_second >> (index - split)) & 1
        )
    )
    return selected, best_score


def exact_pairwise_subset_lower(
    outside: dict,
    anchors: tuple[dict, ...],
) -> Fraction:
    outside_h = int(outside["h"])
    value = sum(
        (
            Fraction(1, outside_h * int(anchor["h"]))
            for anchor in anchors
        ),
        Fraction(0),
    )
    value -= sum(
        (
            Fraction(
                triple_image_index_minors((outside, first, second)),
                outside_h * int(first["h"]) * int(second["h"]),
            )
            for first, second in itertools.combinations(anchors, 2)
        ),
        Fraction(0),
    )
    return value


def select_pairwise_anchor_subset(
    outside: dict,
    anchors: list[dict],
    *,
    chunk_size: int = 256,
) -> tuple[tuple[dict, ...], Fraction, float]:
    compatible = tuple(
        anchor
        for anchor in anchors
        if joint_index(outside, anchor) == 1
    )
    outside_h = int(outside["h"])
    weights = np.array(
        [
            1.0 / (outside_h * int(anchor["h"]))
            for anchor in compatible
        ],
        dtype=np.float64,
    )
    penalties = np.zeros(
        (len(compatible), len(compatible)),
        dtype=np.float64,
    )
    for first in range(len(compatible)):
        for second in range(first):
            penalty = (
                triple_image_index_minors(
                    (
                        outside,
                        compatible[first],
                        compatible[second],
                    )
                )
                / (
                    outside_h
                    * int(compatible[first]["h"])
                    * int(compatible[second]["h"])
                )
            )
            penalties[first, second] = penalty
            penalties[second, first] = penalty
    selected_indices, score = mitm_quadratic_subset(
        weights,
        penalties,
        chunk_size=chunk_size,
    )
    selected = tuple(compatible[index] for index in selected_indices)
    exact = exact_pairwise_subset_lower(outside, selected)
    return selected, exact, score


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("block_certificate", type=Path)
    parser.add_argument("period_certificate", type=Path)
    parser.add_argument(
        "--outside-primes",
        help="comma-separated primes; omit to scan every period outside row",
    )
    parser.add_argument("--chunk-size", type=int, default=256)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if np is None:
        raise RuntimeError("NumPy is required for MITM selection")
    if args.chunk_size < 1:
        raise SystemExit("--chunk-size must be positive")

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
    outside_primes = (
        [
            int(item)
            for item in args.outside_primes.split(",")
            if item
        ]
        if args.outside_primes
        else list(baselines)
    )
    missing = set(outside_primes) - rows_by_prime.keys()
    missing |= set(outside_primes) - baselines.keys()
    if missing:
        raise RuntimeError(f"outside primes are missing: {sorted(missing)}")

    results = []
    for prime in outside_primes:
        outside = rows_by_prime[prime]
        selected, exact, discovery_score = select_pairwise_anchor_subset(
            outside,
            anchors,
            chunk_size=args.chunk_size,
        )
        baseline = baselines[prime]
        improvement = max(exact - baseline, 0)
        results.append(
            {
                "outside_prime": prime,
                "compatible_anchor_count": sum(
                    joint_index(outside, anchor) == 1
                    for anchor in anchors
                ),
                "selected_anchor_primes": [
                    int(anchor["p"]) for anchor in selected
                ],
                "selected_anchor_count": len(selected),
                "discovery_float_score": discovery_score,
                "exact_pairwise_lower_bound": fraction_payload(exact),
                "baseline": fraction_payload(baseline),
                "improvement": fraction_payload(improvement),
            }
        )
        print(
            f"outside={prime} selected={len(selected)} "
            f"exact={exact} improvement={improvement}",
            flush=True,
        )

    output = {
        "schema": "mitm_pairwise_conditional_search_v1",
        "pool": str(args.pool),
        "block_certificate": str(args.block_certificate),
        "period_certificate": str(args.period_certificate),
        "outside_primes": outside_primes,
        "chunk_size": args.chunk_size,
        "results": results,
        "status": "discovery_only_selected_subsets_require_proof_replay",
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
