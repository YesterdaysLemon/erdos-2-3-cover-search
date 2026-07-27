#!/usr/bin/env python3
"""Discover projected conditional-overlap certificate designs.

The search chooses a normalizing anchor, a small exact base, and a set of
anchors whose moduli split into a base projection and pairwise-coprime
residuals.  It evaluates the same exact finite tensor represented by
``certify_projected_conditional_fibre_overlap.py``.

This is a discovery tool.  Any reported winner must still be emitted by the
dedicated certificate generator and replayed by its independent verifier
before it can be used in a proof.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

try:
    import numpy as np
except ImportError:  # pragma: no cover - explicit runtime guard
    np = None

from certify_block_plus_overlap_forest import recorded_anchor_rows
from certify_forced_pair_overlap import joint_index
from certify_projected_conditional_fibre_overlap import (
    primitive_fibre_points,
)


INT64_SAFE_SCORE = (1 << 62) - 1


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


@dataclass(frozen=True)
class Projected:
    prime: int
    projection: int
    residual: int


@dataclass(frozen=True)
class Design:
    normalizer: int
    base_primes: tuple[int, ...]
    projected: tuple[Projected, ...]
    base_period: int
    base_points: int
    tensor_cells: int
    density_score: Fraction

    @property
    def anchor_primes(self) -> tuple[int, ...]:
        return (
            self.normalizer,
            *self.base_primes,
            *(record.prime for record in self.projected),
        )


def deterministic_jitter(prime: int, seed: int) -> int:
    value = (
        prime * 1_103_515_245
        + seed * 12_345
        + 0x9E3779B9
    ) & 0xFFFFFFFF
    value ^= value >> 16
    return value


def projectable(row: dict, period: int) -> Projected | None:
    modulus = int(row["h"])
    projection = math.gcd(modulus, period)
    residual = modulus // projection
    if (
        math.gcd(period // projection, residual) != 1
        or math.gcd(int(row["a"]), int(row["b"]), residual) != 1
    ):
        return None
    return Projected(int(row["p"]), projection, residual)


def greedy_projected_set(
    candidates: list[Projected],
    rows_by_prime: dict[int, dict],
    initial_cells: int,
    max_tensor_cells: int,
    order_kind: str,
    seed: int,
) -> tuple[tuple[Projected, ...], int]:
    def key(record: Projected) -> tuple:
        modulus = int(rows_by_prime[record.prime]["h"])
        if order_kind == "density":
            primary = Fraction(1, modulus)
        elif order_kind == "projection":
            primary = Fraction(1, max(1, record.projection))
        elif order_kind == "ratio":
            primary = Fraction(
                1,
                modulus * max(1, record.projection),
            )
        elif order_kind == "count":
            primary = Fraction(1)
        else:
            primary = Fraction(
                1_000_000 + deterministic_jitter(record.prime, seed),
                1_000_000 * modulus * max(1, record.projection),
            )
        return primary, -record.projection, -modulus, -record.prime

    selected = []
    residuals = []
    cells = initial_cells
    for record in sorted(candidates, key=key, reverse=True):
        if any(
            math.gcd(record.residual, residual) != 1
            for residual in residuals
        ):
            continue
        if cells * record.projection > max_tensor_cells:
            continue
        selected.append(record)
        residuals.append(record.residual)
        cells *= record.projection
    return tuple(sorted(selected, key=lambda item: item.prime)), cells


def candidate_designs(
    outside: dict,
    anchors: list[dict],
    *,
    max_base_anchors: int,
    max_base_period: int,
    max_base_points: int,
    max_tensor_cells: int,
    jitter_seeds: int,
) -> list[Design]:
    rows_by_prime = {int(row["p"]): row for row in anchors}
    designs: dict[
        tuple[int, tuple[int, ...], tuple[int, ...], int],
        Design,
    ] = {}
    for normalizer in anchors:
        if joint_index(outside, normalizer) != 1:
            continue
        normalizer_prime = int(normalizer["p"])
        remaining = [
            row for row in anchors
            if int(row["p"]) != normalizer_prime
        ]
        for base_count in range(max_base_anchors + 1):
            for base_rows in itertools.combinations(
                remaining,
                base_count,
            ):
                base_primes = tuple(
                    sorted(int(row["p"]) for row in base_rows)
                )
                period = math.lcm(
                    int(outside["h"]),
                    int(normalizer["h"]),
                    *(int(row["h"]) for row in base_rows),
                )
                if period > max_base_period:
                    continue
                base_points = (
                    period * period // int(outside["h"])
                )
                if base_points > max_base_points:
                    continue
                initial_cells = math.prod(
                    int(row["h"]) for row in base_rows
                )
                if initial_cells > max_tensor_cells:
                    continue
                base_set = set(base_primes)
                projected_candidates = [
                    record
                    for row in remaining
                    if int(row["p"]) not in base_set
                    for record in [projectable(row, period)]
                    if record is not None
                ]
                orderings = [
                    ("density", 0),
                    ("projection", 0),
                    ("ratio", 0),
                    ("count", 0),
                    *(
                        ("jitter", seed)
                        for seed in range(jitter_seeds)
                    ),
                ]
                for order_kind, seed in orderings:
                    projected, cells = greedy_projected_set(
                        projected_candidates,
                        rows_by_prime,
                        initial_cells,
                        max_tensor_cells,
                        order_kind,
                        seed,
                    )
                    density_score = sum(
                        (
                            Fraction(
                                1,
                                int(rows_by_prime[prime]["h"]),
                            )
                            for prime in (
                                normalizer_prime,
                                *base_primes,
                                *(item.prime for item in projected),
                            )
                        ),
                        Fraction(0),
                    )
                    design = Design(
                        normalizer=normalizer_prime,
                        base_primes=base_primes,
                        projected=projected,
                        base_period=period,
                        base_points=base_points,
                        tensor_cells=cells,
                        density_score=density_score,
                    )
                    signature = (
                        normalizer_prime,
                        base_primes,
                        tuple(item.prime for item in projected),
                        period,
                    )
                    designs[signature] = design
    return list(designs.values())


def select_designs(
    designs: list[Design],
    max_designs: int,
) -> list[Design]:
    def rank(design: Design) -> tuple:
        return (
            design.density_score,
            len(design.anchor_primes),
            -design.tensor_cells,
            -design.base_points,
            -design.base_period,
        )

    selected = sorted(designs, key=rank, reverse=True)[:max_designs]
    by_normalizer: dict[int, Design] = {}
    for design in designs:
        current = by_normalizer.get(design.normalizer)
        if current is None or rank(design) > rank(current):
            by_normalizer[design.normalizer] = design
    selected.extend(by_normalizer.values())
    unique = {}
    for design in selected:
        signature = (
            design.normalizer,
            design.base_primes,
            tuple(item.prime for item in design.projected),
            design.base_period,
        )
        unique[signature] = design
    return sorted(unique.values(), key=rank, reverse=True)


def evaluate_design(
    outside: dict,
    design: Design,
    rows_by_prime: dict[int, dict],
) -> tuple[Fraction, dict] | None:
    normalizer = rows_by_prime[design.normalizer]
    base_rows = [
        rows_by_prime[prime] for prime in design.base_primes
    ]
    projected_rows = [
        (
            rows_by_prime[record.prime],
            record.projection,
            record.residual,
        )
        for record in design.projected
    ]
    shape = (
        *(int(row["h"]) for row in base_rows),
        *(projection for _row, projection, _residual
          in projected_rows),
    )
    counts = np.zeros(shape, dtype=np.int64)
    outside_points = 0
    normalizer_covered = 0
    for k, l in primitive_fibre_points(outside, design.base_period):
        outside_points += 1
        if (
            int(normalizer["a"]) * k
            + int(normalizer["b"]) * l
        ) % int(normalizer["h"]) == 0:
            normalizer_covered += 1
            continue
        index = (
            *(
                (
                    int(row["a"]) * k + int(row["b"]) * l
                ) % int(row["h"])
                for row in base_rows
            ),
            *(
                (
                    int(row["a"]) * k + int(row["b"]) * l
                ) % projection
                for row, projection, _residual in projected_rows
            ),
        )
        counts[index] += 1

    residual_denominator = math.prod(
        residual for _row, _projection, residual in projected_rows
    )
    if outside_points * residual_denominator > INT64_SAFE_SCORE:
        return None
    transformed = counts
    for axis, (_row, projection, residual) in reversed(
        list(enumerate(projected_rows, start=len(base_rows)))
    ):
        matrix = np.full(
            (projection, projection),
            residual,
            dtype=np.int64,
        )
        np.fill_diagonal(matrix, residual - 1)
        transformed = np.tensordot(
            matrix,
            transformed,
            axes=(1, axis),
        )
        transformed = np.moveaxis(transformed, 0, axis)
    for axis in range(len(base_rows)):
        transformed = (
            np.expand_dims(transformed.sum(axis=axis), axis)
            - transformed
        )
    maximum_uncovered_score = int(transformed.max())
    intersection = Fraction(
        outside_points * residual_denominator
        - maximum_uncovered_score,
        residual_denominator
        * design.base_period
        * design.base_period,
    )
    detail = {
        "outside_base_points": outside_points,
        "normalizer_covered_base_points": normalizer_covered,
        "residual_denominator": residual_denominator,
        "maximum_uncovered_score": maximum_uncovered_score,
    }
    return intersection, detail


def design_payload(
    design: Design,
    value: Fraction,
    baseline: Fraction,
    detail: dict,
) -> dict:
    return {
        "normalizer_anchor_prime": design.normalizer,
        "base_anchor_primes": list(design.base_primes),
        "projected_anchors": [
            {
                "prime": record.prime,
                "projection_modulus": record.projection,
                "residual_modulus": record.residual,
            }
            for record in design.projected
        ],
        "anchor_primes": list(design.anchor_primes),
        "base_period": design.base_period,
        "base_points": design.base_points,
        "tensor_cells": design.tensor_cells,
        "density_score": fraction_payload(design.density_score),
        "intersection": fraction_payload(value),
        "baseline": fraction_payload(baseline),
        "improvement": fraction_payload(value - baseline),
        **detail,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("block_certificate", type=Path)
    parser.add_argument("period_certificate", type=Path)
    parser.add_argument("--outside-primes", required=True)
    parser.add_argument("--max-base-anchors", type=int, default=1)
    parser.add_argument("--max-base-period", type=int, default=200_000)
    parser.add_argument("--max-base-points", type=int, default=500_000)
    parser.add_argument("--max-tensor-cells", type=int, default=3_000_000)
    parser.add_argument("--jitter-seeds", type=int, default=8)
    parser.add_argument("--max-designs", type=int, default=20)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if np is None:
        raise RuntimeError("NumPy is required for projected design search")

    source = json.loads(args.pool.read_text())
    block = json.loads(args.block_certificate.read_text())
    period = json.loads(args.period_certificate.read_text())
    rows = source["choices"]
    rows_by_prime = {int(row["p"]): row for row in rows}
    anchors = recorded_anchor_rows(block)
    outside_primes = [
        int(value) for value in args.outside_primes.split(",") if value
    ]
    baseline_by_prime = {
        int(record["outside_prime"]): read_fraction(
            record["used_intersection_lower_bound"]
        )
        for record in period["outside_row_lowers"]
    }
    missing = (
        set(outside_primes)
        - rows_by_prime.keys()
        | set(outside_primes)
        - baseline_by_prime.keys()
    )
    if missing:
        raise RuntimeError(
            f"outside primes absent from source or period: {sorted(missing)}"
        )

    results = []
    for outside_prime in outside_primes:
        outside = rows_by_prime[outside_prime]
        baseline = baseline_by_prime[outside_prime]
        designs = candidate_designs(
            outside,
            anchors,
            max_base_anchors=args.max_base_anchors,
            max_base_period=args.max_base_period,
            max_base_points=args.max_base_points,
            max_tensor_cells=args.max_tensor_cells,
            jitter_seeds=args.jitter_seeds,
        )
        selected = select_designs(designs, args.max_designs)
        best_value = baseline
        best_record = None
        evaluated = 0
        for design in selected:
            evaluated_result = evaluate_design(
                outside,
                design,
                rows_by_prime,
            )
            if evaluated_result is None:
                continue
            evaluated += 1
            value, detail = evaluated_result
            if value <= best_value:
                continue
            best_value = value
            best_record = design_payload(
                design,
                value,
                baseline,
                detail,
            )
            print(
                f"outside={outside_prime} intersection={value} "
                f"improvement={value - baseline} "
                f"normalizer={design.normalizer} "
                f"anchors={len(design.anchor_primes)}",
                flush=True,
            )
        results.append(
            {
                "outside_prime": outside_prime,
                "candidate_design_count": len(designs),
                "evaluated_design_count": evaluated,
                "baseline": fraction_payload(baseline),
                "best_design": best_record,
            }
        )
        print(
            f"outside={outside_prime} candidates={len(designs)} "
            f"evaluated={evaluated} improved={best_record is not None}",
            flush=True,
        )

    output = {
        "schema": "projected_conditional_design_search_v1",
        "pool": str(args.pool),
        "block_certificate": str(args.block_certificate),
        "period_certificate": str(args.period_certificate),
        "outside_primes": outside_primes,
        "max_base_anchors": args.max_base_anchors,
        "max_base_period": args.max_base_period,
        "max_base_points": args.max_base_points,
        "max_tensor_cells": args.max_tensor_cells,
        "jitter_seeds": args.jitter_seeds,
        "max_designs": args.max_designs,
        "results": results,
        "status": "discovery_only_requires_certificate_and_replay",
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
