#!/usr/bin/env python3
"""Compress a relaxed finite-radius obstruction into a weighted inequality.

For the points missed by a base phase, every legal one-row retarget induces a
gain mask.  We seek nonnegative point weights for which every gain mask has
weight at most ``B`` while all missed points have total weight ``W``.  If
``radius * B < W``, no union of that many gain masks can cover every point.

This is a deliberately relaxed certificate: masks need not have distinct row
owners and losses at already covered points are ignored.  Consequently, a
successful inequality is a short exact proof of the finite obstruction.
Failure to find one says nothing about the exact repair problem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import exact_greedy
from certify_anchor_phase_quotient import load_cells
from local_phase_cegis import build_targets


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_point_file(path: Path) -> list[tuple[int, int]]:
    payload = json.loads(path.read_text())
    if (
        isinstance(payload, dict)
        and "points" in payload
        and all(
            isinstance(point, list) and len(point) == 2
            for point in payload["points"]
        )
    ):
        return [tuple(map(int, point)) for point in payload["points"]]
    return load_cells(path)


def build_gain_mask_support(
    rows: list[dict],
    candidates: list[tuple],
    points: list[tuple[int, int]],
    initial: dict[int, int],
    fixed_primes: set[int],
    np,
) -> tuple[
    list[tuple[int, int]],
    dict[int, list[tuple[int, int]]],
]:
    assignment = []
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        target = int(initial.get(prime, residue)) % h
        if h % modulus or target % modulus != residue:
            raise RuntimeError(f"invalid initial phase for p={prime}")
        assignment.append(target)

    targets, _seconds = build_targets(points, candidates, np)
    assignment_array = np.asarray(assignment, dtype=targets.dtype)
    base_cover = np.count_nonzero(
        targets == assignment_array,
        axis=1,
    )
    miss_indices = [
        int(index) for index in np.flatnonzero(base_cover == 0)
    ]
    support: dict[int, list[tuple[int, int]]] = {}
    for row_index, row in enumerate(rows):
        if int(row["p"]) in fixed_primes:
            continue
        current = int(assignment[row_index])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        by_target: dict[int, int] = {}
        for bit, point_index in enumerate(miss_indices):
            target = int(targets[point_index, row_index])
            if target == current or target % modulus != residue:
                continue
            by_target[target] = (
                by_target.get(target, 0) | (1 << bit)
            )
        for target, mask in by_target.items():
            if mask:
                support.setdefault(mask, []).append((row_index, target))
    return [points[index] for index in miss_indices], support


def build_gain_masks(
    rows: list[dict],
    candidates: list[tuple],
    points: list[tuple[int, int]],
    initial: dict[int, int],
    fixed_primes: set[int],
    np,
) -> tuple[list[tuple[int, int]], list[int]]:
    missed_points, support = build_gain_mask_support(
        rows,
        candidates,
        points,
        initial,
        fixed_primes,
        np,
    )
    return missed_points, sorted(support)


def solve_fractional_weights(
    masks: list[int],
    point_count: int,
    linprog,
) -> list[float]:
    # Variables are point weights followed by the maximum mask weight t.
    objective = [0.0] * point_count + [1.0]
    inequalities = []
    for mask in masks:
        row = [
            1.0 if mask & (1 << bit) else 0.0
            for bit in range(point_count)
        ]
        inequalities.append(row + [-1.0])
    result = linprog(
        objective,
        A_ub=inequalities,
        b_ub=[0.0] * len(inequalities),
        A_eq=[[1.0] * point_count + [0.0]],
        b_eq=[1.0],
        bounds=[(0.0, None)] * (point_count + 1),
        method="highs",
    )
    if not result.success:
        raise RuntimeError(f"dual LP failed: {result.message}")
    return [float(value) for value in result.x[:point_count]]


def integerize_weights(
    raw_weights: list[float],
    max_denominator: int,
) -> list[int]:
    fractions = [
        Fraction(max(0.0, value)).limit_denominator(max_denominator)
        for value in raw_weights
    ]
    denominator = math.lcm(
        *(value.denominator for value in fractions)
    )
    weights = [
        value.numerator * (denominator // value.denominator)
        for value in fractions
    ]
    common = math.gcd(*weights)
    return [weight // common for weight in weights]


def summarize_exact_certificate(
    masks: list[int],
    weights: list[int],
    radius: int,
    point_count: int,
) -> dict:
    mask_weights = [
        sum(
            weight
            for bit, weight in enumerate(weights)
            if mask & (1 << bit)
        )
        for mask in masks
    ]
    total = sum(weights)
    maximum = max(mask_weights, default=0)
    positive_mask = sum(
        1 << bit
        for bit, weight in enumerate(weights)
        if weight > 0
    )
    tight_masks = [
        mask for mask, weight in zip(masks, mask_weights) if weight == maximum
    ]
    tight_positive_masks = sorted(
        {mask & positive_mask for mask in tight_masks}
    )

    def tight_partition_exists(
        search_limit: int = 2_000_000,
    ) -> tuple[bool | None, int]:
        """Search exactly for the only cover possible at dual equality.

        If ``radius * maximum == total``, a cover by at most ``radius`` masks
        must use exactly ``radius`` maximum-weight masks.  Those masks must be
        pairwise disjoint on every positive-weight point.  Otherwise overlap
        would make their union weight strictly smaller than ``total``.

        ``None`` means the deliberately bounded discovery search hit its
        state limit; only a completed negative search is certifying.
        """

        if (
            radius < 1
            or maximum < 1
            or radius * maximum != total
            or not positive_mask
        ):
            return None, 0
        containing: dict[int, list[int]] = {
            bit: [
                mask
                for mask in tight_positive_masks
                if mask & (1 << bit)
            ]
            for bit in range(point_count)
            if positive_mask & (1 << bit)
        }
        visited: set[tuple[int, int]] = set()
        nodes = 0

        def search(union: int, used: int) -> bool | None:
            nonlocal nodes
            state = (union, used)
            if state in visited:
                return False
            visited.add(state)
            nodes += 1
            if nodes > search_limit:
                return None
            if union == positive_mask:
                return used == radius
            if used >= radius:
                return False
            uncovered = positive_mask & ~union
            bit = (uncovered & -uncovered).bit_length() - 1
            for mask in containing[bit]:
                if mask & union:
                    continue
                result = search(union | mask, used + 1)
                if result is not False:
                    return result
            return False

        return search(0, 0), nodes

    tight_cover_exists, tight_search_nodes = tight_partition_exists()
    tight_equality_case = (
        radius > 0
        and maximum > 0
        and radius * maximum == total
    )
    tight_equality_certified = (
        tight_equality_case and tight_cover_exists is False
    )
    maximal_masks = [
        mask
        for mask in masks
        if not any(
            mask != other and mask & ~other == 0
            for other in masks
        )
    ]
    pair_union_sizes = [
        (left | right).bit_count()
        for index, left in enumerate(maximal_masks)
        for right in maximal_masks[index:]
    ]
    maximum_pair_union = max(pair_union_sizes, default=0)
    return {
        "integer_weights": weights,
        "total_point_weight": total,
        "maximum_single_mask_weight": maximum,
        "radius_times_maximum": radius * maximum,
        "strict_weight_gap": total - radius * maximum,
        "certified": radius * maximum < total,
        "positive_weight_point_count": positive_mask.bit_count(),
        "tight_mask_count": len(tight_masks),
        "distinct_tight_positive_mask_count": len(
            tight_positive_masks
        ),
        "tight_positive_masks": [
            [
                bit
                for bit in range(point_count)
                if mask & (1 << bit)
            ]
            for mask in tight_positive_masks
        ],
        "tight_equality_case": tight_equality_case,
        "tight_disjoint_cover_exists": tight_cover_exists,
        "tight_disjoint_search_complete": (
            tight_cover_exists is not None
        ),
        "tight_disjoint_search_nodes": tight_search_nodes,
        "tight_equality_certified": tight_equality_certified,
        "pairwise_union_certified": (
            radius == 2 and maximum_pair_union < point_count
        ),
        "gain_mask_count": len(masks),
        "maximum_mask_cardinality": max(
            (mask.bit_count() for mask in masks),
            default=0,
        ),
        "mask_cardinality_histogram": {
            str(size): count
            for size, count in sorted(
                Counter(mask.bit_count() for mask in masks).items()
            )
        },
        "inclusion_maximal_mask_count": len(maximal_masks),
        "inclusion_maximal_masks": [
            [
                bit
                for bit in range(point_count)
                if mask & (1 << bit)
            ]
            for mask in maximal_masks
        ],
        "maximum_pair_union_cardinality": maximum_pair_union,
        "pair_union_cardinality_histogram": {
            str(size): count
            for size, count in sorted(Counter(pair_union_sizes).items())
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path, nargs="+")
    parser.add_argument("--base-phase", type=Path, required=True)
    parser.add_argument("--fixed-primes", default="")
    parser.add_argument("--radius", type=int, required=True)
    parser.add_argument("--max-denominator", type=int, default=10_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.radius < 0:
        raise SystemExit("--radius must be nonnegative")
    if args.max_denominator < 1:
        raise SystemExit("--max-denominator must be positive")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from scipy.optimize import linprog  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    candidates = exact_greedy.load_candidates(args.pool, True)
    points = []
    seen = set()
    for path in args.points:
        for point in load_point_file(path):
            if point not in seen:
                seen.add(point)
                points.append(point)
    initial = {
        int(prime): int(target)
        for prime, target in json.loads(
            args.base_phase.read_text()
        ).items()
    }
    fixed_primes = {
        int(value) for value in args.fixed_primes.split(",") if value
    }
    missed_points, masks = build_gain_masks(
        rows,
        candidates,
        points,
        initial,
        fixed_primes,
        np,
    )
    if not missed_points:
        raise RuntimeError("base phase misses no supplied point")
    raw_weights = solve_fractional_weights(
        masks,
        len(missed_points),
        linprog,
    )
    integer_weights = integerize_weights(
        raw_weights,
        args.max_denominator,
    )
    summary = summarize_exact_certificate(
        masks,
        integer_weights,
        args.radius,
        len(missed_points),
    )
    result = {
        "problem": "relaxed finite-radius gain-mask dual",
        "scope": (
            "the embedded pool, base phase, fixed rows, supplied points, "
            "and declared radius; every mask is a legal one-row gain pattern"
        ),
        "pool": str(args.pool),
        "pool_sha256": sha256_file(args.pool),
        "base_phase": str(args.base_phase),
        "base_phase_sha256": sha256_file(args.base_phase),
        "point_sources": [
            {"path": str(path), "sha256": sha256_file(path)}
            for path in args.points
        ],
        "fixed_primes": sorted(fixed_primes),
        "radius": args.radius,
        "missed_points": [[k, ell] for k, ell in missed_points],
        "positive_weight_points": [
            [k, ell]
            for (k, ell), weight in zip(
                missed_points,
                integer_weights,
                strict=True,
            )
            if weight > 0
        ],
        **summary,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"points={len(missed_points)} masks={len(masks)} "
        f"weights={integer_weights} "
        f"max={summary['maximum_single_mask_weight']} "
        f"total={summary['total_point_weight']} "
        f"weighted={summary['certified']} "
        f"tight_equality={summary['tight_equality_certified']} "
        f"pairwise={summary['pairwise_union_certified']} "
        f"output={args.output}",
        flush=True,
    )
    return 0 if (
        summary["certified"]
        or summary["tight_equality_certified"]
        or summary["pairwise_union_certified"]
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
