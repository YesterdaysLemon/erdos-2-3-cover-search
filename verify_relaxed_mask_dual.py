#!/usr/bin/env python3
"""Independently scalar-replay a relaxed gain-mask dual certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reconstruct_gain_masks(
    rows: list[dict],
    points: list[tuple[int, int]],
    phases: dict[int, int],
    fixed_primes: set[int],
) -> list[int]:
    base_targets = {}
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        target = int(phases.get(prime, residue)) % h
        if h % modulus or target % modulus != residue:
            raise RuntimeError(f"invalid base phase for p={prime}")
        base_targets[prime] = target

    for point_index, (k, ell) in enumerate(points):
        for row in rows:
            prime = int(row["p"])
            if (
                int(row["a"]) * k
                + int(row["b"]) * ell
                - base_targets[prime]
            ) % int(row["h"]) == 0:
                raise RuntimeError(
                    f"point {point_index} is covered at the base phase "
                    f"by p={prime}"
                )

    masks = set()
    for row in rows:
        prime = int(row["p"])
        if prime in fixed_primes:
            continue
        h = int(row["h"])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        current = base_targets[prime]
        by_target: dict[int, int] = {}
        for bit, (k, ell) in enumerate(points):
            target = (
                int(row["a"]) * k + int(row["b"]) * ell
            ) % h
            if target == current or target % modulus != residue:
                continue
            by_target[target] = (
                by_target.get(target, 0) | (1 << bit)
            )
        masks.update(mask for mask in by_target.values() if mask)
    return sorted(masks)


def inclusion_maximal_masks(masks: list[int]) -> list[int]:
    return [
        mask
        for mask in masks
        if not any(
            mask != other and mask & ~other == 0
            for other in masks
        )
    ]


def tight_partition_exists_dp(
    tight_masks: list[int],
    positive_mask: int,
    radius: int,
) -> bool:
    """Independently enumerate disjoint tight-mask unions."""

    if radius < 1 or not positive_mask:
        return False
    reachable = {0}
    for _used in range(radius):
        next_reachable = set()
        for union in reachable:
            for mask in tight_masks:
                if mask & union:
                    continue
                next_reachable.add(union | mask)
        reachable = next_reachable
        if not reachable:
            return False
    return positive_mask in reachable


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    pool_path = args.certificate.parent / certificate["pool"]
    phase_path = args.certificate.parent / certificate["base_phase"]
    if sha256_file(pool_path) != certificate["pool_sha256"]:
        raise RuntimeError("pool SHA-256 mismatch")
    if sha256_file(phase_path) != certificate["base_phase_sha256"]:
        raise RuntimeError("base-phase SHA-256 mismatch")
    rows = json.loads(pool_path.read_text())["choices"]
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(phase_path.read_text()).items()
    }
    points = [
        tuple(map(int, point))
        for point in certificate["missed_points"]
    ]
    fixed_primes = {
        int(prime) for prime in certificate["fixed_primes"]
    }
    masks = reconstruct_gain_masks(
        rows,
        points,
        phases,
        fixed_primes,
    )
    maximal = inclusion_maximal_masks(masks)
    encoded_maximal = [
        sum(1 << int(bit) for bit in bits)
        for bits in certificate["inclusion_maximal_masks"]
    ]
    maximum_pair_union = max(
        (
            (left | right).bit_count()
            for index, left in enumerate(maximal)
            for right in maximal[index:]
        ),
        default=0,
    )
    pair_histogram = Counter(
        (left | right).bit_count()
        for index, left in enumerate(maximal)
        for right in maximal[index:]
    )
    weights = [
        int(weight) for weight in certificate["integer_weights"]
    ]
    mask_weights = [
        sum(
            weight
            for bit, weight in enumerate(weights)
            if mask & (1 << bit)
        )
        for mask in masks
    ]
    total_weight = sum(weights)
    maximum_mask_weight = max(mask_weights, default=0)
    radius = int(certificate["radius"])
    positive_mask = sum(
        1 << bit
        for bit, weight in enumerate(weights)
        if weight > 0
    )
    tight_masks = [
        mask
        for mask, weight in zip(masks, mask_weights)
        if weight == maximum_mask_weight
    ]
    tight_positive_masks = sorted(
        {mask & positive_mask for mask in tight_masks}
    )
    tight_equality_case = (
        radius > 0
        and maximum_mask_weight > 0
        and radius * maximum_mask_weight == total_weight
    )
    tight_cover_exists = (
        tight_partition_exists_dp(
            tight_positive_masks,
            positive_mask,
            radius,
        )
        if tight_equality_case
        else None
    )
    tight_equality_certified = (
        tight_equality_case and tight_cover_exists is False
    )
    tight_fields_present = "tight_equality_certified" in certificate
    if tight_fields_present:
        encoded_tight_masks = [
            sum(1 << int(bit) for bit in bits)
            for bits in certificate["tight_positive_masks"]
        ]
        tight_metadata_matches = (
            positive_mask.bit_count()
            == int(certificate["positive_weight_point_count"])
            and len(tight_masks)
            == int(certificate["tight_mask_count"])
            and len(tight_positive_masks)
            == int(certificate["distinct_tight_positive_mask_count"])
            and tight_positive_masks == encoded_tight_masks
            and tight_equality_case
            == bool(certificate["tight_equality_case"])
            and tight_cover_exists
            == certificate["tight_disjoint_cover_exists"]
            and tight_equality_certified
            == bool(certificate["tight_equality_certified"])
        )
    else:
        tight_metadata_matches = True
    checks = {
        "gain_mask_count": len(masks)
        == int(certificate["gain_mask_count"]),
        "maximal_masks": maximal == encoded_maximal,
        "maximal_mask_count": len(maximal)
        == int(certificate["inclusion_maximal_mask_count"]),
        "maximum_pair_union": maximum_pair_union
        == int(certificate["maximum_pair_union_cardinality"]),
        "pair_union_histogram": {
            str(size): count
            for size, count in sorted(pair_histogram.items())
        }
        == certificate["pair_union_cardinality_histogram"],
        "total_weight": total_weight
        == int(certificate["total_point_weight"]),
        "maximum_mask_weight": maximum_mask_weight
        == int(certificate["maximum_single_mask_weight"]),
        "weighted_flag": (
            radius * maximum_mask_weight < total_weight
        )
        == bool(certificate["certified"]),
        "pairwise_flag": (
            radius == 2 and maximum_pair_union < len(points)
        )
        == bool(certificate["pairwise_union_certified"]),
        "tight_equality_metadata": tight_metadata_matches,
    }
    verified = all(checks.values()) and (
        radius * maximum_mask_weight < total_weight
        or (
            tight_fields_present
            and tight_equality_certified
        )
        or (radius == 2 and maximum_pair_union < len(points))
    )
    result = {
        "certificate": str(args.certificate),
        "verified": verified,
        "checks": checks,
        "point_count": len(points),
        "gain_mask_count": len(masks),
        "inclusion_maximal_mask_count": len(maximal),
        "maximum_pair_union_cardinality": maximum_pair_union,
        "positive_weight_point_count": positive_mask.bit_count(),
        "tight_mask_count": len(tight_masks),
        "distinct_tight_positive_mask_count": len(
            tight_positive_masks
        ),
        "tight_disjoint_cover_exists": tight_cover_exists,
        "tight_equality_certified": tight_equality_certified,
        "scope": (
            "independent scalar reconstruction of every legal one-row "
            "gain mask on the embedded finite point set"
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"verified={verified} points={len(points)} masks={len(masks)} "
        f"maximal={len(maximal)} max_pair_union={maximum_pair_union} "
        f"tight_equality={tight_equality_certified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
