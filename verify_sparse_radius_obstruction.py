#!/usr/bin/env python3
"""Independent scalar verifier for sparse-radius obstruction certificates.

This verifier deliberately does not import the discovery engine, NumPy,
PySAT, Z3, or SciPy.  It reconstructs every finite affine target with Python
integers and exhausts gain-mask skeletons, distinct row owners, and all
zero-gain or duplicate-mask compensator chains.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import defaultdict
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_certificate(
    certificate_path: Path | None,
    *,
    certificate_data: dict | None = None,
    pool_data: dict | None = None,
    phase_data: dict[int, int] | None = None,
) -> dict:
    """Verify a standalone certificate or an embedded partition-tree leaf.

    A partition-tree verifier may supply already authenticated pool data and
    a phase map derived from the tree path.  Standalone callers retain the
    original hash-checked file behavior.
    """
    if certificate_data is None:
        if certificate_path is None:
            raise ValueError("certificate_path is required")
        certificate = json.loads(certificate_path.read_text())
    else:
        certificate = certificate_data

    if pool_data is None:
        pool_path = Path(certificate["pool"])
        if sha256_file(pool_path) != certificate["pool_sha256"]:
            raise RuntimeError("pool hash mismatch")
        pool = json.loads(pool_path.read_text())
    else:
        pool = pool_data

    if phase_data is None:
        phase_path = Path(certificate["base_phase"])
        if sha256_file(phase_path) != certificate["base_phase_sha256"]:
            raise RuntimeError("base-phase hash mismatch")
        phases = {
            int(prime): int(target)
            for prime, target in json.loads(phase_path.read_text()).items()
        }
    else:
        phases = {
            int(prime): int(target)
            for prime, target in phase_data.items()
        }

    rows = pool["choices"]
    points = [tuple(map(int, point)) for point in certificate["points"]]
    if len(points) != len(set(points)):
        raise RuntimeError("certificate repeats a point")
    if any(k < 0 or ell < 0 for k, ell in points):
        raise RuntimeError("certificate contains a negative exponent")
    fixed_primes = set(map(int, certificate["fixed_primes"]))
    max_changes = int(certificate["max_changes"])
    if max_changes < 0:
        raise RuntimeError("negative change radius")

    assignment = []
    target_columns = []
    base_counts = [0] * len(points)
    seen_primes = set()
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        a = int(row["a"])
        b = int(row["b"])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        if prime in seen_primes or h < 1 or modulus < 1 or h % modulus:
            raise RuntimeError("pool row is malformed")
        seen_primes.add(prime)
        current = int(phases.get(prime, residue)) % h
        if current % modulus != residue:
            raise RuntimeError(f"base phase is illegal for p={prime}")
        targets = [
            (a * (k % h) + b * (ell % h)) % h
            for k, ell in points
        ]
        assignment.append(current)
        target_columns.append(targets)
        for point_index, target in enumerate(targets):
            if target == current:
                base_counts[point_index] += 1
    if fixed_primes - seen_primes:
        raise RuntimeError("a fixed prime is absent from the pool")

    miss_indices = [
        index for index, count in enumerate(base_counts) if count == 0
    ]
    miss_position = {
        point_index: bit
        for bit, point_index in enumerate(miss_indices)
    }
    support: dict[int, list[tuple[int, int]]] = defaultdict(list)
    move_gain_mask = {}
    moves_covering_point = [[] for _point in points]
    for row_index, row in enumerate(rows):
        if int(row["p"]) in fixed_primes:
            continue
        current = assignment[row_index]
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        observed: dict[int, list[int]] = defaultdict(list)
        for point_index, target in enumerate(
            target_columns[row_index]
        ):
            if target != current and target % modulus == residue:
                observed[target].append(point_index)
        for target, point_indices in observed.items():
            gain_mask = 0
            for point_index in point_indices:
                if point_index in miss_position:
                    gain_mask |= 1 << miss_position[point_index]
                moves_covering_point[point_index].append(
                    (row_index, target)
                )
            move_gain_mask[row_index, target] = gain_mask
            if gain_mask:
                support[gain_mask].append((row_index, target))
    masks = tuple(sorted(support, key=lambda value: (-value.bit_count(), value)))
    masks_by_bit = [
        tuple(mask for mask in masks if mask & (1 << bit))
        for bit in range(len(miss_indices))
    ]
    full_mask = (1 << len(miss_indices)) - 1

    started = time.monotonic()
    visited = set()
    full_skeletons = 0
    matching_failures = 0
    replay_failures = 0
    minimum_relaxed = None
    repair = None

    def replay_skeleton(selected: frozenset[int]):
        nonlocal matching_failures, replay_failures
        ordered = sorted(
            selected,
            key=lambda mask: len(
                {row_index for row_index, _target in support[mask]}
            ),
        )
        selected_set = set(selected)
        used_rows = set()
        chosen = []
        matched_leaf = False

        def assign_masks(depth: int):
            nonlocal matched_leaf
            if depth == len(ordered):
                matched_leaf = True
                counts = list(base_counts)
                for _mask, row_index, target in chosen:
                    column = target_columns[row_index]
                    current = assignment[row_index]
                    for point_index, observed in enumerate(column):
                        counts[point_index] += (
                            int(observed == target)
                            - int(observed == current)
                        )

                def compensate(remaining: int):
                    deficits = [
                        index
                        for index, count in enumerate(counts)
                        if count == 0
                    ]
                    if not deficits:
                        return list(chosen)
                    if remaining == 0:
                        return None
                    candidates_by_deficit = []
                    for point_index in deficits:
                        candidates = []
                        seen_rows = set()
                        for row_index, target in moves_covering_point[
                            point_index
                        ]:
                            gain_mask = move_gain_mask[row_index, target]
                            if (
                                row_index in used_rows
                                or row_index in seen_rows
                                or (
                                    gain_mask
                                    and gain_mask not in selected_set
                                )
                            ):
                                continue
                            seen_rows.add(row_index)
                            candidates.append(
                                (gain_mask, row_index, target)
                            )
                        candidates_by_deficit.append(candidates)
                    candidates = min(
                        candidates_by_deficit,
                        key=len,
                    )
                    for gain_mask, row_index, target in candidates:
                        used_rows.add(row_index)
                        chosen.append((gain_mask, row_index, target))
                        column = target_columns[row_index]
                        current = assignment[row_index]
                        for point_index, observed in enumerate(column):
                            counts[point_index] += (
                                int(observed == target)
                                - int(observed == current)
                            )
                        result = compensate(remaining - 1)
                        if result is not None:
                            return result
                        for point_index, observed in enumerate(column):
                            counts[point_index] -= (
                                int(observed == target)
                                - int(observed == current)
                            )
                        chosen.pop()
                        used_rows.remove(row_index)
                    return None

                return compensate(max_changes - len(chosen))
            mask = ordered[depth]
            seen_rows = set()
            for row_index, target in support[mask]:
                if row_index in used_rows or row_index in seen_rows:
                    continue
                seen_rows.add(row_index)
                used_rows.add(row_index)
                chosen.append((mask, row_index, target))
                result = assign_masks(depth + 1)
                if result is not None:
                    return result
                chosen.pop()
                used_rows.remove(row_index)
            return None

        result = assign_masks(0)
        if result is not None:
            return result, True
        if not matched_leaf:
            matching_failures += 1
            return None, False
        replay_failures += 1
        return None, True

    def enumerate_skeletons(
        selected: frozenset[int],
        covered: int,
    ):
        nonlocal full_skeletons, minimum_relaxed, repair
        if repair is not None or selected in visited:
            return
        visited.add(selected)
        if covered == full_mask:
            full_skeletons += 1
            if minimum_relaxed is None:
                minimum_relaxed = len(selected)
            else:
                minimum_relaxed = min(
                    minimum_relaxed,
                    len(selected),
                )
            result, matchable = replay_skeleton(selected)
            if result is not None:
                repair = result
                return
            if not matchable or len(selected) == max_changes:
                return
            for mask in masks:
                if mask not in selected:
                    enumerate_skeletons(
                        selected | {mask},
                        covered | mask,
                    )
            return
        if len(selected) == max_changes:
            return
        uncovered_bits = [
            bit
            for bit in range(len(miss_indices))
            if not covered & (1 << bit)
        ]
        bit = min(
            uncovered_bits,
            key=lambda value: sum(
                mask not in selected for mask in masks_by_bit[value]
            ),
        )
        for mask in masks_by_bit[bit]:
            if mask not in selected:
                enumerate_skeletons(
                    selected | {mask},
                    covered | mask,
                )

    if miss_indices and max_changes:
        enumerate_skeletons(frozenset(), 0)
    elif not miss_indices:
        repair = []

    discovery = certificate["discovery"]
    expected = {
        "initial_misses": len(miss_indices),
        "gain_mask_count": len(masks),
        "gain_move_count": sum(map(len, support.values())),
        "minimum_relaxed_mask_count": minimum_relaxed,
    }
    count_match = all(
        discovery.get(key) == value for key, value in expected.items()
    )
    verified = (
        certificate.get("complete") is True
        and discovery.get("status") == "UNSAT"
        and discovery.get("complete_negative") is True
        and repair is None
        and count_match
    )
    return {
        "certificate": (
            str(certificate_path)
            if certificate_path is not None
            else "<embedded partition leaf>"
        ),
        "verified": verified,
        "repair_exists": repair is not None,
        "row_count": len(rows),
        "point_count": len(points),
        "fixed_prime_count": len(fixed_primes),
        "max_changes": max_changes,
        **expected,
        "full_skeleton_count": full_skeletons,
        "matching_failures": matching_failures,
        "replay_failures": replay_failures,
        "discovery_count_match": count_match,
        "elapsed_seconds": time.monotonic() - started,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify_certificate(args.certificate)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"verified={result['verified']} "
        f"repair_exists={result['repair_exists']} "
        f"points={result['point_count']} "
        f"masks={result['gain_mask_count']} "
        f"skeletons={result['full_skeleton_count']} "
        f"elapsed_s={result['elapsed_seconds']:.3f}",
        flush=True,
    )
    return 0 if result["verified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
