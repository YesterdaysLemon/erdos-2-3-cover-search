#!/usr/bin/env python3
"""Greedily phase new derived rows against a saved exact miss set.

Saved phases are immutable.  Each row absent from the saved map may choose
one target satisfying its derived target congruence.  A lazy greedy set-cover
pass chooses targets that cover the largest number of still-uncovered input
points, then assigns the canonical target to every unused new row.
"""

from __future__ import annotations

import argparse
import heapq
import json
from pathlib import Path


def load_points(path: Path) -> list[tuple[int, int]]:
    payload = json.loads(path.read_text())
    if isinstance(payload, dict):
        raw = payload.get("points", payload.get("misses"))
        if raw is None:
            raise RuntimeError("point artifact has neither points nor misses")
    else:
        raw = payload
    points = [tuple(map(int, point)) for point in raw]
    if any(len(point) != 2 for point in points):
        raise RuntimeError("every point must have two coordinates")
    return points


def best_option(
    options: dict[int, int], uncovered: int
) -> tuple[int, int, int]:
    best_gain = -1
    best_target = -1
    best_mask = 0
    for target, mask in options.items():
        active = mask & uncovered
        gain = active.bit_count()
        if gain > best_gain or (
            gain == best_gain and target < best_target
        ):
            best_gain = gain
            best_target = target
            best_mask = active
    return best_gain, best_target, best_mask


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("source_phases", type=Path)
    parser.add_argument("points", type=Path)
    parser.add_argument("--required-coverage", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    if args.required_coverage < 1:
        raise SystemExit("--required-coverage must be positive")

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    source = {
        int(prime): int(target)
        for prime, target in json.loads(args.source_phases.read_text()).items()
    }
    points = load_points(args.points)
    coverage = [0] * len(points)
    phases: dict[str, int] = {}
    missing_rows = []

    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        if h % modulus:
            raise RuntimeError(
                f"target modulus {modulus} does not divide h={h} for p={prime}"
            )
        if prime in source:
            target = source[prime] % h
            if target % modulus != residue:
                raise RuntimeError(
                    f"saved phase violates target restriction for p={prime}"
                )
            phases[str(prime)] = target
            a = int(row["a"])
            b = int(row["b"])
            for index, (k, l) in enumerate(points):
                if (a * k + b * l - target) % h == 0:
                    coverage[index] += 1
        else:
            missing_rows.append(row)

    initially_satisfied = sum(
        count >= args.required_coverage for count in coverage
    )
    option_maps: list[dict[int, int]] = []
    heap: list[tuple[int, int, int]] = []
    uncovered = 0
    for index, count in enumerate(coverage):
        if count < args.required_coverage:
            uncovered |= 1 << index
    for index, row in enumerate(missing_rows):
        h = int(row["h"])
        a = int(row["a"])
        b = int(row["b"])
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        options: dict[int, int] = {}
        for point_index, (k, l) in enumerate(points):
            target = (a * k + b * l) % h
            if target % modulus != residue:
                continue
            options[target] = options.get(target, 0) | (1 << point_index)
        if residue not in options:
            options[residue] = 0
        option_maps.append(options)
        gain, target, _mask = best_option(options, uncovered)
        heapq.heappush(heap, (-gain, index, target))

    selected = []
    assigned: set[int] = set()
    while uncovered and heap:
        negative_gain, index, stale_target = heapq.heappop(heap)
        if index in assigned:
            continue
        gain, target, active_mask = best_option(
            option_maps[index], uncovered
        )
        if gain != -negative_gain or target != stale_target:
            heapq.heappush(heap, (-gain, index, target))
            continue
        if gain <= 0:
            break
        row = missing_rows[index]
        prime = int(row["p"])
        phases[str(prime)] = target
        assigned.add(index)
        full_mask = option_maps[index][target]
        pending = full_mask
        while pending:
            least_bit = pending & -pending
            point_index = least_bit.bit_length() - 1
            coverage[point_index] += 1
            if coverage[point_index] >= args.required_coverage:
                uncovered &= ~least_bit
            pending ^= least_bit
        selected.append(
            {
                "p": prime,
                "h": int(row["h"]),
                "target": target,
                "canonical_target": int(row["target_residue"]),
                "new_points_covered": gain,
            }
        )

    for index, row in enumerate(missing_rows):
        if index in assigned:
            continue
        phases[str(int(row["p"]))] = int(row["target_residue"])

    remaining_indices = [
        index for index in range(len(points)) if uncovered & (1 << index)
    ]
    audit = {
        "pool": str(args.pool),
        "source_phases": str(args.source_phases),
        "points": str(args.points),
        "point_count": len(points),
        "required_coverage": args.required_coverage,
        "source_phase_count": len(source),
        "retained_phase_count": len(rows) - len(missing_rows),
        "new_row_count": len(missing_rows),
        "points_satisfied_by_retained_phases": initially_satisfied,
        "final_minimum_coverage": min(coverage, default=0),
        "selected_new_rows": selected,
        "remaining_point_indices": remaining_indices,
        "remaining_points": [points[index] for index in remaining_indices],
    }
    args.output.write_text(json.dumps(phases) + "\n")
    args.audit.write_text(json.dumps(audit, indent=2) + "\n")
    print(
        f"rows={len(rows)} retained={len(rows) - len(missing_rows)} "
        f"new={len(missing_rows)} points={len(points)} "
        f"required_coverage={args.required_coverage} "
        f"initially_satisfied={initially_satisfied} "
        f"selected={len(selected)} remaining={len(remaining_indices)} "
        f"output={args.output}",
        flush=True,
    )
    return 0 if not remaining_indices else 2


if __name__ == "__main__":
    raise SystemExit(main())
