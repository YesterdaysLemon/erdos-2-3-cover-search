#!/usr/bin/env python3
"""Promote candidate affine rows that hit several adversarial points at once.

Generic pool extension often adds expensive high-component rows that merely
memorize one witness.  This selector instead scores every new row by the
largest set of supplied points sharing one legal phase target and retains
only rows meeting a declared multi-hit threshold.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import exact_uncovered


POINT_KEYS = ("points", "misses", "exact_misses", "cells", "new_cells")


def load_points(paths: list[Path]) -> list[tuple[int, int]]:
    points: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for path in paths:
        payload = json.loads(path.read_text())
        if isinstance(payload, dict):
            key = next((key for key in POINT_KEYS if key in payload), None)
            if key is None:
                raise RuntimeError(f"{path} contains no supported point list")
            payload = payload[key]
        for raw_k, raw_l in payload:
            point = (int(raw_k), int(raw_l))
            if point not in seen:
                seen.add(point)
                points.append(point)
    return points


def best_legal_target(
    row: dict,
    points: list[tuple[int, int]],
) -> tuple[int, int | None, list[int]]:
    h = int(row["h"])
    a = int(row["a"]) % h
    b = int(row["b"]) % h
    modulus = int(row.get("target_modulus", 1))
    residue = int(row.get("target_residue", 0)) % modulus
    if h % modulus:
        raise RuntimeError(
            f"target modulus does not divide h for p={int(row['p'])}"
        )
    by_target: dict[int, list[int]] = defaultdict(list)
    for index, (k, l) in enumerate(points):
        target = (a * k + b * l) % h
        if target % modulus == residue:
            by_target[target].append(index)
    if not by_target:
        return 0, None, []
    target, indices = min(
        by_target.items(),
        key=lambda item: (-len(item[1]), item[0]),
    )
    return len(indices), target, indices


def target_hit_indices(
    row: dict,
    points: list[tuple[int, int]],
    target: int,
) -> list[int]:
    h = int(row["h"])
    normalized = target % h
    return [
        index
        for index, (k, l) in enumerate(points)
        if (int(row["a"]) * k + int(row["b"]) * l) % h == normalized
    ]


def largest_prime_power_component(modulus: int) -> int:
    return max(
        (
            prime**exponent
            for prime, exponent in exact_uncovered.factor(modulus).items()
        ),
        default=1,
    )


def select_rows(
    base_rows: list[dict],
    candidate_rows: list[dict],
    points: list[tuple[int, int]],
    min_hit: int,
    max_rows: int,
    max_component: int = 0,
    validation_points: list[tuple[int, int]] | None = None,
    min_validation_hit: int = 0,
) -> tuple[list[dict], list[dict]]:
    base_by_prime = {int(row["p"]): row for row in base_rows}
    if len(base_by_prime) != len(base_rows):
        raise RuntimeError("base pool contains duplicate primes")
    scored = []
    seen_candidates: set[int] = set()
    for row in candidate_rows:
        prime = int(row["p"])
        if prime in seen_candidates:
            raise RuntimeError(f"candidate pool repeats prime {prime}")
        seen_candidates.add(prime)
        if prime in base_by_prime:
            if row != base_by_prime[prime]:
                raise RuntimeError(f"candidate row conflicts at prime {prime}")
            continue
        component = largest_prime_power_component(int(row["h"]))
        if max_component and component > max_component:
            continue
        hit_count, target, indices = best_legal_target(row, points)
        validation_indices = (
            target_hit_indices(row, validation_points, int(target))
            if target is not None and validation_points
            else []
        )
        if (
            hit_count >= min_hit
            and len(validation_indices) >= min_validation_hit
        ):
            scored.append(
                {
                    "row": row,
                    "p": prime,
                    "h": int(row["h"]),
                    "largest_prime_power_component": component,
                    "hit_count": hit_count,
                    "target": target,
                    "point_indices": indices,
                    "validation_hit_count": len(validation_indices),
                    "validation_point_indices": validation_indices,
                }
            )
    scored.sort(
        key=lambda record: (
            -int(record["validation_hit_count"]),
            -int(record["hit_count"]),
            int(record["largest_prime_power_component"]),
            int(record["h"]),
            int(record["p"]),
            int(record["target"]),
        )
    )
    if max_rows:
        scored = scored[:max_rows]
    rows = list(base_rows) + [record["row"] for record in scored]
    rows.sort(key=lambda row: (int(row["h"]), int(row["p"])))
    audit = [
        {
            key: value
            for key, value in record.items()
            if key != "row"
        }
        for record in scored
    ]
    return rows, audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("base_pool", type=Path)
    parser.add_argument("candidate_pool", type=Path)
    parser.add_argument("points", nargs="+", type=Path)
    parser.add_argument(
        "--validation-points",
        nargs="+",
        type=Path,
        default=[],
        help=(
            "separate point files tested at the training-selected target; "
            "these points do not influence which target is selected"
        ),
    )
    parser.add_argument("--min-hit", type=int, default=2)
    parser.add_argument("--min-validation-hit", type=int, default=0)
    parser.add_argument(
        "--max-rows",
        type=int,
        default=0,
        help="retain at most this many rows; zero retains every qualifying row",
    )
    parser.add_argument(
        "--max-component",
        type=int,
        default=0,
        help=(
            "reject rows with a larger prime-power component; zero disables "
            "the verification-cost guard"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.min_hit < 1:
        raise SystemExit("--min-hit must be positive")
    if args.max_rows < 0:
        raise SystemExit("--max-rows must be nonnegative")
    if args.max_component < 0:
        raise SystemExit("--max-component must be nonnegative")
    if args.min_validation_hit < 0:
        raise SystemExit("--min-validation-hit must be nonnegative")
    if args.min_validation_hit and not args.validation_points:
        raise SystemExit(
            "--min-validation-hit requires --validation-points"
        )

    base = json.loads(args.base_pool.read_text())
    candidates = json.loads(args.candidate_pool.read_text())
    points = load_points(args.points)
    validation_points = load_points(args.validation_points)
    rows, audit = select_rows(
        base["choices"],
        candidates["choices"],
        points,
        args.min_hit,
        args.max_rows,
        args.max_component,
        validation_points,
        args.min_validation_hit,
    )
    result = dict(base)
    result["choices"] = rows
    result["targeted_signature_selection"] = {
        "base_pool": str(args.base_pool),
        "candidate_pool": str(args.candidate_pool),
        "point_files": [str(path) for path in args.points],
        "point_count": len(points),
        "validation_point_files": [
            str(path) for path in args.validation_points
        ],
        "validation_point_count": len(validation_points),
        "min_hit": args.min_hit,
        "min_validation_hit": args.min_validation_hit,
        "max_rows": args.max_rows,
        "max_component": args.max_component,
        "selected": audit,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"points={len(points)} validation={len(validation_points)} "
        f"selected={len(audit)} "
        f"output_rows={len(rows)} output={args.output}",
        flush=True,
    )
    for record in audit:
        print(
            f"p={record['p']} h={record['h']} "
            f"component={record['largest_prime_power_component']} "
            f"hits={record['hit_count']} "
            f"validation_hits={record['validation_hit_count']} "
            f"target={record['target']}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
