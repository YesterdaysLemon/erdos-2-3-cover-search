#!/usr/bin/env python3
"""Generate a family of mutually separated exact finite mask repairs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import exact_greedy
from analyze_relaxed_mask_dual import load_point_file
from finite_sample_mask_repair import search_mask_repair


def phase_assignment(
    rows: list[dict],
    path: Path,
) -> list[int]:
    phase_map = {
        int(prime): int(target)
        for prime, target in json.loads(path.read_text()).items()
    }
    missing = {int(row["p"]) for row in rows} - phase_map.keys()
    if missing:
        raise RuntimeError(
            f"phase {path} omits primes {sorted(missing)[:10]}"
        )
    assignment = []
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        target = int(phase_map[prime]) % h
        if target % modulus != residue:
            raise RuntimeError(f"phase {path} violates p={prime}")
        assignment.append(target)
    return assignment


def hamming_distance(left: list[int], right: list[int]) -> int:
    return sum(
        first != second
        for first, second in zip(left, right, strict=True)
    )


def generate_repairs(
    rows: list[dict],
    candidates: list[tuple],
    points: list[tuple[int, int]],
    initial: dict[int, int],
    fixed_primes: set[int],
    max_changes: int,
    solver_name: str,
    max_mask_models: int,
    time_limit: float,
    avoided: list[list[int]],
    min_hamming_distance: int,
    count: int,
    progress=None,
) -> tuple[list[list[int]], list[dict]]:
    generated = []
    results = []
    all_avoided = list(avoided)
    for _index in range(count):
        result = search_mask_repair(
            rows,
            candidates,
            points,
            initial,
            fixed_primes,
            max_changes,
            solver_name,
            max_mask_models,
            time_limit,
            all_avoided,
            min_hamming_distance,
        )
        results.append(result)
        if result["status"] != "INTEGER_MODEL":
            break
        repaired = [int(target) for target in result["repaired"]]
        if any(
            hamming_distance(repaired, prior) < min_hamming_distance
            for prior in all_avoided
        ):
            raise AssertionError("generated repair violates separation")
        generated.append(repaired)
        all_avoided.append(repaired)
        if progress is not None:
            progress(len(generated), repaired, result)
    return generated, results


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path, nargs="+")
    parser.add_argument("--initial-phases", type=Path, required=True)
    parser.add_argument(
        "--avoid-phase",
        action="append",
        type=Path,
        default=[],
    )
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--min-hamming-distance", type=int, required=True)
    parser.add_argument("--fixed-primes", default="")
    parser.add_argument("--max-changes", type=int, required=True)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--max-mask-models", type=int, default=200_000)
    parser.add_argument("--time-limit", type=float, default=240.0)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    args = parser.parse_args()
    if args.count < 1:
        raise SystemExit("--count must be positive")
    if args.min_hamming_distance < 0:
        raise SystemExit("--min-hamming-distance must be nonnegative")

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
            args.initial_phases.read_text()
        ).items()
    }
    fixed_primes = {
        int(value) for value in args.fixed_primes.split(",") if value
    }
    avoided = [
        phase_assignment(rows, path) for path in args.avoid_phase
    ]
    prefix = str(args.output_prefix)
    records = []

    def save_progress(index, repaired, result):
        phase_path = Path(f"{prefix}_{index:02d}_phase.json")
        result_path = Path(f"{prefix}_{index:02d}_result.json")
        phase_map = {
            str(int(row["p"])): int(target)
            for row, target in zip(rows, repaired, strict=True)
        }
        phase_path.write_text(json.dumps(phase_map) + "\n")
        public_result = {
            key: value
            for key, value in result.items()
            if key != "repaired"
        }
        result_path.write_text(json.dumps(public_result, indent=2) + "\n")
        records.append(
            {
                "phase": str(phase_path),
                "result": str(result_path),
                "changed_phases": result["changed_phases"],
                "mask_models": result["mask_models"],
                "diversity_rejections": result[
                    "diversity_rejections"
                ],
            }
        )
        print(
            f"generated={index}/{args.count} "
            f"mask_models={result['mask_models']} "
            f"diversity_rejections={result['diversity_rejections']}",
            flush=True,
        )

    generated, results = generate_repairs(
        rows,
        candidates,
        points,
        initial,
        fixed_primes,
        args.max_changes,
        args.solver,
        args.max_mask_models,
        args.time_limit,
        avoided,
        args.min_hamming_distance,
        args.count,
        save_progress,
    )
    family = [*avoided, *generated]
    pairwise_distances = [
        hamming_distance(left, right)
        for index, left in enumerate(family)
        for right in family[index + 1 :]
    ]
    manifest = {
        "pool": str(args.pool),
        "points": [str(path) for path in args.points],
        "initial_phases": str(args.initial_phases),
        "seed_phases": [str(path) for path in args.avoid_phase],
        "requested_count": args.count,
        "generated_count": len(generated),
        "minimum_hamming_distance": args.min_hamming_distance,
        "pairwise_distance_minimum": (
            min(pairwise_distances) if pairwise_distances else None
        ),
        "pairwise_distance_maximum": (
            max(pairwise_distances) if pairwise_distances else None
        ),
        "records": records,
        "terminal_status": results[-1]["status"] if results else None,
        "complete": len(generated) == args.count,
        "scope": "exact finite-corpus repairs only",
    }
    args.manifest_output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"requested={args.count} generated={len(generated)} "
        f"family={len(family)} "
        f"distance_min={manifest['pairwise_distance_minimum']} "
        f"status={manifest['terminal_status']}",
        flush=True,
    )
    return 0 if manifest["complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
