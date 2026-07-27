#!/usr/bin/env python3
"""Fast gain-mask search for a bounded finite affine-cover repair.

The current phase's uncovered points define a small bit universe.  Legal row
retargets are projected to gain masks on that universe, and a tiny SAT model
enumerates at most ``R`` masks that cover it.  Every mask model is then given
an exact distinct-row assignment and replayed on *all* supplied points.

Successful output is therefore an exact finite repair.  Exhaustion is only a
one-sided negative result: a genuine repair could contain a zero-initial-gain
move used solely to compensate for coverage lost by another move.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import exact_greedy
from certify_anchor_phase_quotient import load_cells
from local_phase_cegis import build_targets


def search_mask_repair(
    rows: list[dict],
    candidates: list[tuple],
    points: list[tuple[int, int]],
    initial: dict[int, int],
    fixed_primes: set[int],
    max_changes: int,
    solver_name: str,
    max_mask_models: int,
    time_limit: float,
):
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    assignment = []
    for row in rows:
        prime = int(row["p"])
        h = int(row["h"])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        target = int(initial.get(prime, residue)) % h
        if target % modulus != residue:
            raise RuntimeError(f"invalid initial phase for p={prime}")
        assignment.append(target)
    targets, build_seconds = build_targets(points, candidates, np)
    assignment_array = np.asarray(assignment, dtype=targets.dtype)
    base_cover = np.count_nonzero(
        targets == assignment_array,
        axis=1,
    )
    miss_indices = [
        int(index) for index in np.flatnonzero(base_cover == 0)
    ]
    if not miss_indices:
        return {
            "status": "INTEGER_MODEL",
            "repaired": assignment,
            "changed_phases": 0,
            "full_misses": 0,
            "initial_misses": 0,
            "gain_mask_count": 0,
            "gain_move_count": 0,
            "mask_models": 0,
            "matching_failures": 0,
            "replay_failures": 0,
            "complete_negative": False,
            "build_seconds": build_seconds,
            "search_seconds": 0.0,
        }
    if max_changes == 0:
        return {
            "status": "NO_GAIN_MASK_MODEL",
            "repaired": None,
            "changed_phases": None,
            "full_misses": len(miss_indices),
            "initial_misses": len(miss_indices),
            "gain_mask_count": 0,
            "gain_move_count": 0,
            "mask_models": 0,
            "matching_failures": 0,
            "replay_failures": 0,
            "complete_negative": False,
            "build_seconds": build_seconds,
            "search_seconds": 0.0,
        }

    miss_position = {
        point_index: bit
        for bit, point_index in enumerate(miss_indices)
    }
    support: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for row_index, row in enumerate(rows):
        if int(row["p"]) in fixed_primes:
            continue
        current = int(assignment[row_index])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        masks_by_target: dict[int, int] = {}
        for point_index in miss_indices:
            target = int(targets[point_index, row_index])
            if target == current or target % modulus != residue:
                continue
            masks_by_target[target] = (
                masks_by_target.get(target, 0)
                | (1 << miss_position[point_index])
            )
        for target, mask in masks_by_target.items():
            support[mask].append((row_index, target))
    masks = sorted(support, key=lambda mask: (-mask.bit_count(), mask))
    if not masks:
        return {
            "status": "NO_GAIN_MASK_MODEL",
            "repaired": None,
            "changed_phases": None,
            "full_misses": len(miss_indices),
            "initial_misses": len(miss_indices),
            "gain_mask_count": 0,
            "gain_move_count": 0,
            "mask_models": 0,
            "matching_failures": 0,
            "replay_failures": 0,
            "complete_negative": False,
            "build_seconds": build_seconds,
            "search_seconds": 0.0,
        }

    vpool = IDPool()
    mask_variable = {
        mask: vpool.id(("gain_mask", mask)) for mask in masks
    }
    clauses = CardEnc.atmost(
        list(mask_variable.values()),
        bound=max_changes,
        vpool=vpool,
        encoding=EncType.seqcounter,
    ).clauses
    for bit in range(len(miss_indices)):
        clauses.append(
            [
                mask_variable[mask]
                for mask in masks
                if mask & (1 << bit)
            ]
        )
    solver = Solver(name=solver_name, bootstrap_with=clauses)
    started = time.monotonic()
    mask_models = 0
    matching_failures = 0
    replay_failures = 0
    winner = None
    status = "NO_GAIN_MASK_MODEL"
    all_mask_set = set(masks)

    while solver.solve():
        if max_mask_models and mask_models >= max_mask_models:
            status = "SEARCH_LIMIT"
            break
        if time_limit > 0 and time.monotonic() - started >= time_limit:
            status = "SEARCH_LIMIT"
            break
        model = {literal for literal in solver.get_model() if literal > 0}
        selected = [
            mask for mask in masks if mask_variable[mask] in model
        ]
        mask_models += 1
        ordered = sorted(
            selected,
            key=lambda mask: len(
                {row_index for row_index, _target in support[mask]}
            ),
        )
        used_rows = set()
        chosen: list[tuple[int, int, int]] = []
        matched_leaf = False

        def assign_masks(depth: int):
            nonlocal matched_leaf
            if depth == len(ordered):
                matched_leaf = True
                repaired = list(assignment)
                for _mask, row_index, target in chosen:
                    repaired[row_index] = target
                repaired_array = np.asarray(
                    repaired,
                    dtype=targets.dtype,
                )
                cover = np.count_nonzero(
                    targets == repaired_array,
                    axis=1,
                )
                if np.any(cover == 0):
                    return None
                return repaired
            mask = ordered[depth]
            seen_rows = set()

            def loss_score(move):
                row_index, target = move
                return sum(
                    1
                    for point_index, count in enumerate(base_cover)
                    if (
                        count <= 1
                        and int(targets[point_index, row_index])
                        == int(assignment[row_index])
                        and int(targets[point_index, row_index]) != target
                    )
                )

            for row_index, target in sorted(
                support[mask],
                key=loss_score,
            ):
                if row_index in used_rows or row_index in seen_rows:
                    continue
                seen_rows.add(row_index)
                used_rows.add(row_index)
                chosen.append((mask, row_index, target))
                repaired = assign_masks(depth + 1)
                if repaired is not None:
                    return repaired
                chosen.pop()
                used_rows.remove(row_index)
            return None

        repaired = assign_masks(0)
        if repaired is not None:
            winner = repaired
            status = "INTEGER_MODEL"
            break
        if not matched_leaf:
            matching_failures += 1
            # Adding more masks cannot make an unmatchable selected subset
            # matchable, so block every superset of this selection.
            solver.add_clause(
                [-mask_variable[mask] for mask in selected]
            )
        else:
            replay_failures += 1
            # Block only this exact mask set.  A strict superset may repair
            # secondary points that these phase changes uncover.
            selected_set = set(selected)
            solver.add_clause(
                [-mask_variable[mask] for mask in selected]
                + [
                    mask_variable[mask]
                    for mask in all_mask_set - selected_set
                ]
            )
    search_seconds = time.monotonic() - started
    solver.delete()

    changed = None
    full_misses = None
    if winner is not None:
        repaired_array = np.asarray(winner, dtype=targets.dtype)
        cover = np.count_nonzero(targets == repaired_array, axis=1)
        full_misses = int(np.count_nonzero(cover == 0))
        changed = sum(
            old != new for old, new in zip(assignment, winner)
        )
        if full_misses or changed > max_changes:
            raise AssertionError("mask repair failed exact finite replay")
    return {
        "status": status,
        "repaired": winner,
        "changed_phases": changed,
        "full_misses": full_misses,
        "initial_misses": len(miss_indices),
        "gain_mask_count": len(masks),
        "gain_move_count": sum(map(len, support.values())),
        "mask_models": mask_models,
        "matching_failures": matching_failures,
        "replay_failures": replay_failures,
        "complete_negative": False,
        "build_seconds": build_seconds,
        "search_seconds": search_seconds,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path, nargs="+")
    parser.add_argument("--initial-phases", type=Path, required=True)
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    parser.add_argument("--fixed-primes", default="")
    parser.add_argument("--max-changes", type=int, default=4)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--max-mask-models", type=int, default=200_000)
    parser.add_argument("--time-limit", type=float, default=240.0)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    candidates = exact_greedy.load_candidates(args.pool, True)
    points = []
    seen_points = set()
    for path in args.points:
        for point in load_cells(path):
            if point not in seen_points:
                seen_points.add(point)
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
    result = search_mask_repair(
        rows,
        candidates,
        points,
        initial,
        fixed_primes,
        args.max_changes,
        args.solver,
        args.max_mask_models,
        args.time_limit,
    )
    if result["status"] == "INTEGER_MODEL":
        phase_map = {
            str(int(row["p"])): int(target)
            for row, target in zip(rows, result["repaired"])
        }
        args.phase_output.write_text(json.dumps(phase_map) + "\n")
        return_code = 0
    elif result["status"] == "SEARCH_LIMIT":
        return_code = 1
    else:
        return_code = 2
    output = {
        "pool": str(args.pool),
        "points": [str(path) for path in args.points],
        "initial_phases": str(args.initial_phases),
        "result": result["status"],
        "engine": "gain-mask-sat-with-exact-full-corpus-replay",
        "scope": (
            "successful phases are exact finite repairs; negative results "
            "do not exclude zero-initial-gain compensator moves"
        ),
        "fixed_primes": sorted(fixed_primes),
        "max_changes": args.max_changes,
        "point_count": len(points),
        **{
            key: value
            for key, value in result.items()
            if key not in {"status", "repaired"}
        },
    }
    args.result_output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"rows={len(rows)} points={len(points)} "
        f"initial_misses={result['initial_misses']} "
        f"masks={result['gain_mask_count']} "
        f"moves={result['gain_move_count']} "
        f"models={result['mask_models']} "
        f"result={result['status']} "
        f"changed={result['changed_phases']} "
        f"search_s={result['search_seconds']:.4f}",
        flush=True,
    )
    return return_code


if __name__ == "__main__":
    raise SystemExit(main())
