#!/usr/bin/env python3
"""Fast gain-mask search for a bounded finite affine-cover repair.

The current phase's uncovered points define a small bit universe.  Legal row
retargets are projected to gain masks on that universe, and a tiny SAT model
enumerates at most ``R`` masks that cover it.  Every mask model is then given
an exact distinct-row assignment and replayed on *all* supplied points.

Successful output is therefore an exact finite repair.  For each selected
gain-mask skeleton, the search also adds any remaining zero-gain or
duplicate-mask moves needed to cover secondary losses.  Exhaustion is
therefore complete for the declared finite Hamming ball.
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
            "minimum_relaxed_mask_count": 0,
            "mask_models": 0,
            "matching_failures": 0,
            "replay_failures": 0,
            "complete_negative": False,
            "build_seconds": build_seconds,
            "search_seconds": 0.0,
        }
    if max_changes == 0:
        return {
            "status": "UNSAT",
            "repaired": None,
            "changed_phases": None,
            "full_misses": len(miss_indices),
            "initial_misses": len(miss_indices),
            "gain_mask_count": 0,
            "gain_move_count": 0,
            "minimum_relaxed_mask_count": None,
            "mask_models": 0,
            "matching_failures": 0,
            "replay_failures": 0,
            "complete_negative": True,
            "build_seconds": build_seconds,
            "search_seconds": 0.0,
        }

    miss_position = {
        point_index: bit
        for bit, point_index in enumerate(miss_indices)
    }
    support: dict[int, list[tuple[int, int]]] = defaultdict(list)
    move_gain_mask: dict[tuple[int, int], int] = {}
    moves_covering_point: list[list[tuple[int, int]]] = [
        [] for _point in points
    ]
    for row_index, row in enumerate(rows):
        if int(row["p"]) in fixed_primes:
            continue
        current = int(assignment[row_index])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        masks_by_target: dict[int, int] = {}
        observed_points_by_target: dict[int, list[int]] = defaultdict(list)
        for point_index, raw_target in enumerate(targets[:, row_index]):
            target = int(raw_target)
            if target == current or target % modulus != residue:
                continue
            observed_points_by_target[target].append(point_index)
            if point_index in miss_position:
                masks_by_target[target] = (
                    masks_by_target.get(target, 0)
                    | (1 << miss_position[point_index])
                )
        for target, observed_points in observed_points_by_target.items():
            move = (row_index, target)
            gain_mask = masks_by_target.get(target, 0)
            move_gain_mask[move] = gain_mask
            for point_index in observed_points:
                moves_covering_point[point_index].append(move)
        for target, mask in masks_by_target.items():
            support[mask].append((row_index, target))
    masks = sorted(support, key=lambda mask: (-mask.bit_count(), mask))
    if not masks:
        return {
            "status": "UNSAT",
            "repaired": None,
            "changed_phases": None,
            "full_misses": len(miss_indices),
            "initial_misses": len(miss_indices),
            "gain_mask_count": 0,
            "gain_move_count": 0,
            "minimum_relaxed_mask_count": None,
            "mask_models": 0,
            "matching_failures": 0,
            "replay_failures": 0,
            "complete_negative": True,
            "build_seconds": build_seconds,
            "search_seconds": 0.0,
        }

    def relaxed_mask_count() -> int | None:
        for bound in range(1, max_changes + 1):
            local_pool = IDPool()
            local_variables = {
                mask: local_pool.id(("gain_mask", mask))
                for mask in masks
            }
            local_clauses = CardEnc.atmost(
                list(local_variables.values()),
                bound=bound,
                vpool=local_pool,
                encoding=EncType.seqcounter,
            ).clauses
            for bit in range(len(miss_indices)):
                local_clauses.append(
                    [
                        local_variables[mask]
                        for mask in masks
                        if mask & (1 << bit)
                    ]
                )
            local_solver = Solver(
                name=solver_name,
                bootstrap_with=local_clauses,
            )
            sat = local_solver.solve()
            local_solver.delete()
            if sat:
                return bound
        return None

    minimum_relaxed_mask_count = relaxed_mask_count()
    if minimum_relaxed_mask_count is None:
        return {
            "status": "UNSAT",
            "repaired": None,
            "changed_phases": None,
            "full_misses": len(miss_indices),
            "initial_misses": len(miss_indices),
            "gain_mask_count": len(masks),
            "gain_move_count": sum(map(len, support.values())),
            "minimum_relaxed_mask_count": None,
            "mask_models": 0,
            "matching_failures": 0,
            "replay_failures": 0,
            "complete_negative": True,
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
        selected_mask_set = set(selected)
        used_rows = set()
        chosen: list[tuple[int, int, int]] = []
        matched_leaf = False

        def assign_masks(depth: int):
            nonlocal matched_leaf
            if depth == len(ordered):
                matched_leaf = True
                counts = base_cover.astype(np.int32, copy=True)
                for _mask, row_index, target in chosen:
                    counts -= (
                        targets[:, row_index]
                        == int(assignment[row_index])
                    )
                    counts += targets[:, row_index] == target

                def add_compensators(remaining: int):
                    deficits = [
                        int(index)
                        for index in np.flatnonzero(counts == 0)
                    ]
                    if not deficits:
                        repaired = list(assignment)
                        for _mask, row_index, target in chosen:
                            repaired[row_index] = target
                        return repaired
                    if remaining == 0:
                        return None
                    eligible_by_point = []
                    for point_index in deficits:
                        eligible = []
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
                                    and gain_mask not in selected_mask_set
                                )
                            ):
                                continue
                            seen_rows.add(row_index)
                            eligible.append(
                                (gain_mask, row_index, target)
                            )
                        eligible_by_point.append((len(eligible), eligible))
                    _count, eligible = min(eligible_by_point)
                    for gain_mask, row_index, target in eligible:
                        used_rows.add(row_index)
                        chosen.append((gain_mask, row_index, target))
                        counts[:] -= (
                            targets[:, row_index]
                            == int(assignment[row_index])
                        )
                        counts[:] += targets[:, row_index] == target
                        repaired = add_compensators(remaining - 1)
                        if repaired is not None:
                            return repaired
                        counts[:] -= targets[:, row_index] == target
                        counts[:] += (
                            targets[:, row_index]
                            == int(assignment[row_index])
                        )
                        chosen.pop()
                        used_rows.remove(row_index)
                    return None

                return add_compensators(max_changes - len(chosen))
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
    complete_negative = (
        winner is None and status == "NO_GAIN_MASK_MODEL"
    )
    if complete_negative:
        status = "UNSAT"
    return {
        "status": status,
        "repaired": winner,
        "changed_phases": changed,
        "full_misses": full_misses,
        "initial_misses": len(miss_indices),
        "gain_mask_count": len(masks),
        "gain_move_count": sum(map(len, support.values())),
        "minimum_relaxed_mask_count": minimum_relaxed_mask_count,
        "mask_models": mask_models,
        "matching_failures": matching_failures,
        "replay_failures": replay_failures,
        "complete_negative": complete_negative,
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
            "successful phases are exact finite repairs; UNSAT exhausts "
            "all gain-mask skeletons, distinct row owners, and remaining "
            "zero-gain or duplicate-mask compensator chains"
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
        f"minimum_masks={result['minimum_relaxed_mask_count']} "
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
