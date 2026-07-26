#!/usr/bin/env python3
"""Optimize finite-sample coverage by one CRT component at a time.

This is a heuristic synthesis helper, not a proof checker.  Unlike a
whole-phase move, a component move changes one prime-power digit of a row's
target while preserving every coprime digit.  The smoother neighborhood is
useful for building coverage margin before returning to the exact checker.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from pathlib import Path

import exact_greedy
import exact_uncovered
from local_phase_cegis import build_targets


def crt_pair(
    left_residue: int,
    left_modulus: int,
    right_residue: int,
    right_modulus: int,
) -> int:
    """Combine residues for coprime positive moduli."""
    if left_modulus == 1:
        return right_residue % right_modulus
    if right_modulus == 1:
        return left_residue % left_modulus
    multiplier = (
        (right_residue - left_residue)
        * pow(left_modulus, -1, right_modulus)
    ) % right_modulus
    return (
        left_residue + left_modulus * multiplier
    ) % (left_modulus * right_modulus)


def objective(cover, required_coverage: int, np) -> tuple[int, int]:
    return (
        int(np.count_nonzero(cover < required_coverage)),
        int(np.maximum(required_coverage - cover, 0).sum()),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path)
    parser.add_argument("--initial-phases", type=Path, required=True)
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--status-output", type=Path)
    parser.add_argument("--required-coverage", type=int, default=2)
    parser.add_argument("--sweeps", type=int, default=100)
    parser.add_argument("--seed", type=int, default=203_1616615)
    parser.add_argument(
        "--neutral-rate",
        type=float,
        default=0.01,
        help="probability of taking a different equal-score component target",
    )
    parser.add_argument("--progress-every", type=int, default=1)
    args = parser.parse_args()
    if args.required_coverage < 1:
        raise SystemExit("--required-coverage must be positive")
    if args.sweeps < 1:
        raise SystemExit("--sweeps must be positive")
    if not 0 <= args.neutral_rate <= 1:
        raise SystemExit("--neutral-rate must lie in [0,1]")
    if args.progress_every < 1:
        raise SystemExit("--progress-every must be positive")

    dependency_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dependency_path))
    import numpy as np  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    if any(int(row["target_modulus"]) != 1 for row in rows):
        raise RuntimeError(
            "component margin search currently requires target_modulus=1"
        )
    candidates = exact_greedy.load_candidates(args.pool, True)
    if len(candidates) != len(rows):
        raise AssertionError("candidate/metadata row count mismatch")
    points = [
        (int(k), int(l))
        for k, l in json.loads(args.points.read_text())
    ]
    if len(set(points)) != len(points):
        raise RuntimeError("points file contains duplicates")
    initial = {
        int(prime): int(value)
        for prime, value in json.loads(
            args.initial_phases.read_text()
        ).items()
    }
    assignment = np.asarray(
        [
            initial.get(int(row["p"]), int(row["target_residue"]))
            % int(row["h"])
            for row in rows
        ],
        dtype=np.uint32,
    )
    targets, matrix_seconds = build_targets(points, candidates, np)
    cover = np.count_nonzero(
        targets == assignment, axis=1
    ).astype(np.int32)
    current_objective = objective(cover, args.required_coverage, np)
    best_objective = current_objective
    best_assignment = assignment.copy()
    rng = random.Random(args.seed)
    component_powers = []
    for row in rows:
        h = int(row["h"])
        factors = exact_uncovered.factor(h)
        component_powers.append(
            [prime**exponent for prime, exponent in factors.items()]
        )

    print(
        f"rows={len(rows)} points={len(points)} "
        f"components={sum(map(len, component_powers))} "
        f"required={args.required_coverage} "
        f"below={current_objective[0]} deficit={current_objective[1]} "
        f"matrix_s={matrix_seconds:.3f}",
        flush=True,
    )
    started = time.monotonic()
    total_moves = 0
    improving_moves = 0
    neutral_moves = 0
    for sweep in range(1, args.sweeps + 1):
        row_order = list(range(len(rows)))
        rng.shuffle(row_order)
        sweep_moves = 0
        for row_index in row_order:
            h = int(rows[row_index]["h"])
            column = targets[:, row_index]
            powers = component_powers[row_index][:]
            rng.shuffle(powers)
            for component in powers:
                other = h // component
                old_phase = int(assignment[row_index])
                old_hits = column == old_phase
                base_cover = cover - old_hits
                needed = base_cover < args.required_coverage
                other_residue = old_phase % other if other > 1 else 0
                compatible = needed
                if other > 1:
                    compatible = compatible & (
                        column % other == other_residue
                    )
                values = (
                    column[compatible] % component
                ).astype(np.int64, copy=False)
                if not len(values):
                    continue
                unique, counts = np.unique(values, return_counts=True)
                best_score = int(counts.max())
                best_values = unique[counts == best_score]
                old_component = old_phase % component
                old_score = int(
                    np.count_nonzero(
                        compatible
                        & (column % component == old_component)
                    )
                )
                if best_score < old_score:
                    raise AssertionError("best component score regressed")
                alternatives = [
                    int(value)
                    for value in best_values
                    if int(value) != old_component
                ]
                if not alternatives:
                    continue
                improving = best_score > old_score
                if not improving and rng.random() >= args.neutral_rate:
                    continue
                new_component = rng.choice(alternatives)
                new_phase = crt_pair(
                    other_residue,
                    other,
                    new_component,
                    component,
                )
                new_hits = column == new_phase
                cover = base_cover + new_hits
                assignment[row_index] = new_phase
                sweep_moves += 1
                total_moves += 1
                if improving:
                    improving_moves += 1
                else:
                    neutral_moves += 1

        current_objective = objective(
            cover, args.required_coverage, np
        )
        if current_objective < best_objective:
            best_objective = current_objective
            best_assignment = assignment.copy()
            args.phase_output.write_text(
                json.dumps(
                    {
                        str(int(row["p"])): int(best_assignment[index])
                        for index, row in enumerate(rows)
                    }
                )
                + "\n"
            )
        if sweep % args.progress_every == 0 or not current_objective[0]:
            print(
                f"sweep={sweep} moves={sweep_moves} "
                f"below={current_objective[0]} "
                f"deficit={current_objective[1]} "
                f"best_below={best_objective[0]} "
                f"best_deficit={best_objective[1]} "
                f"elapsed_s={time.monotonic() - started:.3f}",
                flush=True,
            )
        if not current_objective[0]:
            best_objective = current_objective
            best_assignment = assignment.copy()
            break

    args.phase_output.write_text(
        json.dumps(
            {
                str(int(row["p"])): int(best_assignment[index])
                for index, row in enumerate(rows)
            }
        )
        + "\n"
    )
    status = {
        "pool": str(args.pool),
        "points": str(args.points),
        "initial_phases": str(args.initial_phases),
        "phase_output": str(args.phase_output),
        "row_count": len(rows),
        "point_count": len(points),
        "required_coverage": args.required_coverage,
        "sweeps_requested": args.sweeps,
        "best_below_required": best_objective[0],
        "best_total_deficit": best_objective[1],
        "total_moves": total_moves,
        "improving_moves": improving_moves,
        "neutral_moves": neutral_moves,
        "matrix_seconds": matrix_seconds,
        "search_seconds": time.monotonic() - started,
        "complete_finite_margin": best_objective[0] == 0,
    }
    if args.status_output:
        args.status_output.write_text(json.dumps(status, indent=2) + "\n")
    print(
        f"result={'MARGIN' if best_objective[0] == 0 else 'INCOMPLETE'} "
        f"best_below={best_objective[0]} "
        f"best_deficit={best_objective[1]}",
        flush=True,
    )
    return 0 if best_objective[0] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
