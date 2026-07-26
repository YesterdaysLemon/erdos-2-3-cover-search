#!/usr/bin/env python3
"""Persistent coordinate-repair CEGIS on a fixed finite point universe."""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from pathlib import Path

from local_phase_cegis import build_targets


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("universe", type=Path)
    parser.add_argument("--initial-phases", type=Path, required=True)
    parser.add_argument("--learned-points", type=Path, required=True)
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--status-output", type=Path, required=True)
    parser.add_argument(
        "--learned-output",
        type=Path,
        help="checkpoint the accumulated learned point values",
    )
    parser.add_argument("--fixed-primes", default="")
    parser.add_argument("--rounds", type=int, default=100)
    parser.add_argument("--sweeps", type=int, default=100)
    parser.add_argument("--zero-move-rate", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=203)
    args = parser.parse_args()
    if args.rounds < 1 or args.sweeps < 1:
        raise SystemExit("--rounds and --sweeps must be positive")
    if not 0 <= args.zero_move_rate <= 1:
        raise SystemExit("--zero-move-rate must lie in [0,1]")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    candidates = [
        (
            int(row["h"]),
            int(row["p"]),
            int(row["a"]),
            int(row["b"]),
            int(row.get("ord2", row["h"])),
            int(row.get("ord3", row["h"])),
        )
        for row in rows
    ]
    points = [
        (int(k), int(l)) for k, l in json.loads(args.universe.read_text())
    ]
    if len(set(points)) != len(points):
        raise RuntimeError("universe contains duplicate points")
    learned_values = {
        (int(k), int(l))
        for k, l in json.loads(args.learned_points.read_text())
    }
    learned_indices = [
        index for index, point in enumerate(points) if point in learned_values
    ]
    if len(learned_indices) != len(learned_values):
        raise RuntimeError("a learned point is absent from the universe")

    saved = {
        int(prime): int(target)
        for prime, target in json.loads(args.initial_phases.read_text()).items()
    }
    fixed_primes = {
        int(value) for value in args.fixed_primes.split(",") if value
    }
    assignment = np.empty(len(rows), dtype=np.uint32)
    mutable = []
    for index, row in enumerate(rows):
        prime = int(row["p"])
        h = int(row["h"])
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        target = saved.get(prime, residue) % h
        if target % modulus != residue:
            raise RuntimeError(f"invalid saved phase for prime {prime}")
        assignment[index] = target
        if prime not in fixed_primes:
            mutable.append(index)
    missing_fixed = fixed_primes - {int(row["p"]) for row in rows}
    if missing_fixed:
        raise RuntimeError(f"fixed primes absent from pool: {missing_fixed}")

    started = time.monotonic()
    targets, build_seconds = build_targets(points, candidates, np)
    full_cover = np.zeros(len(points), dtype=np.int16)
    for index in range(len(rows)):
        full_cover += targets[:, index] == assignment[index]

    rng = random.Random(args.seed)
    status = {
        "pool": str(args.pool),
        "universe": str(args.universe),
        "initial_phases": str(args.initial_phases),
        "initial_learned_points": len(learned_indices),
        "point_count": len(points),
        "row_count": len(rows),
        "build_seconds": build_seconds,
        "rounds": [],
        "complete": False,
        "sample_cover": False,
        "stuck": False,
    }

    def save() -> None:
        args.phase_output.write_text(
            json.dumps(
                {
                    str(int(row["p"])): int(assignment[index])
                    for index, row in enumerate(rows)
                }
            )
            + "\n"
        )
        args.status_output.write_text(json.dumps(status, indent=2) + "\n")
        if args.learned_output:
            args.learned_output.write_text(
                json.dumps(
                    [
                        [points[index][0], points[index][1]]
                        for index in learned_indices
                    ]
                )
                + "\n"
            )

    print(
        f"rows={len(rows)} universe={len(points)} "
        f"learned={len(learned_indices)} matrix_s={build_seconds:.3f}",
        flush=True,
    )
    for round_no in range(1, args.rounds + 1):
        learned_array = np.asarray(learned_indices, dtype=np.int64)
        learned_cover = full_cover[learned_array].copy()
        moves = 0
        zero_moves = 0
        sweeps_used = 0
        for sweep in range(1, args.sweeps + 1):
            holes_before = int(np.count_nonzero(learned_cover == 0))
            if not holes_before:
                break
            order = list(mutable)
            rng.shuffle(order)
            sweep_moves = 0
            for row_index in order:
                values = targets[learned_array, row_index]
                current = int(assignment[row_index])
                current_mask = values == current
                learned_cover -= current_mask
                uncovered = learned_cover == 0
                observed = values[uncovered]
                row = rows[row_index]
                residue = int(row["target_residue"])
                modulus = int(row["target_modulus"])
                if modulus != 1:
                    observed = observed[observed % modulus == residue]
                if not len(observed):
                    learned_cover += current_mask
                    continue
                observed_targets, observed_counts = np.unique(
                    observed, return_counts=True
                )
                maximum = int(observed_counts.max())
                best_targets = observed_targets[observed_counts == maximum]
                target = int(best_targets[rng.randrange(len(best_targets))])
                current_gain = int(np.count_nonzero(observed == current))
                delta = maximum - current_gain
                accept = delta > 0 or (
                    delta == 0
                    and target != current
                    and rng.random() < args.zero_move_rate
                )
                if not accept:
                    learned_cover += current_mask
                    continue

                old_full = targets[:, row_index] == current
                new_full = targets[:, row_index] == target
                full_cover -= old_full
                full_cover += new_full
                learned_cover += values == target
                assignment[row_index] = target
                moves += 1
                sweep_moves += 1
                zero_moves += int(delta == 0)

            sweeps_used = sweep
            holes_after = int(np.count_nonzero(learned_cover == 0))
            print(
                f"round={round_no} sweep={sweep} "
                f"learned_holes={holes_after} moves={sweep_moves}",
                flush=True,
            )
            if holes_after == 0:
                break
            if sweep_moves == 0 or holes_after >= holes_before and not zero_moves:
                break

        learned_holes = int(np.count_nonzero(learned_cover == 0))
        miss_indices = np.flatnonzero(full_cover == 0)
        learned_set = set(learned_indices)
        new_indices = [
            int(index)
            for index in miss_indices
            if int(index) not in learned_set
        ]
        record = {
            "round": round_no,
            "learned_points": len(learned_indices),
            "learned_holes": learned_holes,
            "full_misses": int(len(miss_indices)),
            "new_misses": len(new_indices),
            "moves": moves,
            "zero_moves": zero_moves,
            "sweeps": sweeps_used,
            "elapsed_seconds": time.monotonic() - started,
        }
        status["rounds"].append(record)
        print(
            f"round={round_no} learned={len(learned_indices)} "
            f"learned_holes={learned_holes} full_misses={len(miss_indices)} "
            f"new={len(new_indices)} moves={moves}",
            flush=True,
        )
        save()
        if learned_holes:
            status["complete"] = True
            status["stuck"] = True
            save()
            return 3
        if not len(miss_indices):
            status["complete"] = True
            status["sample_cover"] = True
            save()
            return 0
        if not new_indices:
            raise AssertionError("all full misses were already learned")
        learned_indices.extend(new_indices)

    save()
    return 4


if __name__ == "__main__":
    raise SystemExit(main())
