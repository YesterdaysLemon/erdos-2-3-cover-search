#!/usr/bin/env python3
"""Min-conflicts repair of affine phases on an explicit finite point set."""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from collections import deque
from pathlib import Path

from local_phase_cegis import build_targets


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path)
    parser.add_argument("--initial-phases", type=Path, required=True)
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument(
        "--period",
        type=int,
        default=0,
        help="if nonzero, retain only rows whose modulus divides this period",
    )
    parser.add_argument("--fixed-primes", default="")
    parser.add_argument("--max-steps", type=int, default=100_000)
    parser.add_argument("--noise", type=float, default=0.03)
    parser.add_argument("--tabu", type=int, default=7)
    parser.add_argument(
        "--breakout-every",
        type=int,
        default=0,
        help=(
            "after this many steps without a new best, increment the weight "
            "of every currently uncovered point; zero disables breakout"
        ),
    )
    parser.add_argument(
        "--breakout-increment",
        type=int,
        default=1,
        help="positive clause-weight increment used by breakout",
    )
    parser.add_argument(
        "--restart-every",
        type=int,
        default=0,
        help=(
            "after this many steps without a new best, restore the best "
            "phases and reset clause weights; zero disables restarts"
        ),
    )
    parser.add_argument(
        "--random-phase-rate",
        type=float,
        default=0.0,
        help=(
            "probability of a global random row/phase move instead of a "
            "move chosen to cover a current hole"
        ),
    )
    parser.add_argument("--seed", type=int, default=203)
    parser.add_argument("--progress-every", type=int, default=1_000)
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=10_000,
        help="write the best phase assignment every N steps",
    )
    args = parser.parse_args()
    if args.max_steps < 1:
        raise SystemExit("--max-steps must be positive")
    if not 0 <= args.noise <= 1:
        raise SystemExit("--noise must lie in [0,1]")
    if args.tabu < 0:
        raise SystemExit("--tabu must be nonnegative")
    if args.breakout_every < 0:
        raise SystemExit("--breakout-every must be nonnegative")
    if args.breakout_increment < 1:
        raise SystemExit("--breakout-increment must be positive")
    if args.restart_every < 0:
        raise SystemExit("--restart-every must be nonnegative")
    if not 0 <= args.random_phase_rate <= 1:
        raise SystemExit("--random-phase-rate must lie in [0,1]")
    if args.checkpoint_every < 1:
        raise SystemExit("--checkpoint-every must be positive")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    if args.period:
        rows = [row for row in rows if args.period % int(row["h"]) == 0]
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
        (int(k), int(l)) for k, l in json.loads(args.points.read_text())
    ]
    initial = {
        int(prime): int(target)
        for prime, target in json.loads(args.initial_phases.read_text()).items()
    }
    fixed_primes = {
        int(value) for value in args.fixed_primes.split(",") if value
    }
    phases = np.empty(len(rows), dtype=np.uint32)
    mutable = []
    for index, row in enumerate(rows):
        prime = int(row["p"])
        h = int(row["h"])
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        target = initial.get(prime, residue) % h
        if target % modulus != residue:
            raise RuntimeError(f"invalid initial phase for prime {prime}")
        phases[index] = target
        if prime not in fixed_primes:
            mutable.append(index)
    missing_fixed = fixed_primes - {int(row["p"]) for row in rows}
    if missing_fixed:
        raise RuntimeError(f"fixed primes absent from pool: {missing_fixed}")
    mutable_array = np.asarray(mutable, dtype=np.int32)

    targets, build_seconds = build_targets(points, candidates, np)
    matches = targets == phases
    cover = np.count_nonzero(matches, axis=1).astype(np.int16)
    point_weights = np.ones(len(points), dtype=np.int64)
    singleton_indices = np.flatnonzero(cover == 1)
    singleton_owners = np.argmax(matches[singleton_indices], axis=1)
    singleton_losses = np.bincount(
        singleton_owners,
        weights=point_weights[singleton_indices],
        minlength=len(rows),
    ).astype(np.int64)

    rng = random.Random(args.seed)
    tabu = deque(maxlen=args.tabu)
    best_phases = phases.copy()
    best_uncovered = int(np.count_nonzero(cover == 0))
    last_best_step = 0
    last_breakout_step = 0

    def save_best() -> None:
        phase_map = {
            str(int(row["p"])): int(best_phases[index])
            for index, row in enumerate(rows)
        }
        args.phase_output.write_text(json.dumps(phase_map) + "\n")

    print(
        f"rows={len(rows)} points={len(points)} mutable={len(mutable)} "
        f"matrix_s={build_seconds:.3f} initial_uncovered={best_uncovered} "
        f"initial_singletons={len(singleton_indices)}",
        flush=True,
    )

    started = time.monotonic()
    for step in range(1, args.max_steps + 1):
        uncovered = np.flatnonzero(cover == 0)
        current_uncovered = len(uncovered)
        if not current_uncovered:
            best_phases = phases.copy()
            best_uncovered = 0
            break

        if rng.random() < args.random_phase_rate:
            available_rows = [
                int(index)
                for index in mutable
                if int(index) not in tabu
            ]
            if not available_rows:
                tabu.clear()
                available_rows = list(mutable)
            row_index = available_rows[rng.randrange(len(available_rows))]
            row = rows[row_index]
            h = int(row["h"])
            residue = int(row["target_residue"])
            modulus = int(row["target_modulus"])
            target_count = h // modulus
            old_target = int(phases[row_index])
            if target_count <= 1:
                continue
            new_target = old_target
            while new_target == old_target:
                new_target = (
                    residue + modulus * rng.randrange(target_count)
                )
        else:
            point = int(uncovered[rng.randrange(current_uncovered)])
            proposed = targets[point, mutable_array]
            current = phases[mutable_array]
            changed = proposed != current
            if tabu:
                tabu_mask = np.isin(
                    mutable_array,
                    np.fromiter(tabu, dtype=np.int32),
                )
                if np.any(changed & ~tabu_mask):
                    changed &= ~tabu_mask
            choices = np.flatnonzero(changed)
            if not len(choices):
                tabu.clear()
                continue

            if rng.random() < args.noise:
                # Random-walk moves are biased toward low break-count phases.
                choice_losses = singleton_losses[mutable_array[choices]]
                cutoff = int(np.quantile(choice_losses, 0.25))
                low_break = choices[choice_losses <= cutoff]
                choice = int(low_break[rng.randrange(len(low_break))])
            else:
                uncovered_targets = targets[uncovered][
                    :, mutable_array[choices]
                ]
                gains = np.sum(
                    (uncovered_targets == proposed[choices])
                    * point_weights[uncovered, None],
                    axis=0,
                    dtype=np.int64,
                )
                losses = singleton_losses[mutable_array[choices]]
                scores = gains - losses
                maximum = int(scores.max())
                best_choices = choices[np.flatnonzero(scores == maximum)]
                choice = int(best_choices[rng.randrange(len(best_choices))])

            row_index = int(mutable_array[choice])
            old_target = int(phases[row_index])
            new_target = int(proposed[choice])
        if old_target == new_target:
            raise AssertionError("selected a null move")

        old_indices = np.flatnonzero(matches[:, row_index])
        new_indices = np.flatnonzero(targets[:, row_index] == new_target)
        old_cover = cover[old_indices].copy()
        new_cover = cover[new_indices].copy()

        # Remove the old line. A former double-covered point becomes a
        # singleton owned by its remaining matching row.
        matches[old_indices, row_index] = False
        singleton_losses[row_index] -= int(
            point_weights[old_indices[old_cover == 1]].sum()
        )
        became_singleton = old_indices[old_cover == 2]
        if len(became_singleton):
            owners = np.argmax(matches[became_singleton], axis=1)
            singleton_losses += np.bincount(
                owners,
                weights=point_weights[became_singleton],
                minlength=len(rows),
            ).astype(np.int64)
        cover[old_indices] -= 1

        # Add the new line. A former singleton ceases to be uniquely owned;
        # a former hole becomes a singleton owned by the moved row.
        lost_singleton = new_indices[new_cover == 1]
        if len(lost_singleton):
            owners = np.argmax(matches[lost_singleton], axis=1)
            singleton_losses -= np.bincount(
                owners,
                weights=point_weights[lost_singleton],
                minlength=len(rows),
            ).astype(np.int64)
        singleton_losses[row_index] += int(
            point_weights[new_indices[new_cover == 0]].sum()
        )
        matches[new_indices, row_index] = True
        cover[new_indices] += 1
        phases[row_index] = new_target
        tabu.append(row_index)

        current_uncovered = int(np.count_nonzero(cover == 0))
        if current_uncovered < best_uncovered:
            best_uncovered = current_uncovered
            best_phases = phases.copy()
            last_best_step = step
            print(
                f"step={step} best_uncovered={best_uncovered} "
                f"singletons={int(singleton_losses.sum())} "
                f"elapsed_s={time.monotonic() - started:.3f}",
                flush=True,
            )
        if (
            args.breakout_every
            and current_uncovered
            and step - last_best_step >= args.breakout_every
            and step - last_breakout_step >= args.breakout_every
        ):
            breakout_holes = np.flatnonzero(cover == 0)
            point_weights[breakout_holes] += args.breakout_increment
            last_breakout_step = step
            print(
                f"step={step} breakout={len(breakout_holes)} "
                f"max_weight={int(point_weights.max())}",
                flush=True,
            )
        if (
            args.restart_every
            and step - last_best_step >= args.restart_every
        ):
            phases = best_phases.copy()
            matches = targets == phases
            cover = np.count_nonzero(matches, axis=1).astype(np.int16)
            point_weights.fill(1)
            singleton_indices = np.flatnonzero(cover == 1)
            singleton_owners = np.argmax(
                matches[singleton_indices], axis=1
            )
            singleton_losses = np.bincount(
                singleton_owners,
                weights=point_weights[singleton_indices],
                minlength=len(rows),
            ).astype(np.int64)
            tabu.clear()
            current_uncovered = best_uncovered
            last_best_step = step
            last_breakout_step = step
            print(
                f"step={step} restart best_uncovered={best_uncovered}",
                flush=True,
            )
        elif step % args.progress_every == 0:
            print(
                f"step={step} current_uncovered={current_uncovered} "
                f"best_uncovered={best_uncovered} "
                f"elapsed_s={time.monotonic() - started:.3f}",
                flush=True,
            )
        if step % args.checkpoint_every == 0 or best_uncovered == 0:
            save_best()

        if np.any(singleton_losses < 0):
            raise AssertionError("negative singleton loss count")
        if int(singleton_losses.sum()) != int(
            point_weights[cover == 1].sum()
        ):
            raise AssertionError("singleton accounting drifted")

    save_best()

    reproduced = np.count_nonzero(targets == best_phases, axis=1)
    reproduced_uncovered = int(np.count_nonzero(reproduced == 0))
    if reproduced_uncovered != best_uncovered:
        raise AssertionError("saved best phase count did not reproduce")
    print(
        f"result={'SAT' if best_uncovered == 0 else 'INCOMPLETE'} "
        f"best_uncovered={best_uncovered} output={args.phase_output}",
        flush=True,
    )
    return 0 if best_uncovered == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
