#!/usr/bin/env python3
"""Cross-validated greedy assignment of residual affine-line fibres."""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from pathlib import Path

import exact_greedy

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_pool", type=Path)
    parser.add_argument("--derived-pool", action="store_true")
    parser.add_argument("--base-count", type=int, required=True)
    parser.add_argument("--initial-phase-file", type=Path, required=True)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument("--sample-size", type=int, default=100000)
    parser.add_argument("--min-split-gain", type=int, default=1)
    parser.add_argument("--seed", type=int, default=203407)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    candidates = exact_greedy.load_candidates(args.candidate_pool, args.derived_pool)
    if not 0 < args.base_count <= len(candidates):
        raise SystemExit("invalid --base-count")
    saved = json.loads(args.initial_phase_file.read_text())
    assignment = np.asarray(
        [int(saved.get(str(p), 0)) % h for h, p, *_rest in candidates],
        dtype=np.uint32,
    )

    rng = random.Random(args.seed)
    # A 63-bit uniform box is indistinguishable from uniform modulo every
    # h<=300000 here, while allowing fast vectorized modular arithmetic.
    ks = np.asarray([rng.getrandbits(63) for _ in range(args.sample_size)], dtype=np.uint64)
    ls = np.asarray([rng.getrandbits(63) for _ in range(args.sample_size)], dtype=np.uint64)
    split = args.sample_size // 2
    covered = np.zeros(args.sample_size, dtype=np.bool_)
    started = time.monotonic()
    for index, (h, _p, a, b, _o2, _o3) in enumerate(candidates[: args.base_count]):
        values = (a * (ks % h) + b * (ls % h)) % h
        covered |= values == assignment[index]
    print(
        f"candidates={len(candidates)} base={args.base_count} sample={args.sample_size} "
        f"base_missed={int(np.count_nonzero(~covered))} base_s={time.monotonic()-started:.3f}",
        flush=True,
    )

    selected = 0
    selected_gain = 0
    for index in range(args.base_count, len(candidates)):
        uncovered_indices = np.flatnonzero(~covered)
        if not len(uncovered_indices):
            break
        h, _p, a, b, _o2, _o3 = candidates[index]
        values = (
            a * (ks[uncovered_indices] % h)
            + b * (ls[uncovered_indices] % h)
        ) % h
        train_values = values[uncovered_indices < split]
        valid_values = values[uncovered_indices >= split]
        train_counts = np.bincount(train_values, minlength=h)
        valid_counts = np.bincount(valid_values, minlength=h)
        eligible = (
            (train_counts >= args.min_split_gain)
            & (valid_counts >= args.min_split_gain)
        )
        if not np.any(eligible):
            assignment[index] = 0
            continue
        scores = train_counts + valid_counts
        scores[~eligible] = 0
        target = int(np.argmax(scores))
        gain = int(scores[target])
        assignment[index] = target
        covered[uncovered_indices[values == target]] = True
        selected += 1
        selected_gain += gain
        if selected % 250 == 0:
            print(
                f"selected={selected} scanned={index+1-args.base_count} "
                f"missed={int(np.count_nonzero(~covered))} gain={selected_gain}",
                flush=True,
            )

    args.phase_file.write_text(
        json.dumps(
            {
                str(item[1]): int(assignment[index])
                for index, item in enumerate(candidates)
            }
        )
        + "\n"
    )

    # Include the mandatory target-zero lines that were not selected by the
    # validated greedy pass, then report both halves separately.
    final_covered = np.zeros(args.sample_size, dtype=np.bool_)
    for index, (h, _p, a, b, _o2, _o3) in enumerate(candidates):
        values = (a * (ks % h) + b * (ls % h)) % h
        final_covered |= values == assignment[index]
    train_missed = int(np.count_nonzero(~final_covered[:split]))
    valid_missed = int(np.count_nonzero(~final_covered[split:]))
    print(
        f"DONE selected={selected} selected_gain={selected_gain} "
        f"greedy_missed={int(np.count_nonzero(~covered))} "
        f"final_train_missed={train_missed} final_valid_missed={valid_missed} "
        f"output={args.phase_file}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
