#!/usr/bin/env python3
"""Fast labelled-phase local search for one affine-plane cover."""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--seconds", type=float, default=120.0)
    parser.add_argument("--seed", type=int, default=203)
    parser.add_argument("--kick", type=int, default=4)
    parser.add_argument("--checkpoint-output", type=Path)
    parser.add_argument("--phase-output", type=Path)
    parser.add_argument("--result-output", type=Path)
    args = parser.parse_args()
    q = args.prime
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = [row for row in payload["choices"] if int(row["h"]) == q]
    if any(int(row["target_modulus"]) != 1 for row in rows):
        raise RuntimeError("local search requires target_modulus=1")
    incidence = []
    for row in rows:
        matrix = np.zeros((q, q * q), dtype=np.int16)
        a = int(row["a"]) % q
        b = int(row["b"]) % q
        for k in range(q):
            for l in range(q):
                matrix[(a * k + b * l) % q, k * q + l] = 1
        incidence.append(matrix)

    rng = random.Random(args.seed)
    phases = [int(row.get("c", 0)) % q for row in rows]
    counts = np.zeros(q * q, dtype=np.int16)

    def rebuild() -> None:
        counts.fill(0)
        for matrix, target in zip(incidence, phases):
            counts[matrix[target].astype(bool)] += 1

    def objective() -> int:
        return int(np.count_nonzero(counts == 0))

    def save(best: list[int], score: int, sweeps: int, kicks: int) -> None:
        if not args.checkpoint_output:
            return
        args.checkpoint_output.write_text(
            json.dumps(
                {
                    "pool": str(args.pool),
                    "prime": q,
                    "row_count": len(rows),
                    "best_uncovered": score,
                    "sweeps": sweeps,
                    "kicks": kicks,
                    "phases": {
                        str(int(row["p"])): best[index]
                        for index, row in enumerate(rows)
                    },
                },
                indent=2,
            )
            + "\n"
        )

    rebuild()
    best_score = objective()
    best_phases = list(phases)
    deadline = time.monotonic() + args.seconds
    sweeps = 0
    kicks = 0
    print(
        f"q={q} rows={len(rows)} initial_uncovered={best_score}",
        flush=True,
    )
    while time.monotonic() < deadline and best_score:
        before = objective()
        order = list(range(len(rows)))
        rng.shuffle(order)
        for index in order:
            matrix = incidence[index]
            counts[matrix[phases[index]].astype(bool)] -= 1
            scores = matrix @ (counts == 0).astype(np.int16)
            maximum = int(scores.max())
            choices = np.flatnonzero(scores == maximum).tolist()
            phases[index] = rng.choice(choices)
            counts[matrix[phases[index]].astype(bool)] += 1
        sweeps += 1
        current = objective()
        if current < best_score:
            best_score = current
            best_phases = list(phases)
            save(best_phases, best_score, sweeps, kicks)
            print(
                f"sweep={sweeps} kicks={kicks} "
                f"best_uncovered={best_score}",
                flush=True,
            )
        if current >= before:
            kicks += 1
            for index in rng.sample(
                range(len(rows)),
                min(max(1, args.kick), len(rows)),
            ):
                matrix = incidence[index]
                counts[matrix[phases[index]].astype(bool)] -= 1
                phases[index] = rng.randrange(q)
                counts[matrix[phases[index]].astype(bool)] += 1

    phases[:] = best_phases
    rebuild()
    uncovered = [
        [index // q, index % q]
        for index in np.flatnonzero(counts == 0).tolist()
    ]
    if len(uncovered) != best_score:
        raise AssertionError("best score did not reproduce")
    phase_map = {
        str(int(row["p"])): phases[index]
        for index, row in enumerate(rows)
    }
    result = {
        "pool": str(args.pool),
        "prime": q,
        "row_count": len(rows),
        "sat": not uncovered,
        "best_uncovered": len(uncovered),
        "uncovered_points": uncovered,
        "sweeps": sweeps,
        "kicks": kicks,
        "seed": args.seed,
    }
    if not uncovered and args.phase_output:
        args.phase_output.write_text(json.dumps(phase_map) + "\n")
    if args.result_output:
        args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    save(best_phases, best_score, sweeps, kicks)
    print(
        f"result={'SAT' if not uncovered else 'INCOMPLETE'} "
        f"best_uncovered={len(uncovered)} sweeps={sweeps} kicks={kicks}",
        flush=True,
    )
    return 0 if not uncovered else 2


if __name__ == "__main__":
    raise SystemExit(main())
