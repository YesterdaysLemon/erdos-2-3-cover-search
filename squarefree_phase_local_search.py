#!/usr/bin/env python3
"""Stochastic coordinate search for a two-component fibre cover.

This complements the exact SAT encoding in finite_squarefree_component_cover:
it searches labelled prime-fibre phases directly while exploiting the
q1^2-by-q2^2 product structure for fast exact objective updates.  A reported
zero is exhaustively checked on the whole finite grid before it is written.
"""

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
    parser.add_argument("--primes", required=True)
    parser.add_argument("--seconds", type=float, default=600.0)
    parser.add_argument("--seed", type=int, default=203)
    parser.add_argument("--kick", type=int, default=5)
    parser.add_argument("--checkpoint-output", type=Path)
    parser.add_argument("--phase-output", type=Path)
    parser.add_argument("--result-output", type=Path)
    args = parser.parse_args()

    primes = tuple(int(value) for value in args.primes.split(",") if value)
    if len(primes) != 2 or len(set(primes)) != 2:
        raise SystemExit("--primes must contain exactly two distinct primes")
    q1, q2 = primes
    period = q1 * q2
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = [
        row
        for row in payload["choices"]
        if int(row["h"]) in (q1, q2, period)
    ]
    if not rows:
        raise RuntimeError("no rows retained")
    for row in rows:
        if int(row["target_modulus"]) != 1:
            raise RuntimeError("local search requires target_modulus=1")

    def line_matrix(row: dict, q: int):
        matrix = np.zeros((q, q * q), dtype=np.int16)
        a = int(row["a"]) % q
        b = int(row["b"]) % q
        for k in range(q):
            for l in range(q):
                matrix[(a * k + b * l) % q, k * q + l] = 1
        return matrix

    matrices: dict[tuple[int, int, int], object] = {}
    for row in rows:
        h = int(row["h"])
        for q in primes:
            if h % q == 0:
                key = (q, int(row["a"]) % q, int(row["b"]) % q)
                if key not in matrices:
                    matrices[key] = line_matrix(row, q)

    def matrix_for(row: dict, q: int):
        return matrices[(q, int(row["a"]) % q, int(row["b"]) % q)]

    rng = random.Random(args.seed)
    phases: list[tuple[int, int | None]] = []
    for row in rows:
        h = int(row["h"])
        saved = int(row.get("c", 0))
        if h == q1:
            phases.append((saved % q1, None))
        elif h == q2:
            phases.append((saved % q2, None))
        else:
            phases.append((saved % q1, saved % q2))

    counts = np.zeros((q1 * q1, q2 * q2), dtype=np.int16)

    def change_count(
        row: dict,
        phase: tuple[int, int | None],
        delta: int,
    ) -> None:
        h = int(row["h"])
        if h == q1:
            mask1 = matrix_for(row, q1)[phase[0]].astype(bool)
            counts[mask1, :] += delta
        elif h == q2:
            assert phase[1] is None
            mask2 = matrix_for(row, q2)[phase[0]].astype(bool)
            counts[:, mask2] += delta
        else:
            assert phase[1] is not None
            mask1 = matrix_for(row, q1)[phase[0]].astype(bool)
            mask2 = matrix_for(row, q2)[phase[1]].astype(bool)
            counts[np.ix_(mask1, mask2)] += delta

    def rebuild() -> None:
        counts.fill(0)
        for row, phase in zip(rows, phases):
            change_count(row, phase, 1)

    def objective() -> int:
        return int(np.count_nonzero(counts == 0))

    def best_phase(index: int) -> tuple[int, int | None]:
        row = rows[index]
        old = phases[index]
        change_count(row, old, -1)
        zeros = counts == 0
        h = int(row["h"])
        if h == q1:
            per_point = zeros.sum(axis=1, dtype=np.int64)
            scores = matrix_for(row, q1) @ per_point
            maximum = int(scores.max())
            targets = np.flatnonzero(scores == maximum).tolist()
            chosen = (rng.choice(targets), None)
        elif h == q2:
            per_point = zeros.sum(axis=0, dtype=np.int64)
            scores = matrix_for(row, q2) @ per_point
            maximum = int(scores.max())
            targets = np.flatnonzero(scores == maximum).tolist()
            chosen = (rng.choice(targets), None)
        else:
            scores = (
                matrix_for(row, q1)
                @ zeros.astype(np.int16)
                @ matrix_for(row, q2).T
            )
            maximum = int(scores.max())
            targets = np.argwhere(scores == maximum)
            selected = targets[rng.randrange(len(targets))]
            chosen = (int(selected[0]), int(selected[1]))
        change_count(row, chosen, 1)
        return chosen

    def random_phase(row: dict) -> tuple[int, int | None]:
        h = int(row["h"])
        if h == q1:
            return (rng.randrange(q1), None)
        if h == q2:
            return (rng.randrange(q2), None)
        return (rng.randrange(q1), rng.randrange(q2))

    def crt_pair(t1: int, t2: int) -> int:
        return (
            t1
            + q1
            * (((t2 - t1) * pow(q1, -1, q2)) % q2)
        ) % period

    def phase_payload(values: list[tuple[int, int | None]]) -> dict[str, int]:
        result: dict[str, int] = {}
        for row, phase in zip(rows, values):
            h = int(row["h"])
            if h == period:
                assert phase[1] is not None
                c = crt_pair(phase[0], phase[1])
            else:
                c = phase[0]
            result[str(int(row["p"]))] = c
        return result

    def save_checkpoint(
        best_values: list[tuple[int, int | None]],
        best_objective: int,
        sweeps: int,
        kicks: int,
    ) -> None:
        if not args.checkpoint_output:
            return
        checkpoint = {
            "pool": str(args.pool),
            "primes": list(primes),
            "period": period,
            "row_count": len(rows),
            "best_uncovered": best_objective,
            "sweeps": sweeps,
            "kicks": kicks,
            "phases": phase_payload(best_values),
        }
        args.checkpoint_output.write_text(
            json.dumps(checkpoint, indent=2) + "\n"
        )

    rebuild()
    best_objective = objective()
    best_values = list(phases)
    started = time.monotonic()
    deadline = started + args.seconds
    sweeps = 0
    kicks = 0
    print(
        f"primes={primes} period={period} rows={len(rows)} "
        f"initial_uncovered={best_objective}",
        flush=True,
    )
    while time.monotonic() < deadline and best_objective:
        before = objective()
        order = list(range(len(rows)))
        rng.shuffle(order)
        for index in order:
            phases[index] = best_phase(index)
        sweeps += 1
        current = objective()
        if current < best_objective:
            best_objective = current
            best_values = list(phases)
            save_checkpoint(best_values, best_objective, sweeps, kicks)
            print(
                f"sweep={sweeps} kicks={kicks} "
                f"best_uncovered={best_objective}",
                flush=True,
            )
        if current >= before:
            kicks += 1
            kick_count = min(max(1, args.kick), len(rows))
            for index in rng.sample(range(len(rows)), kick_count):
                change_count(rows[index], phases[index], -1)
                phases[index] = random_phase(rows[index])
                change_count(rows[index], phases[index], 1)

    elapsed = time.monotonic() - started
    phases[:] = best_values
    rebuild()
    exhaustive_uncovered = [
        [index1, index2]
        for index1 in range(q1 * q1)
        for index2 in range(q2 * q2)
        if counts[index1, index2] == 0
    ]
    if len(exhaustive_uncovered) != best_objective:
        raise AssertionError("best checkpoint did not reproduce")
    phase_map = phase_payload(best_values)
    result = {
        "pool": str(args.pool),
        "primes": list(primes),
        "period": period,
        "row_count": len(rows),
        "sat": best_objective == 0,
        "best_uncovered": best_objective,
        "uncovered_points": exhaustive_uncovered,
        "sweeps": sweeps,
        "kicks": kicks,
        "seed": args.seed,
        "elapsed_seconds": elapsed,
    }
    if best_objective == 0 and args.phase_output:
        args.phase_output.write_text(json.dumps(phase_map) + "\n")
    if args.result_output:
        args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    save_checkpoint(best_values, best_objective, sweeps, kicks)
    print(
        f"result={'SAT' if best_objective == 0 else 'INCOMPLETE'} "
        f"best_uncovered={best_objective} sweeps={sweeps} "
        f"kicks={kicks} elapsed_s={elapsed:.3f}",
        flush=True,
    )
    return 0 if best_objective == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
