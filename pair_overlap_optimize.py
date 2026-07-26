#!/usr/bin/env python3
"""Reduce phase-dependent pair intersections of affine prime fibres.

For two fibres, the joint target is either unattainable or their intersection
has density d/(h1*h2), where d is the index of the joint image.  The baseline
1/(h1*h2) is phase-independent; this optimizer minimizes the excess
(d-1)/(h1*h2) over correlated pairs.  Unlike sampled point repair, the
objective is exact on the entire exponent lattice.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import time
from pathlib import Path

from power_anchor_capacity_lp import power_target_congruence


def joint_index(left: dict, right: dict) -> int:
    h1, a1, b1 = left["h"], left["a"], left["b"]
    h2, a2, b2 = right["h"], right["a"], right["b"]
    return math.gcd(
        a1 * b2 - a2 * b1,
        a2 * h1,
        b2 * h1,
        a1 * h2,
        b1 * h2,
        h1 * h2,
    )


def compatible(left: dict, c1: int, right: dict, c2: int, d: int) -> bool:
    h1, a1, b1 = left["h"], left["a"], left["b"]
    h2, a2, b2 = right["h"], right["a"], right["b"]
    return (
        (c1 * a2 - c2 * a1) % d == 0
        and (c1 * b2 - c2 * b1) % d == 0
        and (c2 * h1) % d == 0
        and (c1 * h2) % d == 0
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=2000)
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--candidate-sample", type=int, default=16)
    parser.add_argument("--option-sample", type=int, default=32)
    parser.add_argument("--seed", type=int, default=203)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    saved = json.loads(args.phase_file.read_text())
    rows = []
    for raw in source["choices"]:
        row = {
            key: int(raw[key]) for key in ("h", "p", "a", "b")
        }
        if args.period % row["h"]:
            continue
        try:
            residue, modulus = power_target_congruence(
                row["h"], row["p"], args.power
            )
        except RuntimeError:
            continue
        row["residue"] = residue
        row["modulus"] = modulus
        row["target_count"] = (
            (row["h"] - 1 - residue) // modulus + 1
        )
        row["c"] = int(saved.get(str(row["p"]), residue)) % row["h"]
        if row["c"] % modulus != residue:
            raise RuntimeError(f"saved phase for {row['p']} is incompatible")
        rows.append(row)
    rows.sort(key=lambda row: (row["h"], row["p"]))
    active = rows[: min(args.limit, len(rows))]

    started = time.monotonic()
    adjacency: list[list[tuple[int, int, float]]] = [
        [] for _ in active
    ]
    correlated = 0
    for i, left in enumerate(active):
        for j in range(i):
            right = active[j]
            if math.gcd(left["h"], right["h"]) == 1:
                continue
            d = joint_index(left, right)
            if d == 1:
                continue
            weight = (d - 1) / (left["h"] * right["h"])
            adjacency[i].append((j, d, weight))
            adjacency[j].append((i, d, weight))
            correlated += 1

    phases = [row["c"] for row in active]
    costs = [0.0] * len(active)
    objective = 0.0
    for i, edges in enumerate(adjacency):
        for j, d, weight in edges:
            if j >= i:
                continue
            if compatible(active[i], phases[i], active[j], phases[j], d):
                costs[i] += weight
                costs[j] += weight
                objective += weight
    initial_objective = objective
    print(
        f"rows={len(rows)} active={len(active)} correlated={correlated} "
        f"initial_excess={initial_objective:.15f} "
        f"build_s={time.monotonic() - started:.3f}",
        flush=True,
    )

    rng = random.Random(args.seed)
    improvements = 0
    for step in range(1, args.steps + 1):
        sample = rng.sample(
            range(len(active)),
            min(args.candidate_sample, len(active)),
        )
        index = max(sample, key=costs.__getitem__)
        row = active[index]
        old_target = phases[index]

        options = {old_target}
        for _ in range(args.option_sample):
            options.add(
                row["residue"]
                + row["modulus"] * rng.randrange(row["target_count"])
            )
        best_target = old_target
        best_cost = costs[index]
        for target in options:
            candidate_cost = 0.0
            for neighbor, d, weight in adjacency[index]:
                if compatible(
                    row,
                    target,
                    active[neighbor],
                    phases[neighbor],
                    d,
                ):
                    candidate_cost += weight
            if candidate_cost < best_cost - 1e-18:
                best_cost = candidate_cost
                best_target = target
        if best_target != old_target:
            for neighbor, d, weight in adjacency[index]:
                before = compatible(
                    row,
                    old_target,
                    active[neighbor],
                    phases[neighbor],
                    d,
                )
                after = compatible(
                    row,
                    best_target,
                    active[neighbor],
                    phases[neighbor],
                    d,
                )
                if before == after:
                    continue
                delta = weight if after else -weight
                objective += delta
                costs[neighbor] += delta
            phases[index] = best_target
            costs[index] = best_cost
            improvements += 1
        if step % 1000 == 0:
            print(
                f"step={step} excess={objective:.15f} "
                f"improvements={improvements}",
                flush=True,
            )

    output = dict(saved)
    for row, target in zip(active, phases):
        output[str(row["p"])] = target
    args.output.write_text(json.dumps(output) + "\n")
    print(
        f"final_excess={objective:.15f} "
        f"reduction={initial_objective - objective:.15f} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
