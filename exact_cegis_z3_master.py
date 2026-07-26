#!/usr/bin/env python3
"""Exact CEGIS with an integer-modular Z3 master for large power pools.

Each candidate prime receives one integer target.  A retained exponent pair
adds the exact constraint that at least one target equals the affine-line
value at that pair.  Hence master UNSAT is a valid finite-pool obstruction,
while checker UNSAT is a valid cover.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
import time
from pathlib import Path

import exact_uncovered
import exact_uncovered_z3
from power_anchor_capacity_lp import power_target_congruence

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_pool", type=Path)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument(
        "--period",
        type=int,
        default=0,
        help="if nonzero, retain only candidates whose h divides this period",
    )
    parser.add_argument("--fixed-witness", type=Path)
    parser.add_argument("--seed-points", type=int, default=20)
    parser.add_argument("--random-checks", type=int, default=20000)
    parser.add_argument("--random-batch", type=int, default=10)
    parser.add_argument("--checker-batch", type=int, default=10)
    parser.add_argument("--rounds", type=int, default=10000)
    parser.add_argument("--max-component", type=int, default=200000)
    parser.add_argument("--points-file", type=Path)
    parser.add_argument("--phase-file", type=Path)
    parser.add_argument("--progress-file", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.power < 1:
        raise SystemExit("--power must be positive")
    if args.seed_points < 0 or args.random_checks < 0:
        raise SystemExit("point counts must be nonnegative")
    if args.random_batch < 1 or args.checker_batch < 1:
        raise SystemExit("batch sizes must be positive")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import z3  # type: ignore

    source = json.loads(args.candidate_pool.read_text())
    candidates = []
    incompatible = 0
    seen_primes = set()
    for raw in source["choices"]:
        row = {
            key: int(raw[key])
            for key in ("h", "p", "a", "b", "ord2", "ord3")
        }
        h, p = row["h"], row["p"]
        if args.period and args.period % h:
            continue
        if p in seen_primes:
            raise RuntimeError(f"repeated prime {p}")
        seen_primes.add(p)
        if math.gcd(row["a"], row["b"], h) != 1:
            raise RuntimeError(f"non-surjective row for prime {p}")
        try:
            residue, modulus = power_target_congruence(h, p, args.power)
        except RuntimeError:
            incompatible += 1
            continue
        row["power_residue"] = residue
        row["power_modulus"] = modulus
        candidates.append(row)
    if not candidates:
        raise RuntimeError("no compatible candidates")

    period = math.lcm(*(row["h"] for row in candidates))
    algebraic_primes = tuple(
        prime for prime in exact_uncovered.factor(args.power) if prime % 2
    )
    sophie_germain = args.power % 4 == 0

    fixed_targets = {}
    if args.fixed_witness:
        witness = json.loads(args.fixed_witness.read_text())
        fixed_targets = {
            int(row["p"]): int(row["line"][2]) for row in witness["rows"]
        }
    by_prime = {row["p"]: row for row in candidates}
    for p, target in fixed_targets.items():
        if p not in by_prime:
            raise RuntimeError(f"fixed prime {p} is missing or incompatible")
        row = by_prime[p]
        if not 0 <= target < row["h"]:
            raise RuntimeError(f"fixed target for prime {p} is out of range")
        if target % row["power_modulus"] != row["power_residue"]:
            raise RuntimeError(f"fixed target for prime {p} is incompatible")

    print(
        f"period={period} candidates={len(candidates)} "
        f"incompatible={incompatible} "
        f"density={sum(1 / row['h'] for row in candidates):.12f} "
        f"algebraic_primes={algebraic_primes} "
        f"sophie_germain={sophie_germain}",
        flush=True,
    )

    solver = z3.Solver()
    solver.set(random_seed=203)
    target_vars = {}
    for row in candidates:
        p = row["p"]
        target_index = z3.Int(f"target_index_{p}")
        target_vars[p] = target_index
        target_count = (
            (row["h"] - 1 - row["power_residue"])
            // row["power_modulus"]
            + 1
        )
        solver.add(target_index >= 0, target_index < target_count)
        if p in fixed_targets:
            solver.add(
                target_index
                == (fixed_targets[p] - row["power_residue"])
                // row["power_modulus"]
            )

    def algebraically_covered(k: int, l: int) -> bool:
        if any(k % prime == 0 and l % prime == 0 for prime in algebraic_primes):
            return True
        return sophie_germain and k % 4 == 2 and l % 4 == 0

    points: set[tuple[int, int]] = set()

    def add_point(k: int, l: int) -> bool:
        k %= period
        l %= period
        if algebraically_covered(k, l) or (k, l) in points:
            return False
        points.add((k, l))
        covering = []
        for row in candidates:
            value = (row["a"] * k + row["b"] * l) % row["h"]
            if value % row["power_modulus"] != row["power_residue"]:
                continue
            covering.append(
                target_vars[row["p"]]
                == (value - row["power_residue"]) // row["power_modulus"]
            )
        if not covering:
            raise RuntimeError(f"no compatible line can cover point {(k, l)}")
        solver.add(z3.Or(*covering))
        return True

    if args.points_file and args.points_file.exists():
        for k, l in json.loads(args.points_file.read_text()):
            add_point(int(k), int(l))

    rng = random.Random(203)
    while len(points) < args.seed_points:
        add_point(rng.randrange(period), rng.randrange(period))

    for round_no in range(1, args.rounds + 1):
        master_started = time.monotonic()
        status = solver.check()
        master_seconds = time.monotonic() - master_started
        if status == z3.unsat:
            print(
                f"MASTER_UNSAT round={round_no} points={len(points)} "
                f"master_s={master_seconds:.3f}",
                flush=True,
            )
            return 2
        if status != z3.sat:
            raise RuntimeError(f"master returned {status}")

        model = solver.model()
        rows = []
        compact_rows = []
        phases = {}
        for row in candidates:
            target_index = model.eval(
                target_vars[row["p"]], model_completion=True
            ).as_long()
            c = (
                row["power_residue"]
                + row["power_modulus"] * target_index
            )
            assert 0 <= c < row["h"]
            assert c % row["power_modulus"] == row["power_residue"]
            output_row = {
                key: row[key]
                for key in ("h", "p", "a", "b", "ord2", "ord3")
            }
            output_row["c"] = c
            rows.append(output_row)
            compact_rows.append(
                (row["h"], row["a"], row["b"], c)
            )
            phases[str(row["p"])] = c
        if args.phase_file:
            args.phase_file.write_text(json.dumps(phases) + "\n")

        checker_started = time.monotonic()
        missed = []
        checked = 0
        while checked < args.random_checks and len(missed) < args.random_batch:
            point = (rng.randrange(period), rng.randrange(period))
            checked += 1
            if algebraically_covered(*point):
                continue
            if all(
                (a * point[0] + b * point[1] - c) % h
                for h, a, b, c in compact_rows
            ):
                missed.append(point)
        if missed:
            meta = {
                "engine": "random",
                "sat": True,
                "checks": checked,
                "rows": len(rows),
            }
        else:
            missed, meta = exact_uncovered_z3.find_uncovered(
                rows,
                max_component=args.max_component,
                limit=args.checker_batch,
                algebraic_primes=algebraic_primes,
                sophie_germain=sophie_germain,
            )
        checker_seconds = time.monotonic() - checker_started
        print(
            f"round={round_no} master_points={len(points)} "
            f"checker={'SAT' if missed else 'UNSAT'} new={len(missed)} "
            f"master_s={master_seconds:.3f} checker_s={checker_seconds:.3f}",
            flush=True,
        )

        if args.progress_file:
            args.progress_file.write_text(
                json.dumps(
                    {
                        "round": round_no,
                        "master_points": len(points),
                        "missed": missed,
                        "master_seconds": master_seconds,
                        "checker_seconds": checker_seconds,
                        "checker": meta,
                    },
                    indent=2,
                )
                + "\n"
            )
        if not missed:
            payload = {
                "period": period,
                "candidate_pool": str(args.candidate_pool),
                "choices": rows,
                "power": args.power,
                "algebraic_primes": algebraic_primes,
                "sophie_germain": sophie_germain,
                "fixed_targets": fixed_targets,
                "checker": meta,
            }
            if args.output:
                args.output.write_text(json.dumps(payload, indent=2) + "\n")
            else:
                print(json.dumps(payload, indent=2))
            return 0

        before = len(points)
        for k, l in missed:
            add_point(k, l)
        if len(points) == before:
            raise RuntimeError("checker returned no new counterexample")
        if args.points_file:
            args.points_file.write_text(json.dumps(sorted(points)) + "\n")

    print(f"ROUND_LIMIT points={len(points)}", flush=True)
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
