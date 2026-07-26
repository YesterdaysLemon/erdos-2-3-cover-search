#!/usr/bin/env python3
"""Incrementally decide finite affine-cover point sets with native cardinality."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path)
    parser.add_argument("--batch-size", type=int, default=25_000)
    parser.add_argument(
        "--warm-prefix",
        type=int,
        default=0,
        help=(
            "add and solve this many leading points as one initial batch; "
            "use only with a phase hint already known to satisfy the prefix"
        ),
    )
    parser.add_argument(
        "--assume-hint-prefix",
        action="store_true",
        help=(
            "solve the warm prefix once under the complete initial phase "
            "assignment, then release those assumptions for later batches"
        ),
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=1,
        help="write status and phases every N solved batches",
    )
    parser.add_argument(
        "--print-every",
        type=int,
        default=1,
        help="print progress every N solved batches",
    )
    parser.add_argument("--solver", default="mc")
    parser.add_argument("--fixed-targets", default="")
    parser.add_argument("--initial-phases", type=Path)
    parser.add_argument("--phase-output", type=Path)
    parser.add_argument("--status-output", type=Path)
    args = parser.parse_args()
    if args.batch_size < 1:
        raise SystemExit("--batch-size must be positive")
    if args.warm_prefix < 0:
        raise SystemExit("--warm-prefix must be nonnegative")
    if args.checkpoint_every < 1:
        raise SystemExit("--checkpoint-every must be positive")
    if args.print_every < 1:
        raise SystemExit("--print-every must be positive")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.solvers import Solver  # type: ignore

    load_started = time.monotonic()
    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    points = [(int(k), int(l)) for k, l in json.loads(args.points.read_text())]
    print(
        f"stage=input_loaded rows={len(rows)} points={len(points)} "
        f"seconds={time.monotonic() - load_started:.3f}",
        flush=True,
    )
    if len(set(points)) != len(points):
        raise RuntimeError("points file contains duplicates")
    print(
        f"stage=points_unique seconds={time.monotonic() - load_started:.3f}",
        flush=True,
    )
    if args.warm_prefix > len(points):
        raise SystemExit("--warm-prefix exceeds the point count")
    fixed_targets = {}
    for item in args.fixed_targets.split(","):
        if not item:
            continue
        prime, target = item.split(":", 1)
        fixed_targets[int(prime)] = int(target)

    # For each row store (h, p, a, b, residue, modulus, base variable, count).
    specs = []
    next_variable = 1
    for row in rows:
        h = int(row["h"])
        prime = int(row["p"])
        a = int(row["a"])
        b = int(row["b"])
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        fixed = fixed_targets.get(prime)
        if fixed is not None:
            if not 0 <= fixed < h or fixed % modulus != residue:
                raise RuntimeError(f"invalid fixed target {prime}:{fixed}")
            count = 1
        else:
            count = h // modulus
        specs.append((h, prime, a, b, residue, modulus, next_variable, count, fixed))
        next_variable += count
    option_count = next_variable - 1
    print(
        f"stage=specs_built options={option_count} "
        f"seconds={time.monotonic() - load_started:.3f}",
        flush=True,
    )

    solver = Solver(name=args.solver)
    build_started = time.monotonic()
    for h, prime, a, b, residue, modulus, base, count, fixed in specs:
        variables = list(range(base, base + count))
        if count > 1:
            solver.add_atmost(variables, 1)
        solver.add_clause(variables)

    phase_hints = []
    phase_hint_count = 0
    if args.initial_phases:
        initial = {
            int(prime): int(target)
            for prime, target in json.loads(args.initial_phases.read_text()).items()
        }
        for h, prime, a, b, residue, modulus, base, count, fixed in specs:
            target = fixed if fixed is not None else initial.get(prime)
            if target is None or target % modulus != residue:
                continue
            slot = 0 if fixed is not None else (target - residue) // modulus
            if 0 <= slot < count:
                phase_hints.append(base + slot)
        try:
            solver.set_phases(phase_hints)
            phase_hint_count = len(phase_hints)
        except NotImplementedError:
            phase_hint_count = 0
    build_seconds = time.monotonic() - build_started
    status = {
        "pool": str(args.pool),
        "points": str(args.points),
        "solver": args.solver,
        "fixed_targets": {
            str(prime): target for prime, target in fixed_targets.items()
        },
        "row_count": len(rows),
        "option_count": option_count,
        "point_count": len(points),
        "batch_size": args.batch_size,
        "warm_prefix": args.warm_prefix,
        "assume_hint_prefix": args.assume_hint_prefix,
        "checkpoint_every": args.checkpoint_every,
        "print_every": args.print_every,
        "phase_hint_count": phase_hint_count,
        "base_build_seconds": build_seconds,
        "batches": [],
        "complete": False,
        "sat": None,
    }

    def write_status() -> None:
        if args.status_output:
            args.status_output.write_text(json.dumps(status, indent=2) + "\n")

    print(
        f"rows={len(rows)} options={option_count} points={len(points)} "
        f"native_atmost={sum(spec[7] > 1 for spec in specs)} "
        f"phase_hints={phase_hint_count} base_build_s={build_seconds:.3f}",
        flush=True,
    )
    boundaries = []
    if args.warm_prefix:
        boundaries.append((0, args.warm_prefix))
    boundaries.extend(
        (start, min(start + args.batch_size, len(points)))
        for start in range(args.warm_prefix, len(points), args.batch_size)
    )
    for batch_index, (start, stop) in enumerate(boundaries, start=1):
        add_started = time.monotonic()
        for k, l in points[start:stop]:
            clause = []
            for h, prime, a, b, residue, modulus, base, count, fixed in specs:
                target = (a * k + b * l) % h
                if fixed is not None:
                    if target == fixed:
                        clause.append(base)
                elif target % modulus == residue:
                    clause.append(base + (target - residue) // modulus)
            if not clause:
                raise RuntimeError(f"point {(k, l)} has no valid row")
            solver.add_clause(clause)
        add_seconds = time.monotonic() - add_started
        solve_started = time.monotonic()
        used_hint_assumptions = (
            args.assume_hint_prefix
            and bool(args.warm_prefix)
            and batch_index == 1
            and phase_hint_count == len(specs)
        )
        sat = solver.solve(
            assumptions=phase_hints if used_hint_assumptions else []
        )
        hint_assumptions_satisfied = bool(sat) if used_hint_assumptions else None
        if used_hint_assumptions and not sat:
            # A stale hint is not an obstruction to the unrestricted master.
            sat = solver.solve()
        solve_seconds = time.monotonic() - solve_started
        record = {
            "points": stop,
            "added_points": stop - start,
            "add_seconds": add_seconds,
            "sat": bool(sat),
            "solve_seconds": solve_seconds,
            "used_hint_assumptions": used_hint_assumptions,
            "hint_assumptions_satisfied": hint_assumptions_satisfied,
        }
        status["batches"].append(record)
        status["sat"] = bool(sat)
        terminal = (not sat) or stop == len(points)
        checkpoint = terminal or batch_index % args.checkpoint_every == 0
        if checkpoint:
            write_status()
        if terminal or batch_index % args.print_every == 0:
            print(
                f"points={stop} added={stop-start} add_s={add_seconds:.3f} "
                f"result={'SAT' if sat else 'UNSAT'} solve_s={solve_seconds:.3f}",
                flush=True,
            )
        if not sat:
            status["complete"] = True
            status["unsat_prefix_points"] = stop
            write_status()
            solver.delete()
            return 2

        if args.phase_output and checkpoint:
            positive = {literal for literal in solver.get_model() if literal > 0}
            phases = {}
            for h, prime, a, b, residue, modulus, base, count, fixed in specs:
                slot = next(
                    slot for slot in range(count) if base + slot in positive
                )
                phases[str(prime)] = (
                    fixed if fixed is not None else residue + modulus * slot
                )
            args.phase_output.write_text(json.dumps(phases) + "\n")

    status["complete"] = True
    write_status()
    solver.delete()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
