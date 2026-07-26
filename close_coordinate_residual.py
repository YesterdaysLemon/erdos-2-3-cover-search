#!/usr/bin/env python3
"""Close a sparse ordinary coordinate-grid residual with certified leaf cores.

Closed grid cells are protected by the union of their extracted primary cores.
For the first still-open cell, the script either retargets a mutable row that
covers the whole cell or runs local CEGIS, minimizes the resulting leaf core,
lifts it, and dual-verifies the lift.  The grid is then recertified before the
next cell, so no heuristic result is accepted as progress.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import exact_uncovered
import exact_uncovered_z3_bv


def read_phases(path: Path) -> dict[int, int]:
    return {
        int(prime): int(target)
        for prime, target in json.loads(path.read_text()).items()
    }


def write_phases(path: Path, phases: dict[int, int]) -> None:
    path.write_text(
        json.dumps(
            {str(prime): target for prime, target in phases.items()}
        )
        + "\n"
    )


def read_primes(path: Path) -> set[int]:
    payload = json.loads(path.read_text())
    raw = payload.get("primes", ()) if isinstance(payload, dict) else payload
    return {int(prime) for prime in raw}


def rows_with_phases(pool: dict, phases: dict[int, int]) -> list[dict]:
    rows = []
    for raw in pool["choices"]:
        row = dict(raw)
        row["c"] = phases[int(row["p"])] % int(row["h"])
        rows.append(row)
    return rows


def materialize_cover(
    pool_path: Path,
    phase_path: Path,
    output_path: Path,
) -> None:
    """Write the ordinary exact-cover payload represented by a phase file."""
    pool = json.loads(pool_path.read_text())
    phases = read_phases(phase_path)
    choices = []
    for raw in pool["choices"]:
        row = dict(raw)
        prime = int(row["p"])
        if prime not in phases:
            raise RuntimeError(
                f"recursive phase file omits prime {prime}"
            )
        row["c"] = phases[prime] % int(row["h"])
        choices.append(row)
    output_path.write_text(
        json.dumps(
            {
                "candidate_pool": str(pool_path),
                "candidate_count": len(choices),
                "choices": choices,
                "algebraic_primes": pool.get("algebraic_primes", []),
                "sophie_germain": bool(pool.get("sophie_germain", False)),
            },
            indent=2,
        )
        + "\n"
    )


def classify_grid(
    pool: dict,
    phases: dict[int, int],
    period: int,
    max_component: int,
    solver: str,
) -> tuple[list[list[int]], list[int], dict]:
    rows = rows_with_phases(pool, phases)
    cells = tuple(
        ((period, k, l),)
        for k in range(period)
        for l in range(period)
    )
    _misses, meta = exact_uncovered.find_uncovered(
        rows,
        max_component=max_component,
        limit=1,
        algebraic_primes=tuple(
            int(prime) for prime in pool.get("algebraic_primes", ())
        ),
        sophie_germain=bool(pool.get("sophie_germain", False)),
        solver_name=solver,
        core_coordinate_cells=cells,
    )
    records = meta["core_coordinate_cells"]
    open_cells = [
        [int(record["cell"][0][1]), int(record["cell"][0][2])]
        for record in records
        if not record["closed"]
    ]
    union_indices = sorted(
        {
            int(index)
            for record in records
            if record["closed"]
            for index in record["row_indices"]
        }
    )
    return (
        open_cells,
        [int(rows[index]["p"]) for index in union_indices],
        meta,
    )


def run_logged(command: list[str], log_path: Path) -> int:
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
            log.flush()
        return process.wait()


def verify_rows_on_cell(
    rows: list[dict],
    pool: dict,
    period: int,
    k_residue: int,
    l_residue: int,
    max_component: int,
) -> list[dict]:
    fixed = ((period, k_residue, l_residue),)
    checks = []
    for name, engine in (
        ("pysat", exact_uncovered),
        ("z3bv", exact_uncovered_z3_bv),
    ):
        misses, meta = engine.find_uncovered(
            rows,
            max_component=max_component,
            limit=1,
            algebraic_primes=tuple(
                int(prime) for prime in pool.get("algebraic_primes", ())
            ),
            sophie_germain=bool(pool.get("sophie_germain", False)),
            fixed_coordinate_residues=fixed,
        )
        if misses or meta["sat"]:
            raise RuntimeError(
                f"{name} failed lifted cell {(k_residue, l_residue)}"
            )
        checks.append({"checker": name, "sat": False})
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("parent_pool", type=Path)
    parser.add_argument("initial_phases", type=Path)
    parser.add_argument("immutable_primes", type=Path)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--max-component", type=int, default=1024)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--attempts", type=int, default=4)
    parser.add_argument("--rounds", type=int, default=16)
    parser.add_argument("--seed", type=int, default=270000)
    parser.add_argument(
        "--diversity-schedule",
        default="16;23;29;31;none",
    )
    parser.add_argument(
        "--recursive-periods",
        default="",
        help=(
            "comma-separated coordinate periods used recursively when an "
            "open cell has no mutable full-cell row"
        ),
    )
    parser.add_argument(
        "--density-guard-primes",
        default="",
        help=(
            "comma-separated residual primes whose exact density-at-least-one "
            "condition must survive every direct full-cell retarget"
        ),
    )
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--locks-output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--certificate-output", type=Path, required=True)
    args = parser.parse_args()
    if args.period < 1:
        raise SystemExit("--period must be positive")
    recursive_periods = tuple(
        int(value)
        for value in args.recursive_periods.split(",")
        if value
    )
    if any(period < 2 for period in recursive_periods):
        raise SystemExit("--recursive-periods values must be at least two")
    density_guard_primes = tuple(
        int(value)
        for value in args.density_guard_primes.split(",")
        if value
    )
    if any(prime < 2 for prime in density_guard_primes):
        raise SystemExit("--density-guard-primes values must be at least two")

    pool = json.loads(args.parent_pool.read_text())
    phases = read_phases(args.initial_phases)
    immutable = read_primes(args.immutable_primes)
    initial_immutable_phases = {
        prime: phases[prime]
        for prime in immutable
        if prime in phases
    }
    manifest = {
        "parent_pool": str(args.parent_pool),
        "initial_phases": str(args.initial_phases),
        "immutable_primes": str(args.immutable_primes),
        "period": args.period,
        "iterations": [],
        "complete": False,
    }
    write_phases(args.phase_output, phases)
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")

    iteration = 0
    while True:
        iteration += 1
        grid_started = time.monotonic()
        open_cells, core_primes, grid_meta = classify_grid(
            pool,
            phases,
            args.period,
            args.max_component,
            args.solver,
        )
        locks = immutable | set(core_primes)
        args.locks_output.write_text(
            json.dumps(sorted(locks), indent=2) + "\n"
        )
        grid_path = Path(f"{args.prefix}_grid_i{iteration:02d}.json")
        grid_path.write_text(
            json.dumps(
                {
                    "phase_file": str(args.phase_output),
                    "period": args.period,
                    "open_cells": open_cells,
                    "closed_cell_count": args.period**2 - len(open_cells),
                    "union_primes": core_primes,
                    "checker": grid_meta,
                    "elapsed_seconds": time.monotonic() - grid_started,
                },
                indent=2,
            )
            + "\n"
        )
        print(
            f"GRID iteration={iteration} open={len(open_cells)} "
            f"closed={args.period**2 - len(open_cells)} "
            f"locks={len(locks)}",
            flush=True,
        )
        if not open_cells:
            break

        k_residue, l_residue = open_cells[0]
        tag = (
            f"{args.prefix}_i{iteration:02d}_"
            f"{k_residue}_{l_residue}"
        )
        child_pool_path = Path(f"{tag}_pool.json")
        child_initial_path = Path(f"{tag}_initial_phases.json")
        seed_path = Path(f"{tag}_parent_seed.json")
        points_path = Path(f"{tag}_points.json")
        seed_path.write_text(
            json.dumps([[k_residue, l_residue]]) + "\n"
        )
        condition_command = [
            sys.executable,
            "condition_derived_cell.py",
            str(args.parent_pool),
            "--cell-period",
            str(args.period),
            "--k-residue",
            str(k_residue),
            "--l-residue",
            str(l_residue),
            "--phase-file",
            str(args.phase_output),
            "--phase-output",
            str(child_initial_path),
            "--fixed-primes",
            ",".join(map(str, sorted(locks))),
            "--output",
            str(child_pool_path),
        ]
        subprocess.run(condition_command, check=True)
        subprocess.run(
            [
                sys.executable,
                "transform_cell_points.py",
                str(child_pool_path),
                str(seed_path),
                "--output",
                str(points_path),
            ],
            check=True,
        )
        child_pool = json.loads(child_pool_path.read_text())
        mutable_full_cell = [
            int(prime)
            for prime in child_pool.get("full_cell_primes", ())
            if int(prime) not in locks
        ]
        record = {
            "iteration": iteration,
            "cell": [k_residue, l_residue],
            "grid_before": str(grid_path),
        }

        full_cell_choice = None
        for prime in sorted(mutable_full_cell):
            parent_row = next(
                row
                for row in pool["choices"]
                if int(row["p"]) == prime
            )
            h = int(parent_row["h"])
            target = (
                int(parent_row["a"]) * k_residue
                + int(parent_row["b"]) * l_residue
            ) % h
            if (
                target % int(parent_row["target_modulus"])
                != int(parent_row["target_residue"])
            ):
                raise AssertionError("full-cell target violates restriction")
            tentative_phases = dict(phases)
            tentative_phases[prime] = target
            guard_records = {}
            guarded = True
            if density_guard_primes:
                from component_density_cegis import (
                    find_low_density_cells,
                )

                for guard_prime in density_guard_primes:
                    violations, guard_meta = find_low_density_cells(
                        pool["choices"],
                        tentative_phases,
                        guard_prime,
                        1,
                        args.max_component,
                    )
                    guard_records[str(guard_prime)] = guard_meta
                    if violations:
                        guarded = False
                        break
            if guarded:
                full_cell_choice = (
                    prime,
                    target,
                    parent_row,
                    guard_records,
                )
                break
            print(
                f"REJECT_FULL_CELL cell={(k_residue, l_residue)} "
                f"prime={prime} density_guard_failed",
                flush=True,
            )

        if full_cell_choice is not None:
            prime, target, parent_row, guard_records = full_cell_choice
            phases[prime] = target
            lifted_row = dict(parent_row)
            lifted_row["c"] = target
            verification = verify_rows_on_cell(
                [lifted_row],
                pool,
                args.period,
                k_residue,
                l_residue,
                args.max_component,
            )
            core_path = Path(f"{tag}_parent_core.json")
            core_path.write_text(
                json.dumps(
                    {
                        "cell": [args.period, k_residue, l_residue],
                        "rows": [lifted_row],
                        "primes": [prime],
                        "verification": verification,
                    },
                    indent=2,
                )
                + "\n"
            )
            record.update(
                {
                    "method": "full_cell_prime",
                    "prime": prime,
                    "target": target,
                    "parent_core": str(core_path),
                    "density_guards": guard_records,
                }
            )
            print(
                f"FULL_CELL cell={(k_residue, l_residue)} "
                f"prime={prime} target={target}",
                flush=True,
            )
        else:
            present = {
                int(row["p"]) for row in child_pool["choices"]
            }
            active_locks = sorted(locks & present)
            cover_path = None
            working_path = None
            method = "cegis"
            method_details = {}
            if recursive_periods:
                recursive_period = recursive_periods[0]
                remaining_periods = recursive_periods[1:]
                recursive_prefix = f"{tag}_r{recursive_period}"
                recursive_immutable = Path(
                    f"{recursive_prefix}_immutable_primes.json"
                )
                recursive_immutable.write_text(
                    json.dumps(sorted(locks), indent=2) + "\n"
                )
                working_path = Path(
                    f"{recursive_prefix}_phases.json"
                )
                recursive_locks = Path(
                    f"{recursive_prefix}_locks.json"
                )
                recursive_manifest = Path(
                    f"{recursive_prefix}_manifest.json"
                )
                recursive_certificate = Path(
                    f"{recursive_prefix}_certificate.json"
                )
                recursive_log = Path(f"{recursive_prefix}.log")
                command = [
                    sys.executable,
                    "close_coordinate_residual.py",
                    str(child_pool_path),
                    str(child_initial_path),
                    str(recursive_immutable),
                    "--period",
                    str(recursive_period),
                    "--prefix",
                    recursive_prefix,
                    "--max-component",
                    str(args.max_component),
                    "--solver",
                    args.solver,
                    "--attempts",
                    str(args.attempts),
                    "--rounds",
                    str(args.rounds),
                    "--seed",
                    str(args.seed + 1000 * iteration),
                    "--diversity-schedule",
                    args.diversity_schedule,
                    "--recursive-periods",
                    ",".join(map(str, remaining_periods)),
                    "--density-guard-primes",
                    ",".join(map(str, density_guard_primes)),
                    "--phase-output",
                    str(working_path),
                    "--locks-output",
                    str(recursive_locks),
                    "--manifest",
                    str(recursive_manifest),
                    "--certificate-output",
                    str(recursive_certificate),
                ]
                print(
                    f"RECURSE cell={(k_residue, l_residue)} "
                    f"period={recursive_period} "
                    f"remaining={remaining_periods}",
                    flush=True,
                )
                return_code = run_logged(command, recursive_log)
                if return_code != 0 or not recursive_certificate.exists():
                    raise RuntimeError(
                        "recursive coordinate close failed for "
                        f"{(k_residue, l_residue)} at period "
                        f"{recursive_period}"
                    )
                recursive_result = json.loads(
                    recursive_certificate.read_text()
                )
                if not recursive_result.get("complete"):
                    raise RuntimeError(
                        "recursive coordinate certificate is incomplete"
                    )
                cover_path = Path(f"{recursive_prefix}_cover.json")
                materialize_cover(
                    child_pool_path,
                    working_path,
                    cover_path,
                )
                method = "recursive_coordinate"
                method_details = {
                    "recursive_period": recursive_period,
                    "recursive_certificate": str(
                        recursive_certificate
                    ),
                }
            else:
                for attempt in range(1, args.attempts + 1):
                    working_path = Path(
                        f"{tag}_a{attempt}_working_phases.json"
                    )
                    cover_path = Path(f"{tag}_a{attempt}_cover.json")
                    checker_path = Path(f"{tag}_a{attempt}_checker.json")
                    log_path = Path(f"{tag}_a{attempt}.log")
                    command = [
                        sys.executable,
                        "local_phase_cegis.py",
                        str(child_pool_path),
                        "--derived-pool",
                        "--derived-targets",
                        "--points-file",
                        str(points_path),
                        "--phase-file",
                        str(working_path),
                        "--initial-phase-file",
                        str(child_initial_path),
                        "--output",
                        str(cover_path),
                        "--checker-checkpoint-file",
                        str(checker_path),
                        "--rounds",
                        str(args.rounds),
                        "--seed",
                        str(args.seed + 100 * iteration + attempt),
                        "--randomize-mutable-initial",
                        "--required-coverage",
                        "1",
                        "--random-checks",
                        "0",
                        "--exact-batch",
                        "5000",
                        "--diversity-coordinate-schedule",
                        args.diversity_schedule,
                        "--diversity-stage-rounds",
                        "2",
                        "--diversity-quota",
                        "2",
                        "--max-component",
                        str(args.max_component),
                        "--retain-target-matrix-during-exact",
                    ]
                    if active_locks:
                        insert_at = command.index("--output")
                        command[insert_at:insert_at] = [
                            "--fixed-primes",
                            ",".join(map(str, active_locks)),
                        ]
                    print(
                        f"CEGIS cell={(k_residue, l_residue)} "
                        f"attempt={attempt} "
                        f"active_locks={len(active_locks)}",
                        flush=True,
                    )
                    return_code = run_logged(command, log_path)
                    if (
                        return_code == 0
                        and cover_path.exists()
                        and working_path.exists()
                    ):
                        break
                else:
                    raise RuntimeError(
                        "no leaf cover found for "
                        f"{(k_residue, l_residue)}"
                    )
            assert cover_path is not None and working_path is not None

            min_core_path = Path(f"{tag}_min_core.json")
            subprocess.run(
                [
                    sys.executable,
                    "minimize_cover_core.py",
                    str(cover_path),
                    "--max-component",
                    str(args.max_component),
                    "--output",
                    str(min_core_path),
                ],
                check=True,
            )
            next_phase_path = Path(f"{tag}_lifted_phases.json")
            parent_core_path = Path(f"{tag}_parent_core.json")
            subprocess.run(
                [
                    sys.executable,
                    "lift_derived_conditioned_phases.py",
                    str(args.parent_pool),
                    str(child_pool_path),
                    str(cover_path),
                    str(min_core_path),
                    str(args.phase_output),
                    "--phase-output",
                    str(next_phase_path),
                    "--rows-output",
                    str(parent_core_path),
                ],
                check=True,
            )
            next_phases = read_phases(next_phase_path)
            changed_locks = [
                prime
                for prime in locks
                if prime in phases
                and next_phases[prime] != phases[prime]
            ]
            if changed_locks:
                raise RuntimeError(
                    f"protected phases changed: {changed_locks[:5]}"
                )
            parent_core = json.loads(parent_core_path.read_text())
            verification = verify_rows_on_cell(
                parent_core["rows"],
                pool,
                args.period,
                k_residue,
                l_residue,
                args.max_component,
            )
            phases = next_phases
            record.update(
                {
                    "method": method,
                    "cover": str(cover_path),
                    "min_core": str(min_core_path),
                    "parent_core": str(parent_core_path),
                    "verification": verification,
                    **method_details,
                }
            )

        changed_immutable = [
            prime
            for prime, target in initial_immutable_phases.items()
            if phases[prime] != target
        ]
        if changed_immutable:
            raise RuntimeError(
                f"immutable phases changed: {changed_immutable[:5]}"
            )
        write_phases(args.phase_output, phases)
        manifest["iterations"].append(record)
        args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")

    final_rows = rows_with_phases(pool, phases)
    final_checks = {}
    for name, engine in (
        ("pysat", exact_uncovered),
        ("z3bv", exact_uncovered_z3_bv),
    ):
        misses, meta = engine.find_uncovered(
            final_rows,
            max_component=args.max_component,
            limit=1,
            algebraic_primes=tuple(
                int(prime) for prime in pool.get("algebraic_primes", ())
            ),
            sophie_germain=bool(pool.get("sophie_germain", False)),
        )
        if misses or meta["sat"]:
            raise RuntimeError(f"final {name} global verification failed")
        final_checks[name] = meta
    certificate = {
        "parent_pool": str(args.parent_pool),
        "phase_file": str(args.phase_output),
        "period": args.period,
        "manifest": str(args.manifest),
        "grid_iterations": iteration,
        "immutable_primes": sorted(immutable),
        "final_verification": final_checks,
        "complete": True,
    }
    args.certificate_output.write_text(
        json.dumps(certificate, indent=2) + "\n"
    )
    manifest["complete"] = True
    manifest["certificate"] = str(args.certificate_output)
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"COMPLETE parent={args.parent_pool} iterations={iteration}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
