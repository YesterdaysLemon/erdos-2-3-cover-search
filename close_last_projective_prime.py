#!/usr/bin/env python3
"""Close every remaining direction for the last algebraic prime.

The parent pool must have exactly one algebraic prime q.  Conditioning an
open projective q-direction removes that final algebraic component, producing
an ordinary leaf cover problem.  Every accepted leaf is minimized, lifted,
and verified on all q-1 nonzero scalar cells by both PySAT and Z3-BV.
Existing protected phases are never changed.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import exact_uncovered
import exact_uncovered_z3_bv


def normalize_direction(k: int, l: int, q: int) -> tuple[int, int]:
    k %= q
    l %= q
    if k:
        inverse = pow(k, -1, q)
        return 1, l * inverse % q
    if l:
        return 0, 1
    raise ValueError("the origin has no projective direction")


def rows_with_phases(pool: dict, phases: dict[int, int]) -> list[dict]:
    rows = []
    for raw in pool["choices"]:
        row = dict(raw)
        row["c"] = phases[int(row["p"])]
        rows.append(row)
    return rows


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
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
            log.flush()
        return process.wait()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("parent_pool", type=Path)
    parser.add_argument("initial_phases", type=Path)
    parser.add_argument("locked_primes", type=Path)
    parser.add_argument("closed_core_bundle", type=Path)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--max-component", type=int, default=1024)
    parser.add_argument("--attempts", type=int, default=3)
    parser.add_argument("--rounds", type=int, default=12)
    parser.add_argument("--seed", type=int, default=719000)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--locks-output", type=Path, required=True)
    parser.add_argument("--certificate-output", type=Path, required=True)
    args = parser.parse_args()

    pool = json.loads(args.parent_pool.read_text())
    algebraic = tuple(int(p) for p in pool["algebraic_primes"])
    if len(algebraic) != 1:
        raise SystemExit(
            f"expected exactly one algebraic prime, got {algebraic}"
        )
    q = algebraic[0]
    phases = {
        int(prime): int(value)
        for prime, value in json.loads(args.initial_phases.read_text()).items()
    }
    lock_payload = json.loads(args.locked_primes.read_text())
    raw_locked = (
        lock_payload.get("primes", ())
        if isinstance(lock_payload, dict)
        else lock_payload
    )
    locked = {
        int(prime)
        for prime in raw_locked
    }
    closed_bundle = json.loads(args.closed_core_bundle.read_text())
    certified_primes = {
        int(prime) for prime in closed_bundle["union_primes"]
    }
    if not certified_primes <= locked:
        raise RuntimeError("closed-core primes are not all protected")
    initial_locked_phases = {p: phases[p] for p in locked}
    manifest = {
        "parent_pool": str(args.parent_pool),
        "initial_phases": str(args.initial_phases),
        "locked_primes": str(args.locked_primes),
        "closed_core_bundle": str(args.closed_core_bundle),
        "projective_prime": q,
        "directions": [],
        "complete": False,
    }
    args.phase_output.write_text(
        json.dumps(
            {str(prime): value for prime, value in phases.items()}
        )
        + "\n"
    )
    args.locks_output.write_text(
        json.dumps(sorted(locked), indent=2) + "\n"
    )
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")

    iteration = 0
    while True:
        parent_rows = rows_with_phases(pool, phases)
        witnesses, _ = exact_uncovered.find_uncovered(
            parent_rows,
            max_component=args.max_component,
            limit=q * q,
            algebraic_primes=algebraic,
            diversity_coordinate_moduli=(q,),
            diversity_quota=1,
        )
        directions = sorted(
            {
                normalize_direction(k, l, q)
                for k, l in witnesses
            }
        )
        if not directions:
            break
        iteration += 1
        direction = directions[0]
        representative = next(
            (k, l)
            for k, l in witnesses
            if (k % q, l % q) == direction
        )
        tag = (
            f"{args.prefix}_d{iteration:02d}_"
            f"{direction[0]}_{direction[1]}"
        )
        child_pool_path = Path(f"{tag}_pool.json")
        child_initial_path = Path(f"{tag}_initial_phases.json")
        seed_path = Path(f"{tag}_parent_seed.json")
        points_path = Path(f"{tag}_points.json")
        seed_path.write_text(json.dumps([list(representative)]) + "\n")
        subprocess.run(
            [
                sys.executable,
                "condition_derived_cell.py",
                str(args.parent_pool),
                "--cell-period",
                str(q),
                "--k-residue",
                str(direction[0]),
                "--l-residue",
                str(direction[1]),
                "--phase-file",
                str(args.phase_output),
                "--phase-output",
                str(child_initial_path),
                "--fixed-primes",
                ",".join(map(str, sorted(locked))),
                "--canonicalize-algebraic",
                "--output",
                str(child_pool_path),
            ],
            check=True,
        )
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
        if child_pool["algebraic_primes"]:
            raise RuntimeError("conditioning did not remove the last prime")
        present = {int(row["p"]) for row in child_pool["choices"]}
        active_locks = sorted(locked & present)

        cover_path = None
        working_path = None
        for attempt in range(1, args.attempts + 1):
            working_path = Path(f"{tag}_a{attempt}_working_phases.json")
            cover_path = Path(f"{tag}_a{attempt}_cover.json")
            checker_path = Path(f"{tag}_a{attempt}_checker.json")
            log_path = Path(f"{tag}_a{attempt}.log")
            shutil.copyfile(child_initial_path, working_path)
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
                "--fixed-primes",
                ",".join(map(str, active_locks)),
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
                "2",
                "--random-checks",
                "0",
                "--exact-batch",
                "5000",
                "--diversity-coordinate-schedule",
                "16;16,27;none",
                "--diversity-stage-rounds",
                "2",
                "--diversity-quota",
                "2",
                "--max-component",
                str(args.max_component),
                "--retain-target-matrix-during-exact",
            ]
            print(
                f"DIRECTION iteration={iteration} direction={direction} "
                f"attempt={attempt}",
                flush=True,
            )
            return_code = run_logged(command, log_path)
            if return_code == 0 and cover_path.exists():
                break
        else:
            raise RuntimeError(
                f"no leaf cover found for direction {direction}"
            )
        assert cover_path is not None and working_path is not None

        cover = json.loads(cover_path.read_text())
        cover_rows = cover["choices"]
        misses, core_meta = exact_uncovered.find_uncovered(
            cover_rows,
            max_component=args.max_component,
            limit=1,
            algebraic_primes=(),
            core_coordinate_cells=((),),
        )
        if misses or core_meta["sat"]:
            raise RuntimeError("primary leaf core extraction was not UNSAT")
        kept = list(
            core_meta["core_coordinate_cells"][0]["row_indices"]
        )
        raw_core_size = len(kept)
        for index in list(kept):
            trial = [value for value in kept if value != index]
            trial_rows = [cover_rows[value] for value in trial]
            zmisses, _ = exact_uncovered_z3_bv.find_uncovered(
                trial_rows,
                max_component=args.max_component,
                limit=1,
                algebraic_primes=(),
            )
            if not zmisses:
                kept = trial
        minimized_rows = [cover_rows[index] for index in kept]
        pmisses, pmeta = exact_uncovered.find_uncovered(
            minimized_rows,
            max_component=args.max_component,
            limit=1,
            algebraic_primes=(),
        )
        zmisses, zmeta = exact_uncovered_z3_bv.find_uncovered(
            minimized_rows,
            max_component=args.max_component,
            limit=1,
            algebraic_primes=(),
        )
        if pmisses or zmisses:
            raise RuntimeError("minimized leaf core failed verification")
        core_path = Path(f"{tag}_min_core.json")
        core_path.write_text(
            json.dumps(
                {
                    "cover": str(cover_path),
                    "row_indices": kept,
                    "primes": [
                        int(cover_rows[index]["p"]) for index in kept
                    ],
                    "raw_core_size": raw_core_size,
                    "primary_verification": pmeta,
                    "independent_verification": zmeta,
                },
                indent=2,
            )
            + "\n"
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
                str(core_path),
                str(args.phase_output),
                "--phase-output",
                str(next_phase_path),
                "--rows-output",
                str(parent_core_path),
            ],
            check=True,
        )
        next_phases = {
            int(prime): int(value)
            for prime, value in json.loads(
                next_phase_path.read_text()
            ).items()
        }
        changed_locked = [
            prime
            for prime in locked
            if next_phases[prime] != phases[prime]
        ]
        if changed_locked:
            raise RuntimeError(
                f"protected phases changed: {changed_locked[:5]}"
            )
        parent_core = json.loads(parent_core_path.read_text())
        lifted_rows = parent_core["rows"]
        scalar_checks = []
        for scalar in range(1, q):
            kr = scalar * direction[0] % q
            lr = scalar * direction[1] % q
            fixed = ((q, kr, lr),)
            for name, engine in (
                ("pysat", exact_uncovered),
                ("z3bv", exact_uncovered_z3_bv),
            ):
                check_misses, check_meta = engine.find_uncovered(
                    lifted_rows,
                    max_component=args.max_component,
                    limit=1,
                    algebraic_primes=algebraic,
                    fixed_coordinate_residues=fixed,
                )
                if check_misses or check_meta["sat"]:
                    raise RuntimeError(
                        f"{name} failed lifted scalar {(kr, lr)}"
                    )
                scalar_checks.append(
                    {
                        "cell": [kr, lr],
                        "checker": name,
                        "meta": check_meta,
                    }
                )
        scalar_path = Path(f"{tag}_scalar_verification.json")
        scalar_path.write_text(
            json.dumps(
                {
                    "parent_core": str(parent_core_path),
                    "checks": scalar_checks,
                },
                indent=2,
            )
            + "\n"
        )

        new_primes = {
            int(prime) for prime in parent_core["lifted_phases"]
        }
        certified_primes.update(new_primes)
        locked.update(new_primes)
        phases = next_phases
        args.phase_output.write_text(
            json.dumps(
                {str(prime): value for prime, value in phases.items()}
            )
            + "\n"
        )
        args.locks_output.write_text(
            json.dumps(sorted(locked), indent=2) + "\n"
        )
        if any(phases[p] != initial_locked_phases[p] for p in initial_locked_phases):
            raise RuntimeError("an initial protected phase changed")
        direction_record = {
            "iteration": iteration,
            "direction": list(direction),
            "cover": str(cover_path),
            "min_core": str(core_path),
            "parent_core": str(parent_core_path),
            "scalar_verification": str(scalar_path),
            "raw_core_size": raw_core_size,
            "min_core_size": len(kept),
            "new_primes": sorted(new_primes),
            "remaining_before": len(directions),
        }
        manifest["directions"].append(direction_record)
        args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
        print(
            f"DIRECTION_CLOSED direction={direction} "
            f"core={len(kept)} remaining_before={len(directions)}",
            flush=True,
        )

    final_rows = rows_with_phases(pool, phases)
    row_by_prime = {int(row["p"]): row for row in final_rows}
    certificate_rows = [
        row_by_prime[prime] for prime in sorted(certified_primes)
    ]
    pmisses, pmeta = exact_uncovered.find_uncovered(
        certificate_rows,
        max_component=args.max_component,
        limit=1,
        algebraic_primes=algebraic,
    )
    zmisses, zmeta = exact_uncovered_z3_bv.find_uncovered(
        certificate_rows,
        max_component=args.max_component,
        limit=1,
        algebraic_primes=algebraic,
    )
    if pmisses or zmisses:
        raise RuntimeError("final projective certificate is not a cover")
    args.certificate_output.write_text(
        json.dumps(
            {
                "parent_pool": str(args.parent_pool),
                "phase_file": str(args.phase_output),
                "projective_prime": q,
                "primes": sorted(certified_primes),
                "rows": certificate_rows,
                "primary_verification": pmeta,
                "independent_verification": zmeta,
            },
            indent=2,
        )
        + "\n"
    )
    manifest["complete"] = True
    manifest["final_phase_file"] = str(args.phase_output)
    manifest["final_locks_file"] = str(args.locks_output)
    manifest["certificate"] = str(args.certificate_output)
    manifest["certificate_rows"] = len(certificate_rows)
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"PROJECTIVE_COVER q={q} directions_closed={iteration} "
        f"certificate_rows={len(certificate_rows)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
