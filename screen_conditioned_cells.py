#!/usr/bin/env python3
"""Screen ranked open cells with protected, exact conditioned CEGIS runs.

The script is intentionally a proposal generator.  A reported cover must
still pass core extraction, independent Z3-BV verification, lifting to parent
coordinates, and direct parent-coordinate verification before it is merged
into the monotone certificate set.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path


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
    parser.add_argument("parent_phases", type=Path)
    parser.add_argument("locked_primes", type=Path)
    parser.add_argument("open_cells", type=Path)
    parser.add_argument("ranking", type=Path)
    parser.add_argument("--start-rank", type=int, default=1)
    parser.add_argument("--count", type=int, default=10)
    parser.add_argument("--rounds", type=int, default=16)
    parser.add_argument("--seed", type=int, default=203000)
    parser.add_argument("--prefix", default="conditioned_screen")
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    if args.start_rank < 1 or args.count < 1 or args.rounds < 1:
        raise SystemExit("rank, count, and rounds must be positive")

    lock_payload = json.loads(args.locked_primes.read_text())
    raw_locked = (
        lock_payload.get("primes", ())
        if isinstance(lock_payload, dict)
        else lock_payload
    )
    locked = sorted(
        int(prime)
        for prime in raw_locked
    )
    locked_csv = ",".join(map(str, locked))
    open_payload = json.loads(args.open_cells.read_text())
    points = [tuple(map(int, point)) for point in open_payload["points"]]
    point_by_cell = {
        (k % 16, l % 16, k % 27, l % 27): (k, l)
        for k, l in points
    }
    ranked = json.loads(args.ranking.read_text())["ranked"]
    selected = ranked[
        args.start_rank - 1 : args.start_rank - 1 + args.count
    ]
    manifest = {
        "parent_pool": str(args.parent_pool),
        "parent_phases": str(args.parent_phases),
        "locked_primes": str(args.locked_primes),
        "open_cells": str(args.open_cells),
        "ranking": str(args.ranking),
        "start_rank": args.start_rank,
        "requested_count": args.count,
        "rounds": args.rounds,
        "runs": [],
        "cover": None,
    }
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")

    for offset, record in enumerate(selected):
        rank = args.start_rank + offset
        cell = tuple(map(int, record["cell"]))
        tag = (
            f"{args.prefix}_r{rank:03d}_"
            f"{cell[0]}_{cell[1]}_{cell[2]}_{cell[3]}"
        )
        child_pool = Path(f"{tag}_pool.json")
        initial_phases = Path(f"{tag}_initial_phases.json")
        working_phases = Path(f"{tag}_working_phases.json")
        parent_seed = Path(f"{tag}_parent_seed.json")
        child_points = Path(f"{tag}_points.json")
        cover = Path(f"{tag}_cover.json")
        checker = Path(f"{tag}_checker.json")
        log = Path(f"{tag}.log")
        point = point_by_cell.get(cell)
        if point is None:
            raise RuntimeError(f"no exact witness stored for cell {cell}")
        parent_seed.write_text(json.dumps([list(point)]) + "\n")

        print(
            f"SCREEN rank={rank} cell={cell} "
            f"density={record['density']:.12f}",
            flush=True,
        )
        condition_command = [
            sys.executable,
            "condition_derived_cell.py",
            str(args.parent_pool),
            "--cell-period",
            "432",
            "--k-residue",
            str(record["k_residue"]),
            "--l-residue",
            str(record["l_residue"]),
            "--phase-file",
            str(args.parent_phases),
            "--phase-output",
            str(initial_phases),
            "--fixed-primes",
            locked_csv,
            "--canonicalize-algebraic",
            "--output",
            str(child_pool),
        ]
        subprocess.run(condition_command, check=True)
        subprocess.run(
            [
                sys.executable,
                "transform_cell_points.py",
                str(child_pool),
                str(parent_seed),
                "--output",
                str(child_points),
            ],
            check=True,
        )

        child_payload = json.loads(child_pool.read_text())
        present = {int(row["p"]) for row in child_payload["choices"]}
        active_locks = sorted(set(locked) & present)
        shutil.copyfile(initial_phases, working_phases)
        cell_seed = args.seed + rank
        command = [
            sys.executable,
            "local_phase_cegis.py",
            str(child_pool),
            "--derived-pool",
            "--derived-targets",
            "--points-file",
            str(child_points),
            "--phase-file",
            str(working_phases),
            "--fixed-primes",
            ",".join(map(str, active_locks)),
            "--output",
            str(cover),
            "--checker-checkpoint-file",
            str(checker),
            "--rounds",
            str(args.rounds),
            "--seed",
            str(cell_seed),
            "--randomize-mutable-initial",
            "--required-coverage",
            "2",
            "--random-checks",
            "0",
            "--exact-batch",
            "2000",
            "--diversity-coordinate-schedule",
            "7;7,11;7,13;none",
            "--diversity-stage-rounds",
            "2",
            "--diversity-quota",
            "2",
            "--max-component",
            "1024",
            "--retain-target-matrix-during-exact",
        ]
        return_code = run_logged(command, log)
        run_record = {
            "rank": rank,
            "cell": list(cell),
            "density": record["density"],
            "active_locks": active_locks,
            "seed": cell_seed,
            "return_code": return_code,
            "pool": str(child_pool),
            "initial_phases": str(initial_phases),
            "working_phases": str(working_phases),
            "points": str(child_points),
            "cover": str(cover) if cover.exists() else None,
            "checker": str(checker),
            "log": str(log),
        }
        manifest["runs"].append(run_record)
        if return_code == 0 and cover.exists():
            manifest["cover"] = run_record
            args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
            print(f"SCREEN_COVER rank={rank} cover={cover}", flush=True)
            return 0
        args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
        print(
            f"SCREEN_CONTINUE rank={rank} return_code={return_code}",
            flush=True,
        )

    print(
        f"SCREEN_EXHAUSTED tested={len(manifest['runs'])} "
        f"manifest={args.manifest}",
        flush=True,
    )
    return 4


if __name__ == "__main__":
    raise SystemExit(main())
