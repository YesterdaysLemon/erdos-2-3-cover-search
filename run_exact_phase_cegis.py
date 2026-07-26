#!/usr/bin/env python3
"""Drive exact finite phase synthesis against the exact CRT checker."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--points", type=Path, required=True)
    parser.add_argument("--phases", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--start-index", type=int, default=1)
    parser.add_argument("--iterations", type=int, default=20)
    parser.add_argument("--batch", type=int, default=5000)
    parser.add_argument("--min-parent-common", type=int, default=0)
    parser.add_argument("--max-component", type=int, default=300000)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument(
        "--fixed-targets",
        default="",
        help=(
            "comma-separated prime:target phases fixed by an external "
            "symmetry proof; forwarded to every finite master"
        ),
    )
    parser.add_argument("--diversity-coordinate-moduli", default="")
    parser.add_argument("--status-output", type=Path)
    args = parser.parse_args()
    if args.iterations < 1 or args.batch < 1:
        raise SystemExit("iterations and batch must be positive")

    points_path = args.points
    phases_path = args.phases
    points = [
        (int(k), int(l))
        for k, l in json.loads(points_path.read_text())
    ]
    seen = set(points)
    started = time.monotonic()
    status = {
        "pool": str(args.pool),
        "prefix": args.prefix,
        "fixed_targets": args.fixed_targets,
        "complete": False,
        "finite_obstruction": False,
        "iterations": [],
    }

    def write_status() -> None:
        if args.status_output:
            args.status_output.write_text(
                json.dumps(status, indent=2) + "\n"
            )

    for offset in range(args.iterations):
        index = args.start_index + offset
        stem = f"{args.prefix}_i{index:02d}"
        next_phases = Path(f"{stem}_phases.json")
        master_result = Path(f"{stem}_master_result.json")
        master_command = [
            sys.executable,
            "finite_phase_master.py",
            str(args.pool),
            str(points_path),
            "--min-parent-common",
            str(args.min_parent_common),
            "--solver",
            args.solver,
            "--initial-phases",
            str(phases_path),
            "--phase-output",
            str(next_phases),
            "--result-output",
            str(master_result),
        ]
        if args.fixed_targets:
            master_command.extend(
                ["--fixed-targets", args.fixed_targets]
            )
        master_started = time.monotonic()
        master = subprocess.run(master_command, check=False)
        master_seconds = time.monotonic() - master_started
        master_payload = json.loads(master_result.read_text())
        record = {
            "index": index,
            "input_points": len(points),
            "master_returncode": master.returncode,
            "master_sat": bool(master_payload["sat"]),
            "master_seconds": master_seconds,
        }
        status["iterations"].append(record)
        if not master_payload["sat"]:
            status["finite_obstruction"] = True
            status["last_points"] = str(points_path)
            status["elapsed_seconds"] = time.monotonic() - started
            write_status()
            print(
                f"iteration={index} result=FINITE_UNSAT "
                f"points={len(points)}",
                flush=True,
            )
            return 2

        checker_result = Path(f"{stem}_checker_result.json")
        checker_command = [
            sys.executable,
            "exact_derived_phase_misses.py",
            str(args.pool),
            "--phase-file",
            str(next_phases),
            "--max-component",
            str(args.max_component),
            "--limit",
            str(args.batch),
            "--solver",
            args.solver,
            "--output",
            str(checker_result),
        ]
        if args.diversity_coordinate_moduli:
            checker_command.extend(
                [
                    "--diversity-coordinate-moduli",
                    args.diversity_coordinate_moduli,
                ]
            )
        checker_started = time.monotonic()
        checker = subprocess.run(checker_command, check=False)
        checker_seconds = time.monotonic() - checker_started
        checker_payload = json.loads(checker_result.read_text())
        misses = [
            (int(k), int(l))
            for k, l in checker_payload["misses"]
        ]
        record.update(
            {
                "phase_file": str(next_phases),
                "checker_returncode": checker.returncode,
                "checker_misses": len(misses),
                "checker_seconds": checker_seconds,
            }
        )
        if not misses:
            status["complete"] = True
            status["phase_file"] = str(next_phases)
            status["checker_result"] = str(checker_result)
            status["elapsed_seconds"] = time.monotonic() - started
            write_status()
            print(
                f"iteration={index} result=COVER "
                f"points={len(points)}",
                flush=True,
            )
            return 0

        added = 0
        for point in misses:
            if point in seen:
                continue
            seen.add(point)
            points.append(point)
            added += 1
        if not added:
            raise AssertionError("checker returned no new witness")
        next_points = Path(f"{stem}_points.json")
        next_points.write_text(
            json.dumps([[k, l] for k, l in points]) + "\n"
        )
        record["added_points"] = added
        record["output_points"] = len(points)
        points_path = next_points
        phases_path = next_phases
        status["last_points"] = str(points_path)
        status["last_phases"] = str(phases_path)
        status["elapsed_seconds"] = time.monotonic() - started
        write_status()
        print(
            f"iteration={index} result=REFINE added={added} "
            f"points={len(points)} master_s={master_seconds:.3f} "
            f"checker_s={checker_seconds:.3f}",
            flush=True,
        )

    status["elapsed_seconds"] = time.monotonic() - started
    write_status()
    print(
        f"result=INCOMPLETE iterations={args.iterations} "
        f"points={len(points)}",
        flush=True,
    )
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
