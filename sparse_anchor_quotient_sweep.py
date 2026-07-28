#!/usr/bin/env python3
"""Resumable sparse-repair sweep over the exact five-anchor quotient.

This is an experiment driver, not a proof checker.  It enumerates the 162
representatives obtained from the exact p=41/p=17 parity quotient, freezes the
five low anchors, and invokes ``finite_sample_milp_repair.py`` on a common
finite point corpus.  Every successful child solve replays all supplied
points before its result is accepted by the driver.
"""

from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
import time
from pathlib import Path


LOW_PRIMES = (41, 17, 19, 37, 73)


def low_anchor_representatives() -> tuple[tuple[int, ...], ...]:
    """Return the 162 representatives of the exact low-anchor quotient."""
    representatives = []
    for p41 in range(2):
        # On the p=41 residual half-plane, the two p=17 targets with parity
        # p41 are both inactive and hence equivalent.  Keep one canonical
        # target from that pair and both targets of the opposite parity.
        p17_targets = (
            (p41, *(target for target in range(4) if target % 2 != p41))
        )
        for p17, p19, p37, p73 in itertools.product(
            p17_targets,
            range(3),
            range(3),
            range(3),
        ):
            representatives.append((p41, p17, p19, p37, p73))
    result = tuple(representatives)
    if len(result) != 162 or len(set(result)) != 162:
        raise AssertionError("low-anchor quotient enumeration failed")
    return result


def branch_slug(index: int, targets: tuple[int, ...]) -> str:
    encoded = "_".join(map(str, targets))
    return f"branch_{index:03d}_{encoded}"


def write_summary(
    output: Path,
    *,
    pool: Path,
    point_paths: list[Path],
    base_phase: Path,
    max_changes: int,
    backend: str,
    first_branch: int,
    last_branch: int,
    records: list[dict],
    elapsed_seconds: float,
) -> None:
    counts: dict[str, int] = {}
    for record in records:
        status = str(record["status"])
        counts[status] = counts.get(status, 0) + 1
    payload = {
        "pool": str(pool),
        "points": [str(path) for path in point_paths],
        "base_phase": str(base_phase),
        "low_primes": list(LOW_PRIMES),
        "quotient_class_count": 162,
        "max_changes": max_changes,
        "backend": backend,
        "first_branch": first_branch,
        "last_branch": last_branch,
        "status_counts": counts,
        "records": records,
        "elapsed_seconds": elapsed_seconds,
    }
    output.write_text(json.dumps(payload, indent=2) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("points", type=Path, nargs="+")
    parser.add_argument("--base-phase", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    parser.add_argument("--max-changes", type=int, default=4)
    parser.add_argument("--time-limit", type=float, default=240.0)
    parser.add_argument(
        "--backend",
        choices=("milp", "sat", "rc2", "z3", "mask"),
        default="sat",
    )
    parser.add_argument("--first-branch", type=int, default=0)
    parser.add_argument("--last-branch", type=int, default=161)
    parser.add_argument(
        "--rerun",
        action="store_true",
        help="rerun branches that already have a result file",
    )
    args = parser.parse_args()
    if args.max_changes < 0:
        raise SystemExit("--max-changes must be nonnegative")
    if not 0 <= args.first_branch <= args.last_branch < 162:
        raise SystemExit("branch interval must lie in 0..161")

    payload = json.loads(args.pool.read_text())
    row_count = len(payload["choices"])
    row_primes = {int(row["p"]) for row in payload["choices"]}
    missing_low = set(LOW_PRIMES) - row_primes
    if missing_low:
        raise RuntimeError(f"low anchors absent from pool: {missing_low}")
    base = {
        str(prime): int(target)
        for prime, target in json.loads(args.base_phase.read_text()).items()
    }
    args.work_dir.mkdir(parents=True, exist_ok=True)
    representatives = low_anchor_representatives()
    records = []
    started = time.monotonic()

    for index in range(args.first_branch, args.last_branch + 1):
        targets = representatives[index]
        slug = branch_slug(index, targets)
        phase_input = args.work_dir / f"{slug}_input_phase.json"
        phase_output = args.work_dir / f"{slug}_repair_phase.json"
        result_output = args.work_dir / f"{slug}_repair_result.json"
        branch_phase = dict(base)
        branch_phase.update(
            {
                str(prime): target
                for prime, target in zip(LOW_PRIMES, targets)
            }
        )
        phase_input.write_text(json.dumps(branch_phase) + "\n")

        if result_output.exists() and not args.rerun:
            result = json.loads(result_output.read_text())
            status = str(result.get("result", "MALFORMED_RESULT"))
            return_code = None
            reused = True
        else:
            if result_output.exists():
                result_output.unlink()
            if phase_output.exists():
                phase_output.unlink()
            command = [
                sys.executable,
                str(
                    Path(__file__).with_name(
                        "finite_sample_sat_repair.py"
                        if args.backend in {"sat", "rc2"}
                        else (
                            "finite_sample_z3_repair.py"
                            if args.backend == "z3"
                            else (
                                "finite_sample_mask_repair.py"
                                if args.backend == "mask"
                                else "finite_sample_milp_repair.py"
                            )
                        )
                    )
                ),
                str(args.pool),
                *(str(path) for path in args.points),
                "--initial-phases",
                str(phase_input),
                "--phase-output",
                str(phase_output),
                "--result-output",
                str(result_output),
                "--fixed-primes",
                ",".join(map(str, LOW_PRIMES)),
                "--max-changes",
                str(args.max_changes),
                "--time-limit",
                (
                    "0"
                    if args.backend in {"sat", "rc2"}
                    else str(args.time_limit)
                ),
            ]
            if args.backend == "milp":
                command.extend(
                    (
                        "--core-max-coverage",
                        str(row_count),
                        "--feasibility-only",
                    )
                )
            elif args.backend == "rc2":
                command.append("--minimize-changes")
            branch_started = time.monotonic()
            try:
                completed = subprocess.run(
                    command,
                    check=False,
                    text=True,
                    capture_output=True,
                    timeout=(
                        args.time_limit + 15
                        if args.time_limit > 0
                        else None
                    ),
                )
            except subprocess.TimeoutExpired as error:
                branch_seconds = time.monotonic() - branch_started
                stdout = error.stdout or ""
                stderr = error.stderr or ""
                if isinstance(stdout, bytes):
                    stdout = stdout.decode(errors="replace")
                if isinstance(stderr, bytes):
                    stderr = stderr.decode(errors="replace")
                if stdout.strip():
                    print(stdout.strip(), flush=True)
                if stderr.strip():
                    print(stderr.strip(), file=sys.stderr, flush=True)
                result = {
                    "result": "TIME_LIMIT",
                    "driver_seconds": branch_seconds,
                    "driver_timeout_seconds": args.time_limit + 15,
                }
                result_output.write_text(json.dumps(result, indent=2) + "\n")
                status = "TIME_LIMIT"
                return_code = None
                reused = False
                completed = None
            branch_seconds = time.monotonic() - branch_started
            if completed is not None:
                return_code = completed.returncode
                reused = False
                if result_output.exists():
                    result = json.loads(result_output.read_text())
                    status = str(
                        result.get("result", "MALFORMED_RESULT")
                    )
                else:
                    result = {}
                    status = "DRIVER_ERROR"
                if completed.stdout.strip():
                    print(completed.stdout.strip(), flush=True)
                if completed.stderr.strip():
                    print(
                        completed.stderr.strip(),
                        file=sys.stderr,
                        flush=True,
                    )
            result.setdefault("driver_seconds", branch_seconds)

        if status == "INTEGER_MODEL":
            if int(result.get("full_misses", -1)) != 0:
                status = "REPLAY_FAILED"
            elif not phase_output.exists():
                status = "MISSING_PHASE"
        record = {
            "branch_index": index,
            "targets": list(targets),
            "status": status,
            "return_code": return_code,
            "reused": reused,
            "phase_input": str(phase_input),
            "phase_output": str(phase_output),
            "result_output": str(result_output),
            "changed_phases": result.get("changed_phases"),
            "minimum_changes": result.get("minimum_changes"),
            "solve_seconds": result.get("solve_seconds"),
        }
        records.append(record)
        print(
            f"branch={index}/161 targets={targets} status={status} "
            f"changed={record['changed_phases']} "
            f"minimum={record['minimum_changes']} reused={reused}",
            flush=True,
        )
        write_summary(
            args.summary_output,
            pool=args.pool,
            point_paths=args.points,
            base_phase=args.base_phase,
            max_changes=args.max_changes,
            backend=args.backend,
            first_branch=args.first_branch,
            last_branch=args.last_branch,
            records=records,
            elapsed_seconds=time.monotonic() - started,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
