#!/usr/bin/env python3
"""Ask Z3 for an affine cover using one quantified lattice formula.

The phase targets are existential constants and the exponent coordinates are
universal integers.  Because every row predicate is periodic, quantifying all
integers is equivalent to quantifying a common finite period without
materializing its square.

This is an experimental second encoding.  ``unsat`` proves that the supplied
finite row family has no phase cover.  ``sat`` is replayed symbolically before
being reported.  ``unknown`` is explicitly inconclusive.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path


def solve_rows(rows: list[dict], timeout_ms: int):
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import z3  # type: ignore

    solver = z3.Solver()
    if timeout_ms > 0:
        solver.set(timeout=timeout_ms)
    targets = []
    for index, row in enumerate(rows):
        h = int(row["h"])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        if h < 1 or modulus < 1 or h % modulus:
            raise ValueError(f"invalid target restriction at row {index}")
        target = z3.Int(f"target_{index}")
        targets.append(target)
        solver.add(target >= 0, target < h, target % modulus == residue)
    k = z3.Int("universal_k")
    ell = z3.Int("universal_l")
    coverage = z3.Or(
        *[
            (
                int(row["a"]) * k
                + int(row["b"]) * ell
                - target
            )
            % int(row["h"])
            == 0
            for row, target in zip(rows, targets, strict=True)
        ]
    )
    solver.add(z3.ForAll([k, ell], coverage))
    started = time.monotonic()
    check = solver.check()
    elapsed = time.monotonic() - started
    model_targets = None
    if check == z3.sat:
        model = solver.model()
        model_targets = [
            model.eval(target, model_completion=True).as_long()
            for target in targets
        ]
        # Replay the same universally quantified predicate with concrete
        # targets.  Its negation must be UNSAT for a reported cover.
        replay = z3.Solver()
        if timeout_ms > 0:
            replay.set(timeout=timeout_ms)
        replay.add(
            z3.Not(
                z3.Or(
                    *[
                        (
                            int(row["a"]) * k
                            + int(row["b"]) * ell
                            - target
                        )
                        % int(row["h"])
                        == 0
                        for row, target in zip(
                            rows,
                            model_targets,
                            strict=True,
                        )
                    ]
                )
            )
        )
        replay_check = replay.check()
        if replay_check != z3.unsat:
            raise RuntimeError(
                "quantified SAT model failed concrete symbolic replay"
            )
    return {
        "check": str(check),
        "targets": model_targets,
        "elapsed_seconds": elapsed,
        "reason_unknown": (
            solver.reason_unknown() if check == z3.unknown else None
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--max-h", type=int, default=0)
    parser.add_argument("--time-limit", type=float, default=240.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.max_h < 0 or args.time_limit < 0:
        raise SystemExit("limits must be nonnegative")

    payload = json.loads(args.input.read_text())
    source_rows = payload.get("choices", payload.get("rows"))
    if source_rows is None:
        raise RuntimeError("input has neither choices nor rows")
    rows = [
        row
        for row in source_rows
        if not args.max_h or int(row["h"]) <= args.max_h
    ]
    if not rows:
        raise RuntimeError("selected row family is empty")
    result = solve_rows(rows, round(1000 * args.time_limit))
    period = math.lcm(*(int(row["h"]) for row in rows))
    output = {
        "input": str(args.input),
        "row_count": len(rows),
        "max_h": args.max_h,
        "period": period,
        "encoding": (
            "existential integer phase constants and two universally "
            "quantified integer exponent coordinates"
        ),
        "result": result["check"].upper(),
        "proved_no_cover": result["check"] == "unsat",
        "proved_cover": result["check"] == "sat",
        "targets": result["targets"],
        "elapsed_seconds": result["elapsed_seconds"],
        "reason_unknown": result["reason_unknown"],
        "scope": "all phase assignments for exactly the selected finite rows",
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"rows={len(rows)} period={period} result={output['result']} "
        f"seconds={result['elapsed_seconds']:.3f} "
        f"reason={result['reason_unknown']}",
        flush=True,
    )
    return 0 if result["check"] in {"sat", "unsat"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
