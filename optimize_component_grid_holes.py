#!/usr/bin/env python3
"""Exactly minimize holes left by one-prime affine fibres on F_q^2."""

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
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--fixed-targets", default="")
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--result-output", type=Path, required=True)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType
    from pysat.examples.rc2 import RC2
    from pysat.formula import IDPool, WCNF

    payload = json.loads(args.pool.read_text())
    rows = [
        row
        for row in payload["choices"]
        if int(row["h"]) == args.prime
    ]
    fixed_targets = {}
    for item in args.fixed_targets.split(","):
        if not item:
            continue
        row_prime, target = item.split(":", 1)
        fixed_targets[int(row_prime)] = int(target)
    fixed_rows = [
        row for row in rows if int(row["p"]) in fixed_targets
    ]
    mutable_rows = [
        row for row in rows if int(row["p"]) not in fixed_targets
    ]

    vpool = IDPool()
    formula = WCNF()
    for row_index, _row in enumerate(mutable_rows):
        literals = [
            vpool.id((row_index, target))
            for target in range(args.prime)
        ]
        encoding = CardEnc.equals(
            literals,
            bound=1,
            vpool=vpool,
            encoding=EncType.seqcounter,
        )
        for clause in encoding.clauses:
            formula.append(clause)

    initially_open = []
    for k in range(args.prime):
        for l in range(args.prime):
            if any(
                (
                    int(row["a"]) * k
                    + int(row["b"]) * l
                    - fixed_targets[int(row["p"])]
                )
                % args.prime
                == 0
                for row in fixed_rows
            ):
                continue
            clause = [
                vpool.id(
                    (
                        row_index,
                        (
                            int(row["a"]) * k
                            + int(row["b"]) * l
                        )
                        % args.prime,
                    )
                )
                for row_index, row in enumerate(mutable_rows)
            ]
            formula.append(clause, weight=1)
            initially_open.append((k, l))

    print(
        f"prime={args.prime} rows={len(rows)} fixed={len(fixed_rows)} "
        f"mutable={len(mutable_rows)} initially_open={len(initially_open)} "
        f"variables={vpool.top}",
        flush=True,
    )
    started = time.monotonic()
    with RC2(
        formula,
        solver=args.solver,
        adapt=True,
        exhaust=True,
        incr=False,
        verbose=1,
    ) as optimizer:
        model = optimizer.compute()
        optimum = optimizer.cost
    elapsed = time.monotonic() - started
    if model is None:
        raise RuntimeError("hard phase constraints are unsatisfiable")
    positive = {literal for literal in model if literal > 0}
    phases = {}
    for row_index, row in enumerate(mutable_rows):
        phases[int(row["p"])] = next(
            target
            for target in range(args.prime)
            if vpool.id((row_index, target)) in positive
        )
    misses = []
    for k, l in initially_open:
        if not any(
            (
                int(row["a"]) * k
                + int(row["b"]) * l
                - phases[int(row["p"])]
            )
            % args.prime
            == 0
            for row in mutable_rows
        ):
            misses.append((k, l))
    if len(misses) != optimum:
        raise AssertionError("decoded assignment disagrees with RC2 optimum")

    result = {
        "pool": str(args.pool),
        "prime": args.prime,
        "fixed_rows": len(fixed_rows),
        "mutable_rows": len(mutable_rows),
        "initially_open": len(initially_open),
        "optimal_holes": optimum,
        "mutable_phases": {
            str(row_prime): target
            for row_prime, target in phases.items()
        },
        "misses": [[k, l] for k, l in misses],
        "solver": args.solver,
        "elapsed_seconds": elapsed,
    }
    args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"OPTIMUM holes={optimum} elapsed_s={elapsed:.3f}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
