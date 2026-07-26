#!/usr/bin/env python3
"""Optimize a single affine-plane component of a conditioned fibre pool.

Rows with modulus q and the same projective normal direction are
interchangeable.  Select at most the available number of parallel lines in
each direction and minimize the number of uncovered points of F_q^2.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--solver", default="g4")
    parser.add_argument("--result-output", type=Path)
    args = parser.parse_args()
    q = args.prime
    if q < 2:
        raise SystemExit("--prime must be at least two")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.examples.rc2 import RC2  # type: ignore
    from pysat.formula import IDPool, WCNF  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = [row for row in payload["choices"] if int(row["h"]) == q]
    groups: dict[tuple[int, int], list[tuple[dict, int]]] = defaultdict(list)
    for row in rows:
        if int(row["target_modulus"]) != 1:
            raise RuntimeError("optimization requires target_modulus=1")
        a = int(row["a"]) % q
        b = int(row["b"]) % q
        if not a and not b:
            raise RuntimeError("zero affine normal")
        scale = pow(a, -1, q) if a else pow(b, -1, q)
        direction = (a * scale % q, b * scale % q)
        groups[direction].append((row, scale))

    vpool = IDPool()
    formula = WCNF()
    line_var: dict[tuple[tuple[int, int], int], int] = {}
    for direction, members in groups.items():
        variables = []
        for target in range(q):
            variable = vpool.id(("line", direction, target))
            line_var[(direction, target)] = variable
            variables.append(variable)
        if len(members) < q:
            encoding = CardEnc.atmost(
                variables,
                bound=len(members),
                vpool=vpool,
                encoding=EncType.seqcounter,
            )
            for clause in encoding.clauses:
                formula.append(clause)

    for k in range(q):
        for l in range(q):
            formula.append(
                [
                    line_var[
                        (
                            direction,
                            (direction[0] * k + direction[1] * l) % q,
                        )
                    ]
                    for direction in groups
                ],
                weight=1,
            )

    started = time.monotonic()
    with RC2(formula, solver=args.solver, adapt=True, exhaust=True) as rc2:
        model = rc2.compute()
        cost = rc2.cost
    solve_seconds = time.monotonic() - started
    positive = {literal for literal in model if literal > 0}
    selected = {
        direction: [
            target
            for target in range(q)
            if line_var[(direction, target)] in positive
        ]
        for direction in groups
    }
    uncovered = [
        [k, l]
        for k in range(q)
        for l in range(q)
        if not any(
            (
                direction[0] * k + direction[1] * l
            )
            % q
            in selected[direction]
            for direction in groups
        )
    ]
    if len(uncovered) != cost:
        raise AssertionError(
            f"decoded uncovered count {len(uncovered)} != RC2 cost {cost}"
        )

    result = {
        "pool": str(args.pool),
        "prime": q,
        "row_count": len(rows),
        "direction_count": len(groups),
        "direction_capacities": {
            f"{a},{b}": len(members)
            for (a, b), members in groups.items()
        },
        "selected_targets": {
            f"{a},{b}": targets
            for (a, b), targets in selected.items()
        },
        "minimum_uncovered": cost,
        "uncovered_points": uncovered,
        "solver": args.solver,
        "solve_seconds": solve_seconds,
    }
    if args.result_output:
        args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"q={q} rows={len(rows)} directions={len(groups)} "
        f"minimum_uncovered={cost} solve_s={solve_seconds:.3f}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
