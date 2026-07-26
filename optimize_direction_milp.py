#!/usr/bin/env python3
"""Exact MILP minimum-uncovered optimization for one affine plane."""

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
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--result-output", type=Path)
    parser.add_argument("--phase-output", type=Path)
    args = parser.parse_args()
    q = args.prime
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from scipy.optimize import Bounds, LinearConstraint, milp  # type: ignore
    from scipy.sparse import lil_matrix  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = [row for row in payload["choices"] if int(row["h"]) == q]
    groups: dict[tuple[int, int], list[tuple[dict, int]]] = defaultdict(list)
    for row in rows:
        if int(row["target_modulus"]) != 1:
            raise RuntimeError("optimization requires target_modulus=1")
        a = int(row["a"]) % q
        b = int(row["b"]) % q
        scale = pow(a, -1, q) if a else pow(b, -1, q)
        direction = (a * scale % q, b * scale % q)
        groups[direction].append((row, scale))

    directions = sorted(groups)
    line_index = {
        (direction, target): index
        for index, (direction, target) in enumerate(
            (item for direction in directions for item in (
                (direction, target) for target in range(q)
            ))
        )
    }
    line_count = len(line_index)
    point_count = q * q
    variable_count = line_count + point_count
    constraint_count = point_count + len(directions)
    matrix = lil_matrix(
        (constraint_count, variable_count),
        dtype=np.float64,
    )
    lower = np.full(constraint_count, -np.inf)
    upper = np.full(constraint_count, np.inf)

    for k in range(q):
        for l in range(q):
            point = k * q + l
            for direction in directions:
                target = (
                    direction[0] * k + direction[1] * l
                ) % q
                matrix[point, line_index[(direction, target)]] = 1
            matrix[point, line_count + point] = 1
            lower[point] = 1

    for offset, direction in enumerate(directions):
        row_index = point_count + offset
        for target in range(q):
            matrix[row_index, line_index[(direction, target)]] = 1
        upper[row_index] = len(groups[direction])

    objective = np.zeros(variable_count)
    objective[line_count:] = 1
    integrality = np.ones(variable_count, dtype=np.uint8)
    started = time.monotonic()
    solved = milp(
        objective,
        integrality=integrality,
        bounds=Bounds(0, 1),
        constraints=LinearConstraint(matrix.tocsr(), lower, upper),
        options={
            "time_limit": args.time_limit,
            "presolve": True,
            "mip_rel_gap": 0.0,
        },
    )
    solve_seconds = time.monotonic() - started
    selected = defaultdict(list)
    uncovered = []
    phases = {}
    if solved.x is not None:
        for direction in directions:
            for target in range(q):
                if solved.x[line_index[(direction, target)]] > 0.5:
                    selected[direction].append(target)
        uncovered = [
            [point // q, point % q]
            for point in range(point_count)
            if solved.x[line_count + point] > 0.5
        ]
        for direction, members in groups.items():
            targets = selected[direction] or [0]
            for index, (row, scale) in enumerate(members):
                canonical = targets[min(index, len(targets) - 1)]
                phases[str(int(row["p"]))] = (
                    canonical * pow(scale, -1, q) % q
                )
        if not uncovered and args.phase_output:
            args.phase_output.write_text(json.dumps(phases) + "\n")

    result = {
        "pool": str(args.pool),
        "prime": q,
        "row_count": len(rows),
        "direction_count": len(groups),
        "status": int(solved.status),
        "success": bool(solved.success),
        "message": solved.message,
        "minimum_uncovered": (
            int(round(solved.fun)) if solved.fun is not None else None
        ),
        "mip_gap": getattr(solved, "mip_gap", None),
        "mip_node_count": getattr(solved, "mip_node_count", None),
        "uncovered_points": uncovered,
        "selected_targets": {
            f"{a},{b}": values
            for (a, b), values in selected.items()
        },
        "solve_seconds": solve_seconds,
    }
    if args.result_output:
        args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"q={q} status={solved.status} success={solved.success} "
        f"minimum_uncovered={result['minimum_uncovered']} "
        f"gap={result['mip_gap']} solve_s={solve_seconds:.3f}",
        flush=True,
    )
    if solved.success and result["minimum_uncovered"] == 0:
        return 0
    if solved.success:
        return 2
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
