#!/usr/bin/env python3
"""Lift a phase-independent conditioned-cell obstruction to its parent.

If every component cell missed by the parent's pure component lines has the
same obstructed residual geometry, those pure lines would need to cover the
entire affine component plane.  The Jamison--Brouwer--Schrijver threshold
then rules this out when there is no full parallel class.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path


def direction(row: dict, prime: int) -> tuple[int, int]:
    a = int(row["a"]) % prime
    b = int(row["b"]) % prime
    if a:
        return 1, b * pow(a, -1, prime) % prime
    if b:
        return 0, 1
    raise RuntimeError("degenerate component line")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("parent_pool", type=Path)
    parser.add_argument("conditioned_pool", type=Path)
    parser.add_argument("conditioned_obstruction", type=Path)
    parser.add_argument("--component-prime", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    parent = json.loads(args.parent_pool.read_text())
    child = json.loads(args.conditioned_pool.read_text())
    obstruction = json.loads(args.conditioned_obstruction.read_text())
    component = args.component_prime
    cell = child["cell_condition"]
    if int(cell["period"]) != component:
        raise RuntimeError("child was conditioned on another period")
    if not obstruction.get("proved_no_cover"):
        raise RuntimeError("conditioned obstruction is not proved")
    if (
        Path(obstruction["grand_pool"]).name
        != args.conditioned_pool.name
    ):
        raise RuntimeError("obstruction uses another conditioned pool")

    pure_rows = [
        row
        for row in parent["choices"]
        if int(row["h"]) == component
    ]
    residual_parent_rows = [
        row
        for row in parent["choices"]
        if int(row["h"]) != component
    ]
    child_by_prime = {
        int(row["p"]): row for row in child["choices"]
    }
    if sorted(child_by_prime) != sorted(
        int(row["p"]) for row in residual_parent_rows
    ):
        raise RuntimeError(
            "conditioned pool is not exactly the non-pure parent rows"
        )

    structural_mismatches = []
    for row in residual_parent_rows:
        h = int(row["h"])
        a = int(row["a"]) % h
        b = int(row["b"]) % h
        common = math.gcd(component * a, component * b, h)
        expected = (
            h // common,
            (component * a // common) % (h // common),
            (component * b // common) % (h // common),
        )
        conditioned = child_by_prime[int(row["p"])]
        actual = (
            int(conditioned["h"]),
            int(conditioned["a"]),
            int(conditioned["b"]),
        )
        if expected != actual:
            structural_mismatches.append(
                {
                    "p": int(row["p"]),
                    "expected": list(expected),
                    "actual": list(actual),
                }
            )
    if structural_mismatches:
        raise RuntimeError("conditioned row geometry is not cell-invariant")

    directions = Counter(
        direction(row, component) for row in pure_rows
    )
    cover_threshold = 2 * component - 1
    maximum_direction_count = max(directions.values(), default=0)
    proved = (
        len(pure_rows) < cover_threshold
        and maximum_direction_count < component
    )
    if not proved:
        raise RuntimeError("pure component rows are outside cover theorem")

    result = {
        "parent_pool": str(args.parent_pool),
        "conditioned_pool": str(args.conditioned_pool),
        "conditioned_obstruction": str(
            args.conditioned_obstruction
        ),
        "component_prime": component,
        "parent_row_count": len(parent["choices"]),
        "pure_component_rows": len(pure_rows),
        "residual_rows": len(residual_parent_rows),
        "affine_cover_threshold": cover_threshold,
        "direction_counts": {
            f"{first},{second}": count
            for (first, second), count in sorted(directions.items())
        },
        "maximum_direction_count": maximum_direction_count,
        "conditioned_geometry_cell_invariant": True,
        "inherited_fixed_inactive_primes": [
            int(prime)
            for prime in parent.get("fixed_inactive_primes", ())
        ],
        "proved_no_cover": True,
        "scope": (
            "all phase assignments of rows present in the parent conditioned "
            "pool; its inherited fixed-inactive rows remain omitted"
        ),
        "argument": (
            "Any component cell missed by the pure component rows leaves the "
            "same residual row geometry as the certified impossible child, "
            "up to phase translation. Thus the pure rows would have to cover "
            "the whole affine plane. They contain no parallel class and are "
            "fewer than the nonparallel affine-cover threshold."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"parent={len(parent['choices'])} pure={len(pure_rows)} "
        f"residual={len(residual_parent_rows)} "
        f"threshold={cover_threshold} max_direction="
        f"{maximum_direction_count}",
        flush=True,
    )
    print("PROVED parent component-cell obstruction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
