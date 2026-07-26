#!/usr/bin/env python3
"""Lift a phase-independent conditioned obstruction through a finite period.

Conditioning a congruence row

    a*k + b*l == c (mod h)

on a cell (k,l) == (r,s) (mod P) either removes the row or produces a
congruence with modulus h/gcd(P*a, P*b, h).  The resulting modulus and
coefficients do not depend on (r,s); only the phase does.

Rows for which gcd(P*a, P*b, h) == h are "pure period" rows: on each
P-by-P residue cell they are either identically active or inactive.  If the
remaining conditioned geometry is impossible for every phase assignment,
then the pure rows would have to occupy every P-cell.  A union-capacity bound
rules that out when the sum of their exact per-row cell capacities is below
P**2.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def conditioned_geometry(row: dict, period: int) -> tuple[int, int, int]:
    """Return the cell-independent active geometry after conditioning."""
    h = int(row["h"])
    a = int(row["a"]) % h
    b = int(row["b"]) % h
    common = math.gcd(period * a, period * b, h)
    modulus = h // common
    return (
        modulus,
        (period * a // common) % modulus,
        (period * b // common) % modulus,
    )


def pure_cell_capacity(row: dict, period: int) -> int:
    """Maximum P-cells occupied by one pure row, over all phases."""
    h = int(row["h"])
    a = int(row["a"]) % h
    b = int(row["b"]) % h
    if math.gcd(period * a, period * b, h) != h:
        raise ValueError("row is not pure for this period")
    if period % h:
        raise ValueError("pure row modulus does not divide period")
    coefficient_gcd = math.gcd(a, b, h)
    return period * period * coefficient_gcd // h


def obstruction_uses_child(obstruction: dict, child_path: Path) -> bool:
    child_name = child_path.name
    for key in ("parent_pool", "grand_pool", "conditioned_pool", "source_pool"):
        value = obstruction.get(key)
        if value and Path(value).name == child_name:
            return True
    return False


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("parent_pool", type=Path)
    parser.add_argument("conditioned_pool", type=Path)
    parser.add_argument("conditioned_obstruction", type=Path)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.period <= 0:
        raise RuntimeError("period must be positive")

    parent = json.loads(args.parent_pool.read_text())
    child = json.loads(args.conditioned_pool.read_text())
    obstruction = json.loads(args.conditioned_obstruction.read_text())
    period = args.period

    cell = child.get("cell_condition", {})
    if int(cell.get("period", 0)) != period:
        raise RuntimeError("conditioned pool uses another period")
    if not obstruction.get("proved_no_cover"):
        raise RuntimeError("conditioned obstruction is not proved")
    if not obstruction_uses_child(obstruction, args.conditioned_pool):
        raise RuntimeError("conditioned obstruction does not use child pool")
    if child.get("full_cell_primes"):
        raise RuntimeError("conditioned pool contains full-cell rows")

    pure_rows: list[dict] = []
    residual_rows: list[dict] = []
    for row in parent["choices"]:
        geometry = conditioned_geometry(row, period)
        if geometry[0] == 1:
            pure_rows.append(row)
        else:
            residual_rows.append(row)

    child_by_prime = {int(row["p"]): row for row in child["choices"]}
    residual_primes = sorted(int(row["p"]) for row in residual_rows)
    if sorted(child_by_prime) != residual_primes:
        raise RuntimeError(
            "conditioned pool is not exactly the residual parent rows"
        )

    structural_mismatches = []
    for row in residual_rows:
        expected = conditioned_geometry(row, period)
        conditioned = child_by_prime[int(row["p"])]
        actual_h = int(conditioned["h"])
        actual = (
            actual_h,
            int(conditioned["a"]) % actual_h,
            int(conditioned["b"]) % actual_h,
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

    capacities = []
    for row in pure_rows:
        capacity = pure_cell_capacity(row, period)
        capacities.append(
            {
                "p": int(row["p"]),
                "h": int(row["h"]),
                "coefficient_gcd": math.gcd(
                    int(row["a"]), int(row["b"]), int(row["h"])
                ),
                "cell_capacity_upper_bound": capacity,
            }
        )
    capacity_sum = sum(item["cell_capacity_upper_bound"] for item in capacities)
    cell_count = period * period
    if capacity_sum >= cell_count:
        raise RuntimeError("pure-row union capacity does not force a missed cell")

    fixed_inactive = sorted(
        int(prime) for prime in child.get("fixed_inactive_primes", ())
    )
    pure_primes = sorted(int(row["p"]) for row in pure_rows)
    if fixed_inactive != pure_primes:
        raise RuntimeError(
            "child fixed-inactive list does not exactly record pure rows"
        )

    result = {
        "parent_pool": str(args.parent_pool),
        "conditioned_pool": str(args.conditioned_pool),
        "conditioned_obstruction": str(args.conditioned_obstruction),
        "period": period,
        "parent_row_count": len(parent["choices"]),
        "pure_period_row_count": len(pure_rows),
        "residual_row_count": len(residual_rows),
        "pure_period_rows": capacities,
        "period_cell_count": cell_count,
        "pure_union_capacity_upper_bound": capacity_sum,
        "forced_missed_cells_lower_bound": cell_count - capacity_sum,
        "conditioned_geometry_cell_invariant": True,
        "child_records_exact_pure_inactive_set": True,
        "inherited_fixed_inactive_primes": [
            int(prime)
            for prime in parent.get("fixed_inactive_primes", ())
        ],
        "proved_no_cover": True,
        "scope": (
            "all phase assignments of rows present in the parent pool; its "
            "inherited fixed-inactive rows remain omitted"
        ),
        "argument": (
            "Every period cell missed by the pure-period rows has at most the "
            "same residual geometry as the certified phase-independent "
            "impossible child. Hence the pure rows would have to occupy every "
            "period cell. Their union occupies at most the sum of the exact "
            "per-row cell capacities, which is strictly below P^2."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"parent={len(parent['choices'])} pure={len(pure_rows)} "
        f"residual={len(residual_rows)} capacity={capacity_sum}/{cell_count} "
        f"missed>={cell_count - capacity_sum}",
        flush=True,
    )
    print("PROVED parent period-cell obstruction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
