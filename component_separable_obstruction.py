#!/usr/bin/env python3
"""Certify impossibility for a pool split across two independent components.

If every row has modulus p or q, then the uncovered set is the Cartesian
product of the points missed by the p-rows and the points missed by the
q-rows.  Fewer than p affine lines cannot cover F_p^2, since each line has
exactly p points; likewise for q.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 1 if divisor == 2 else 2
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--first-prime", type=int, required=True)
    parser.add_argument("--second-prime", type=int, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    args = parser.parse_args()

    p = args.first_prime
    q = args.second_prime
    if p == q or not is_prime(p) or not is_prime(q):
        raise SystemExit("components must be distinct primes")
    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    unsupported = sorted(
        {
            int(row["h"])
            for row in rows
            if int(row["h"]) not in (p, q)
        }
    )
    if unsupported:
        raise RuntimeError(f"unsupported row moduli: {unsupported}")
    if any(int(row["target_modulus"]) != 1 for row in rows):
        raise RuntimeError("all row targets must remain freely selectable")
    for row in rows:
        modulus = int(row["h"])
        if math.gcd(
            int(row["a"]),
            int(row["b"]),
            modulus,
        ) != 1:
            raise RuntimeError(
                f"row p={row['p']} is not an affine line modulo {modulus}"
            )

    first_rows = [row for row in rows if int(row["h"]) == p]
    second_rows = [row for row in rows if int(row["h"]) == q]
    first_covered_upper = len(first_rows) * p
    second_covered_upper = len(second_rows) * q
    first_holes_lower = p * p - first_covered_upper
    second_holes_lower = q * q - second_covered_upper
    proved = first_holes_lower > 0 and second_holes_lower > 0

    result = {
        "pool": str(args.pool),
        "components": [p, q],
        "row_counts": {str(p): len(first_rows), str(q): len(second_rows)},
        "plane_sizes": {str(p): p * p, str(q): q * q},
        "covered_point_upper_bounds": {
            str(p): first_covered_upper,
            str(q): second_covered_upper,
        },
        "hole_lower_bounds": {
            str(p): first_holes_lower,
            str(q): second_holes_lower,
        },
        "product_hole_lower_bound": first_holes_lower * second_holes_lower,
        "proved_no_cover": proved,
        "scope": "all phase assignments for every row in the supplied pool",
        "argument": (
            "The p-rows and q-rows depend on disjoint CRT components. "
            "Their global uncovered set is U_p x U_q. Each affine line in "
            "F_r^2 has r points, so the stated row counts leave both factors "
            "nonempty."
        ),
    }
    args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"p={p} rows={len(first_rows)} holes_lb={first_holes_lower} "
        f"q={q} rows={len(second_rows)} holes_lb={second_holes_lower} "
        f"product_lb={first_holes_lower * second_holes_lower}",
        flush=True,
    )
    print(
        "PROVED separable-component obstruction"
        if proved
        else "NO obstruction from separable line counts",
        flush=True,
    )
    return 0 if proved else 1


if __name__ == "__main__":
    raise SystemExit(main())
