#!/usr/bin/env python3
"""Certify reduction of a two-high-digit cover to its coarse rows.

The supported pool has a squarefree coarse period Q and every noncoarse row
has exactly one additional factor p or q.  On a coarse cell missed by every
coarse row, the surviving refiners split into affine lines on independent
F_p^2 and F_q^2 digit planes.  If there are fewer than p and q respectively,
their Cartesian product contains an uncovered refinement.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--coarse-period", type=int, required=True)
    parser.add_argument("--first-high-prime", type=int, required=True)
    parser.add_argument("--second-high-prime", type=int, required=True)
    parser.add_argument("--coarse-pool-output", type=Path, required=True)
    parser.add_argument("--certificate-output", type=Path, required=True)
    args = parser.parse_args()

    first = args.first_high_prime
    second = args.second_high_prime
    if first == second:
        raise SystemExit("high-digit primes must be distinct")
    if args.coarse_period % first or args.coarse_period % second:
        raise SystemExit("both high-digit primes must divide the coarse period")

    payload = json.loads(args.pool.read_text())
    coarse_rows = []
    first_rows = []
    second_rows = []
    unsupported = []
    for row in payload["choices"]:
        h = int(row["h"])
        coarse_part = math.gcd(h, args.coarse_period)
        refinement = h // coarse_part
        if refinement == 1:
            if args.coarse_period % h:
                unsupported.append(
                    {
                        "p": int(row["p"]),
                        "h": h,
                        "reason": "coarse modulus does not divide period",
                    }
                )
            else:
                coarse_rows.append(row)
        elif refinement == first:
            first_rows.append(row)
        elif refinement == second:
            second_rows.append(row)
        else:
            unsupported.append(
                {
                    "p": int(row["p"]),
                    "h": h,
                    "refinement": refinement,
                }
            )

    first_holes = (
        first * first - len(first_rows) * first
    )
    second_holes = (
        second * second - len(second_rows) * second
    )
    proved = (
        not unsupported
        and first_holes > 0
        and second_holes > 0
    )
    if not proved:
        raise RuntimeError("pool does not satisfy the coarse-reduction bounds")

    coarse_payload = dict(payload)
    coarse_payload["source_pool"] = str(args.pool)
    coarse_payload["coarse_period"] = args.coarse_period
    coarse_payload["source_choice_count"] = len(payload["choices"])
    coarse_payload["retained_choice_count"] = len(coarse_rows)
    coarse_payload["choices"] = coarse_rows
    args.coarse_pool_output.write_text(
        json.dumps(coarse_payload, indent=2) + "\n"
    )

    certificate = {
        "pool": str(args.pool),
        "coarse_pool": str(args.coarse_pool_output),
        "coarse_period": args.coarse_period,
        "high_digit_primes": [first, second],
        "row_counts": {
            "all": len(payload["choices"]),
            "coarse": len(coarse_rows),
            str(first): len(first_rows),
            str(second): len(second_rows),
            "unsupported": len(unsupported),
        },
        "refinement_hole_lower_bounds": {
            str(first): first_holes,
            str(second): second_holes,
        },
        "product_hole_lower_bound": first_holes * second_holes,
        "unsupported_rows": unsupported,
        "proved_equivalence": True,
        "equivalence": (
            "The full pool covers Z^2 for a phase assignment iff its coarse "
            "rows cover every coordinate pair modulo the coarse period."
        ),
        "argument": (
            "In a coarse cell missed by all coarse rows, each remaining row "
            "depends on exactly one of two independent higher-digit planes. "
            "At most n affine lines cover n*r points of F_r^2. Both residual "
            "hole factors are nonempty, so their Cartesian product contains "
            "an uncovered refinement."
        ),
    }
    args.certificate_output.write_text(
        json.dumps(certificate, indent=2) + "\n"
    )
    print(
        f"all={len(payload['choices'])} coarse={len(coarse_rows)} "
        f"high_{first}={len(first_rows)} holes_lb={first_holes} "
        f"high_{second}={len(second_rows)} holes_lb={second_holes}",
        flush=True,
    )
    print("PROVED coarse-cover equivalence", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
