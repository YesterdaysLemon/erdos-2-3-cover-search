#!/usr/bin/env python3
"""Independently verify a separable-component obstruction certificate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    pool_path = Path(certificate["pool"])
    pool = json.loads(pool_path.read_text())
    p, q = map(int, certificate["components"])
    rows = pool["choices"]
    checks = {
        "only_two_component_moduli": all(
            int(row["h"]) in (p, q) for row in rows
        ),
        "unconstrained_targets": all(
            int(row["target_modulus"]) == 1 for row in rows
        ),
        "all_rows_nondegenerate": all(
            math.gcd(
                int(row["a"]),
                int(row["b"]),
                int(row["h"]),
            )
            == 1
            for row in rows
        ),
    }
    counts = {
        p: sum(int(row["h"]) == p for row in rows),
        q: sum(int(row["h"]) == q for row in rows),
    }
    holes = {
        p: p * p - counts[p] * p,
        q: q * q - counts[q] * q,
    }
    checks.update(
        {
            "first_factor_has_hole": holes[p] > 0,
            "second_factor_has_hole": holes[q] > 0,
            "certificate_counts_match": (
                counts[p]
                == int(certificate["row_counts"][str(p)])
                and counts[q]
                == int(certificate["row_counts"][str(q)])
            ),
            "certificate_hole_bounds_match": (
                holes[p]
                == int(certificate["hole_lower_bounds"][str(p)])
                and holes[q]
                == int(certificate["hole_lower_bounds"][str(q)])
            ),
            "certificate_product_matches": (
                holes[p] * holes[q]
                == int(certificate["product_hole_lower_bound"])
            ),
        }
    )
    passed = all(checks.values()) and bool(certificate["proved_no_cover"])
    result = {
        "certificate": str(args.certificate),
        "pool": str(pool_path),
        "checks": checks,
        "recomputed_row_counts": {
            str(p): counts[p],
            str(q): counts[q],
        },
        "recomputed_hole_lower_bounds": {
            str(p): holes[p],
            str(q): holes[q],
        },
        "verified": passed,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS" if passed else "FAIL", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
