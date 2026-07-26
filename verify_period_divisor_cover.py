#!/usr/bin/env python3
"""Brute-force verifier for a finite divisor-period cover.

The cover itself is replayed by direct enumeration.  For small mutable sets,
the verifier also exhausts every proper subset and every allowed target tuple
to prove that all named mutable rows are necessary relative to the fixed
targets.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path


def row_covers(row: dict, target: int, x: int, y: int) -> bool:
    return (
        int(row["a"]) * x + int(row["b"]) * y - target
    ) % int(row["h"]) == 0


def target_options(row: dict) -> range:
    return range(
        int(row["target_residue"]),
        int(row["h"]),
        int(row["target_modulus"]),
    )


def assignment_covers(
    period: int,
    rows: dict[int, dict],
    phases: dict[int, int],
) -> bool:
    return all(
        any(row_covers(rows[prime], target, x, y)
            for prime, target in phases.items())
        for x in range(period)
        for y in range(period)
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    pool = json.loads(Path(certificate["pool"]).read_text())
    period = int(certificate["period"])
    all_rows = {int(row["p"]): row for row in pool["choices"]}
    fixed = {
        int(prime): int(target)
        for prime, target in certificate["fixed_targets"].items()
    }
    mutable_targets = {
        int(prime): int(target)
        for prime, target in certificate["mutable_targets"].items()
    }
    mutable = [int(prime) for prime in certificate["mutable_primes"]]
    selected = dict(fixed)
    selected.update(mutable_targets)

    selected_exist = set(selected) <= set(all_rows)
    divisor_moduli = selected_exist and all(
        period % int(all_rows[prime]["h"]) == 0 for prime in selected
    )
    targets_allowed = selected_exist and all(
        target % int(all_rows[prime]["target_modulus"])
        == int(all_rows[prime]["target_residue"])
        for prime, target in selected.items()
    )
    full_cover = (
        selected_exist
        and assignment_covers(period, all_rows, selected)
    )

    core_phases = {
        int(row["p"]): int(row["c"])
        for row in certificate["choices"]
    }
    core_rows_match = all(
        prime in all_rows
        and all(
            int(row[field]) == int(all_rows[prime][field])
            for field in ("h", "a", "b")
        )
        for row in certificate["choices"]
        for prime in [int(row["p"])]
    )
    core_cover = (
        core_rows_match
        and assignment_covers(period, all_rows, core_phases)
    )

    smaller_subset_cover = None
    exhaustive_assignments = 0
    if len(mutable) <= 8:
        smaller_subset_cover = False
        for count in range(len(mutable)):
            for subset in itertools.combinations(mutable, count):
                option_lists = [
                    list(target_options(all_rows[prime]))
                    for prime in subset
                ]
                for targets in itertools.product(*option_lists):
                    exhaustive_assignments += 1
                    phases = dict(fixed)
                    phases.update(zip(subset, targets))
                    if assignment_covers(period, all_rows, phases):
                        smaller_subset_cover = True
                        break
                if smaller_subset_cover:
                    break
            if smaller_subset_cover:
                break

    checks = {
        "selected_rows_exist": selected_exist,
        "selected_moduli_divide_period": divisor_moduli,
        "targets_respect_restrictions": targets_allowed,
        "direct_full_cover": full_cover,
        "core_rows_match_pool": core_rows_match,
        "direct_core_cover": core_cover,
        "no_proper_mutable_subset_covers": smaller_subset_cover is False,
        "recorded_zero_misses": int(certificate["verified_misses"]) == 0,
        "certificate_claims_cover": bool(certificate["proved_cover"]),
    }
    passed = all(checks.values())
    result = {
        "certificate": str(args.certificate),
        "checks": checks,
        "recomputed": {
            "period_cells": period * period,
            "selected_rows": len(selected),
            "core_rows": len(core_phases),
            "mutable_rows": len(mutable),
            "proper_subset_assignments_exhausted": exhaustive_assignments,
        },
        "verified": passed,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS" if passed else "FAIL", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
