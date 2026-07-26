#!/usr/bin/env python3
"""Independently replay a composed locked-branch obstruction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    grand = json.loads(Path(certificate["grand_pool"]).read_text())
    independent = json.loads(
        Path(certificate["independent_pool"]).read_text()
    )
    reduction = json.loads(
        Path(certificate["coarse_reduction"]).read_text()
    )
    obstruction = json.loads(
        Path(certificate["coarse_obstruction"]).read_text()
    )
    component = int(certificate["component_prime"])
    expected_independent = sorted(
        int(row["p"])
        for row in grand["choices"]
        if int(row["h"]) % component != 0
    )
    actual_independent = sorted(
        int(row["p"]) for row in independent["choices"]
    )
    incident_count = sum(
        int(row["h"]) % component == 0
        for row in grand["choices"]
    )
    hole_lower = (
        component * component - incident_count * component
    )
    checks = {
        "independent_rows_exactly_match": (
            expected_independent == actual_independent
        ),
        "incident_count_below_component": (
            incident_count < component
        ),
        "component_hole_bound_positive": hole_lower > 0,
        "component_counts_match": (
            len(actual_independent)
            == int(certificate["component_independent_rows"])
            and incident_count
            == int(certificate["component_incident_rows"])
        ),
        "component_hole_bound_matches": (
            hole_lower
            == int(certificate["component_hole_lower_bound"])
        ),
        "reduction_source_matches": (
            Path(reduction["pool"]).name
            == Path(certificate["independent_pool"]).name
        ),
        "reduction_proved": bool(
            reduction["proved_equivalence"]
        ),
        "obstruction_source_matches": (
            Path(obstruction["pool"]).name
            == Path(reduction["coarse_pool"]).name
        ),
        "coarse_obstruction_proved": bool(
            obstruction["proved_no_cover"]
        ),
        "certificate_claims_no_cover": bool(
            certificate["proved_no_cover"]
        ),
    }
    passed = all(checks.values())
    result = {
        "certificate": str(args.certificate),
        "checks": checks,
        "recomputed": {
            "independent_rows": len(actual_independent),
            "incident_rows": incident_count,
            "component_hole_lower_bound": hole_lower,
        },
        "verified": passed,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS" if passed else "FAIL", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
