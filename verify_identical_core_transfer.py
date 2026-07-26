#!/usr/bin/env python3
"""Transfer an exact UNSAT result to independently verified identical cores."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def row_key(row: dict[str, object]) -> tuple[int, ...]:
    return tuple(
        int(row[field])
        for field in (
            "h",
            "p",
            "a",
            "b",
            "target_residue",
            "target_modulus",
        )
    )


def parse_case(text: str) -> tuple[int, Path, Path]:
    period_text, separator, rest = text.partition(":")
    if not separator:
        raise argparse.ArgumentTypeError(
            "case must be PERIOD:CORE:CORE_VERIFICATION"
        )
    core_text, separator, verification_text = rest.partition(":")
    if not separator:
        raise argparse.ArgumentTypeError(
            "case must be PERIOD:CORE:CORE_VERIFICATION"
        )
    return int(period_text), Path(core_text), Path(verification_text)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("base_core", type=Path)
    parser.add_argument("base_unsat_audit", type=Path)
    parser.add_argument("--case", type=parse_case, action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.case:
        raise SystemExit("at least one --case is required")

    base_core_payload = json.loads(args.base_core.read_text())
    base_rows = tuple(
        sorted(row_key(row) for row in base_core_payload["choices"])
    )
    unsat = json.loads(args.base_unsat_audit.read_text())
    if not bool(unsat["all_phase_assignments_excluded"]):
        raise AssertionError("base audit does not exclude every phase assignment")
    if int(unsat["selected_rows"]) != len(base_rows):
        raise AssertionError("base audit row count does not match base core")

    transferred = []
    seen_periods = set()
    for period, core_path, verification_path in args.case:
        if period in seen_periods:
            raise AssertionError(f"duplicate transferred period {period}")
        seen_periods.add(period)
        core = json.loads(core_path.read_text())
        rows = tuple(sorted(row_key(row) for row in core["choices"]))
        if rows != base_rows:
            raise AssertionError(f"period {period} core is not identical")
        verification = json.loads(verification_path.read_text())
        if not bool(verification["verified"]):
            raise AssertionError(
                f"period {period} component reduction is not verified"
            )
        if Path(verification["result"]) != core_path:
            raise AssertionError(
                f"period {period} verification points to another core"
            )
        source_path = Path(verification["source"])
        source = json.loads(source_path.read_text())
        if int(source["period_filter"]) != period:
            raise AssertionError(
                f"period {period} source has another period filter"
            )
        if int(verification["survivors"]) != len(base_rows):
            raise AssertionError(
                f"period {period} survivor count mismatch"
            )
        transferred.append(
            {
                "period": period,
                "source": str(source_path),
                "core": str(core_path),
                "component_core_verification": str(verification_path),
                "rows": len(rows),
                "no_cover": True,
            }
        )

    result = {
        "base_core": str(args.base_core),
        "base_unsat_audit": str(args.base_unsat_audit),
        "base_rows": len(base_rows),
        "transferred_periods": transferred,
        "all_transfers_verified": True,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"PASS base_rows={len(base_rows)} "
        f"transferred={len(transferred)} "
        f"periods={sorted(seen_periods)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
