#!/usr/bin/env python3
"""Find component-core artifacts that used a bound outside its valid domain."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--directory", type=Path, default=Path("."), help="artifact directory"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    affected = []
    inspected = 0
    component_artifacts = 0
    for path in sorted(args.directory.glob("*.json")):
        inspected += 1
        try:
            payload = json.loads(path.read_text())
        except (OSError, UnicodeDecodeError, ValueError):
            continue
        if not isinstance(payload, dict):
            continue
        core = payload.get("component_core")
        if not isinstance(core, dict):
            continue
        component_artifacts += 1
        invalid_records = [
            {
                "round": int(record["round"]),
                "prime": int(record["prime"]),
                "exponent": int(record["exponent"]),
                "rows": int(record["rows"]),
                "max_parallel": int(record["max_parallel"]),
            }
            for record in core.get("audit", ())
            if record.get("reason") == "parallel_class_union_bound"
            and int(record["max_parallel"]) >= int(record["prime"])
        ]
        if invalid_records:
            affected.append(
                {
                    "path": path.name,
                    "invalid_records": invalid_records,
                }
            )

    result = {
        "criterion": (
            "parallel_class_union_bound with max_parallel >= prime"
        ),
        "inspected_json_files": inspected,
        "component_core_artifacts": component_artifacts,
        "affected_artifacts": len(affected),
        "artifacts": affected,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"inspected={inspected} component_core={component_artifacts} "
        f"affected={len(affected)} output={args.output}"
    )
    return 1 if affected else 0


if __name__ == "__main__":
    raise SystemExit(main())
