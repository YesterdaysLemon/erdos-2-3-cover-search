#!/usr/bin/env python3
"""Merge compatible affine-fibre pool artifacts by prime."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("base", type=Path)
    parser.add_argument("additions", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    base = json.loads(args.base.read_text())
    rows = {int(row["p"]): row for row in base["choices"]}
    base_count = len(rows)
    if base_count != len(base["choices"]):
        raise RuntimeError("base pool contains duplicate primes")

    addition_counts = []
    duplicates = 0
    incomplete_sources = []
    for path in args.additions:
        payload = json.loads(path.read_text())
        source_scan_status = payload.get("source_scan_status", {})
        if int(
            payload.get(
                "unresolved",
                source_scan_status.get("unresolved", 0),
            )
        ):
            incomplete_sources.append(str(path))
        added_here = 0
        for row in payload["choices"]:
            prime = int(row["p"])
            if prime in rows:
                if rows[prime] != row:
                    raise RuntimeError(
                        f"conflicting rows for prime {prime} from {path}"
                    )
                duplicates += 1
                continue
            rows[prime] = row
            added_here += 1
        addition_counts.append(
            {
                "path": str(path),
                "input_rows": len(payload["choices"]),
                "new_rows": added_here,
            }
        )

    choices = sorted(
        rows.values(), key=lambda row: (int(row["h"]), int(row["p"]))
    )
    result = dict(base)
    result["choices"] = choices
    result["source_union"] = [str(args.base)] + [
        str(path) for path in args.additions
    ]
    result["union_metadata"] = {
        "base_rows": base_count,
        "additions": addition_counts,
        "duplicate_rows": duplicates,
        "output_rows": len(choices),
        "incomplete_sources": incomplete_sources,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base={base_count} added={len(choices) - base_count} "
        f"duplicates={duplicates} output_rows={len(choices)} "
        f"incomplete_sources={len(incomplete_sources)} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
