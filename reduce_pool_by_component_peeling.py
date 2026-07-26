#!/usr/bin/env python3
"""Materialize the residual pool from a verified peeling certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.pool.read_text())
    certificate = json.loads(args.certificate.read_text())
    threshold = int(certificate["min_parent_common"])
    retained = [
        row
        for row in payload["choices"]
        if int(row["h"]) > 1
        and int(row.get("parent_common", 1)) >= threshold
    ]
    for step in certificate["steps"]:
        prime = int(step["prime"])
        incident = [
            row for row in retained if int(row["h"]) % prime == 0
        ]
        expected = [int(value) for value in step["removed_fibre_primes"]]
        if sorted(int(row["p"]) for row in incident) != expected:
            raise AssertionError(
                f"certificate row mismatch at component {prime}"
            )
        retained = [
            row for row in retained if int(row["h"]) % prime != 0
        ]
    expected_remaining = [
        int(value) for value in certificate["remaining_fibre_primes"]
    ]
    if sorted(int(row["p"]) for row in retained) != expected_remaining:
        raise AssertionError("certificate residual mismatch")

    result = dict(payload)
    result["choices"] = retained
    result["source_pool"] = str(args.pool)
    result["component_peeling_certificate"] = str(args.certificate)
    result["source_choice_count"] = len(payload["choices"])
    result["retained_choice_count"] = len(retained)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"source={len(payload['choices'])} retained={len(retained)} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
