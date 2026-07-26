#!/usr/bin/env python3
"""Select a certified base pool plus rows chosen by a greedy audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("base_pool", type=Path)
    parser.add_argument("selection_audit", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    by_prime = {int(row["p"]): row for row in payload["choices"]}
    if len(by_prime) != len(payload["choices"]):
        raise RuntimeError("augmented pool contains duplicate primes")

    base = json.loads(args.base_pool.read_text())
    base_primes = {int(row["p"]) for row in base["choices"]}
    for row in base["choices"]:
        prime = int(row["p"])
        if prime not in by_prime:
            raise RuntimeError(f"base prime {prime} is absent from augmented pool")
        if by_prime[prime] != row:
            raise RuntimeError(f"base row differs for prime {prime}")

    audit = json.loads(args.selection_audit.read_text())
    selected_primes = {
        int(record["p"]) for record in audit["selected_new_rows"]
    }
    if selected_primes & base_primes:
        raise RuntimeError("selection audit unexpectedly contains a base prime")
    missing = selected_primes - set(by_prime)
    if missing:
        raise RuntimeError(
            f"{len(missing)} selected primes are absent from augmented pool"
        )

    keep = base_primes | selected_primes
    rows = sorted(
        (by_prime[prime] for prime in keep),
        key=lambda row: (int(row["h"]), int(row["p"])),
    )
    result = dict(payload)
    result["choices"] = rows
    result["support_selection"] = {
        "source_pool": str(args.pool),
        "base_pool": str(args.base_pool),
        "selection_audit": str(args.selection_audit),
        "base_rows": len(base_primes),
        "selected_rows": len(selected_primes),
        "output_rows": len(rows),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"base={len(base_primes)} selected={len(selected_primes)} "
        f"output_rows={len(rows)} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
