#!/usr/bin/env python3
"""Restrict a prime-to-target phase map to the rows present in a pool."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("phases", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = json.loads(args.pool.read_text())["choices"]
    source = {
        int(prime): int(target)
        for prime, target in json.loads(args.phases.read_text()).items()
    }
    pool_primes = {int(row["p"]) for row in rows}
    if len(pool_primes) != len(rows):
        raise RuntimeError("pool contains duplicate primes")

    missing = sorted(pool_primes - set(source))
    if missing:
        raise RuntimeError(
            f"phase map is missing {len(missing)} pool primes; first={missing[0]}"
        )

    restricted = {
        str(prime): source[prime]
        for prime in sorted(pool_primes)
    }
    args.output.write_text(json.dumps(restricted) + "\n")
    print(
        f"pool_rows={len(rows)} source_phases={len(source)} "
        f"output_phases={len(restricted)} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
