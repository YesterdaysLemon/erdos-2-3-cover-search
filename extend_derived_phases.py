#!/usr/bin/env python3
"""Extend a derived-pool phase map without changing any saved phase.

Rows already present in the source phase map retain their exact target.
New rows receive the canonical allowed target_residue.  All targets are
checked against the derived row's target congruence before the result is
written.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("source_phases", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    source = {
        int(prime): int(target)
        for prime, target in json.loads(args.source_phases.read_text()).items()
    }
    extended: dict[str, int] = {}
    retained = 0
    added = 0
    for row in payload["choices"]:
        prime = int(row["p"])
        h = int(row["h"])
        residue = int(row["target_residue"])
        modulus = int(row["target_modulus"])
        if h % modulus:
            raise RuntimeError(
                f"target modulus {modulus} does not divide h={h} for p={prime}"
            )
        if prime in source:
            target = source[prime] % h
            retained += 1
        else:
            target = residue
            added += 1
        if target % modulus != residue:
            raise RuntimeError(
                f"phase {target} violates target {residue} mod {modulus} "
                f"for p={prime}"
            )
        extended[str(prime)] = target

    args.output.write_text(json.dumps(extended) + "\n")
    print(
        f"rows={len(extended)} retained={retained} added={added} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
