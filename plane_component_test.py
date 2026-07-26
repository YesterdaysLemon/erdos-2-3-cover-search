#!/usr/bin/env python3
"""Run one exact maximal-component plane relaxation with a logged verdict."""

from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path

import component_core
import exact_uncovered


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("prime", type=int)
    parser.add_argument("--solver", default="maplesat")
    args = parser.parse_args()

    rows = json.loads(args.pool.read_text())["choices"]
    factors = [
        {
            prime: component_core.valuation(int(row["h"]), prime)
            for prime in exact_uncovered.factor(int(row["h"]))
        }
        for row in rows
    ]
    exponent = max(item.get(args.prime, 0) for item in factors)
    indices = [
        index
        for index, item in enumerate(factors)
        if item.get(args.prime, 0) == exponent
    ]
    counts = collections.Counter(
        component_core.direction(rows[index], args.prime) for index in indices
    )
    print(
        json.dumps(
            {
                "prime": args.prime,
                "exponent": exponent,
                "rows": len(indices),
                "directions": len(counts),
                "max_parallel": max(counts.values()),
                "capacities": sorted(counts.values(), reverse=True),
            }
        ),
        flush=True,
    )
    sat, metadata = component_core.relaxed_plane_cover_sat(
        args.prime, counts, args.solver
    )
    print(json.dumps({"sat": sat, "metadata": metadata}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
