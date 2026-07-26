#!/usr/bin/env python3
"""Extract, greedily minimize, and dual-verify an exact cover core."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import exact_uncovered
import exact_uncovered_z3_bv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("cover", type=Path)
    parser.add_argument("--max-component", type=int, default=16384)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cover = json.loads(args.cover.read_text())
    rows = cover["choices"]
    algebraic_primes = tuple(
        int(prime) for prime in cover.get("algebraic_primes", ())
    )
    sophie_germain = bool(cover.get("sophie_germain", False))
    misses, core_meta = exact_uncovered.find_uncovered(
        rows,
        max_component=args.max_component,
        limit=1,
        algebraic_primes=algebraic_primes,
        sophie_germain=sophie_germain,
        core_coordinate_cells=((),),
    )
    if misses or core_meta["sat"]:
        raise RuntimeError("the supplied phase assignment is not a cover")
    kept = [
        int(index)
        for index in core_meta["core_coordinate_cells"][0]["row_indices"]
    ]
    raw_core_size = len(kept)
    print(f"PRIMARY core={raw_core_size}", flush=True)

    for number, index in enumerate(list(kept), 1):
        trial = [value for value in kept if value != index]
        trial_rows = [rows[value] for value in trial]
        trial_misses, trial_meta = exact_uncovered_z3_bv.find_uncovered(
            trial_rows,
            max_component=args.max_component,
            limit=1,
            algebraic_primes=algebraic_primes,
            sophie_germain=sophie_germain,
        )
        if not trial_misses and not trial_meta["sat"]:
            kept = trial
            print(
                f"DELETE {number}/{raw_core_size} kept={len(kept)}",
                flush=True,
            )

    minimized_rows = [rows[index] for index in kept]
    primary_misses, primary_meta = exact_uncovered.find_uncovered(
        minimized_rows,
        max_component=args.max_component,
        limit=1,
        algebraic_primes=algebraic_primes,
        sophie_germain=sophie_germain,
    )
    z3_misses, z3_meta = exact_uncovered_z3_bv.find_uncovered(
        minimized_rows,
        max_component=args.max_component,
        limit=1,
        algebraic_primes=algebraic_primes,
        sophie_germain=sophie_germain,
    )
    if (
        primary_misses
        or primary_meta["sat"]
        or z3_misses
        or z3_meta["sat"]
    ):
        raise RuntimeError("minimized core failed final verification")

    payload = {
        "cover": str(args.cover),
        "row_indices": kept,
        "primes": [int(rows[index]["p"]) for index in kept],
        "raw_core_size": raw_core_size,
        "minimized_core_size": len(kept),
        "primary_verification": primary_meta,
        "independent_verification": z3_meta,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"PASS raw={raw_core_size} minimized={len(kept)} "
        f"primes={payload['primes']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
