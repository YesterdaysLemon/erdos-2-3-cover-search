#!/usr/bin/env python3
"""Find exact exponent pairs missed by every supplied phase assignment.

For several phase maps on one affine pool, a common miss must avoid the union
of every distinct fibre selected by any map.  Recent sparse repairs differ on
few rows, so this union is usually only slightly larger than one phase even
when many historical repairs are included.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import exact_greedy
import exact_uncovered


def build_union_rows(
    candidates: list[tuple],
    row_by_prime: dict[int, dict],
    phase_maps: list[dict[int, int]],
) -> list[dict]:
    union: dict[tuple[int, int], dict] = {}
    for phases in phase_maps:
        for candidate in candidates:
            h, prime, *_rest = candidate
            source = row_by_prime[int(prime)]
            residue = int(source.get("target_residue", 0))
            modulus = int(source.get("target_modulus", 1))
            target = int(phases.get(int(prime), residue)) % int(h)
            if int(h) % modulus or target % modulus != residue:
                raise RuntimeError(
                    f"phase for p={int(prime)} violates target restriction"
                )
            row = exact_greedy.as_row(candidate, target)
            row["target_modulus"] = modulus
            row["target_residue"] = residue
            union[int(prime), target] = row
    return sorted(
        union.values(),
        key=lambda row: (
            int(row["h"]),
            int(row["p"]),
            int(row["c"]),
        ),
    )


def replay_common_misses(
    rows: list[dict],
    points: list[tuple[int, int]],
    algebraic_primes: tuple[int, ...],
    sophie_germain: bool,
) -> bool:
    for k, l in points:
        if any(k % prime == 0 and l % prime == 0 for prime in algebraic_primes):
            return False
        if sophie_germain and k % 4 == 2 and l % 4 == 0:
            return False
        if any(
            (
                int(row["a"]) * k
                + int(row["b"]) * l
                - int(row["c"])
            )
            % int(row["h"])
            == 0
            for row in rows
        ):
            return False
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("phase_files", nargs="+", type=Path)
    parser.add_argument("--max-component", type=int, default=256)
    parser.add_argument("--limit", type=int, default=1)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument(
        "--diversity-coordinate-moduli",
        default="",
        help=(
            "comma-separated prime-power coordinate moduli whose complete "
            "residue fingerprints are blocked between common misses"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    candidates = exact_greedy.load_candidates(args.pool, True)
    row_by_prime = {
        int(row["p"]): row for row in payload["choices"]
    }
    phase_maps = [
        {
            int(prime): int(target)
            for prime, target in json.loads(path.read_text()).items()
        }
        for path in args.phase_files
    ]
    union_rows = build_union_rows(
        candidates,
        row_by_prime,
        phase_maps,
    )
    algebraic_primes = tuple(
        int(value) for value in payload.get("algebraic_primes", ())
    )
    sophie_germain = bool(payload.get("sophie_germain", False))
    diversity_coordinate_moduli = tuple(
        int(value)
        for value in args.diversity_coordinate_moduli.split(",")
        if value
    )
    misses, meta = exact_uncovered.find_uncovered(
        union_rows,
        max_component=args.max_component,
        limit=args.limit,
        algebraic_primes=algebraic_primes,
        solver_name=args.solver,
        diversity_coordinate_moduli=diversity_coordinate_moduli,
        sophie_germain=sophie_germain,
    )
    replayed = replay_common_misses(
        union_rows,
        misses,
        algebraic_primes,
        sophie_germain,
    )
    if not replayed:
        raise AssertionError("common misses failed scalar replay")
    result = {
        "pool": str(args.pool),
        "phase_files": [str(path) for path in args.phase_files],
        "phase_count": len(phase_maps),
        "base_row_count": len(candidates),
        "union_row_count": len(union_rows),
        "additional_distinct_fibres": len(union_rows) - len(candidates),
        "diversity_coordinate_moduli": list(
            diversity_coordinate_moduli
        ),
        "checker": meta,
        "misses": [[k, l] for k, l in misses],
        "scalar_replay": replayed,
        "scope": (
            "each point is missed by every supplied phase; UNSAT only means "
            "their selected fibre union covers the declared exact domain"
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"phases={len(phase_maps)} base_rows={len(candidates)} "
        f"union_rows={len(union_rows)} misses={len(misses)} "
        f"sat={meta['sat']} replay={replayed} output={args.output}",
        flush=True,
    )
    return 1 if misses else 0


if __name__ == "__main__":
    raise SystemExit(main())
