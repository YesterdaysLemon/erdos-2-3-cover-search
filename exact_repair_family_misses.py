#!/usr/bin/env python3
"""Find exact holes common to every repair in an enumerated family."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import exact_greedy
import exact_uncovered
from exact_common_phase_misses import (
    build_union_rows,
    replay_common_misses,
)


def build_phase_maps(
    base: dict[int, int],
    repairs: list[dict],
) -> list[dict[int, int]]:
    phase_maps = []
    for repair in repairs:
        phases = dict(base)
        seen = set()
        for move in repair["moves"]:
            prime = int(move["p"])
            target = int(move["target"])
            if prime not in phases:
                raise KeyError(f"base phase omits prime {prime}")
            if prime in seen:
                raise ValueError(
                    f"repair updates prime {prime} more than once"
                )
            seen.add(prime)
            phases[prime] = target
        phase_maps.append(phases)
    return phase_maps


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("repair_report", type=Path)
    parser.add_argument("--max-component", type=int, default=256)
    parser.add_argument("--limit", type=int, default=1)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument(
        "--diversity-coordinate-moduli",
        default="",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    report = json.loads(args.repair_report.read_text())
    candidates = exact_greedy.load_candidates(args.pool, True)
    row_by_prime = {
        int(row["p"]): row for row in payload["choices"]
    }
    base_path = Path(report["base_phase"])
    base = {
        int(prime): int(target)
        for prime, target in json.loads(base_path.read_text()).items()
    }
    repairs = report["exact_repairs"]
    if len(repairs) != int(report["exact_repair_count"]):
        raise RuntimeError("repair report count is inconsistent")
    phase_maps = build_phase_maps(base, repairs)
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
        raise AssertionError("repair-family misses failed scalar replay")
    result = {
        "pool": str(args.pool),
        "repair_report": str(args.repair_report),
        "base_phase": str(base_path),
        "repair_count": len(repairs),
        "base_row_count": len(candidates),
        "union_row_count": len(union_rows),
        "additional_distinct_fibres": len(union_rows) - len(candidates),
        "diversity_coordinate_moduli": list(
            diversity_coordinate_moduli
        ),
        "checker": meta,
        "misses": [
            [exponent_k, exponent_l]
            for exponent_k, exponent_l in misses
        ],
        "scalar_replay": replayed,
        "scope": (
            "each point is missed by every exact finite repair in the "
            "declared complete owner-replay family"
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"repairs={len(repairs)} union_rows={len(union_rows)} "
        f"misses={len(misses)} sat={meta['sat']} replay={replayed}",
        flush=True,
    )
    return 1 if misses else 0


if __name__ == "__main__":
    raise SystemExit(main())
