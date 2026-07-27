#!/usr/bin/env python3
"""Independently filter base-pool holes through a targeted pool extension.

This avoids putting a costly new prime-power component into the exact SAT
checker.  Every retained point is replayed with scalar modular arithmetic
against every row of the augmented pool and against the declared algebraic
exclusions.  It is a sound witness generator, not a completeness proof when
all supplied base holes are intercepted by the extension.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def row_target(row: dict, point: tuple[int, int]) -> int:
    k, l = point
    h = int(row["h"])
    return (int(row["a"]) * k + int(row["b"]) * l) % h


def verify_layered_misses(
    base: dict,
    augmented: dict,
    phases: dict[int, int],
    points: list[tuple[int, int]],
) -> dict:
    base_by_prime = {int(row["p"]): row for row in base["choices"]}
    augmented_by_prime = {
        int(row["p"]): row for row in augmented["choices"]
    }
    if len(base_by_prime) != len(base["choices"]):
        raise RuntimeError("base pool contains duplicate primes")
    if len(augmented_by_prime) != len(augmented["choices"]):
        raise RuntimeError("augmented pool contains duplicate primes")
    if base_by_prime.keys() - augmented_by_prime.keys():
        raise RuntimeError("augmented pool omits a base prime")
    for prime, row in base_by_prime.items():
        if augmented_by_prime[prime] != row:
            raise RuntimeError(f"augmented pool changes base row p={prime}")
    if augmented_by_prime.keys() - phases.keys():
        raise RuntimeError("phase file omits an augmented-pool prime")

    normalized_phases = {}
    for prime, row in augmented_by_prime.items():
        h = int(row["h"])
        modulus = int(row.get("target_modulus", 1))
        residue = int(row.get("target_residue", 0)) % modulus
        target = int(phases[prime]) % h
        if h % modulus or target % modulus != residue:
            raise RuntimeError(f"phase restriction fails for p={prime}")
        normalized_phases[prime] = target

    algebraic_primes = tuple(
        int(value) for value in augmented.get("algebraic_primes", ())
    )
    sophie_germain = bool(augmented.get("sophie_germain", False))
    extras = sorted(augmented_by_prime.keys() - base_by_prime.keys())
    retained = []
    intercepted = []
    for point in points:
        k, l = point
        if any(k % prime == 0 and l % prime == 0 for prime in algebraic_primes):
            raise RuntimeError("source point violates an algebraic exclusion")
        if sophie_germain and k % 4 == 2 and l % 4 == 0:
            raise RuntimeError("source point violates Sophie Germain exclusion")
        base_covers = [
            prime
            for prime, row in base_by_prime.items()
            if row_target(row, point) == normalized_phases[prime]
        ]
        if base_covers:
            raise RuntimeError(
                f"source point is covered by base primes {base_covers[:5]}"
            )
        extra_covers = [
            prime
            for prime in extras
            if row_target(augmented_by_prime[prime], point)
            == normalized_phases[prime]
        ]
        record = {
            "point": [k, l],
            "covering_extra_primes": extra_covers,
        }
        if extra_covers:
            intercepted.append(record)
        else:
            retained.append(record)
    return {
        "base_row_count": len(base_by_prime),
        "augmented_row_count": len(augmented_by_prime),
        "extra_primes": extras,
        "input_miss_count": len(points),
        "exact_augmented_miss_count": len(retained),
        "intercepted_count": len(intercepted),
        "misses": [record["point"] for record in retained],
        "exact_augmented_misses": [
            record["point"] for record in retained
        ],
        "intercepted": intercepted,
        "verified": True,
        "scope": (
            "each retained point is an exact scalar-replayed miss; an empty "
            "retained list does not prove the augmented phase covers"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("base_pool", type=Path)
    parser.add_argument("augmented_pool", type=Path)
    parser.add_argument("base_checker_result", type=Path)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    base_bytes = args.base_pool.read_bytes()
    augmented_bytes = args.augmented_pool.read_bytes()
    checker_bytes = args.base_checker_result.read_bytes()
    phase_bytes = args.phase_file.read_bytes()
    base = json.loads(base_bytes)
    augmented = json.loads(augmented_bytes)
    checker = json.loads(checker_bytes)
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(phase_bytes).items()
    }
    raw_points = checker.get("misses", checker.get("exact_misses"))
    if raw_points is None:
        raise RuntimeError("checker result contains no supported miss list")
    points = [tuple(map(int, point)) for point in raw_points]
    result = verify_layered_misses(base, augmented, phases, points)
    result.update(
        {
            "base_pool": str(args.base_pool),
            "augmented_pool": str(args.augmented_pool),
            "base_checker_result": str(args.base_checker_result),
            "phase_file": str(args.phase_file),
            "sha256": {
                "base_pool": hashlib.sha256(base_bytes).hexdigest(),
                "augmented_pool": hashlib.sha256(augmented_bytes).hexdigest(),
                "base_checker_result": hashlib.sha256(
                    checker_bytes
                ).hexdigest(),
                "phase_file": hashlib.sha256(phase_bytes).hexdigest(),
            },
        }
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"input={result['input_miss_count']} "
        f"exact_augmented={result['exact_augmented_miss_count']} "
        f"intercepted={result['intercepted_count']} "
        f"output={args.output}",
        flush=True,
    )
    return 1 if result["exact_augmented_miss_count"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
