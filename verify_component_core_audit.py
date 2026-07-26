#!/usr/bin/env python3
"""Independently replay the recorded finite-plane component-core eliminations.

This verifier intentionally does not import ``component_core``.  It rebuilds
the maximal prime-power groups from the original pool, recomputes their
direction capacities, checks the stated finite-plane inequality, removes the
recorded groups round by round, and compares the survivors with the output
pool.
"""

from __future__ import annotations

import argparse
import collections
import itertools
import json
from fractions import Fraction
from pathlib import Path


def smallest_prime_factors(limit: int) -> list[int]:
    spf = list(range(limit + 1))
    if limit >= 1:
        spf[1] = 1
    for value in range(2, int(limit**0.5) + 1):
        if spf[value] != value:
            continue
        for multiple in range(value * value, limit + 1, value):
            if spf[multiple] == multiple:
                spf[multiple] = value
    return spf


def factor(value: int, spf: list[int]) -> dict[int, int]:
    result: dict[int, int] = {}
    while value > 1:
        prime = spf[value]
        result[prime] = result.get(prime, 0) + 1
        value //= prime
    return result


def factor_by_trial_division(value: int) -> dict[int, int]:
    result: dict[int, int] = {}
    divisor = 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            result[divisor] = result.get(divisor, 0) + 1
            value //= divisor
        divisor += 1 if divisor == 2 else 2
    if value > 1:
        result[value] = result.get(value, 0) + 1
    return result


def direction(row: dict, prime: int) -> tuple[int, int]:
    a = int(row["a"]) % prime
    b = int(row["b"]) % prime
    if a:
        return 1, b * pow(a, -1, prime) % prime
    if b:
        return 0, 1
    raise AssertionError(
        f"row p={row['p']} h={row['h']} is not surjective mod {prime}"
    )


def tangent_bound(counts: collections.Counter, prime: int) -> bool:
    total = sum(counts.values())
    if prime == 2:
        return False
    if max(counts.values()) >= prime or len(counts) >= prime + 1:
        return False
    capacities = sorted(counts.values(), reverse=True)
    prefix = [0]
    for capacity in capacities:
        prefix.append(prefix[-1] + capacity)
    minimum_size = 3 * (prime + 1) // 2
    for size in range(minimum_size, total + 2):
        occupied = min(size - prime, len(capacities))
        if prefix[occupied] >= size - 1:
            return False
    return True


def parallel_bound(counts: collections.Counter, prime: int) -> bool:
    total = sum(counts.values())
    if total < prime:
        return True
    # q parallel lines cover the whole affine plane.  The inequality below
    # is valid only for a proper parallel class; raw multiplicities above q
    # merely force duplicate offsets and cannot prove impossibility.
    if max(counts.values()) >= prime:
        return False
    return any(
        total * prime - capacity * (total - capacity) < prime * prime
        for capacity in counts.values()
    )


def brute_essential_rows(
    prime: int,
    row_directions: list[tuple[int, int]],
    excluded_points: frozenset[tuple[int, int]],
) -> set[int]:
    """Independently enumerate rows essential in some top-plane cover."""
    points = [
        (x, y)
        for x in range(prime)
        for y in range(prime)
        if (x, y) not in excluded_points
    ]
    line_points: dict[tuple[int, int], tuple[int, ...]] = {}
    for row_index, (kind, slope) in enumerate(row_directions):
        for offset in range(prime):
            line_points[(row_index, offset)] = tuple(
                point_index
                for point_index, (x, y) in enumerate(points)
                if (
                    y if kind == 0 else (x + slope * y) % prime
                )
                == offset
            )

    essential: set[int] = set()
    for offsets in itertools.product(
        range(prime), repeat=len(row_directions)
    ):
        covers = [0] * len(points)
        for row_index, offset in enumerate(offsets):
            for point_index in line_points[(row_index, offset)]:
                covers[point_index] += 1
        if 0 in covers:
            continue
        for row_index, offset in enumerate(offsets):
            if any(
                covers[point_index] == 1
                for point_index in line_points[(row_index, offset)]
            ):
                essential.add(row_index)
        if len(essential) == len(row_directions):
            break
    return essential


def row_key(row: dict) -> tuple[int, ...]:
    return tuple(
        int(row[name])
        for name in ("h", "p", "a", "b", "ord2", "ord3")
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("result", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--require-density-below-one",
        action="store_true",
        help=(
            "after replaying every proof-safe elimination, also require "
            "the exact reciprocal density of the survivors to be below one"
        ),
    )
    args = parser.parse_args()

    source = json.loads(args.source.read_text())
    result = json.loads(args.result.read_text())
    rows = source["choices"]
    audit = result["component_core"]["audit"]
    max_h = max((int(row["h"]) for row in rows), default=1)
    if max_h <= 5_000_000:
        spf = smallest_prime_factors(max_h)
        factors = [factor(int(row["h"]), spf) for row in rows]
    else:
        factors = [
            factor_by_trial_division(int(row["h"]))
            for row in rows
        ]
    alive = set(range(len(rows)))

    by_round: dict[int, list[dict]] = collections.defaultdict(list)
    for record in audit:
        by_round[int(record["round"])].append(record)
    expected_rounds = list(range(1, max(by_round, default=0) + 1))
    if sorted(by_round) != expected_rounds:
        raise AssertionError("audit round numbers are not contiguous")

    reason_counts: collections.Counter = collections.Counter()
    algebraic_primes = {
        int(prime)
        for prime in result["component_core"].get(
            "algebraic_primes", ()
        )
    }
    for round_no in expected_rounds:
        groups: dict[int, dict[int, list[int]]] = collections.defaultdict(
            lambda: collections.defaultdict(list)
        )
        for index in alive:
            for prime, exponent in factors[index].items():
                groups[prime][exponent].append(index)

        removed: set[int] = set()
        seen_components: set[int] = set()
        for record in by_round[round_no]:
            prime = int(record["prime"])
            exponent = int(record["exponent"])
            if prime in seen_components:
                raise AssertionError(
                    f"round {round_no} records component {prime} twice"
                )
            seen_components.add(prime)
            if prime not in groups or exponent != max(groups[prime]):
                raise AssertionError(
                    f"round {round_no} component {prime}^{exponent} "
                    "is not the maximal live exponent"
                )
            indices = groups[prime][exponent]
            counts = collections.Counter(direction(rows[index], prime) for index in indices)
            if len(indices) != int(record["rows"]):
                raise AssertionError(f"round {round_no} row count mismatch")
            if len(counts) != int(record["directions"]):
                raise AssertionError(f"round {round_no} direction count mismatch")
            if max(counts.values()) != int(record["max_parallel"]):
                raise AssertionError(f"round {round_no} parallel count mismatch")

            reason = str(record["reason"])
            if reason == "blokhuis_brouwer_tangent_capacity_bound":
                valid = tangent_bound(counts, prime)
                removed_indices = set(indices)
            elif reason == "parallel_class_union_bound":
                valid = parallel_bound(counts, prime)
                removed_indices = set(indices)
            elif reason == "never_essential_in_any_top_plane_cover":
                if int(record.get("removed_rows", 0)) != 1:
                    raise AssertionError(
                        "inessential-row record must remove exactly one row"
                    )
                removed_prime = int(record["removed_prime"])
                matching = [
                    index
                    for index in indices
                    if int(rows[index]["p"]) == removed_prime
                ]
                if len(matching) != 1:
                    raise AssertionError(
                        "recorded inessential prime is not unique in group"
                    )
                row_directions = [
                    direction(rows[index], prime) for index in indices
                ]
                essential = brute_essential_rows(
                    prime,
                    row_directions,
                    (
                        frozenset({(0, 0)})
                        if prime in algebraic_primes
                        else frozenset()
                    ),
                )
                position = indices.index(matching[0])
                valid = position not in essential
                removed_indices = {matching[0]}
            else:
                raise AssertionError(f"unsupported audit reason: {reason}")
            if not valid:
                raise AssertionError(
                    f"round {round_no} reason fails for {prime}^{exponent}"
                )
            reason_counts[reason] += 1
            removed.update(removed_indices)

        if not removed:
            raise AssertionError(f"round {round_no} records no eliminations")
        alive.difference_update(removed)

    survivor_keys = sorted(row_key(rows[index]) for index in alive)
    result_keys = sorted(row_key(row) for row in result["choices"])
    if survivor_keys != result_keys:
        raise AssertionError("replayed survivors differ from result choices")
    metadata = result["component_core"]
    if len(rows) != int(metadata["input_rows"]):
        raise AssertionError("input row metadata mismatch")
    if len(alive) != int(metadata["output_rows"]):
        raise AssertionError("output row metadata mismatch")
    survivor_density = sum(
        (
            Fraction(1, int(rows[index]["h"]))
            for index in alive
        ),
        Fraction(0),
    )
    density_below_one = survivor_density < 1
    if args.require_density_below_one and not density_below_one:
        raise AssertionError(
            "replayed survivor reciprocal density is not below one"
        )

    summary = {
        "source": str(args.source),
        "result": str(args.result),
        "input_rows": len(rows),
        "rounds": len(expected_rounds),
        "records": len(audit),
        "survivors": len(alive),
        "reason_counts": dict(reason_counts),
        "survivors_match": True,
        "survivor_density_numerator": survivor_density.numerator,
        "survivor_density_denominator": survivor_density.denominator,
        "survivor_density_decimal": float(survivor_density),
        "survivor_density_below_one": density_below_one,
        "proved_no_cover": (
            args.require_density_below_one and density_below_one
        ),
        "verified": True,
    }
    if args.output:
        args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"PASS input={len(rows)} rounds={len(expected_rounds)} "
        f"records={len(audit)} survivors={len(alive)} "
        f"tangent={reason_counts['blokhuis_brouwer_tangent_capacity_bound']} "
        f"parallel={reason_counts['parallel_class_union_bound']} "
        f"inessential={reason_counts['never_essential_in_any_top_plane_cover']} "
        f"density={survivor_density.numerator}/"
        f"{survivor_density.denominator} "
        f"proved_no_cover={summary['proved_no_cover']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
