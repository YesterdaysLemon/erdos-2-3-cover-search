#!/usr/bin/env python3
"""Independent verifier for derived_lower_digit_core.py certificates."""

from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path


def valuation(value: int, prime: int) -> int:
    exponent = 0
    while value % prime == 0:
        exponent += 1
        value //= prime
    return exponent


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_pool", type=Path)
    parser.add_argument("result", type=Path)
    args = parser.parse_args()

    source = json.loads(args.input_pool.read_text())
    result = json.loads(args.result.read_text())
    rows = source["choices"]
    certificate = result["derived_lower_digit_core"]
    algebraic_primes = {
        int(value) for value in certificate["algebraic_primes"]
    }
    audit_by_round = collections.defaultdict(list)
    for record in certificate["audit"]:
        if record["reason"] != "derived_lower_digit_all_states_impossible":
            raise RuntimeError(f"unsupported reason {record['reason']}")
        audit_by_round[int(record["round"])].append(record)

    alive = set(range(len(rows)))
    verified_states = 0
    removed_total = 0
    for round_no in sorted(audit_by_round):
        round_remove = set()
        for record in audit_by_round[round_no]:
            prime = int(record["prime"])
            exponent = int(record["exponent"])
            maximal_exponent = max(
                (
                    valuation(int(rows[index]["h"]), prime)
                    for index in alive
                ),
                default=0,
            )
            if exponent != maximal_exponent:
                raise AssertionError(
                    f"round {round_no} q={prime}: non-maximal exponent"
                )
            indices = [
                index
                for index in alive
                if valuation(int(rows[index]["h"]), prime) == exponent
            ]
            if len(indices) != int(record["rows"]):
                raise AssertionError(
                    f"round {round_no} q={prime}: row-count mismatch"
                )
            if exponent < 2:
                raise AssertionError("lower-digit peel requires exponent >= 2")
            if any(
                valuation(int(rows[index]["target_modulus"]), prime)
                >= exponent
                for index in indices
            ):
                raise AssertionError("a removed row fixes its top digit")

            state_count = 0
            for x in range(prime):
                for y in range(prime):
                    if prime in algebraic_primes and x == 0 and y == 0:
                        continue
                    available = []
                    for index in indices:
                        row = rows[index]
                        target_valuation = valuation(
                            int(row["target_modulus"]),
                            prime,
                        )
                        if not target_valuation or (
                            int(row["a"]) * x
                            + int(row["b"]) * y
                            - int(row["target_residue"])
                        ) % prime == 0:
                            available.append(index)
                    # Fewer than q affine lines contain fewer than q^2
                    # points, regardless of directions or offsets.
                    if len(available) >= prime:
                        raise AssertionError(
                            f"round {round_no} q={prime} state=({x},{y}) "
                            f"has {len(available)} available rows"
                        )
                    state_count += 1
            expected_states = prime * prime - int(
                prime in algebraic_primes
            )
            if state_count != expected_states:
                raise AssertionError("required-state count mismatch")
            verified_states += state_count
            round_remove.update(indices)
        alive.difference_update(round_remove)
        removed_total += len(round_remove)

    kept = [row for index, row in enumerate(rows) if index in alive]
    if kept != result["choices"]:
        raise AssertionError("final survivor list does not match replay")
    if len(rows) != int(certificate["input_rows"]):
        raise AssertionError("input count mismatch")
    if len(kept) != int(certificate["output_rows"]):
        raise AssertionError("output count mismatch")
    print(
        f"PASS input={len(rows)} removed={removed_total} "
        f"survivors={len(kept)} states={verified_states}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
