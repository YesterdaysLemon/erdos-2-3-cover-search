#!/usr/bin/env python3
"""Independently replay a homogeneous-refinement sibling obstruction."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path


def prime_power_factors(n: int) -> list[tuple[int, int]]:
    result = []
    prime = 2
    while prime * prime <= n:
        exponent = 0
        while n % prime == 0:
            exponent += 1
            n //= prime
        if exponent:
            result.append((prime, exponent))
        prime += 1 if prime == 2 else 2
    if n > 1:
        result.append((n, 1))
    return result


def normal_class(
    a: int,
    b: int,
    modulus: int,
    factors: list[tuple[int, int]],
) -> tuple[tuple[int, int, int, int], ...]:
    """Canonicalize the equation normal, independently of direction code."""
    pieces = []
    for prime, exponent in factors:
        prime_power = prime**exponent
        left = a % prime_power
        right = b % prime_power
        if left % prime:
            pieces.append(
                (
                    prime,
                    exponent,
                    0,
                    right * pow(left, -1, prime_power) % prime_power,
                )
            )
        elif right % prime:
            pieces.append(
                (
                    prime,
                    exponent,
                    1,
                    left * pow(right, -1, prime_power) % prime_power,
                )
            )
        else:
            raise ValueError(
                f"normal {(a, b)} is not primitive modulo {modulus}"
            )
    return tuple(pieces)


def replay(rows: list[dict]) -> dict:
    factors_by_index: dict[int, list[tuple[int, int]]] = {}
    leaf_classes = set()
    child_classes = defaultdict(set)
    required = {}

    for row in rows:
        index = int(row["h"])
        factors = factors_by_index.setdefault(
            index,
            prime_power_factors(index),
        )
        child = normal_class(
            int(row["a"]),
            int(row["b"]),
            index,
            factors,
        )
        leaf_classes.add((index, child))
        for refinement_prime, exponent in factors:
            parent_index = index // refinement_prime
            parent_factors = [
                (
                    prime,
                    local_exponent - (
                        1 if prime == refinement_prime else 0
                    ),
                )
                for prime, local_exponent in factors
                if local_exponent
                - (1 if prime == refinement_prime else 0)
                > 0
            ]
            parent = normal_class(
                int(row["a"]),
                int(row["b"]),
                parent_index,
                parent_factors,
            )
            bucket = (index, refinement_prime, parent)
            child_classes[bucket].add(child)
            required[bucket] = (
                refinement_prime
                if exponent > 1
                else refinement_prime + 1
            )

    complete = []
    best = Fraction(0, 1)
    best_record = None
    refinement_primes = set()
    multi_child_buckets = 0
    max_observed_children = 0
    for bucket, children in child_classes.items():
        index, refinement_prime, parent = bucket
        need = required[bucket]
        observed = len(children)
        ratio = Fraction(observed, need)
        if ratio > best:
            best = ratio
            # The discovery certificate uses lattice directions (b,-a).
            # Rotating the normal (a,b) to that direction maps a normalized
            # normal chart to the opposite direction chart.  Recompute only
            # the scalar summary here; the parent key is diagnostic.
            best_record = {
                "child_index": index,
                "parent_index": index // refinement_prime,
                "refinement_prime": refinement_prime,
                "observed_children": observed,
                "required_children": need,
            }
        if observed == need:
            complete.append(
                {
                    "child_index": index,
                    "parent_index": index // refinement_prime,
                    "refinement_prime": refinement_prime,
                    "children": observed,
                }
            )
        refinement_primes.add(refinement_prime)
        multi_child_buckets += observed > 1
        max_observed_children = max(max_observed_children, observed)

    return {
        "rows": len(rows),
        "unique_leaf_lattices": len(leaf_classes),
        "duplicate_leaf_rows": len(rows) - len(leaf_classes),
        "parent_buckets": len(child_classes),
        "complete_sibling_groups": complete,
        "best_completion": best_record,
        "best_completion_fraction": {
            "numerator": best.numerator,
            "denominator": best.denominator,
        },
        "refinement_primes": len(refinement_primes),
        "multi_child_buckets": multi_child_buckets,
        "max_observed_children": max_observed_children,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--pool", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    if (
        certificate.get("schema")
        != "erdos203-homogeneous-refinement-obstruction-v1"
    ):
        raise RuntimeError("unsupported certificate schema")
    pool_path = args.pool
    if pool_path is None:
        pool_path = args.certificate.parent / Path(
            certificate["source_pool"]
        ).name
    raw_bytes = pool_path.read_bytes()
    source_hash = hashlib.sha256(raw_bytes).hexdigest()
    if source_hash != certificate["source_sha256"]:
        raise RuntimeError("source pool hash does not match certificate")
    source = json.loads(raw_bytes)
    if int(source.get("unresolved", -1)) != 0:
        raise RuntimeError("source does not report a complete factor scan")
    if int(source["max_h"]) != int(certificate["source_max_h"]):
        raise RuntimeError("source max_h does not match certificate")

    result = replay(source["choices"])
    expected = certificate["analysis"]
    comparable_keys = (
        "rows",
        "unique_leaf_lattices",
        "duplicate_leaf_rows",
        "parent_buckets",
        "best_completion_fraction",
        "refinement_primes",
        "multi_child_buckets",
        "max_observed_children",
    )
    mismatches = {
        key: {"certificate": expected[key], "replay": result[key]}
        for key in comparable_keys
        if expected[key] != result[key]
    }
    expected_complete_count = len(expected["complete_sibling_groups"])
    replay_complete_count = len(result["complete_sibling_groups"])
    verified = (
        not mismatches
        and expected_complete_count == replay_complete_count == 0
        and bool(certificate.get("obstruction"))
    )
    report = {
        "schema": "erdos203-homogeneous-refinement-verification-v1",
        "certificate": str(args.certificate),
        "source_pool": str(pool_path),
        "source_sha256": source_hash,
        "replay": result,
        "mismatches": mismatches,
        "verified": verified,
        "scope": certificate["scope"],
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"verified={verified} rows={result['rows']} "
        f"buckets={result['parent_buckets']} "
        f"complete={replay_complete_count} output={args.output}",
        flush=True,
    )
    return 0 if verified else 2


if __name__ == "__main__":
    raise SystemExit(main())
