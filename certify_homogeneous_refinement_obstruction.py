#!/usr/bin/env python3
"""Certify the absence of a bottom sibling group for homogeneous refinement.

A nontrivial finite refinement of the trivial homogeneous lattice cover has
an internal node of greatest depth.  Every child of that node is a leaf, so
the available leaf lattices must contain all q-descendants of one parent for
some prime q.  This program checks that necessary condition in a finite prime
fibre pool.

The emitted obstruction is finite only.  It says nothing about fibres beyond
the supplied pool and nothing about general affine (nonhomogeneous) covers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path


def factor(n: int) -> dict[int, int]:
    factors: dict[int, int] = {}
    divisor = 2
    while divisor * divisor <= n:
        while n % divisor == 0:
            factors[divisor] = factors.get(divisor, 0) + 1
            n //= divisor
        divisor += 1 if divisor == 2 else 2
    if n > 1:
        factors[n] = 1
    return factors


def projective_direction_key(
    c: int,
    d: int,
    modulus: int,
    factorization: dict[int, int] | None = None,
) -> tuple[tuple[int, int, int, int], ...]:
    """Canonicalize a primitive direction in P^1(Z/modulus Z).

    The fibre equation ``a*k+b*l=0`` has lattice direction ``(b,-a)``.
    Projective CRT components are normalized by whichever coordinate is a
    unit at the component prime.
    """
    if modulus == 1:
        return ()
    local = []
    for prime, exponent in (
        factorization if factorization is not None else factor(modulus)
    ).items():
        prime_power = prime**exponent
        first = c % prime_power
        second = d % prime_power
        if math.gcd(first, prime) == 1:
            local.append(
                (
                    prime,
                    exponent,
                    0,
                    second * pow(first, -1, prime_power) % prime_power,
                )
            )
        else:
            if math.gcd(second, prime) != 1:
                raise ValueError(
                    f"direction {(c, d)} is not primitive modulo {modulus}"
                )
            local.append(
                (
                    prime,
                    exponent,
                    1,
                    first * pow(second, -1, prime_power) % prime_power,
                )
            )
    return tuple(local)


def reduce_direction_key(
    key: tuple[tuple[int, int, int, int], ...],
    prime: int,
) -> tuple[tuple[int, int, int, int], ...]:
    """Reduce a projective class from index M to M/prime."""
    reduced = []
    for local_prime, exponent, chart, value in key:
        if local_prime != prime:
            reduced.append((local_prime, exponent, chart, value))
        elif exponent > 1:
            smaller_power = local_prime ** (exponent - 1)
            reduced.append(
                (
                    local_prime,
                    exponent - 1,
                    chart,
                    value % smaller_power,
                )
            )
    return tuple(reduced)


def analyze_rows(rows: list[dict]) -> dict:
    factor_cache: dict[int, dict[int, int]] = {}
    leaf_states: set[
        tuple[int, tuple[tuple[int, int, int, int], ...]]
    ] = set()
    buckets: dict[
        tuple[
            int,
            int,
            tuple[tuple[int, int, int, int], ...],
        ],
        set[tuple[tuple[int, int, int, int], ...]],
    ] = defaultdict(set)
    required_by_bucket = {}

    for row in rows:
        index = int(row["h"])
        factors = factor_cache.setdefault(index, factor(index))
        child_key = projective_direction_key(
            int(row["b"]),
            -int(row["a"]),
            index,
            factors,
        )
        leaf_states.add((index, child_key))
        for refinement_prime, exponent in factors.items():
            parent_index = index // refinement_prime
            parent_key = reduce_direction_key(
                child_key,
                refinement_prime,
            )
            bucket = (
                index,
                refinement_prime,
                parent_key,
            )
            buckets[bucket].add(child_key)
            required_by_bucket[bucket] = (
                refinement_prime
                if exponent >= 2
                else refinement_prime + 1
            )

    complete_groups = []
    best_fraction = Fraction(0, 1)
    best_bucket = None
    refinement_primes = set()
    multi_child_buckets = 0
    max_observed_children = 0
    for bucket, children in buckets.items():
        child_index, refinement_prime, parent_key = bucket
        required = required_by_bucket[bucket]
        observed = len(children)
        ratio = Fraction(observed, required)
        if ratio > best_fraction:
            best_fraction = ratio
            best_bucket = {
                "child_index": child_index,
                "parent_index": child_index // refinement_prime,
                "refinement_prime": refinement_prime,
                "observed_children": observed,
                "required_children": required,
                "parent_key": [list(item) for item in parent_key],
            }
        if observed == required:
            complete_groups.append(
                {
                    "child_index": child_index,
                    "parent_index": child_index // refinement_prime,
                    "refinement_prime": refinement_prime,
                    "children": observed,
                    "parent_key": [list(item) for item in parent_key],
                }
            )
        refinement_primes.add(refinement_prime)
        multi_child_buckets += observed > 1
        max_observed_children = max(max_observed_children, observed)

    return {
        "rows": len(rows),
        "unique_leaf_lattices": len(leaf_states),
        "duplicate_leaf_rows": len(rows) - len(leaf_states),
        "parent_buckets": len(buckets),
        "complete_sibling_groups": complete_groups,
        "best_completion": best_bucket,
        "best_completion_fraction": {
            "numerator": best_fraction.numerator,
            "denominator": best_fraction.denominator,
        },
        "refinement_primes": len(refinement_primes),
        "multi_child_buckets": multi_child_buckets,
        "max_observed_children": max_observed_children,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    raw_bytes = args.pool.read_bytes()
    source = json.loads(raw_bytes)
    if int(source.get("unresolved", -1)) != 0:
        raise RuntimeError(
            "the source must report zero unresolved factorization cofactors"
        )
    rows = source["choices"]
    analysis = analyze_rows(rows)
    certificate = {
        "schema": "erdos203-homogeneous-refinement-obstruction-v1",
        "source_pool": str(args.pool),
        "source_sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "source_max_h": int(source["max_h"]),
        "source_unresolved": int(source["unresolved"]),
        "analysis": analysis,
        "obstruction": not analysis["complete_sibling_groups"],
        "argument": (
            "A nontrivial finite refinement tree has a deepest internal "
            "node, all of whose q-descendants are leaves. The source leaf "
            "lattices contain no complete q-sibling group, so no homogeneous "
            "refinement of the trivial cover can use only this finite pool."
        ),
        "scope": (
            "Finite pool only; this does not obstruct general affine covers, "
            "non-refinement homogeneous covers, or higher-index fibres."
        ),
    }
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    print(
        f"rows={analysis['rows']} "
        f"unique_leaves={analysis['unique_leaf_lattices']} "
        f"buckets={analysis['parent_buckets']} "
        f"complete={len(analysis['complete_sibling_groups'])} "
        f"obstruction={certificate['obstruction']} "
        f"output={args.output}",
        flush=True,
    )
    return 0 if certificate["obstruction"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
