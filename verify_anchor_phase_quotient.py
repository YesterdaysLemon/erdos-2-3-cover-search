#!/usr/bin/env python3
"""Independently replay a frozen-remainder anchor quotient certificate."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path


EXPECTED_COLUMNS = [
    "p",
    "h",
    "a",
    "b",
    "base_phase",
    "target_modulus",
    "target_residue",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_certificate(certificate: dict) -> dict:
    if certificate.get("schema_version") != 1:
        raise ValueError("unsupported certificate schema")
    if certificate.get("row_columns") != EXPECTED_COLUMNS:
        raise ValueError("unexpected compact-row schema")

    rows = []
    row_by_prime = {}
    for raw_row in certificate["rows"]:
        if not isinstance(raw_row, list) or len(raw_row) != len(EXPECTED_COLUMNS):
            raise ValueError("malformed compact row")
        prime, h, a, b, phase, modulus, residue = map(int, raw_row)
        if prime in row_by_prime:
            raise ValueError(f"duplicate row prime {prime}")
        if h < 1 or modulus < 1 or h % modulus:
            raise ValueError(f"invalid modulus data for p={prime}")
        if math.gcd(math.gcd(a, b), h) != 1:
            raise ValueError(f"nonsurjective affine row for p={prime}")
        a %= h
        b %= h
        phase %= h
        residue %= modulus
        if phase % modulus != residue:
            raise ValueError(f"forbidden base phase for p={prime}")
        row = (prime, h, a, b, phase, modulus, residue)
        rows.append(row)
        row_by_prime[prime] = row

    anchors = certificate["anchors"]
    anchor_primes = tuple(int(anchor["p"]) for anchor in anchors)
    if len(set(anchor_primes)) != len(anchor_primes):
        raise ValueError("duplicate anchor prime")
    if set(anchor_primes) - row_by_prime.keys():
        raise ValueError("anchor row is absent")

    legal_spaces = []
    for anchor, prime in zip(anchors, anchor_primes):
        _p, h, a, b, _phase, modulus, residue = row_by_prime[prime]
        legal = tuple(
            target for target in range(h) if target % modulus == residue
        )
        expected_anchor = {
            "p": prime,
            "h": h,
            "a": a,
            "b": b,
            "target_modulus": modulus,
            "target_residue": residue,
            "legal_targets": list(legal),
        }
        if anchor != expected_anchor:
            raise ValueError(f"anchor metadata mismatch for p={prime}")
        legal_spaces.append(legal)

    expected_branches = list(itertools.product(*legal_spaces))
    branches = certificate["branches"]
    if len(branches) != len(expected_branches):
        raise ValueError("certificate does not contain every legal branch")

    algebraic_primes = tuple(
        int(prime) for prime in certificate.get("algebraic_primes", ())
    )
    if len(set(algebraic_primes)) != len(algebraic_primes):
        raise ValueError("duplicate algebraic prime")
    verified_witnesses = set()
    minimum_claimed_multiplicity = None
    maximum_claimed_multiplicity = None
    for expected_targets, branch in zip(expected_branches, branches):
        targets = tuple(map(int, branch["targets"]))
        if targets != expected_targets:
            raise ValueError("branches are missing, repeated, or out of order")
        raw_witness = branch["witness"]
        if not isinstance(raw_witness, list) or len(raw_witness) != 2:
            raise ValueError("malformed branch witness")
        k, l = map(int, raw_witness)
        if k < 0 or l < 0:
            raise ValueError("branch witness has a negative exponent")
        if any(
            k % prime == 0 and l % prime == 0
            for prime in algebraic_primes
        ):
            raise ValueError("branch witness is algebraically covered")

        branch_target = dict(zip(anchor_primes, targets))
        for prime, h, a, b, base_phase, _modulus, _residue in rows:
            target = branch_target.get(prime, base_phase)
            if (a * k + b * l - target) % h == 0:
                raise ValueError(
                    f"branch witness is covered by row p={prime}"
                )
        multiplicity = int(branch["available_witnesses"])
        if multiplicity < 1:
            raise ValueError("claimed witness multiplicity is not positive")
        minimum_claimed_multiplicity = (
            multiplicity
            if minimum_claimed_multiplicity is None
            else min(minimum_claimed_multiplicity, multiplicity)
        )
        maximum_claimed_multiplicity = (
            multiplicity
            if maximum_claimed_multiplicity is None
            else max(maximum_claimed_multiplicity, multiplicity)
        )
        verified_witnesses.add((k, l))

    summary = certificate["summary"]
    if int(summary["rows"]) != len(rows):
        raise ValueError("summary row count mismatch")
    if int(summary["anchor_rows"]) != len(anchor_primes):
        raise ValueError("summary anchor count mismatch")
    if int(summary["nonanchor_rows"]) != len(rows) - len(anchor_primes):
        raise ValueError("summary nonanchor count mismatch")
    if int(summary["legal_anchor_branches"]) != len(expected_branches):
        raise ValueError("summary branch count mismatch")
    if int(summary["open_anchor_branches"]) != 0:
        raise ValueError("certificate claims an open branch")
    if int(summary["minimum_available_witnesses"]) != (
        minimum_claimed_multiplicity
    ):
        raise ValueError("summary minimum multiplicity mismatch")
    if int(summary["maximum_available_witnesses"]) != (
        maximum_claimed_multiplicity
    ):
        raise ValueError("summary maximum multiplicity mismatch")

    return {
        "verified": True,
        "rows": len(rows),
        "anchor_rows": len(anchor_primes),
        "nonanchor_rows": len(rows) - len(anchor_primes),
        "legal_anchor_branches": len(expected_branches),
        "verified_branch_witnesses": len(branches),
        "distinct_recorded_witnesses": len(verified_witnesses),
        "algebraic_primes": list(algebraic_primes),
        "engine": "independent-python-modular-replay",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    result = verify_certificate(certificate)
    result["certificate"] = args.certificate.name
    result["certificate_sha256"] = sha256_file(args.certificate)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"VERIFIED rows={result['rows']} "
        f"branches={result['legal_anchor_branches']} "
        f"witnesses={result['distinct_recorded_witnesses']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
