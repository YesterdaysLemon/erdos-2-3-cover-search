#!/usr/bin/env python3
"""Independently verify a bounded homogeneous non-cover witness.

The discovery program checks stored signature equations.  This verifier
instead performs modular exponentiation from the source primes themselves:
the homogeneous fibre for ``p`` contains ``(k,l)`` exactly when
``2^k * 3^l == 1 (mod p)``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def modular_hits(
    rows: list[dict],
    k: int,
    l: int,
) -> list[int]:
    """Return source primes whose homogeneous fibres contain ``(k,l)``."""
    hits = []
    for row in rows:
        prime = int(row["p"])
        if pow(2, k, prime) * pow(3, l, prime) % prime == 1:
            hits.append(prime)
    return hits


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--pool", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    if (
        certificate.get("schema")
        != "erdos203-bounded-homogeneous-noncover-v1"
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
    rows = source["choices"]
    k = int(certificate["witness"]["k"])
    l = int(certificate["witness"]["l"])
    hits = modular_hits(rows, k, l)
    checks = {
        "source_complete": int(source.get("unresolved", -1)) == 0,
        "source_max_h_matches": (
            int(source["max_h"]) == int(certificate["source_max_h"])
        ),
        "source_row_count_matches": (
            len(rows) == int(certificate["source_rows"])
        ),
        "certificate_claims_zero_signature_hits": (
            int(certificate["signature_hit_count"]) == 0
        ),
        "independent_modular_hit_count_zero": not hits,
        "certificate_claims_obstruction": bool(
            certificate.get("obstruction")
        ),
    }
    verified = all(checks.values())
    report = {
        "schema": "erdos203-bounded-homogeneous-noncover-verification-v1",
        "certificate": str(args.certificate),
        "source_pool": str(pool_path),
        "source_sha256": source_hash,
        "source_rows": len(rows),
        "witness": {"k": k, "l": l},
        "independent_modular_hit_count": len(hits),
        "independent_modular_hit_primes": hits,
        "checks": checks,
        "verified": verified,
        "scope": certificate["scope"],
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"verified={verified} rows={len(rows)} modular_hits={len(hits)} "
        f"witness=({k},{l}) output={args.output}",
        flush=True,
    )
    return 0 if verified else 2


if __name__ == "__main__":
    raise SystemExit(main())
