#!/usr/bin/env python3
"""Emit a finite witness missed by every homogeneous fibre in a pool.

For a prime-signature row ``(h, p, a, b)``, the homogeneous fibre is

    a*k + b*l == 0 (mod h).

If one exponent pair misses every row in the complete supplied pool, then
the union of all of those homogeneous fibres is not a cover.  Consequently
no homogeneous subfamily of that pool can be a cover either.

This is a bounded statement only.  It says nothing about affine phases,
higher-index fibres, or homogeneous fibres outside the supplied pool.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path


DEFAULT_SEED = 2_031_050_000
DEFAULT_BITS = 63
DEFAULT_ATTEMPTS = 20_000


def signature_hits(
    rows: list[dict],
    k: int,
    l: int,
) -> list[tuple[int, int]]:
    """Return ``(h,p)`` for every signature lattice containing ``(k,l)``."""
    hits = []
    for row in rows:
        modulus = int(row["h"])
        value = (
            int(row["a"]) * (k % modulus)
            + int(row["b"]) * (l % modulus)
        ) % modulus
        if value == 0:
            hits.append((modulus, int(row["p"])))
    return hits


def find_witness(
    rows: list[dict],
    *,
    seed: int = DEFAULT_SEED,
    bits: int = DEFAULT_BITS,
    attempts: int = DEFAULT_ATTEMPTS,
) -> tuple[int, int, int] | None:
    """Deterministically search for a point outside every supplied lattice."""
    rng = random.Random(seed)
    for attempt in range(1, attempts + 1):
        k = rng.getrandbits(bits)
        l = rng.getrandbits(bits)
        if not signature_hits(rows, k, l):
            return k, l, attempt
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--bits", type=int, default=DEFAULT_BITS)
    parser.add_argument("--attempts", type=int, default=DEFAULT_ATTEMPTS)
    args = parser.parse_args()

    if args.bits < 1 or args.attempts < 1:
        raise SystemExit("--bits and --attempts must be positive")
    raw_bytes = args.pool.read_bytes()
    source = json.loads(raw_bytes)
    if int(source.get("unresolved", -1)) != 0:
        raise RuntimeError(
            "the source must report zero unresolved factorization cofactors"
        )
    rows = source["choices"]
    witness = find_witness(
        rows,
        seed=args.seed,
        bits=args.bits,
        attempts=args.attempts,
    )
    if witness is None:
        print(
            f"rows={len(rows)} witness=not-found attempts={args.attempts}",
            flush=True,
        )
        return 2
    k, l, attempt = witness
    certificate = {
        "schema": "erdos203-bounded-homogeneous-noncover-v1",
        "source_pool": str(args.pool),
        "source_sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "source_max_h": int(source["max_h"]),
        "source_unresolved": int(source["unresolved"]),
        "source_rows": len(rows),
        "search": {
            "seed": args.seed,
            "coordinate_bits": args.bits,
            "attempt_limit": args.attempts,
            "witness_attempt": attempt,
        },
        "witness": {"k": k, "l": l},
        "signature_hit_count": 0,
        "obstruction": True,
        "argument": (
            "The displayed exponent pair lies in none of the homogeneous "
            "prime-signature lattices in the complete supplied finite pool. "
            "Their full union is therefore not a cover, so no homogeneous "
            "subfamily of this pool can cover Z^2."
        ),
        "scope": (
            "Finite pool only; this does not obstruct affine phase "
            "assignments, homogeneous fibres outside the source bound, or "
            "higher-index fibres."
        ),
    }
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    print(
        f"rows={len(rows)} witness=({k},{l}) attempt={attempt} "
        f"signature_hits=0 output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
