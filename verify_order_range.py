#!/usr/bin/env python3
"""Independently verify every emitted row in an order-scan artifact.

By default the source must also report no unresolved cofactors.  The explicit
``--allow-unresolved-source`` mode verifies the emitted rows without treating
their incomplete source scan as an exclusion or completeness certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("range_file", type=Path)
    parser.add_argument("--known-pool", type=Path)
    parser.add_argument(
        "--allow-unresolved-source",
        action="store_true",
        help=(
            "verify emitted rows even when the source has unresolved "
            "cofactors; this does not certify source completeness"
        ),
    )
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import sympy  # type: ignore

    raw_bytes = args.range_file.read_bytes()
    payload = json.loads(raw_bytes)
    source_unresolved = int(payload.get("unresolved", -1))
    if source_unresolved and not args.allow_unresolved_source:
        raise RuntimeError("range contains unresolved cofactors")
    start_h = int(payload["start_h"])
    end_h = int(payload["end_h"])
    rows = payload["choices"]
    known_primes: set[int] = set()
    if args.known_pool:
        known = json.loads(args.known_pool.read_text())
        if int(known.get("unresolved", -1)):
            raise RuntimeError("known pool contains unresolved cofactors")
        known_primes = {int(row["p"]) for row in known["choices"]}

    seen: set[int] = set()
    density = 0.0
    for index, raw in enumerate(rows):
        p = int(raw["p"])
        h = int(raw["h"])
        a = int(raw["a"]) % h
        b = int(raw["b"]) % h
        stored_ord2 = int(raw["ord2"])
        stored_ord3 = int(raw["ord3"])
        if p in seen:
            raise RuntimeError(f"duplicate prime {p}")
        if p in known_primes:
            raise RuntimeError(f"prime {p} already occurs in the known pool")
        seen.add(p)
        if not sympy.isprime(p):
            raise RuntimeError(f"row {index} has composite p={p}")
        ord2 = int(sympy.n_order(2, p))
        ord3 = int(sympy.n_order(3, p))
        expected_h = math.lcm(ord2, ord3)
        if (stored_ord2, stored_ord3, h) != (ord2, ord3, expected_h):
            raise RuntimeError(
                f"row {index} order mismatch: "
                f"stored={(stored_ord2, stored_ord3, h)} "
                f"expected={(ord2, ord3, expected_h)}"
            )
        if not start_h <= h <= end_h:
            raise RuntimeError(f"row {index} has out-of-range h={h}")
        if math.gcd(a, b, h) != 1:
            raise RuntimeError(f"row {index} has non-surjective signature")
        if h // math.gcd(a, h) != ord2:
            raise RuntimeError(f"row {index} coefficient a has wrong order")
        if h // math.gcd(b, h) != ord3:
            raise RuntimeError(f"row {index} coefficient b has wrong order")
        if pow(2, b, p) != pow(3, a, p):
            raise RuntimeError(f"row {index} coefficients fail log relation")
        density += 1 / h

    print(
        f"PASS rows={len(rows)} start_h={start_h} end_h={end_h} "
        f"density={density:.12f} known_primes={len(known_primes)} "
        f"source_unresolved={source_unresolved} "
        f"complete_source={source_unresolved == 0} "
        f"sha256={hashlib.sha256(raw_bytes).hexdigest()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
