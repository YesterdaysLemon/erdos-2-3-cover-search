#!/usr/bin/env python3
"""Search for a common projective linear form across many prime fibres.

If (a_p,b_p) is a unit multiple of one fixed (A,B) modulo h_p, then every
prime fibre is a congruence class of n=A*k+B*l.  A one-dimensional covering
system in n would solve the full two-dimensional problem.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

import cegis_cover
import search_cover


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-limit", type=int, default=2_000_000)
    parser.add_argument("--max-h", type=int, default=5_000)
    parser.add_argument("--bound", type=int, default=2_000)
    parser.add_argument("--top", type=int, default=20)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    candidates = []
    for p in search_cover.sieve(args.prime_limit):
        ord2 = search_cover.multiplicative_order(2, p)
        ord3 = search_cover.multiplicative_order(3, p)
        h = math.lcm(ord2, ord3)
        if h <= args.max_h:
            h, a, b, ord2, ord3 = cegis_cover.signature(p, ord2, ord3)
            candidates.append((h, p, a, b, ord2, ord3))
    candidates.sort()

    score = np.zeros((args.bound, args.bound), dtype=np.float64)
    for h, p, a, b, ord2, ord3 in candidates:
        weight = 1.0 / h
        for unit in range(h):
            if math.gcd(unit, h) != 1:
                continue
            A = unit * a % h
            B = unit * b % h
            score[A::h, B::h] += weight

    flat = score.ravel()
    count = min(args.top, flat.size)
    indices = np.argpartition(flat, -count)[-count:]
    indices = indices[np.argsort(flat[indices])[::-1]]
    print(
        f"candidates={len(candidates)} total_density="
        f"{sum(1/item[0] for item in candidates):.12f}"
    )
    for index in indices:
        A, B = divmod(int(index), args.bound)
        compatible = []
        for item in candidates:
            h, p, a, b, ord2, ord3 = item
            if any(
                (unit * a - A) % h == 0 and (unit * b - B) % h == 0
                for unit in range(h)
                if math.gcd(unit, h) == 1
            ):
                compatible.append(item)
        print(
            f"A={A} B={B} density={flat[index]:.12f} "
            f"primes={len(compatible)} rows={compatible}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
