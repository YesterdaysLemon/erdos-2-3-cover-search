#!/usr/bin/env python3
"""Generate a reusable corrected low-subgroup-order prime pool."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import local_cover


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-limit", type=int, required=True)
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    candidates = local_cover.get_candidates(None, args.prime_limit, args.count)
    payload = {
        "prime_limit": args.prime_limit,
        "choices": [
            {
                "h": h,
                "p": p,
                "a": a,
                "b": b,
                "ord2": ord2,
                "ord3": ord3,
                "c": 0,
            }
            for h, p, a, b, ord2, ord3 in candidates
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"wrote={args.output} candidates={len(candidates)} "
        f"density={sum(1 / row[0] for row in candidates):.12f} "
        f"max_h={max(row[0] for row in candidates)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
