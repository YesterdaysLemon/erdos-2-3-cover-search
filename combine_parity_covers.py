#!/usr/bin/env python3
"""Map three exact induced parity covers back to original prime fibres."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cegis_cover
import exact_uncovered
import exact_uncovered_z3
import search_cover


def load_and_map(cover_path: Path, pool_path: Path) -> list[dict]:
    cover = json.loads(cover_path.read_text())
    if cover.get("checker", {}).get("sat", True):
        raise RuntimeError(f"{cover_path} is not an exact checker-UNSAT cover")
    pool = json.loads(pool_path.read_text())
    metadata = {int(row["p"]): row for row in pool["choices"]}
    mapped = []
    for row in cover["choices"]:
        p = int(row["p"])
        source = metadata[p]
        ord2 = search_cover.multiplicative_order(2, p)
        ord3 = search_cover.multiplicative_order(3, p)
        h, a, b, ord2, ord3 = cegis_cover.signature(p, ord2, ord3)
        if (h, a, b) != (
            int(source["original_h"]),
            int(source["original_a"]),
            int(source["original_b"]),
        ):
            raise RuntimeError(f"stale original signature metadata for prime {p}")
        scale = int(source["target_scale"])
        c = scale * int(row["c"]) % h
        mapped.append(
            {
                "h": h,
                "p": p,
                "a": a,
                "b": b,
                "ord2": ord2,
                "ord3": ord3,
                "c": c,
            }
        )
    return mapped


def main() -> int:
    parser = argparse.ArgumentParser()
    for name in ("k", "l", "sum"):
        parser.add_argument(f"--{name}-cover", type=Path, required=True)
        parser.add_argument(f"--{name}-pool", type=Path, required=True)
    parser.add_argument("--max-component", type=int, default=100000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    seen = set()
    for cover_path, pool_path in (
        (args.k_cover, args.k_pool),
        (args.l_cover, args.l_pool),
        (args.sum_cover, args.sum_pool),
    ):
        for row in load_and_map(cover_path, pool_path):
            if row["p"] in seen:
                raise RuntimeError(f"prime {row['p']} is allocated twice")
            seen.add(row["p"])
            rows.append(row)

    missed, primary = exact_uncovered.find_uncovered(
        rows, max_component=args.max_component, limit=1
    )
    if missed:
        raise RuntimeError(f"combined parity rows miss {missed[0]}")
    missed_z3, secondary = exact_uncovered_z3.find_uncovered(
        rows, max_component=args.max_component, limit=1
    )
    if missed_z3:
        raise RuntimeError(f"Z3 checker found a miss {missed_z3[0]}")
    payload = {
        "construction": "three disjoint induced parity covers",
        "choices": rows,
        "checker": primary,
        "independent_checker": secondary,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"wrote={args.output} rows={len(rows)} "
        f"primary=UNSAT secondary=UNSAT"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
