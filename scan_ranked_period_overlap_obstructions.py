#!/usr/bin/env python3
"""Apply the forced-overlap matching bound to ranked period families."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from certify_forced_pair_overlap import joint_index


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("ranking", type=Path)
    parser.add_argument("--ranking-key", default="smallest_period")
    parser.add_argument("--limit", type=int, default=200)
    parser.add_argument(
        "--block-certificate",
        type=Path,
        action="append",
        default=[],
        help=(
            "repeatable verified small-anchor block certificate whose "
            "overlap loss is combined with disjoint pair intersections"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.limit < 1:
        raise SystemExit("--limit must be positive")

    source = json.loads(args.pool.read_text())
    ranking = json.loads(args.ranking.read_text())
    rows = source["choices"]
    source_by_prime = {int(row["p"]): row for row in rows}
    blocks = []
    for block_path in args.block_certificate:
        block = json.loads(block_path.read_text())
        recorded_rows = list(block.get("anchor_rows", ()))
        if not recorded_rows and block.get("base_certificate"):
            base_path = Path(block["base_certificate"])
            base = json.loads(base_path.read_text())
            recorded_rows.extend(base["anchor_rows"])
            recorded_rows.append(block["extra_row"])
        for recorded in recorded_rows:
            prime = int(recorded["p"])
            if prime not in source_by_prime or any(
                int(source_by_prime[prime][key]) != int(recorded[key])
                for key in ("h", "a", "b")
            ):
                raise RuntimeError(
                    "block certificate anchor does not match source pool"
                )
        blocks.append(
            {
                "path": str(block_path),
                "primes": {
                    int(prime) for prime in block["anchor_primes"]
                },
                "loss": Fraction(
                    int(block["forced_overlap_loss"]["numerator"]),
                    int(block["forced_overlap_loss"]["denominator"]),
                ),
            }
        )
    families = ranking[args.ranking_key][: args.limit]
    results = []
    for family in families:
        period = int(family["period"])
        selected_rows = [
            row for row in rows if period % int(row["h"]) == 0
        ]
        total_density = sum(
            (Fraction(1, int(row["h"])) for row in selected_rows),
            Fraction(0),
        )
        edges = sorted(
            (
                (
                    Fraction(1, int(left["h"]) * int(right["h"])),
                    min(int(left["p"]), int(right["p"])),
                    max(int(left["p"]), int(right["p"])),
                    left,
                    right,
                )
                for index, left in enumerate(selected_rows)
                for right in selected_rows[:index]
                if joint_index(left, right) == 1
            ),
            key=lambda edge: (-edge[0], edge[1], edge[2]),
        )
        present_primes = {int(row["p"]) for row in selected_rows}
        block_options = [
            {"path": None, "primes": set(), "loss": Fraction(0)},
            *(
                block
                for block in blocks
                if block["primes"] <= present_primes
            ),
        ]
        best = None
        for block in block_options:
            used = set(block["primes"])
            overlap = block["loss"]
            pairs = []
            for weight, _first, _second, left, right in edges:
                left_prime = int(left["p"])
                right_prime = int(right["p"])
                if left_prime in used or right_prime in used:
                    continue
                used.update((left_prime, right_prime))
                overlap += weight
                pairs.append(
                    {
                        "primes": [left_prime, right_prime],
                        "moduli": [int(left["h"]), int(right["h"])],
                        "overlap": fraction_payload(weight),
                    }
                )
            upper = total_density - overlap
            candidate = {
                "block": block,
                "pairs": pairs,
                "overlap": overlap,
                "upper": upper,
            }
            if best is None or candidate["upper"] < best["upper"]:
                best = candidate
        assert best is not None
        overlap = best["overlap"]
        upper = best["upper"]
        pairs = best["pairs"]
        chosen_block = best["block"]
        results.append(
            {
                "period": period,
                "components": family["components"],
                "rows": len(selected_rows),
                "total_density": fraction_payload(total_density),
                "selected_pairs": pairs,
                "block_certificate_used": (
                    chosen_block["path"]
                ),
                "block_overlap_loss": fraction_payload(
                    chosen_block["loss"]
                ),
                "forced_overlap_sum": fraction_payload(overlap),
                "union_upper_bound": fraction_payload(upper),
                "proved_no_cover": upper < 1,
            }
        )

    proved_count = sum(item["proved_no_cover"] for item in results)
    unresolved = [item for item in results if not item["proved_no_cover"]]
    payload = {
        "pool": str(args.pool),
        "ranking": str(args.ranking),
        "ranking_key": args.ranking_key,
        "block_certificate": (
            [str(path) for path in args.block_certificate]
        ),
        "families_checked": len(results),
        "proved_no_cover": proved_count,
        "unresolved": len(unresolved),
        "first_unresolved": unresolved[0] if unresolved else None,
        "results": results,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"checked={len(results)} proved={proved_count} "
        f"unresolved={len(unresolved)}",
        flush=True,
    )
    for item in unresolved[:10]:
        print(
            f"unresolved period={item['period']} rows={item['rows']} "
            f"density={item['total_density']['decimal']:.12f} "
            f"upper={item['union_upper_bound']['decimal']:.12f}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
