#!/usr/bin/env python3
"""Scan ranked divisor-period families with verified block-star bounds."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

from certify_block_pair_anchor_star import triple_image_index_minors
from certify_block_plus_overlap_forest import (
    read_fraction,
    recorded_anchor_rows,
)
from certify_forced_pair_overlap import joint_index


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def determinant4(matrix: tuple[tuple[int, int, int, int], ...]) -> int:
    """Return a 4-by-4 determinant by cofactor expansion."""
    if len(matrix) != 4 or any(len(row) != 4 for row in matrix):
        raise ValueError("determinant4 requires a 4-by-4 matrix")

    def det3(rows: tuple[tuple[int, int, int], ...]) -> int:
        (a, b, c), (d, e, f), (g, h, i) = rows
        return (
            a * e * i + b * f * g + c * d * h
            - c * e * g - b * d * i - a * f * h
        )

    answer = 0
    for column in range(4):
        minor = tuple(
            tuple(row[index] for index in range(4) if index != column)
            for row in matrix[1:]
        )
        answer += (-1) ** column * matrix[0][column] * det3(minor)
    return answer


def quadruple_image_index_minors(rows: tuple[dict, dict, dict, dict]) -> int:
    """Index of the two-coordinate target map into four cyclic factors."""
    moduli = [int(row["h"]) for row in rows]
    columns = (
        (moduli[0], 0, 0, 0),
        (0, moduli[1], 0, 0),
        (0, 0, moduli[2], 0),
        (0, 0, 0, moduli[3]),
        tuple(int(row["a"]) for row in rows),
        tuple(int(row["b"]) for row in rows),
    )
    minors = [
        abs(
            determinant4(
                tuple(
                    tuple(columns[column][row] for column in selected)
                    for row in range(4)
                )
            )
        )
        for selected in itertools.combinations(range(6), 4)
    ]
    return math.gcd(*minors)


def best_anchor_lower(
    outside: dict,
    anchors: list[dict],
    *,
    fourway_triples: bool = False,
) -> Fraction:
    compatible = [
        anchor
        for anchor in anchors
        if joint_index(outside, anchor) == 1
    ]
    if not compatible:
        return Fraction(0)
    outside_h = int(outside["h"])
    singles = [
        Fraction(1, outside_h * int(anchor["h"]))
        for anchor in compatible
    ]
    triple_uppers = {
        (first, second): Fraction(
            triple_image_index_minors(
                (
                    outside,
                    compatible[first],
                    compatible[second],
                )
            ),
            outside_h
            * int(compatible[first]["h"])
            * int(compatible[second]["h"]),
        )
        for first in range(len(compatible))
        for second in range(first)
    }
    best = Fraction(0)
    subset_values = [Fraction(0)] * (1 << len(compatible))
    for mask in range(1, 1 << len(compatible)):
        bit = mask & -mask
        added = bit.bit_length() - 1
        previous = mask ^ bit
        lower = subset_values[previous] + singles[added]
        remaining = previous
        while remaining:
            other_bit = remaining & -remaining
            other = other_bit.bit_length() - 1
            lower -= triple_uppers[
                (max(added, other), min(added, other))
            ]
            remaining ^= other_bit
        subset_values[mask] = lower
        best = max(best, lower)

    if fourway_triples:
        for selected in itertools.combinations(range(len(compatible)), 3):
            first, second, third = selected
            lower = singles[first] + singles[second] + singles[third]
            lower -= triple_uppers[(max(first, second), min(first, second))]
            lower -= triple_uppers[(max(first, third), min(first, third))]
            lower -= triple_uppers[(max(second, third), min(second, third))]
            fixed_term = Fraction(
                1,
                outside_h
                * int(compatible[first]["h"])
                * int(compatible[second]["h"])
                * int(compatible[third]["h"]),
            )
            if lower + fixed_term <= best:
                continue
            if (
                quadruple_image_index_minors(
                    (
                        outside,
                        compatible[first],
                        compatible[second],
                        compatible[third],
                    )
                )
                == 1
            ):
                best = lower + fixed_term
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("ranking", type=Path)
    parser.add_argument(
        "--ranking-key",
        default="smallest_period",
        choices=("smallest_period", "fewest_rows", "highest_density"),
    )
    parser.add_argument(
        "--block-certificate",
        action="append",
        type=Path,
        default=[],
    )
    parser.add_argument("--block-catalog", type=Path)
    parser.add_argument("--lower-cache", type=Path)
    parser.add_argument(
        "--fourway-triples",
        action="store_true",
        help=(
            "restore the fixed four-fold term for jointly-surjective "
            "outside-plus-three-anchor target maps"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    ranking = json.loads(args.ranking.read_text())
    pool_sha256 = hashlib.sha256(args.pool.read_bytes()).hexdigest()
    lower_method = (
        "pair_bonferroni_plus_surjective_fourway"
        if args.fourway_triples
        else "pair_bonferroni"
    )
    cache = {
        "pool_sha256": pool_sha256,
        "lower_method": lower_method,
        "blocks": {},
    }
    if args.lower_cache and args.lower_cache.exists():
        loaded_cache = json.loads(args.lower_cache.read_text())
        if (
            loaded_cache.get("pool_sha256") == pool_sha256
            and loaded_cache.get("lower_method", "pair_bonferroni")
            == lower_method
        ):
            cache = loaded_cache
    block_paths = list(args.block_certificate)
    if args.block_catalog:
        catalog = json.loads(args.block_catalog.read_text())
        block_paths.extend(
            Path(
                entry["certificate"]
                if isinstance(entry, dict)
                else entry
            )
            for entry in catalog["blocks"]
        )
    if not block_paths:
        raise SystemExit("provide --block-certificate or --block-catalog")
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    if len(by_prime) != len(rows):
        raise RuntimeError("pool contains repeated fibre primes")

    blocks = []
    seen_block_signatures = set()
    for path in block_paths:
        payload = json.loads(path.read_text())
        anchors = recorded_anchor_rows(payload)
        anchor_primes = tuple(int(row["p"]) for row in anchors)
        for recorded in anchors:
            prime = int(recorded["p"])
            if prime not in by_prime or any(
                int(by_prime[prime][key]) != int(recorded[key])
                for key in ("h", "a", "b")
            ):
                raise RuntimeError(f"block {path} differs from source pool")
        loss = read_fraction(payload["forced_overlap_loss"])
        signature = (frozenset(anchor_primes), loss)
        if signature in seen_block_signatures:
            continue
        seen_block_signatures.add(signature)
        blocks.append(
            {
                "path": str(path),
                "anchor_primes": frozenset(anchor_primes),
                "anchors": [by_prime[p] for p in anchor_primes],
                "loss": loss,
            }
        )

    for block in blocks:
        anchor_primes = block["anchor_primes"]
        certificate_sha256 = hashlib.sha256(
            Path(block["path"]).read_bytes()
        ).hexdigest()
        cached = cache["blocks"].get(block["path"])
        if (
            cached
            and cached.get("certificate_sha256") == certificate_sha256
            and len(cached.get("row_lowers", {})) == len(rows)
        ):
            block["row_lowers"] = {
                int(prime): Fraction(int(value[0]), int(value[1]))
                for prime, value in cached["row_lowers"].items()
            }
        else:
            block["row_lowers"] = {
                int(row["p"]): (
                    Fraction(0)
                    if int(row["p"]) in anchor_primes
                    else best_anchor_lower(
                        row,
                        block["anchors"],
                        fourway_triples=args.fourway_triples,
                    )
                )
                for row in rows
            }
            cache["blocks"][block["path"]] = {
                "certificate_sha256": certificate_sha256,
                "row_lowers": {
                    str(prime): [value.numerator, value.denominator]
                    for prime, value in block["row_lowers"].items()
                },
            }
    if args.lower_cache:
        args.lower_cache.write_text(json.dumps(cache) + "\n")

    results = []
    resolved = 0
    for family in ranking[args.ranking_key]:
        period = int(family["period"])
        selected = [
            row for row in rows if period % int(row["h"]) == 0
        ]
        selected_primes = {int(row["p"]) for row in selected}
        total_density = sum(
            (Fraction(1, int(row["h"])) for row in selected),
            Fraction(0),
        )
        choices = []
        for block in blocks:
            if not block["anchor_primes"] <= selected_primes:
                continue
            star = sum(
                (
                    block["row_lowers"][int(row["p"])]
                    for row in selected
                    if int(row["p"]) not in block["anchor_primes"]
                ),
                Fraction(0),
            )
            upper = total_density - block["loss"] - star
            choices.append((upper, -block["loss"], block["path"], star))
        if choices:
            upper, negative_loss, block_path, star = min(choices)
            block_loss = -negative_loss
            proved = upper < 1
        else:
            upper = total_density
            block_loss = Fraction(0)
            star = Fraction(0)
            block_path = None
            proved = False
        resolved += int(proved)
        results.append(
            {
                "period": period,
                "components": family["components"],
                "rows": len(selected),
                "total_density": fraction_payload(total_density),
                "block_certificate": block_path,
                "block_overlap_loss": fraction_payload(block_loss),
                "pair_anchor_star_loss": fraction_payload(star),
                "union_upper_bound": fraction_payload(upper),
                "proved_no_cover": proved,
            }
        )

    unresolved = [
        result for result in results if not result["proved_no_cover"]
    ]
    unresolved.sort(
        key=lambda item: (
            item["union_upper_bound"]["decimal"],
            item["period"],
        )
    )
    output = {
        "pool": str(args.pool),
        "ranking": str(args.ranking),
        "ranking_key": args.ranking_key,
        "fourway_triples": args.fourway_triples,
        "lower_method": lower_method,
        "block_certificates": [block["path"] for block in blocks],
        "lower_cache": str(args.lower_cache) if args.lower_cache else None,
        "families_checked": len(results),
        "proved_no_cover": resolved,
        "unresolved_count": len(unresolved),
        "unresolved": unresolved,
        "results": results,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"families={len(results)} proved={resolved} "
        f"unresolved={len(unresolved)}",
        flush=True,
    )
    for result in unresolved[:20]:
        print(
            f"period={result['period']} rows={result['rows']} "
            f"density={result['total_density']['decimal']:.12f} "
            f"upper={result['union_upper_bound']['decimal']:.12f}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
