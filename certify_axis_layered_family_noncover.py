#!/usr/bin/env python3
"""Certify that no active-class placement of a layered row family can cover.

The certificate has two exact stages.

First, a prime-deficit argument iteratively removes residual moduli.  If
fewer than ``q`` current moduli are divisible by prime ``q``, the
``q``-divisible residue classes are redundant in every one-dimensional
cover, regardless of their phases.

Second, every remaining row chooses one active layer class.  A necessary
condition for a cover is reciprocal residual capacity at least one in every
layer column.  A reduced multi-valued decision diagram (MDD) exhausts all
active-class choices against a supplied finite column core.  Terminal zero
is an exact obstruction, not a floating-point MILP result.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def prime_factors(value: int) -> tuple[int, ...]:
    if value < 1:
        raise ValueError("modulus must be positive")
    factors = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1 if divisor == 2 else 2
    if value > 1:
        factors.append(value)
    return tuple(factors)


def layer_rows(payload: dict) -> list[dict]:
    axis = payload["layer_axis"]
    if axis not in {"k", "l"}:
        raise ValueError("layer axis must be k or l")
    rows = []
    seen = set()
    for raw in payload["choices"]:
        prime = int(raw["p"])
        if prime in seen:
            raise ValueError("duplicate source prime")
        seen.add(prime)
        h = int(raw["h"])
        a = int(raw["a"]) % h
        b = int(raw["b"]) % h
        if h < 1 or math.gcd(a, b, h) != 1:
            raise ValueError(f"invalid affine row for p={prime}")
        if axis == "k":
            active_modulus = math.gcd(b, h)
            active_coefficient = a
        else:
            active_modulus = math.gcd(a, h)
            active_coefficient = b
        if math.gcd(active_coefficient, active_modulus) != 1:
            raise AssertionError("active-coordinate map is not surjective")
        rows.append(
            {
                "p": prime,
                "h": h,
                "active_modulus": active_modulus,
                "residual_modulus": h // active_modulus,
            }
        )
    return rows


def prune_prime_deficits(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    current = list(rows)
    rounds = []
    while True:
        counts = Counter(
            prime
            for row in current
            for prime in prime_factors(row["residual_modulus"])
        )
        deficient = sorted(
            prime for prime, count in counts.items() if count < prime
        )
        if not deficient:
            return current, rounds
        prime = deficient[0]
        removed = [
            row
            for row in current
            if row["residual_modulus"] % prime == 0
        ]
        if not removed or len(removed) != counts[prime]:
            raise AssertionError("prime-deficit removal count mismatch")
        rounds.append(
            {
                "prime": prime,
                "row_count_before": len(current),
                "divisible_row_count": len(removed),
                "divisible_count_below_prime": len(removed) < prime,
                "removed_primes": [row["p"] for row in removed],
                "row_count_after": len(current) - len(removed),
            }
        )
        removed_primes = {row["p"] for row in removed}
        current = [
            row for row in current if row["p"] not in removed_primes
        ]


def exact_mdd_capacity_obstruction(
    rows: list[dict],
    columns: list[int],
) -> dict:
    if not rows or not columns:
        raise ValueError("MDD rows and columns must be nonempty")
    period = math.lcm(*(row["active_modulus"] for row in rows))
    normalized_columns = [column % period for column in columns]
    if len(set(normalized_columns)) != len(normalized_columns):
        raise ValueError("column core contains duplicates modulo its period")
    scale = math.lcm(*(row["residual_modulus"] for row in rows))
    fixed_weight = sum(
        scale // row["residual_modulus"]
        for row in rows
        if row["active_modulus"] == 1
    )
    flexible = [
        {
            **row,
            "scaled_weight": scale // row["residual_modulus"],
        }
        for row in rows
        if row["active_modulus"] > 1
    ]
    flexible.sort(
        key=lambda row: (
            -row["scaled_weight"],
            -row["active_modulus"],
            row["p"],
        )
    )
    required = scale - fixed_weight
    if required <= 0:
        raise RuntimeError("fixed rows already meet every capacity constraint")

    domains = [row["active_modulus"] for row in flexible]
    weights = [row["scaled_weight"] for row in flexible]
    variable_count = len(flexible)
    suffix = [0] * (variable_count + 1)
    for index in range(variable_count - 1, -1, -1):
        suffix[index] = suffix[index + 1] + weights[index]

    # Terminals are 0=false and 1=true.  Every other node is reduced by
    # merging equal children and hash-consed by (variable level, children).
    levels = [variable_count, variable_count]
    children: list[tuple[int, ...]] = [(), ()]
    unique: dict[tuple[int, tuple[int, ...]], int] = {}

    def make_node(level: int, node_children: tuple[int, ...]) -> int:
        if all(child == node_children[0] for child in node_children):
            return node_children[0]
        key = (level, node_children)
        known = unique.get(key)
        if known is not None:
            return known
        node = len(levels)
        unique[key] = node
        levels.append(level)
        children.append(node_children)
        return node

    def build_column(column: int) -> int:
        @functools.lru_cache(maxsize=None)
        def build(level: int, needed: int) -> int:
            if needed <= 0:
                return 1
            if level == variable_count or suffix[level] < needed:
                return 0
            miss = build(level + 1, needed)
            hit = build(level + 1, needed - weights[level])
            hit_class = column % domains[level]
            return make_node(
                level,
                tuple(
                    hit if active_class == hit_class else miss
                    for active_class in range(domains[level])
                ),
            )

        return build(0, required)

    @functools.lru_cache(maxsize=None)
    def conjunction(left: int, right: int) -> int:
        if left == 0 or right == 0:
            return 0
        if left == 1:
            return right
        if right == 1 or left == right:
            return left
        if left > right:
            left, right = right, left
        left_level = levels[left]
        right_level = levels[right]
        level = min(left_level, right_level)
        domain = domains[level]
        left_children = (
            children[left] if left_level == level else (left,) * domain
        )
        right_children = (
            children[right] if right_level == level else (right,) * domain
        )
        return make_node(
            level,
            tuple(
                conjunction(a, b)
                for a, b in zip(
                    left_children,
                    right_children,
                    strict=True,
                )
            ),
        )

    root = 1
    processed_columns = 0
    for column in normalized_columns:
        root = conjunction(root, build_column(column))
        processed_columns += 1
        if root == 0:
            break
    if root != 0:
        raise RuntimeError("supplied column core is MDD-feasible")
    return {
        "engine": "exact-reduced-multivalued-decision-diagram",
        "capacity_period": period,
        "capacity_scale": scale,
        "fixed_scaled_weight": fixed_weight,
        "required_flexible_scaled_weight": required,
        "row_count": len(rows),
        "fixed_row_count": len(rows) - len(flexible),
        "flexible_row_count": len(flexible),
        "flexible_order": [
            {
                "p": row["p"],
                "active_modulus": row["active_modulus"],
                "residual_modulus": row["residual_modulus"],
                "scaled_weight": row["scaled_weight"],
            }
            for row in flexible
        ],
        "core_columns": normalized_columns,
        "processed_columns": processed_columns,
        "terminal_root": root,
        "mdd_node_count_including_terminals": len(levels),
        "unique_nonterminal_node_count": len(unique),
        "apply_cache_entries": conjunction.cache_info().currsize,
        "proved_no_capacity_placement": True,
    }


def fraction_record(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def build_certificate(
    source_path: Path,
    payload: dict,
    columns: list[int],
) -> dict:
    source_rows = layer_rows(payload)
    surviving, pruning = prune_prime_deficits(source_rows)
    raw_density = sum(
        (Fraction(1, row["h"]) for row in surviving),
        Fraction(),
    )
    mdd = exact_mdd_capacity_obstruction(surviving, columns)
    return {
        "schema": "axis_layered_family_prime_pruned_mdd_obstruction_v1",
        "source": str(source_path),
        "source_sha256": sha256(source_path),
        "layer_axis": payload["layer_axis"],
        "coordinate_basis": payload.get("coordinate_basis"),
        "layer_period": int(payload["layer_period"]),
        "initial_row_count": len(source_rows),
        "prime_deficit_pruning": pruning,
        "surviving_row_count": len(surviving),
        "surviving_rows": surviving,
        "surviving_raw_reciprocal_density": fraction_record(raw_density),
        "capacity_obstruction": mdd,
        "proved_no_active_class_placement": True,
        "proved_no_declared_layered_family_cover": True,
        "proof_chain": [
            (
                "at each recorded prime q, fewer than q residual moduli "
                "contain q, so those rows are redundant in every column cover"
            ),
            (
                "after all proof-safe removals, any cover would still need "
                "one active class per surviving row with residual capacity "
                "at least one in every column"
            ),
            (
                "the exact reduced MDD exhausts every such active-class "
                "assignment and reaches terminal false on the column core"
            ),
        ],
        "scope": (
            "all active-class placements of exactly the declared finite "
            "layered row family; this does not rule out other directions, "
            "periods, source rows beyond the bounded pool, or the original "
            "infinite problem"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument(
        "--core-columns",
        required=True,
        help="comma-separated layer columns forming an exact MDD obstruction",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    columns = [
        int(value) for value in args.core_columns.split(",") if value
    ]
    payload = json.loads(args.input.read_text())
    if payload.get("schema") != "axis_layered_pool_v1":
        raise RuntimeError("input is not an axis-layered pool artifact")
    certificate = build_certificate(args.input, payload, columns)
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    mdd = certificate["capacity_obstruction"]
    print(
        f"initial_rows={certificate['initial_row_count']} "
        f"surviving_rows={certificate['surviving_row_count']} "
        f"pruning_rounds={len(certificate['prime_deficit_pruning'])} "
        f"core_columns={len(mdd['core_columns'])} "
        f"mdd_nodes={mdd['mdd_node_count_including_terminals']} "
        f"terminal_root={mdd['terminal_root']} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
