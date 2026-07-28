#!/usr/bin/env python3
"""Independently replay a prime-pruned layered-family MDD obstruction."""

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


def factor_set(value: int) -> set[int]:
    if value < 1:
        raise AssertionError("modulus must be positive")
    result = set()
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            result.add(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1 if divisor == 2 else 2
    if value > 1:
        result.add(value)
    return result


def reconstruct_rows(source: dict) -> list[dict]:
    axis = source["layer_axis"]
    if axis not in {"k", "l"}:
        raise AssertionError("invalid layer axis")
    seen = set()
    rows = []
    for raw in source["choices"]:
        prime = int(raw["p"])
        if prime in seen:
            raise AssertionError("duplicate row prime")
        seen.add(prime)
        h = int(raw["h"])
        if h < 1:
            raise AssertionError("row modulus must be positive")
        a = int(raw["a"]) % h
        b = int(raw["b"]) % h
        if math.gcd(a, b, h) != 1:
            raise AssertionError("row map is not surjective")
        active_modulus = math.gcd(b if axis == "k" else a, h)
        active_coefficient = a if axis == "k" else b
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


def replay_pruning(
    rows: list[dict],
    claimed_rounds: list[dict],
) -> list[dict]:
    current = list(rows)
    for expected_round in claimed_rounds:
        counts = Counter(
            prime
            for row in current
            for prime in factor_set(row["residual_modulus"])
        )
        deficient = sorted(
            prime for prime, count in counts.items() if count < prime
        )
        if not deficient:
            raise AssertionError("certificate contains an extra pruning round")
        prime = deficient[0]
        removed = [
            row
            for row in current
            if row["residual_modulus"] % prime == 0
        ]
        actual = {
            "prime": prime,
            "row_count_before": len(current),
            "divisible_row_count": len(removed),
            "divisible_count_below_prime": len(removed) < prime,
            "removed_primes": [row["p"] for row in removed],
            "row_count_after": len(current) - len(removed),
        }
        if expected_round != actual:
            raise AssertionError("prime-deficit pruning round mismatch")
        removed_primes = {row["p"] for row in removed}
        current = [
            row for row in current if row["p"] not in removed_primes
        ]
    counts = Counter(
        prime
        for row in current
        for prime in factor_set(row["residual_modulus"])
    )
    if any(count < prime for prime, count in counts.items()):
        raise AssertionError("certificate stopped pruning too early")
    return current


def replay_mdd(rows: list[dict], expected: dict) -> dict:
    period = math.lcm(*(row["active_modulus"] for row in rows))
    scale = math.lcm(*(row["residual_modulus"] for row in rows))
    fixed_weight = sum(
        scale // row["residual_modulus"]
        for row in rows
        if row["active_modulus"] == 1
    )
    flexible = []
    for row in rows:
        if row["active_modulus"] == 1:
            continue
        flexible.append(
            {
                **row,
                "scaled_weight": scale // row["residual_modulus"],
            }
        )
    flexible.sort(
        key=lambda row: (
            -row["scaled_weight"],
            -row["active_modulus"],
            row["p"],
        )
    )
    required = scale - fixed_weight
    if required <= 0:
        raise AssertionError("fixed rows already meet the capacity threshold")
    domains = [row["active_modulus"] for row in flexible]
    weights = [row["scaled_weight"] for row in flexible]
    variable_count = len(flexible)
    suffix = [0] * (variable_count + 1)
    for index in range(variable_count - 1, -1, -1):
        suffix[index] = suffix[index + 1] + weights[index]
    columns = [int(value) % period for value in expected["core_columns"]]
    if len(columns) != len(set(columns)):
        raise AssertionError("column core contains duplicates")

    levels = [variable_count, variable_count]
    children: list[tuple[int, ...]] = [(), ()]
    unique: dict[tuple[int, tuple[int, ...]], int] = {}

    def node(level: int, branches: tuple[int, ...]) -> int:
        if all(branch == branches[0] for branch in branches):
            return branches[0]
        key = (level, branches)
        known = unique.get(key)
        if known is not None:
            return known
        identifier = len(levels)
        unique[key] = identifier
        levels.append(level)
        children.append(branches)
        return identifier

    def threshold_root(column: int) -> int:
        @functools.lru_cache(maxsize=None)
        def visit(level: int, needed: int) -> int:
            if needed <= 0:
                return 1
            if level == variable_count or suffix[level] < needed:
                return 0
            miss = visit(level + 1, needed)
            hit = visit(level + 1, needed - weights[level])
            selected = column % domains[level]
            return node(
                level,
                tuple(
                    hit if value == selected else miss
                    for value in range(domains[level])
                ),
            )

        return visit(0, required)

    @functools.lru_cache(maxsize=None)
    def intersect(left: int, right: int) -> int:
        if left == 0 or right == 0:
            return 0
        if left == 1:
            return right
        if right == 1 or left == right:
            return left
        if left > right:
            left, right = right, left
        level = min(levels[left], levels[right])
        domain = domains[level]
        left_branches = (
            children[left]
            if levels[left] == level
            else (left,) * domain
        )
        right_branches = (
            children[right]
            if levels[right] == level
            else (right,) * domain
        )
        return node(
            level,
            tuple(
                intersect(a, b)
                for a, b in zip(
                    left_branches,
                    right_branches,
                    strict=True,
                )
            ),
        )

    root = 1
    processed = 0
    for column in columns:
        root = intersect(root, threshold_root(column))
        processed += 1
        if root == 0:
            break
    actual = {
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
        "core_columns": columns,
        "processed_columns": processed,
        "terminal_root": root,
        "mdd_node_count_including_terminals": len(levels),
        "unique_nonterminal_node_count": len(unique),
        "apply_cache_entries": intersect.cache_info().currsize,
    }
    for key, value in actual.items():
        if expected.get(key) != value:
            raise AssertionError(f"MDD replay mismatch for {key}")
    if (
        expected.get("engine")
        != "exact-reduced-multivalued-decision-diagram"
        or expected.get("proved_no_capacity_placement") is not True
        or root != 0
    ):
        raise AssertionError("MDD did not prove capacity infeasibility")
    return actual


def replay(source: dict, certificate: dict) -> dict:
    if source.get("schema") != "axis_layered_pool_v1":
        raise AssertionError("source schema mismatch")
    if (
        certificate.get("schema")
        != "axis_layered_family_prime_pruned_mdd_obstruction_v1"
    ):
        raise AssertionError("certificate schema mismatch")
    if (
        certificate["layer_axis"] != source["layer_axis"]
        or certificate.get("coordinate_basis")
        != source.get("coordinate_basis")
        or int(certificate["layer_period"])
        != int(source["layer_period"])
    ):
        raise AssertionError("certificate metadata does not match source")
    rows = reconstruct_rows(source)
    if int(certificate["initial_row_count"]) != len(rows):
        raise AssertionError("initial row count mismatch")
    surviving = replay_pruning(
        rows,
        certificate["prime_deficit_pruning"],
    )
    if (
        int(certificate["surviving_row_count"]) != len(surviving)
        or certificate["surviving_rows"] != surviving
    ):
        raise AssertionError("surviving row family mismatch")
    density = sum(
        (Fraction(1, row["h"]) for row in surviving),
        Fraction(),
    )
    density_record = certificate["surviving_raw_reciprocal_density"]
    stored_density = Fraction(
        int(density_record["numerator"]),
        int(density_record["denominator"]),
    )
    if stored_density != density or float(density_record["decimal"]) != float(
        density
    ):
        raise AssertionError("surviving raw density mismatch")
    mdd = replay_mdd(surviving, certificate["capacity_obstruction"])
    if not (
        certificate.get("proved_no_active_class_placement") is True
        and certificate.get("proved_no_declared_layered_family_cover")
        is True
    ):
        raise AssertionError("certificate does not declare its proved claims")
    return {
        "verified": True,
        "initial_row_count": len(rows),
        "prime_deficit_pruning_rounds": len(
            certificate["prime_deficit_pruning"]
        ),
        "surviving_row_count": len(surviving),
        "capacity_period": mdd["capacity_period"],
        "capacity_scale": mdd["capacity_scale"],
        "core_column_count": len(mdd["core_columns"]),
        "mdd_node_count_including_terminals": mdd[
            "mdd_node_count_including_terminals"
        ],
        "terminal_root": mdd["terminal_root"],
        "proved_no_active_class_placement": True,
        "proved_no_declared_layered_family_cover": True,
        "engine": "independent-prime-pruning-and-exact-mdd-replay",
        "scope": certificate["scope"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    certificate = json.loads(args.certificate.read_text())
    source_path = Path(certificate["source"])
    if sha256(source_path) != certificate["source_sha256"]:
        raise AssertionError("source SHA-256 mismatch")
    source = json.loads(source_path.read_text())
    report = replay(source, certificate)
    report.update(
        {
            "certificate": str(args.certificate),
            "certificate_sha256": sha256(args.certificate),
            "source": str(source_path),
            "source_sha256": certificate["source_sha256"],
        }
    )
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"verified={report['verified']} "
        f"initial_rows={report['initial_row_count']} "
        f"surviving_rows={report['surviving_row_count']} "
        f"core_columns={report['core_column_count']} "
        f"mdd_nodes={report['mdd_node_count_including_terminals']} "
        f"terminal_root={report['terminal_root']} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
