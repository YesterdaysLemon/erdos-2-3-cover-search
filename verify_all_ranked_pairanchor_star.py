#!/usr/bin/env python3
"""Independently verify an aggregate ranked pair-anchor-star scan."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def recorded_anchor_rows(block: dict) -> list[dict]:
    if block.get("anchor_rows"):
        return list(block["anchor_rows"])
    if block.get("base_certificate"):
        base = json.loads(Path(block["base_certificate"]).read_text())
        return [*recorded_anchor_rows(base), block["extra_row"]]
    raise RuntimeError("block certificate does not expose its anchor rows")


def pair_index(first: dict, second: dict) -> int:
    h1, a1, b1 = (int(first[key]) for key in ("h", "a", "b"))
    h2, a2, b2 = (int(second[key]) for key in ("h", "a", "b"))
    return math.gcd(
        abs(h1 * h2),
        abs(h1 * a2),
        abs(h1 * b2),
        abs(h2 * a1),
        abs(h2 * b1),
        abs(a1 * b2 - a2 * b1),
    )


def determinant(
    first: tuple[int, int, int],
    second: tuple[int, int, int],
    third: tuple[int, int, int],
) -> int:
    a, b, c = first
    d, e, f = second
    g, h, i = third
    return (
        a * e * i + b * f * g + c * d * h
        - c * e * g - b * d * i - a * f * h
    )


def triple_index(rows: tuple[dict, dict, dict]) -> int:
    moduli = [int(row["h"]) for row in rows]
    columns = (
        (moduli[0], 0, 0),
        (0, moduli[1], 0),
        (0, 0, moduli[2]),
        tuple(int(row["a"]) for row in rows),
        tuple(int(row["b"]) for row in rows),
    )
    values = [
        abs(determinant(columns[i], columns[j], columns[k]))
        for i, j, k in itertools.combinations(range(5), 3)
    ]
    return math.gcd(*values)


def determinant4_explicit(matrix: tuple[tuple[int, int, int, int], ...]) -> int:
    """Independent Leibniz-form 4-by-4 determinant."""
    if len(matrix) != 4 or any(len(row) != 4 for row in matrix):
        raise ValueError("determinant4_explicit requires a 4-by-4 matrix")
    total = 0
    for permutation in itertools.permutations(range(4)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(4)
            for j in range(i + 1, 4)
        )
        product = math.prod(
            matrix[row][permutation[row]] for row in range(4)
        )
        total += (-1 if inversions % 2 else 1) * product
    return total


def quadruple_index_explicit(rows: tuple[dict, dict, dict, dict]) -> int:
    moduli = [int(row["h"]) for row in rows]
    generators = (
        (moduli[0], 0, 0, 0),
        (0, moduli[1], 0, 0),
        (0, 0, moduli[2], 0),
        (0, 0, 0, moduli[3]),
        tuple(int(row["a"]) for row in rows),
        tuple(int(row["b"]) for row in rows),
    )
    determinants = []
    for selected in itertools.combinations(range(6), 4):
        matrix = tuple(
            tuple(generators[column][row] for column in selected)
            for row in range(4)
        )
        determinants.append(abs(determinant4_explicit(matrix)))
    return math.gcd(*determinants)


def best_lower(
    outside: dict,
    anchors: list[dict],
    *,
    fourway_triples: bool = False,
) -> Fraction:
    compatible = [
        anchor for anchor in anchors if pair_index(outside, anchor) == 1
    ]
    count = len(compatible)
    if not count:
        return Fraction(0)
    outside_h = int(outside["h"])
    singles = [
        Fraction(1, outside_h * int(anchor["h"]))
        for anchor in compatible
    ]
    penalties = {
        (first, second): Fraction(
            triple_index(
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
        for first in range(count)
        for second in range(first)
    }
    answer = Fraction(0)

    def visit(
        next_index: int,
        selected: tuple[int, ...],
        value: Fraction,
    ) -> None:
        nonlocal answer
        if next_index == count:
            answer = max(answer, value)
            return
        visit(next_index + 1, selected, value)
        included = value + singles[next_index]
        included -= sum(
            (
                penalties[
                    (max(next_index, previous), min(next_index, previous))
                ]
                for previous in selected
            ),
            Fraction(0),
        )
        visit(next_index + 1, (*selected, next_index), included)

    visit(0, (), Fraction(0))

    if fourway_triples:
        for selected in itertools.combinations(range(count), 3):
            first, second, third = selected
            value = singles[first] + singles[second] + singles[third]
            value -= penalties[(max(first, second), min(first, second))]
            value -= penalties[(max(first, third), min(first, third))]
            value -= penalties[(max(second, third), min(second, third))]
            rows = (
                outside,
                compatible[first],
                compatible[second],
                compatible[third],
            )
            fixed_term = Fraction(
                1,
                math.prod(int(row["h"]) for row in rows),
            )
            if value + fixed_term <= answer:
                continue
            if quadruple_index_explicit(rows) == 1:
                answer = value + fixed_term
    return answer


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("certificate", type=Path)
    parser.add_argument(
        "--block-verification",
        action="append",
        type=Path,
        default=[],
    )
    parser.add_argument("--block-catalog", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    certificate = json.loads(args.certificate.read_text())
    fourway_triples = bool(certificate.get("fourway_triples", False))
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    verification_paths = list(args.block_verification)
    if args.block_catalog:
        catalog = json.loads(args.block_catalog.read_text())
        verification_paths.extend(
            Path(entry["verification"])
            for entry in catalog["blocks"]
            if isinstance(entry, dict) and entry.get("verification")
        )
    if not verification_paths:
        raise SystemExit("provide --block-verification or --block-catalog")
    verification_by_certificate = {}
    for path in verification_paths:
        payload = json.loads(path.read_text())
        if payload.get("verified") and payload.get("certificate"):
            verification_by_certificate[payload["certificate"]] = payload

    blocks = {}
    block_errors = []
    for block_path in certificate["block_certificates"]:
        payload = json.loads(Path(block_path).read_text())
        anchors = recorded_anchor_rows(payload)
        anchor_primes = frozenset(int(row["p"]) for row in anchors)
        anchors_match = all(
            prime in by_prime
            and all(
                int(by_prime[prime][key]) == int(recorded[key])
                for key in ("h", "a", "b")
            )
            for prime, recorded in (
                (int(row["p"]), row) for row in anchors
            )
        )
        verification = verification_by_certificate.get(block_path)
        loss = read_fraction(payload["forced_overlap_loss"])
        verification_valid = (
            verification is not None
            and bool(verification.get("verified"))
            and read_fraction(verification["forced_overlap_loss"]) == loss
        )
        if not anchors_match or not verification_valid:
            block_errors.append(block_path)
            continue
        block_rows = [by_prime[prime] for prime in anchor_primes]
        blocks[block_path] = {
            "anchor_primes": anchor_primes,
            "loss": loss,
            "row_lowers": {
                int(row["p"]): (
                    Fraction(0)
                    if int(row["p"]) in anchor_primes
                    else best_lower(
                        row,
                        block_rows,
                        fourway_triples=fourway_triples,
                    )
                )
                for row in rows
            },
        }

    invalid_records = []
    proved = 0
    unresolved_periods = []
    for record in certificate["results"]:
        period = int(record["period"])
        selected = [
            row for row in rows if period % int(row["h"]) == 0
        ]
        selected_primes = {int(row["p"]) for row in selected}
        total_density = sum(
            (Fraction(1, int(row["h"])) for row in selected),
            Fraction(0),
        )
        block_path = record["block_certificate"]
        if block_path is None:
            block_loss = Fraction(0)
            star_loss = Fraction(0)
            block_applicable = True
        elif block_path in blocks:
            block = blocks[block_path]
            block_applicable = block["anchor_primes"] <= selected_primes
            block_loss = block["loss"]
            star_loss = sum(
                (
                    block["row_lowers"][int(row["p"])]
                    for row in selected
                    if int(row["p"]) not in block["anchor_primes"]
                ),
                Fraction(0),
            )
        else:
            block_applicable = False
            block_loss = Fraction(0)
            star_loss = Fraction(0)
        upper = total_density - block_loss - star_loss
        no_cover = upper < 1
        valid = (
            int(record["rows"]) == len(selected)
            and read_fraction(record["total_density"]) == total_density
            and block_applicable
            and read_fraction(record["block_overlap_loss"]) == block_loss
            and read_fraction(record["pair_anchor_star_loss"]) == star_loss
            and read_fraction(record["union_upper_bound"]) == upper
            and bool(record["proved_no_cover"]) == no_cover
        )
        if not valid:
            invalid_records.append(period)
        proved += int(no_cover)
        if not no_cover:
            unresolved_periods.append(period)

    verified = (
        str(args.pool) == certificate["pool"]
        and not block_errors
        and not invalid_records
        and int(certificate["families_checked"])
        == len(certificate["results"])
        and int(certificate["proved_no_cover"]) == proved
        and int(certificate["unresolved_count"]) == len(unresolved_periods)
        and {
            int(record["period"]) for record in certificate["unresolved"]
        }
        == set(unresolved_periods)
    )
    output = {
        "source": str(args.pool),
        "certificate": str(args.certificate),
        "fourway_triples": fourway_triples,
        "blocks_checked": len(certificate["block_certificates"]),
        "block_errors": block_errors,
        "families_checked": len(certificate["results"]),
        "proved_no_cover": proved,
        "unresolved_count": len(unresolved_periods),
        "unresolved_periods": sorted(unresolved_periods),
        "invalid_records": invalid_records,
        "verified": verified,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"blocks={len(certificate['block_certificates'])} "
        f"families={len(certificate['results'])} proved={proved} "
        f"unresolved={len(unresolved_periods)} verified={verified}",
        flush=True,
    )
    return 0 if verified else 1


if __name__ == "__main__":
    raise SystemExit(main())
