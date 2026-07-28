#!/usr/bin/env python3
"""Inventory arithmetically realizable line slots on a finite quotient.

An abstract affine cover of a small finite abelian group is useful only when
its line slots can be assigned to distinct primes from the arithmetic pool.
This preflight transforms each source signature into a declared unimodular
basis, keeps exactly the rows whose affine predicate descends to
Z/U x Z/V, and groups geometrically equivalent directions with their actual
prime capacities.

The output is an inventory, not a cover and not evidence that a cover exists.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

from build_axis_layered_pool import sha256


def prime_power_components(value: int) -> list[int]:
    components = []
    divisor = 2
    while divisor * divisor <= value:
        component = 1
        while value % divisor == 0:
            value //= divisor
            component *= divisor
        if component > 1:
            components.append(component)
        divisor += 1 if divisor == 2 else 2
    if value > 1:
        components.append(value)
    return components


def axis_rectangle_moduli(
    h: int,
    a: int,
    b: int,
) -> tuple[int, int] | None:
    """Return the independent k/l moduli when a fibre is a CRT rectangle."""

    k_modulus = 1
    l_modulus = 1
    for component in prime_power_components(h):
        if b % component == 0 and math.gcd(a, component) == 1:
            k_modulus *= component
        elif a % component == 0 and math.gcd(b, component) == 1:
            l_modulus *= component
        else:
            return None
    if k_modulus * l_modulus != h:
        raise AssertionError("rectangle components do not reconstruct h")
    return k_modulus, l_modulus


def transformed_coefficients(
    row: dict,
    direction: tuple[int, int],
    transverse: tuple[int, int],
) -> tuple[int, int, int]:
    h = int(row["h"])
    if h < 1:
        raise ValueError("row modulus must be positive")
    source_a = int(row["a"]) % h
    source_b = int(row["b"]) % h
    if math.gcd(source_a, source_b, h) != 1:
        raise ValueError(f"non-surjective row p={row['p']}")
    a = (
        source_a * direction[0] + source_b * direction[1]
    ) % h
    b = (
        source_a * transverse[0] + source_b * transverse[1]
    ) % h
    if math.gcd(a, b, h) != 1:
        raise AssertionError("unimodular transform lost surjectivity")
    return h, a, b


def canonical_projective_direction(
    h: int,
    a: int,
    b: int,
) -> tuple[int, int]:
    """Canonicalize a primitive pair under multiplication by units mod h."""

    return min(
        (
            multiplier * a % h,
            multiplier * b % h,
        )
        for multiplier in range(1, h + 1)
        if math.gcd(multiplier, h) == 1
    )


def transformed_target_restriction(
    row: dict,
    h: int,
    a: int,
    b: int,
    canonical: tuple[int, int],
) -> tuple[int, int]:
    """Transport a target congruence through direction normalization."""

    modulus = int(row.get("target_modulus", 1))
    if modulus < 1 or h % modulus:
        raise ValueError(f"invalid target restriction p={row['p']}")
    residue = int(row.get("target_residue", 0)) % modulus
    multiplier = next(
        value
        for value in range(1, h + 1)
        if math.gcd(value, h) == 1
        and (value * a % h, value * b % h) == canonical
    )
    return modulus, multiplier * residue % modulus


def pair_map_cokernel_index(
    left: dict,
    right: dict,
) -> int:
    """Return the index of the joint affine image in its codomain."""

    left_h = int(left["h"])
    right_h = int(right["h"])
    left_a = int(left["a"]) % left_h
    left_b = int(left["b"]) % left_h
    right_a = int(right["a"]) % right_h
    right_b = int(right["b"]) % right_h
    return math.gcd(
        left_a * right_b - left_b * right_a,
        left_h * right_a,
        left_h * right_b,
        right_h * left_a,
        right_h * left_b,
        left_h * right_h,
    )


def maximum_forced_overlap_forest(
    line_types: list[dict],
    cell_count: int,
) -> tuple[list[dict], Fraction]:
    """Return a maximum-weight Hunter forest of forced intersections."""

    slots = sorted(
        (
            {
                "p": int(prime),
                "h": int(record["h"]),
                "a": int(record["a"]) % int(record["h"]),
                "b": int(record["b"]) % int(record["h"]),
            }
            for record in line_types
            for prime in record["primes"]
        ),
        key=lambda slot: slot["p"],
    )
    remaining = set(range(len(slots)))
    selected_pairs: list[tuple[dict, dict, Fraction]] = []
    total = Fraction()
    while remaining:
        root = min(remaining, key=lambda index: slots[index]["p"])
        remaining.remove(root)
        best: dict[int, tuple[Fraction, int]] = {}
        for index in remaining:
            if pair_map_cokernel_index(slots[root], slots[index]) == 1:
                best[index] = (
                    Fraction(1, slots[root]["h"] * slots[index]["h"]),
                    root,
                )
        while best:
            child = min(
                best,
                key=lambda index: (
                    -best[index][0],
                    min(
                        slots[best[index][1]]["p"],
                        slots[index]["p"],
                    ),
                    max(
                        slots[best[index][1]]["p"],
                        slots[index]["p"],
                    ),
                ),
            )
            weight, parent = best.pop(child)
            if child not in remaining:
                continue
            remaining.remove(child)
            left, right = sorted(
                (slots[parent], slots[child]),
                key=lambda slot: slot["p"],
            )
            selected_pairs.append((left, right, weight))
            total += weight
            for index in remaining:
                if (
                    pair_map_cokernel_index(
                        slots[child],
                        slots[index],
                    )
                    != 1
                ):
                    continue
                candidate = Fraction(
                    1,
                    slots[child]["h"] * slots[index]["h"],
                )
                previous = best.get(index)
                candidate_key = (
                    min(slots[child]["p"], slots[index]["p"]),
                    max(slots[child]["p"], slots[index]["p"]),
                )
                previous_key = (
                    (
                        min(
                            slots[previous[1]]["p"],
                            slots[index]["p"],
                        ),
                        max(
                            slots[previous[1]]["p"],
                            slots[index]["p"],
                        ),
                    )
                    if previous is not None
                    else None
                )
                if (
                    previous is None
                    or candidate > previous[0]
                    or (
                        candidate == previous[0]
                        and candidate_key < previous_key
                    )
                ):
                    best[index] = (candidate, child)

    selected = []
    for left, right, weight in selected_pairs:
        weight = Fraction(1, left["h"] * right["h"])
        if cell_count % weight.denominator:
            raise AssertionError(
                "surjective joint affine map has fractional fibre size"
            )
        selected.append(
            {
                "left_prime": left["p"],
                "left_modulus": left["h"],
                "right_prime": right["p"],
                "right_modulus": right["h"],
                "intersection_density_numerator": weight.numerator,
                "intersection_density_denominator": weight.denominator,
                "intersection_cell_count": (
                    cell_count * weight.numerator // weight.denominator
                ),
                "joint_map_cokernel_index": 1,
                "reason": (
                    "joint affine map is surjective, forcing this "
                    "intersection for every pair of phases"
                ),
            }
        )
    return selected, total


def inventory(
    payload: dict,
    source_path: Path,
    width: int,
    height: int,
    direction: tuple[int, int],
    transverse: tuple[int, int],
) -> dict:
    authenticated_payload = json.loads(source_path.read_text())
    if authenticated_payload != payload:
        raise ValueError("source payload differs from authenticated file")
    if width < 1 or height < 1:
        raise ValueError("quotient dimensions must be positive")
    determinant = (
        direction[0] * transverse[1]
        - direction[1] * transverse[0]
    )
    if abs(determinant) != 1:
        raise ValueError("basis must be unimodular")
    rows = payload["choices"]
    seen_primes = set()
    groups: dict[tuple[int, int, int, int, int], dict] = {}
    excluded_not_descending = 0
    for row in rows:
        prime = int(row["p"])
        if prime in seen_primes:
            raise ValueError(f"duplicate source prime {prime}")
        seen_primes.add(prime)
        h, a, b = transformed_coefficients(
            row,
            direction,
            transverse,
        )
        if (a * width) % h or (b * height) % h:
            excluded_not_descending += 1
            continue
        if (width * height) % h:
            raise AssertionError("descending surjection has nonintegral fibre")
        canonical = canonical_projective_direction(h, a, b)
        target_modulus, target_residue = transformed_target_restriction(
            row,
            h,
            a,
            b,
            canonical,
        )
        rectangle = axis_rectangle_moduli(h, *canonical)
        key = (
            h,
            canonical[0],
            canonical[1],
            target_modulus,
            target_residue,
        )
        record = groups.setdefault(
            key,
            {
                "h": h,
                "a": canonical[0],
                "b": canonical[1],
                "target_modulus": target_modulus,
                "target_residue": target_residue,
                "fibre_size": width * height // h,
                "legal_target_count_per_prime": h // target_modulus,
                "axis_rectangle": (
                    {
                        "k_modulus": rectangle[0],
                        "l_modulus": rectangle[1],
                        "proper_two_axis_rectangle": (
                            rectangle[0] > 1 and rectangle[1] > 1
                        ),
                    }
                    if rectangle is not None
                    else None
                ),
                "primes": [],
                "source_signatures": [],
            },
        )
        record["primes"].append(prime)
        record["source_signatures"].append(
            {
                "p": prime,
                "a": a,
                "b": b,
            }
        )

    ordered_groups = []
    for key in sorted(groups):
        record = groups[key]
        record["primes"].sort()
        record["source_signatures"].sort(key=lambda item: item["p"])
        record["prime_capacity"] = len(record["primes"])
        ordered_groups.append(record)
    descending_count = sum(
        record["prime_capacity"] for record in ordered_groups
    )
    axis_separable_count = sum(
        record["prime_capacity"]
        for record in ordered_groups
        if record["axis_rectangle"] is not None
    )
    proper_rectangle_count = sum(
        record["prime_capacity"]
        for record in ordered_groups
        if record["axis_rectangle"] is not None
        and record["axis_rectangle"]["proper_two_axis_rectangle"]
    )
    density = sum(
        (
            Fraction(record["prime_capacity"], record["h"])
            for record in ordered_groups
        ),
        Fraction(),
    )
    overlap_forest, forced_overlap = maximum_forced_overlap_forest(
        ordered_groups,
        width * height,
    )
    union_upper_bound = density - forced_overlap
    phase_independent_noncover = union_upper_bound < 1
    unrestricted_groups = [
        record
        for record in ordered_groups
        if record["target_modulus"] == 1
    ]
    unrestricted_count = sum(
        record["prime_capacity"] for record in unrestricted_groups
    )
    avoidance_probability = Fraction(1)
    for record in unrestricted_groups:
        avoidance_probability *= Fraction(
            record["h"] - 1,
            record["h"],
        ) ** int(record["prime_capacity"])
    expected_uncovered = width * height * avoidance_probability
    first_moment_cover_exists = expected_uncovered < 1
    if first_moment_cover_exists and phase_independent_noncover:
        raise AssertionError(
            "first-moment existence conflicts with Hunter noncover bound"
        )
    return {
        "schema": "reverse_group_arithmetic_line_inventory_v3",
        "source": str(source_path),
        "source_sha256": sha256(source_path),
        "basis": {
            "direction": list(direction),
            "transverse": list(transverse),
            "determinant": determinant,
        },
        "group": {
            "width": width,
            "height": height,
            "cell_count": width * height,
        },
        "source_row_count": len(rows),
        "descending_row_count": descending_count,
        "excluded_not_descending_count": excluded_not_descending,
        "geometric_type_count": len(ordered_groups),
        "axis_separable_row_count": axis_separable_count,
        "proper_two_axis_rectangle_row_count": proper_rectangle_count,
        "raw_density_numerator": density.numerator,
        "raw_density_denominator": density.denominator,
        "raw_density_decimal": float(density),
        "density_at_least_one": density >= 1,
        "forced_overlap_forest": overlap_forest,
        "forced_overlap_density_numerator": (
            forced_overlap.numerator
        ),
        "forced_overlap_density_denominator": (
            forced_overlap.denominator
        ),
        "phase_independent_union_upper_bound_numerator": (
            union_upper_bound.numerator
        ),
        "phase_independent_union_upper_bound_denominator": (
            union_upper_bound.denominator
        ),
        "phase_independent_union_upper_bound_decimal": float(
            union_upper_bound
        ),
        "finite_group_cover_impossible_by_density_overlap": (
            phase_independent_noncover
        ),
        "unrestricted_descending_row_count": unrestricted_count,
        "random_phase_avoidance_probability_numerator": (
            avoidance_probability.numerator
        ),
        "random_phase_avoidance_probability_denominator": (
            avoidance_probability.denominator
        ),
        "expected_uncovered_cell_count_numerator": (
            expected_uncovered.numerator
        ),
        "expected_uncovered_cell_count_denominator": (
            expected_uncovered.denominator
        ),
        "expected_uncovered_cell_count_decimal": float(
            expected_uncovered
        ),
        "first_moment_cover_exists": first_moment_cover_exists,
        "line_types": ordered_groups,
        "claim": {
            "inventory_exact": True,
            "finite_group_cover_exists_by_first_moment": (
                first_moment_cover_exists
            ),
            "finite_group_cover_impossible_by_density_overlap": (
                phase_independent_noncover
            ),
            "explicit_finite_group_cover_found": False,
            "distinct_prime_matching_found": False,
            "integer_m_found": False,
        },
        "scope": (
            "exact inventory of source-prime line slots whose predicates "
            "descend to the declared finite group, including exact CRT "
            "rectangle signatures and a maximum-weight, phase-independent "
            "Hunter forest; any noncover claim applies only to this "
            "declared finite quotient family"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--width", type=int, required=True)
    parser.add_argument("--height", type=int, required=True)
    parser.add_argument(
        "--direction",
        nargs=2,
        type=int,
        default=(1, 0),
    )
    parser.add_argument(
        "--transverse",
        nargs=2,
        type=int,
        default=(0, 1),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.pool.read_text())
    result = inventory(
        payload,
        args.pool,
        args.width,
        args.height,
        tuple(args.direction),
        tuple(args.transverse),
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"descending={result['descending_row_count']} "
        f"types={result['geometric_type_count']} "
        f"rectangles={result['proper_two_axis_rectangle_row_count']} "
        f"density={result['raw_density_numerator']}/"
        f"{result['raw_density_denominator']} "
        f"at_least_one={result['density_at_least_one']} "
        f"first_moment={result['first_moment_cover_exists']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
