#!/usr/bin/env python3
"""Regression tests for the surjective four-fold star correction."""

from __future__ import annotations

import itertools
import math
import random
import unittest
from fractions import Fraction

from scan_all_ranked_pairanchor_star import (
    best_anchor_lower,
    quintuple_image_index_minors,
    quadruple_image_index_minors,
)
from certify_ranked_period_conditional_star import select_star_anchors
from verify_all_ranked_pairanchor_star import (
    best_lower,
    quintuple_index_explicit,
    quadruple_index_explicit,
)
from verify_ranked_period_conditional_star import (
    select_star_anchors_independent,
)


def brute_image_index(rows: tuple[dict, ...]) -> int:
    moduli = tuple(int(row["h"]) for row in rows)
    period = math.lcm(*moduli)
    image = {
        tuple(
            (
                int(row["a"]) * first
                + int(row["b"]) * second
            )
            % modulus
            for row, modulus in zip(rows, moduli)
        )
        for first in range(period)
        for second in range(period)
    }
    return math.prod(moduli) // len(image)


def brute_outside_anchor_union(
    outside: dict,
    anchors: list[dict],
) -> Fraction:
    rows = [outside, *anchors]
    period = math.lcm(*(int(row["h"]) for row in rows))

    def holds(row: dict, first: int, second: int) -> bool:
        return (
            int(row["a"]) * first
            + int(row["b"]) * second
            - int(row["c"])
        ) % int(row["h"]) == 0

    covered = sum(
        holds(outside, first, second)
        and any(holds(anchor, first, second) for anchor in anchors)
        for first in range(period)
        for second in range(period)
    )
    return Fraction(covered, period * period)


class FourwayIndexTests(unittest.TestCase):
    def test_limited_star_subset_is_compatible_and_smallest_first(
        self,
    ) -> None:
        outside = {"p": 2, "h": 2, "a": 1, "b": 0}
        anchors = [
            {"p": 11, "h": 11, "a": 1, "b": 0},
            {"p": 4, "h": 4, "a": 1, "b": 0},
            {"p": 3, "h": 3, "a": 1, "b": 0},
            {"p": 7, "h": 7, "a": 1, "b": 0},
            {"p": 5, "h": 5, "a": 1, "b": 0},
        ]
        expected = [anchors[2], anchors[4], anchors[3]]
        self.assertEqual(
            select_star_anchors(outside, anchors, 3),
            expected,
        )
        self.assertEqual(
            select_star_anchors_independent(outside, anchors, 3),
            expected,
        )

    def test_random_small_maps_match_brute_image(self) -> None:
        rng = random.Random(203)
        for _ in range(40):
            moduli = tuple(rng.randint(2, 6) for _ in range(4))
            rows = tuple(
                {
                    "h": modulus,
                    "a": rng.randrange(modulus),
                    "b": rng.randrange(modulus),
                }
                for modulus in moduli
            )
            expected = brute_image_index(rows)
            self.assertEqual(quadruple_image_index_minors(rows), expected)
            self.assertEqual(quadruple_index_explicit(rows), expected)

    def test_random_five_factor_maps_match_brute_image(self) -> None:
        rng = random.Random(211)
        for _ in range(30):
            moduli = tuple(rng.randint(2, 5) for _ in range(5))
            rows = tuple(
                {
                    "h": modulus,
                    "a": rng.randrange(modulus),
                    "b": rng.randrange(modulus),
                }
                for modulus in moduli
            )
            expected = brute_image_index(rows)
            self.assertEqual(quintuple_image_index_minors(rows), expected)
            self.assertEqual(quintuple_index_explicit(rows), expected)

    def test_surjective_term_strictly_improves_three_anchor_bound(self) -> None:
        outside = {"h": 2, "a": 1, "b": 0}
        anchors = [
            {"h": modulus, "a": 1, "b": 0}
            for modulus in (3, 5, 7)
        ]
        old = best_anchor_lower(outside, anchors)
        strengthened = best_anchor_lower(
            outside,
            anchors,
            fourway_triples=True,
        )
        self.assertEqual(old, Fraction(56, 210))
        self.assertEqual(strengthened, Fraction(57, 210))
        self.assertEqual(
            strengthened,
            best_lower(outside, anchors, fourway_triples=True),
        )

    def test_fourth_order_bound_is_exact_for_independent_anchors(self) -> None:
        outside = {"h": 2, "a": 1, "b": 0}
        anchors = [
            {"h": modulus, "a": 1, "b": 0}
            for modulus in (3, 5, 7, 11)
        ]
        old = best_anchor_lower(
            outside,
            anchors,
            fourway_triples=True,
        )
        strengthened = best_anchor_lower(
            outside,
            anchors,
            fourway_triples=True,
            fourterm_quads=True,
        )
        expected = Fraction(1, 2) * (
            1
            - math.prod(
                Fraction(modulus - 1, modulus)
                for modulus in (3, 5, 7, 11)
            )
        )
        self.assertGreater(strengthened, old)
        self.assertEqual(strengthened, expected)
        self.assertEqual(
            strengthened,
            best_lower(
                outside,
                anchors,
                fourway_triples=True,
                fourterm_quads=True,
            ),
        )

    def test_random_fourth_order_bounds_are_below_exact_unions(self) -> None:
        rng = random.Random(223)
        for _ in range(50):
            rows = [
                {
                    "h": modulus,
                    "a": rng.randrange(modulus),
                    "b": rng.randrange(modulus),
                    "c": rng.randrange(modulus),
                }
                for modulus in (
                    rng.randint(2, 5),
                    rng.randint(2, 5),
                    rng.randint(2, 5),
                    rng.randint(2, 5),
                    rng.randint(2, 5),
                )
            ]
            outside, anchors = rows[0], rows[1:]
            generated = best_anchor_lower(
                outside,
                anchors,
                fourway_triples=True,
                fourterm_quads=True,
            )
            independently_replayed = best_lower(
                outside,
                anchors,
                fourway_triples=True,
                fourterm_quads=True,
            )
            self.assertEqual(generated, independently_replayed)
            self.assertLessEqual(
                generated,
                brute_outside_anchor_union(outside, anchors),
            )


if __name__ == "__main__":
    unittest.main()
