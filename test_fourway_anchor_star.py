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
    quadruple_image_index_minors,
)
from verify_all_ranked_pairanchor_star import (
    best_lower,
    quadruple_index_explicit,
)


def brute_image_index(rows: tuple[dict, dict, dict, dict]) -> int:
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


class FourwayIndexTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
