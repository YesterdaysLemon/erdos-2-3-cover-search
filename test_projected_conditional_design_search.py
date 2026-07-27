#!/usr/bin/env python3
"""Tests for projected conditional-design discovery."""

from __future__ import annotations

import math
import unittest
from fractions import Fraction

from search_projected_conditional_designs import (
    Design,
    Projected,
    evaluate_design,
    projectable,
)


class ProjectedConditionalDesignTests(unittest.TestCase):
    def test_projectable_uses_exact_coprime_residual_split(self) -> None:
        row = {"p": 13, "h": 12, "a": 1, "b": 4}
        self.assertEqual(projectable(row, 20), Projected(13, 4, 3))

        # This row's residual equation is not primitive modulo three.
        nonprimitive = {"p": 13, "h": 12, "a": 3, "b": 6}
        self.assertIsNone(projectable(nonprimitive, 4))

    def test_small_design_matches_complete_phase_enumeration(self) -> None:
        outside = {"p": 2, "h": 2, "a": 1, "b": 0}
        rows = {
            3: {"p": 3, "h": 3, "a": 0, "b": 1},
            5: {"p": 5, "h": 5, "a": 1, "b": 1},
            7: {"p": 7, "h": 7, "a": 2, "b": 1},
        }
        design = Design(
            normalizer=3,
            base_primes=(5,),
            projected=(Projected(7, 1, 7),),
            base_period=30,
            base_points=450,
            tensor_cells=5,
            density_score=Fraction(1, 3)
            + Fraction(1, 5)
            + Fraction(1, 7),
        )
        evaluated = evaluate_design(outside, design, rows)
        self.assertIsNotNone(evaluated)
        actual, _detail = evaluated

        period = math.lcm(2, 3, 5, 7)
        expected_count = period * period
        for target5 in range(5):
            for target7 in range(7):
                covered = 0
                for first in range(period):
                    if first % 2:
                        continue
                    for second in range(period):
                        if (
                            second % 3 == 0
                            or (first + second - target5) % 5 == 0
                            or (
                                2 * first + second - target7
                            ) % 7 == 0
                        ):
                            covered += 1
                expected_count = min(expected_count, covered)
        self.assertEqual(
            actual,
            Fraction(expected_count, period * period),
        )


if __name__ == "__main__":
    unittest.main()
