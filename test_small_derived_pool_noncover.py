#!/usr/bin/env python3
"""Tests for the independent small-pool phase enumerator."""

from __future__ import annotations

import unittest

import verify_small_derived_pool_noncover


class SmallDerivedPoolNoncoverTests(unittest.TestCase):
    def test_multiplicative_order_and_power_congruence(self):
        self.assertEqual(
            verify_small_derived_pool_noncover.multiplicative_order(2, 41),
            20,
        )
        self.assertEqual(
            verify_small_derived_pool_noncover.multiplicative_order(3, 41),
            8,
        )
        self.assertEqual(
            verify_small_derived_pool_noncover.power_target_congruence(
                40,
                41,
                1616615,
            ),
            (0, 5),
        )

    def test_one_binary_row_cannot_cover_both_points(self):
        rows = [{"h": 2, "p": 5, "a": 1, "b": 0}]
        result = verify_small_derived_pool_noncover.enumerate_phase_cover(
            rows,
            [(0, 0), (1, 0)],
        )
        cover_exists, enumerated, assignment_count, phases = result
        self.assertFalse(cover_exists)
        self.assertEqual(enumerated, 2)
        self.assertEqual(assignment_count, 2)
        self.assertIsNone(phases)

    def test_two_binary_rows_can_split_the_points(self):
        rows = [
            {"h": 2, "p": 5, "a": 1, "b": 0},
            {"h": 2, "p": 7, "a": 1, "b": 0},
        ]
        cover_exists, _enumerated, assignment_count, phases = (
            verify_small_derived_pool_noncover.enumerate_phase_cover(
                rows,
                [(0, 0), (1, 0)],
            )
        )
        self.assertTrue(cover_exists)
        self.assertEqual(assignment_count, 4)
        self.assertEqual(set(phases), {0, 1})


if __name__ == "__main__":
    unittest.main()
