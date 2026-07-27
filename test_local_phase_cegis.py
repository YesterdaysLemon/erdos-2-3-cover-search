#!/usr/bin/env python3
"""Regression tests for structural exact-witness lessons."""

from __future__ import annotations

import unittest

import local_phase_cegis


class TopComponentTileTests(unittest.TestCase):
    def test_two_independent_top_digits_form_complete_tile(self) -> None:
        # Period 12 = 2^2 * 3.  The top binary and ternary digits give
        # 2*3 coordinate values independently in k and l.
        points = local_phase_cegis.expand_top_component_tiles(
            [(1, 2)],
            [4, 3, 12],
            [2, 3],
        )
        self.assertEqual(len(points), 36)
        self.assertEqual(len(set(points)), 36)
        self.assertTrue(
            all(k % 2 == 1 and 0 <= k < 12 for k, _l in points)
        )
        self.assertTrue(all(0 <= l < 12 for _k, l in points))

    def test_duplicate_input_points_are_deduplicated(self) -> None:
        points = local_phase_cegis.expand_top_component_tiles(
            [(0, 0), (0, 0)],
            [8],
            [2],
        )
        self.assertEqual(len(points), 4)

    def test_low_digit_changes_every_divisible_row_component(self) -> None:
        # In period 12, changing the lowest binary digit uses step 3.  It
        # preserves the ternary coordinate but changes parity.
        points = local_phase_cegis.expand_component_digit_tiles(
            [(0, 0)],
            [12],
            [(2, 0)],
        )
        self.assertEqual(set(points), {(0, 0), (0, 3), (3, 0), (3, 3)})

    def test_missing_component_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            local_phase_cegis.expand_top_component_tiles(
                [(0, 0)],
                [8],
                [3],
            )


if __name__ == "__main__":
    unittest.main()
