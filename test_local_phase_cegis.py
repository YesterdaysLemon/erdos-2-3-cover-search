#!/usr/bin/env python3
"""Regression tests for structural exact-witness lessons."""

from __future__ import annotations

import random
import unittest

import local_phase_cegis

try:
    import numpy as np
except ImportError:  # pragma: no cover - dependency-gated test environment
    np = None


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

    def test_multiple_tile_groups_are_unioned_without_duplicates(self) -> None:
        points = local_phase_cegis.expand_component_digit_tile_groups(
            [(0, 0)],
            [420],
            (
                ((2, 0), (3, 0), (5, 0)),
                ((7, 0),),
            ),
        )
        # The 900-point 2*3*5 tile and 49-point 7 tile share only their
        # common base point.
        self.assertEqual(len(points), 900 + 49 - 1)
        self.assertEqual(len(points), len(set(points)))

    def test_missing_component_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            local_phase_cegis.expand_top_component_tiles(
                [(0, 0)],
                [8],
                [3],
            )


@unittest.skipIf(np is None, "NumPy is unavailable")
class StreamingRepairTests(unittest.TestCase):
    candidates = [
        (2, 5, 1, 0, 2, 1),
        (2, 7, 1, 0, 2, 1),
    ]

    def test_streamed_columns_match_dense_target_matrix(self) -> None:
        points = [(0, 0), (1, 0), (2, 3), (7, 8)]
        dense, _seconds = local_phase_cegis.build_targets(
            points,
            self.candidates,
            np,
        )
        ks, ls = local_phase_cegis.streaming_coordinate_arrays(
            points,
            self.candidates,
            np,
        )
        for column, candidate in enumerate(self.candidates):
            streamed = local_phase_cegis.streaming_target_values(
                ks,
                ls,
                candidate,
                np,
            )
            np.testing.assert_array_equal(streamed, dense[:, column])

    def test_component_residue_state_matches_big_integer_targets(self) -> None:
        points = [
            ((1 << 100) + 12345, (1 << 96) + 6789),
            ((1 << 120) + 77, (1 << 111) + 91),
        ]
        candidates = [
            (12, 13, 5, 7, 12, 12),
            (35, 71, 11, 19, 35, 35),
        ]
        state = local_phase_cegis.streaming_coordinate_state(
            points,
            candidates,
            np,
        )
        self.assertEqual(state["mode"], "components")
        for candidate in candidates:
            h, _p, a, b, _ord2, _ord3 = candidate
            expected = np.asarray(
                [
                    (a * k + b * l) % h
                    for k, l in points
                ],
                dtype=np.uint32,
            )
            actual = local_phase_cegis.streaming_state_target_values(
                state,
                candidate,
                np,
            )
            np.testing.assert_array_equal(actual, expected)

    def test_streaming_repair_can_split_two_identical_rows(self) -> None:
        for coordinate_descent in (False, True):
            with self.subTest(coordinate_descent=coordinate_descent):
                assignment = np.asarray([0, 0], dtype=np.uint32)
                solved, steps, before, after = (
                    local_phase_cegis.repair_streaming(
                        [(0, 0), (1, 0)],
                        self.candidates,
                        assignment,
                        [0, 1],
                        random.Random(203),
                        np,
                        max_steps=10,
                        sample_size=2,
                        valid_moduli=np.asarray([1, 1], dtype=np.uint32),
                        valid_residues=np.asarray([0, 0], dtype=np.uint32),
                        candidate_moduli=np.asarray(
                            [2, 2],
                            dtype=np.uint32,
                        ),
                        coordinate_descent=coordinate_descent,
                    )
                )
                self.assertTrue(solved)
                self.assertEqual(steps, 1)
                self.assertEqual((before, after), (1, 0))
                self.assertEqual(set(map(int, assignment)), {0, 1})

    def test_streaming_cover_state_can_be_reused(self) -> None:
        points = [(0, 0), (1, 0)]
        assignment = np.asarray([0, 1], dtype=np.uint32)
        state = local_phase_cegis.streaming_coordinate_state(
            points,
            self.candidates,
            np,
        )
        cover = local_phase_cegis.streaming_cover_counts(
            state,
            self.candidates,
            assignment,
            np,
        )
        result = local_phase_cegis.repair_streaming(
            points,
            self.candidates,
            assignment,
            [0, 1],
            random.Random(203),
            np,
            max_steps=10,
            sample_size=2,
            valid_moduli=np.asarray([1, 1], dtype=np.uint32),
            valid_residues=np.asarray([0, 0], dtype=np.uint32),
            candidate_moduli=np.asarray([2, 2], dtype=np.uint32),
            coordinate_state=state,
            initial_cover=cover,
            return_state=True,
        )
        solved, steps, before, after, returned_state, returned_cover = result
        self.assertTrue(solved)
        self.assertEqual((steps, before, after), (0, 0, 0))
        self.assertIs(returned_state, state)
        np.testing.assert_array_equal(returned_cover, cover)


if __name__ == "__main__":
    unittest.main()
