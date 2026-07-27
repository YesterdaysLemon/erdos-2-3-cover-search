#!/usr/bin/env python3

from __future__ import annotations

import unittest

import anchor_phase_master


class AnchorPhaseMasterTests(unittest.TestCase):
    def test_master_finds_coordinated_anchor_targets(self):
        rows = [
            {"p": 5, "h": 2, "a": 1, "b": 0},
            {"p": 7, "h": 2, "a": 0, "b": 1},
        ]
        phases, metadata = anchor_phase_master.solve_anchor_master(
            rows,
            {5: 0, 7: 0},
            (5, 7),
            [(0, 0), (1, 1)],
        )
        self.assertTrue(metadata["sat"])
        self.assertIsNotNone(phases)
        for k, l in [(0, 0), (1, 1)]:
            self.assertTrue(
                any(
                    (row["a"] * k + row["b"] * l - phases[row["p"]])
                    % row["h"]
                    == 0
                    for row in rows
                )
            )

    def test_master_proves_one_binary_row_cannot_cover_both_classes(self):
        rows = [{"p": 5, "h": 2, "a": 1, "b": 0}]
        phases, metadata = anchor_phase_master.solve_anchor_master(
            rows,
            {5: 0},
            (5,),
            [(0, 0), (1, 0)],
        )
        self.assertFalse(metadata["sat"])
        self.assertIsNone(phases)

    def test_frozen_and_algebraic_points_are_discharged(self):
        rows = [
            {"p": 5, "h": 2, "a": 1, "b": 0},
            {"p": 7, "h": 2, "a": 0, "b": 1},
        ]
        _phases, metadata = anchor_phase_master.solve_anchor_master(
            rows,
            {5: 0, 7: 0},
            (5,),
            [(1, 0), (3, 3)],
            (3,),
        )
        self.assertEqual(metadata["frozen_discharged_points"], 1)
        self.assertEqual(metadata["algebraically_discharged_points"], 1)
        self.assertEqual(metadata["eligible_points"], 0)

    def test_minimum_change_master_preserves_preferred_targets(self):
        rows = [
            {"p": 5, "h": 2, "a": 1, "b": 0},
            {"p": 7, "h": 2, "a": 0, "b": 1},
        ]
        phases, metadata = anchor_phase_master.solve_anchor_master(
            rows,
            {5: 0, 7: 0},
            (5, 7),
            [(0, 0)],
            preferred_anchor_phases={5: 0, 7: 0},
        )
        self.assertEqual(phases[5], 0)
        self.assertEqual(phases[7], 0)
        self.assertEqual(metadata["minimum_anchor_changes"], 0)

    def test_minimum_change_master_makes_only_a_required_change(self):
        rows = [
            {"p": 5, "h": 2, "a": 1, "b": 0},
            {"p": 7, "h": 2, "a": 0, "b": 1},
        ]
        phases, metadata = anchor_phase_master.solve_anchor_master(
            rows,
            {5: 0, 7: 0},
            (5, 7),
            [(1, 1)],
            preferred_anchor_phases={5: 0, 7: 0},
        )
        self.assertEqual(metadata["minimum_anchor_changes"], 1)
        self.assertEqual(sum(phases[p] != 0 for p in (5, 7)), 1)


if __name__ == "__main__":
    unittest.main()
