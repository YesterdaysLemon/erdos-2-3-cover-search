#!/usr/bin/env python3
"""Regression tests for gain-mask finite repair search."""

from __future__ import annotations

import unittest

import finite_sample_mask_repair

try:
    import pysat  # noqa: F401
except ImportError:  # pragma: no cover - dependency-gated environment
    pysat = None


@unittest.skipIf(pysat is None, "PySAT is unavailable")
class FiniteSampleMaskRepairTests(unittest.TestCase):
    rows = [
        {
            "h": 2,
            "p": 5,
            "a": 1,
            "b": 0,
            "ord2": 2,
            "ord3": 1,
            "target_modulus": 1,
            "target_residue": 0,
        },
        {
            "h": 2,
            "p": 7,
            "a": 1,
            "b": 0,
            "ord2": 2,
            "ord3": 1,
            "target_modulus": 1,
            "target_residue": 0,
        },
    ]
    candidates = [
        (
            int(row["h"]),
            int(row["p"]),
            int(row["a"]),
            int(row["b"]),
            int(row["ord2"]),
            int(row["ord3"]),
        )
        for row in rows
    ]
    points = [(0, 0), (1, 0)]

    def test_one_change_repair_replays_every_point(self) -> None:
        result = finite_sample_mask_repair.search_mask_repair(
            self.rows,
            self.candidates,
            self.points,
            {5: 0, 7: 0},
            set(),
            1,
            "cadical195",
            100,
            10.0,
        )
        self.assertEqual(result["status"], "INTEGER_MODEL")
        self.assertEqual(result["changed_phases"], 1)
        self.assertEqual(result["full_misses"], 0)
        self.assertEqual(set(result["repaired"]), {0, 1})

    def test_negative_result_is_explicitly_incomplete(self) -> None:
        result = finite_sample_mask_repair.search_mask_repair(
            self.rows,
            self.candidates,
            self.points,
            {5: 0, 7: 0},
            {5, 7},
            1,
            "cadical195",
            100,
            10.0,
        )
        self.assertEqual(result["status"], "NO_GAIN_MASK_MODEL")
        self.assertFalse(result["complete_negative"])


if __name__ == "__main__":
    unittest.main()
