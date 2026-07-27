#!/usr/bin/env python3
"""Regression tests for adversarial exact-witness diversity."""

from __future__ import annotations

import unittest

import exact_uncovered


class PerTargetDiversityTests(unittest.TestCase):
    rows = [
        {"h": 2, "p": 5, "a": 1, "b": 0, "c": 0},
        {"h": 3, "p": 7, "a": 0, "b": 1, "c": 0},
    ]

    def test_cap_limits_reuse_of_each_selected_row_target(self) -> None:
        witnesses, meta = exact_uncovered.find_uncovered(
            self.rows,
            max_component=3,
            limit=10,
            diversity_primes=(5, 7),
            diversity_target_cap=1,
        )
        # Every hole has k=1 mod 2, so the p=5 target value is identical.
        # Once that value supplies one witness, the cap exhausts this
        # deliberately diversified checker batch.
        self.assertEqual(len(witnesses), 1)
        self.assertEqual(meta["diversity_target_cap"], 1)

    def test_positive_cap_requires_selected_rows(self) -> None:
        with self.assertRaises(ValueError):
            exact_uncovered.find_uncovered(
                self.rows,
                max_component=3,
                limit=1,
                diversity_target_cap=1,
            )


if __name__ == "__main__":
    unittest.main()
