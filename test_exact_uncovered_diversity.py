#!/usr/bin/env python3
"""Regression tests for adversarial exact-witness diversity."""

from __future__ import annotations

import unittest

import exact_derived_phase_misses
import exact_uncovered


class PerTargetDiversityTests(unittest.TestCase):
    def test_phase_overrides_are_explicit_and_unique(self):
        self.assertEqual(
            exact_derived_phase_misses.parse_phase_overrides(
                ["41:1", "73:0"],
                {41, 73},
            ),
            {41: 1, 73: 0},
        )
        with self.assertRaises(ValueError):
            exact_derived_phase_misses.parse_phase_overrides(
                ["41:0", "41:1"],
                {41},
            )
        with self.assertRaises(ValueError):
            exact_derived_phase_misses.parse_phase_overrides(
                ["19:0"],
                {41},
            )

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

    def test_cap_ignores_targets_forbidden_by_row_restriction(self) -> None:
        rows = [
            {
                "h": 2,
                "p": 5,
                "a": 1,
                "b": 0,
                "c": 0,
                "target_modulus": 2,
                "target_residue": 0,
            }
        ]
        witnesses, _meta = exact_uncovered.find_uncovered(
            rows,
            max_component=2,
            limit=10,
            diversity_primes=(5,),
            diversity_target_cap=1,
        )
        # Both k=1 points have synthesized row target 1, but target 1 is
        # illegal for this restricted row and therefore cannot be used to
        # absorb either witness by changing its phase.
        self.assertEqual(len(witnesses), 2)
        self.assertEqual({k for k, _l in witnesses}, {1})
        self.assertEqual({l for _k, l in witnesses}, {0, 1})

    def test_coprime_selector_adds_every_full_weight_row(self) -> None:
        rows = [
            {"h": 6, "p": 5},
            {"h": 10, "p": 7},
            {"h": 7, "p": 11},
        ]
        selected = exact_derived_phase_misses.expanded_diversity_primes(
            rows,
            (7,),
            5,
        )
        self.assertEqual(selected, (7, 5, 11))

    def test_coprime_selector_rejects_one(self) -> None:
        with self.assertRaisesRegex(ValueError, "at least two"):
            exact_derived_phase_misses.expanded_diversity_primes(
                [],
                (),
                1,
            )


if __name__ == "__main__":
    unittest.main()
