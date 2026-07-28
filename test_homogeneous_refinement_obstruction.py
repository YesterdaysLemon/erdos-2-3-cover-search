#!/usr/bin/env python3
"""Regression tests for the homogeneous-refinement sibling obstruction."""

from __future__ import annotations

import unittest

import certify_homogeneous_refinement_obstruction as certify
import verify_homogeneous_refinement_obstruction as verify


def row(index: int, normal: tuple[int, int], prime: int) -> dict:
    return {
        "h": index,
        "p": prime,
        "a": normal[0],
        "b": normal[1],
        "ord2": index,
        "ord3": index,
        "c": 0,
    }


class HomogeneousRefinementObstructionTests(unittest.TestCase):
    def assert_replays_same_counts(self, rows: list[dict]) -> None:
        first = certify.analyze_rows(rows)
        second = verify.replay(rows)
        for key in (
            "rows",
            "unique_leaf_lattices",
            "duplicate_leaf_rows",
            "parent_buckets",
            "best_completion_fraction",
            "refinement_primes",
            "multi_child_buckets",
            "max_observed_children",
        ):
            self.assertEqual(first[key], second[key])
        self.assertEqual(
            len(first["complete_sibling_groups"]),
            len(second["complete_sibling_groups"]),
        )

    def test_full_index_two_cover_is_detected(self) -> None:
        rows = [
            row(2, (1, 0), 5),
            row(2, (0, 1), 7),
            row(2, (1, 1), 11),
        ]
        self.assert_replays_same_counts(rows)
        result = certify.analyze_rows(rows)
        self.assertEqual(len(result["complete_sibling_groups"]), 1)
        group = result["complete_sibling_groups"][0]
        self.assertEqual(
            (
                group["parent_index"],
                group["child_index"],
                group["refinement_prime"],
                group["children"],
            ),
            (1, 2, 2, 3),
        )

    def test_two_binary_lifts_over_even_parent_are_detected(self) -> None:
        rows = [
            row(4, (1, 0), 5),
            row(4, (1, 2), 13),
        ]
        self.assert_replays_same_counts(rows)
        result = certify.analyze_rows(rows)
        self.assertEqual(len(result["complete_sibling_groups"]), 1)
        self.assertEqual(
            result["complete_sibling_groups"][0]["parent_index"],
            2,
        )

    def test_missing_sibling_is_an_obstruction(self) -> None:
        rows = [row(4, (1, 0), 5)]
        self.assert_replays_same_counts(rows)
        result = certify.analyze_rows(rows)
        self.assertFalse(result["complete_sibling_groups"])
        self.assertEqual(
            result["best_completion_fraction"],
            {"numerator": 1, "denominator": 2},
        )


if __name__ == "__main__":
    unittest.main()
