#!/usr/bin/env python3
"""Tests for the resumable sparse low-anchor quotient sweep."""

from __future__ import annotations

import unittest

import sparse_anchor_quotient_sweep


class LowAnchorRepresentativeTests(unittest.TestCase):
    def test_exact_quotient_has_162_unique_representatives(self) -> None:
        representatives = (
            sparse_anchor_quotient_sweep.low_anchor_representatives()
        )
        self.assertEqual(len(representatives), 162)
        self.assertEqual(len(set(representatives)), 162)

    def test_inactive_p17_parity_has_one_canonical_target(self) -> None:
        representatives = (
            sparse_anchor_quotient_sweep.low_anchor_representatives()
        )
        for p41 in range(2):
            observed = {
                targets[1]
                for targets in representatives
                if targets[0] == p41
            }
            self.assertEqual(
                observed,
                {
                    p41,
                    *(target for target in range(4) if target % 2 != p41),
                },
            )
            self.assertNotIn(p41 + 2, observed)

    def test_branch_slugs_are_stable_and_distinct(self) -> None:
        representatives = (
            sparse_anchor_quotient_sweep.low_anchor_representatives()
        )
        slugs = {
            sparse_anchor_quotient_sweep.branch_slug(index, targets)
            for index, targets in enumerate(representatives)
        }
        self.assertEqual(len(slugs), 162)
        self.assertEqual(
            sparse_anchor_quotient_sweep.branch_slug(
                0,
                representatives[0],
            ),
            "branch_000_0_0_0_0_0",
        )


if __name__ == "__main__":
    unittest.main()
