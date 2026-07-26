#!/usr/bin/env python3
"""Regression tests for proof-safe affine-plane component bounds."""

from __future__ import annotations

import collections
import unittest

import component_core
import verify_component_core_audit


class ParallelClassBoundTests(unittest.TestCase):
    def test_full_parallel_class_is_feasible(self) -> None:
        # Three distinct parallel lines already cover F_3^2.  Extra rows in
        # that direction may duplicate offsets but can never turn this into
        # an impossibility certificate.
        counts = collections.Counter(
            {
                (1, 0): 6,
                (1, 1): 2,
                (1, 2): 2,
                (0, 1): 1,
            }
        )
        self.assertFalse(component_core.cheap_impossible(3, counts))
        self.assertFalse(verify_component_core_audit.parallel_bound(counts, 3))
        self.assertTrue(
            component_core.normalized_plane_cover_enumeration(
                3, counts, 1_000_000
            )
        )

    def test_proper_parallel_classes_can_be_impossible(self) -> None:
        counts = collections.Counter({(1, 0): 2, (1, 1): 2})
        self.assertTrue(component_core.cheap_impossible(3, counts))
        self.assertTrue(verify_component_core_audit.parallel_bound(counts, 3))
        self.assertFalse(
            component_core.normalized_plane_cover_enumeration(
                3, counts, 1_000_000
            )
        )


if __name__ == "__main__":
    unittest.main()
