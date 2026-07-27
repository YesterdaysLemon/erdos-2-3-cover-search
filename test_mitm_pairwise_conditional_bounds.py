#!/usr/bin/env python3
"""Tests for meet-in-the-middle pairwise subset selection."""

from __future__ import annotations

import itertools
import random
import unittest

import numpy as np

from search_mitm_pairwise_conditional_bounds import (
    mitm_quadratic_subset,
)


def brute_quadratic_subset(
    weights: object,
    penalties: object,
) -> tuple[tuple[int, ...], float]:
    best_subset: tuple[int, ...] = ()
    best = 0.0
    for mask in range(1 << len(weights)):
        selected = tuple(
            index
            for index in range(len(weights))
            if (mask >> index) & 1
        )
        value = sum(float(weights[index]) for index in selected)
        value -= sum(
            float(penalties[first, second])
            for first, second in itertools.combinations(selected, 2)
        )
        if value > best:
            best = value
            best_subset = selected
    return best_subset, best


class MitmPairwiseConditionalTests(unittest.TestCase):
    def test_random_small_instances_match_complete_enumeration(self) -> None:
        rng = random.Random(901)
        for size in range(1, 10):
            for _ in range(8):
                weights = np.array(
                    [rng.random() for _ in range(size)],
                    dtype=np.float64,
                )
                penalties = np.zeros((size, size), dtype=np.float64)
                for first in range(size):
                    for second in range(first):
                        value = rng.random() / 3
                        penalties[first, second] = value
                        penalties[second, first] = value
                selected, value = mitm_quadratic_subset(
                    weights,
                    penalties,
                    chunk_size=3,
                )
                expected_selected, expected = brute_quadratic_subset(
                    weights,
                    penalties,
                )
                self.assertAlmostEqual(value, expected, places=12)
                selected_value = sum(weights[index] for index in selected)
                selected_value -= sum(
                    penalties[first, second]
                    for first, second in itertools.combinations(selected, 2)
                )
                self.assertAlmostEqual(selected_value, expected, places=12)
                if selected != expected_selected:
                    self.assertAlmostEqual(selected_value, expected, places=12)


if __name__ == "__main__":
    unittest.main()
