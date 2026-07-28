#!/usr/bin/env python3
"""Tests for phase-uniform sixth-order conditional bounds."""

from __future__ import annotations

import itertools
import math
import unittest
from fractions import Fraction

from search_sixterm_conditional_bounds import (
    sixterm_subset_lower,
    target_image_index,
)


def actual_intersection(
    outside: dict,
    anchors: tuple[dict, ...],
    targets: tuple[int, ...],
) -> Fraction:
    period = math.lcm(
        int(outside["h"]),
        *(int(row["h"]) for row in anchors),
    )
    covered = 0
    for k in range(period):
        for l in range(period):
            if (
                int(outside["a"]) * k + int(outside["b"]) * l
            ) % int(outside["h"]):
                continue
            if any(
                (
                    int(row["a"]) * k + int(row["b"]) * l - target
                ) % int(row["h"])
                == 0
                for row, target in zip(anchors, targets)
            ):
                covered += 1
    return Fraction(covered, period * period)


class SixtermConditionalTests(unittest.TestCase):
    def test_target_image_index_matches_small_surjective_map(self) -> None:
        rows = (
            {"h": 2, "a": 1, "b": 0},
            {"h": 3, "a": 0, "b": 1},
        )
        self.assertEqual(target_image_index(rows), 1)
        repeated = rows + ({"h": 2, "a": 1, "b": 0},)
        self.assertEqual(target_image_index(repeated), 2)

    def test_sixterm_bound_holds_for_every_small_target_choice(self) -> None:
        outside = {"p": 2, "h": 2, "a": 1, "b": 0}
        anchors = (
            {"p": 3, "h": 2, "a": 0, "b": 1},
            {"p": 5, "h": 2, "a": 1, "b": 1},
            {"p": 7, "h": 2, "a": 0, "b": 1},
            {"p": 11, "h": 2, "a": 1, "b": 1},
            {"p": 13, "h": 2, "a": 0, "b": 1},
            {"p": 17, "h": 2, "a": 1, "b": 1},
        )
        lower = sixterm_subset_lower(outside, anchors)
        for targets in itertools.product(range(2), repeat=len(anchors)):
            self.assertLessEqual(
                lower,
                actual_intersection(outside, anchors, targets),
            )


if __name__ == "__main__":
    unittest.main()
