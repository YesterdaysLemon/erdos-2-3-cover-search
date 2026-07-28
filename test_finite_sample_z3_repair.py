#!/usr/bin/env python3
"""Regression tests for sparse Z3 pseudo-Boolean repair."""

from __future__ import annotations

import unittest

import finite_sample_z3_repair

try:
    import z3  # noqa: F401
except ImportError:  # pragma: no cover - dependency-gated environment
    z3 = None


@unittest.skipIf(z3 is None, "Z3 is unavailable")
class FiniteSampleZ3RepairTests(unittest.TestCase):
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

    def test_one_change_splits_two_identical_rows(self) -> None:
        result = finite_sample_z3_repair.solve_repair(
            self.rows,
            self.candidates,
            self.points,
            {5: 0, 7: 0},
            set(),
            1,
            10.0,
        )
        self.assertEqual(result["check"], "sat")
        self.assertEqual(result["changed_phases"], 1)
        self.assertEqual(result["full_misses"], 0)
        self.assertEqual(set(result["repaired"]), {0, 1})

    def test_fixed_rows_return_point_core(self) -> None:
        result = finite_sample_z3_repair.solve_repair(
            self.rows,
            self.candidates,
            self.points,
            {5: 0, 7: 0},
            {5, 7},
            1,
            10.0,
        )
        self.assertEqual(result["check"], "unsat")
        self.assertEqual(result["core_indices"], [1])


if __name__ == "__main__":
    unittest.main()
