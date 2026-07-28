#!/usr/bin/env python3
"""Tests for projected conditional-design discovery."""

from __future__ import annotations

import argparse
import math
import unittest
from fractions import Fraction
from pathlib import Path

from search_projected_conditional_designs import (
    Design,
    Projected,
    candidate_designs,
    evaluate_design,
    projectable,
    search_output_payload,
)
from verify_projected_conditional_fibre_overlap import (
    primitive_kernel_residues,
)


class ProjectedConditionalDesignTests(unittest.TestCase):
    def test_cyclic_kernel_enumeration_matches_brute_force(self) -> None:
        for h in range(2, 18):
            for a in range(h):
                for b in range(h):
                    if math.gcd(a, b, h) != 1:
                        continue
                    row = {"h": h, "a": a, "b": b}
                    expected = {
                        (k, ell)
                        for k in range(h)
                        for ell in range(h)
                        if (a * k + b * ell) % h == 0
                    }
                    actual = set(primitive_kernel_residues(row))
                    self.assertEqual(actual, expected)
                    self.assertEqual(len(actual), h)

    def test_checkpoint_payload_records_completed_primes(self) -> None:
        args = argparse.Namespace(
            pool=Path("pool.json"),
            block_certificate=Path("block.json"),
            period_certificate=Path("period.json"),
            max_base_anchors=1,
            max_base_period=100,
            max_base_points=200,
            max_tensor_cells=300,
            jitter_seeds=4,
            max_designs=5,
        )
        results = [
            {"outside_prime": 37},
            {"outside_prime": 61},
        ]
        payload = search_output_payload(
            args,
            [37, 61, 89],
            results,
            "discovery_in_progress",
        )
        self.assertEqual(payload["completed_outside_primes"], [37, 61])
        self.assertEqual(payload["results"], results)

    def test_projectable_uses_exact_coprime_residual_split(self) -> None:
        row = {"p": 13, "h": 12, "a": 1, "b": 4}
        self.assertEqual(projectable(row, 20), Projected(13, 4, 3))

        # This row's residual equation is not primitive modulo three.
        nonprimitive = {"p": 13, "h": 12, "a": 3, "b": 6}
        self.assertIsNone(projectable(nonprimitive, 4))

    def test_small_design_matches_complete_phase_enumeration(self) -> None:
        outside = {"p": 2, "h": 2, "a": 1, "b": 0}
        rows = {
            3: {"p": 3, "h": 3, "a": 0, "b": 1},
            5: {"p": 5, "h": 5, "a": 1, "b": 1},
            7: {"p": 7, "h": 7, "a": 2, "b": 1},
        }
        design = Design(
            normalizer=3,
            base_primes=(5,),
            projected=(Projected(7, 1, 7),),
            base_period=30,
            base_points=450,
            tensor_cells=5,
            density_score=Fraction(1, 3)
            + Fraction(1, 5)
            + Fraction(1, 7),
        )
        evaluated = evaluate_design(outside, design, rows)
        self.assertIsNotNone(evaluated)
        actual, _detail = evaluated

        period = math.lcm(2, 3, 5, 7)
        expected_count = period * period
        for target5 in range(5):
            for target7 in range(7):
                covered = 0
                for first in range(period):
                    if first % 2:
                        continue
                    for second in range(period):
                        if (
                            second % 3 == 0
                            or (first + second - target5) % 5 == 0
                            or (
                                2 * first + second - target7
                            ) % 7 == 0
                        ):
                            covered += 1
                expected_count = min(expected_count, covered)
        self.assertEqual(
            actual,
            Fraction(expected_count, period * period),
        )

    def test_shared_residual_pair_matches_complete_phase_enumeration(
        self,
    ) -> None:
        outside = {"p": 2, "h": 2, "a": 1, "b": 0}
        rows = {
            3: {"p": 3, "h": 3, "a": 0, "b": 1},
            7: {"p": 7, "h": 14, "a": 1, "b": 1},
            11: {"p": 11, "h": 7, "a": 2, "b": 1},
        }
        design = Design(
            normalizer=3,
            base_primes=(),
            projected=(),
            base_period=6,
            base_points=18,
            tensor_cells=2,
            density_score=Fraction(1, 3)
            + Fraction(1, 14)
            + Fraction(1, 7),
            paired_projected=Projected(7, 2, 7),
            paired_shared=11,
        )
        evaluated = evaluate_design(outside, design, rows)
        self.assertIsNotNone(evaluated)
        actual, _detail = evaluated

        period = math.lcm(2, 3, 14, 7)
        expected_count = period * period
        for target7 in range(14):
            for target11 in range(7):
                covered = 0
                for first in range(period):
                    if first % 2:
                        continue
                    for second in range(period):
                        if (
                            second % 3 == 0
                            or (
                                first + second - target7
                            ) % 14 == 0
                            or (
                                2 * first + second - target11
                            ) % 7 == 0
                        ):
                            covered += 1
                expected_count = min(expected_count, covered)
        self.assertEqual(
            actual,
            Fraction(expected_count, period * period),
        )

    def test_pair_candidate_requires_residual_coprime_to_base_period(
        self,
    ) -> None:
        outside = {"p": 2, "h": 2, "a": 1, "b": 0}
        anchors = [
            {"p": 3, "h": 3, "a": 0, "b": 1},
            {"p": 5, "h": 4, "a": 1, "b": 1},
            {"p": 7, "h": 2, "a": 0, "b": 1},
        ]
        designs = candidate_designs(
            outside,
            anchors,
            max_base_anchors=0,
            max_base_period=100,
            max_base_points=1_000,
            max_tensor_cells=1_000,
            jitter_seeds=1,
        )
        self.assertTrue(designs)
        self.assertTrue(
            all(design.paired_projected is None for design in designs)
        )


if __name__ == "__main__":
    unittest.main()
