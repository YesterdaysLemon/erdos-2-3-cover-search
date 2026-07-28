#!/usr/bin/env python3

from __future__ import annotations

import copy
import unittest

import verify_anchor_phase_quotient


class AnchorPhaseQuotientTests(unittest.TestCase):
    def certificate(self) -> dict:
        return {
            "schema_version": 1,
            "row_columns": [
                "p",
                "h",
                "a",
                "b",
                "base_phase",
                "target_modulus",
                "target_residue",
            ],
            "rows": [
                [5, 2, 1, 0, 0, 1, 0],
                [7, 3, 0, 1, 0, 3, 0],
            ],
            "algebraic_primes": [],
            "anchors": [
                {
                    "p": 5,
                    "h": 2,
                    "a": 1,
                    "b": 0,
                    "target_modulus": 1,
                    "target_residue": 0,
                    "legal_targets": [0, 1],
                }
            ],
            "branches": [
                {
                    "targets": [0],
                    "witness": [1, 1],
                    "available_witnesses": 1,
                },
                {
                    "targets": [1],
                    "witness": [0, 1],
                    "available_witnesses": 1,
                },
            ],
            "summary": {
                "rows": 2,
                "anchor_rows": 1,
                "nonanchor_rows": 1,
                "legal_anchor_branches": 2,
                "open_anchor_branches": 0,
                "minimum_available_witnesses": 1,
                "maximum_available_witnesses": 1,
            },
        }

    def test_independent_verifier_accepts_exact_branch_witnesses(self):
        result = verify_anchor_phase_quotient.verify_certificate(
            self.certificate()
        )
        self.assertTrue(result["verified"])
        self.assertEqual(result["legal_anchor_branches"], 2)
        self.assertEqual(result["distinct_recorded_witnesses"], 2)

    def test_independent_verifier_rejects_a_covered_witness(self):
        certificate = copy.deepcopy(self.certificate())
        certificate["branches"][0]["witness"] = [0, 1]
        with self.assertRaises(ValueError):
            verify_anchor_phase_quotient.verify_certificate(certificate)

    def test_independent_verifier_rejects_an_omitted_branch(self):
        certificate = copy.deepcopy(self.certificate())
        certificate["branches"].pop()
        with self.assertRaises(ValueError):
            verify_anchor_phase_quotient.verify_certificate(certificate)


if __name__ == "__main__":
    unittest.main()
