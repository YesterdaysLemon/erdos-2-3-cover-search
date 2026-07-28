#!/usr/bin/env python3

from __future__ import annotations

import copy
import unittest

import verify_anchor_phase_master_unsat


class AnchorPhaseMasterUnsatTests(unittest.TestCase):
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
            "rows": [[5, 2, 1, 0, 0, 1, 0]],
            "algebraic_primes": [],
            "anchors": [
                {
                    "p": 5,
                    "h": 2,
                    "a": 1,
                    "b": 0,
                    "target_modulus": 1,
                    "target_residue": 0,
                    "legal_target_count": 2,
                }
            ],
            "points": [[0, 0], [1, 0]],
            "master": {
                "sat": False,
                "eligible_points": 2,
                "empty_clause_point": None,
            },
            "summary": {
                "rows": 1,
                "anchor_rows": 1,
                "frozen_rows": 0,
                "core_points": 2,
                "legal_target_space_sizes": [2],
                "joint_legal_target_assignments": 2,
                "pysat_unsat": True,
            },
        }

    def test_independent_integer_master_confirms_unsat(self):
        result = verify_anchor_phase_master_unsat.verify_certificate(
            self.certificate()
        )
        self.assertTrue(result["verified"])
        self.assertEqual(result["z3_verdict"], "unsat")

    def test_verifier_rejects_a_satisfiable_tamper(self):
        certificate = copy.deepcopy(self.certificate())
        certificate["points"].pop()
        certificate["summary"]["core_points"] = 1
        certificate["master"]["eligible_points"] = 1
        with self.assertRaises(ValueError):
            verify_anchor_phase_master_unsat.verify_certificate(certificate)


if __name__ == "__main__":
    unittest.main()
