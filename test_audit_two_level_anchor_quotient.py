#!/usr/bin/env python3

from __future__ import annotations

import unittest

import audit_two_level_anchor_quotient as audit


class TwoLevelAnchorQuotientTests(unittest.TestCase):
    def test_fingerprint_filter_discharges_outer_and_algebraic_points(self):
        rows = [
            {"p": 5, "h": 2, "a": 1, "b": 0},
            {"p": 7, "h": 2, "a": 0, "b": 1},
            {"p": 11, "h": 2, "a": 1, "b": 1},
        ]
        fingerprints, metadata = audit.point_fingerprints(
            rows,
            {5: 0, 7: 0, 11: 0},
            (5,),
            (7,),
            [(1, 0), (1, 1), (3, 3)],
            (3,),
        )
        self.assertEqual(fingerprints, {((1,), (0,))})
        self.assertEqual(metadata["frozen_discharged_points"], 1)
        self.assertEqual(metadata["algebraically_discharged_points"], 1)

    def test_single_incremental_model_audits_all_low_branches(self):
        rows = [
            {"p": 5, "h": 2, "a": 1, "b": 0},
            {"p": 7, "h": 2, "a": 0, "b": 1},
        ]
        fingerprints = {
            ((0,), (0,)),
            ((0,), (1,)),
            ((1,), (0,)),
        }
        result = audit.audit_fingerprints(
            rows,
            {5: 0, 7: 0},
            (5,),
            (7,),
            fingerprints,
        )
        by_branch = {
            tuple(outcome["branch"]): outcome
            for outcome in result["outcomes"]
        }
        self.assertTrue(by_branch[(0,)]["sat"])
        self.assertFalse(by_branch[(1,)]["sat"])
        self.assertEqual(result["summary"]["unsat_branches"], 1)

    def test_restricted_illegal_target_can_make_an_empty_clause(self):
        rows = [
            {
                "p": 5,
                "h": 2,
                "a": 1,
                "b": 0,
                "target_modulus": 2,
                "target_residue": 0,
            },
            {
                "p": 7,
                "h": 2,
                "a": 0,
                "b": 1,
                "target_modulus": 2,
                "target_residue": 0,
            },
        ]
        result = audit.audit_fingerprints(
            rows,
            {5: 0, 7: 0},
            (5,),
            (7,),
            {((1,), (1,))},
        )
        self.assertEqual(result["summary"]["low_branches"], 1)
        self.assertEqual(result["summary"]["unsat_branches"], 1)

    def test_binary_refinement_symmetry_detects_the_41_17_relation(self):
        rows = (
            (41, 2, 0, 1, 1, 0),
            (17, 4, 2, 3, 1, 0),
            (19, 3, 1, 1, 1, 0),
        )
        relations = audit.binary_refinement_symmetries(rows)
        self.assertEqual(len(relations), 1)
        relation = relations[0]
        self.assertEqual(relation["parent_prime"], 41)
        self.assertEqual(relation["child_prime"], 17)
        self.assertEqual(relation["raw_local_branches"], 8)
        self.assertEqual(relation["quotient_local_classes"], 6)
        self.assertEqual(
            relation["inactive_child_targets_by_parent_phase"],
            {"0": [0, 2], "1": [1, 3]},
        )


if __name__ == "__main__":
    unittest.main()
