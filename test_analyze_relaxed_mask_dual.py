import unittest

import analyze_relaxed_mask_dual as dual


class RelaxedMaskDualTests(unittest.TestCase):
    def test_exact_summary_certifies_uniform_three_singletons(self):
        summary = dual.summarize_exact_certificate(
            [0b001, 0b010, 0b100],
            [1, 1, 1],
            2,
            3,
        )
        self.assertTrue(summary["certified"])
        self.assertEqual(summary["strict_weight_gap"], 1)
        self.assertEqual(summary["maximum_mask_cardinality"], 1)

    def test_exact_summary_rejects_nonstrict_bound(self):
        summary = dual.summarize_exact_certificate(
            [0b0011, 0b1100],
            [1, 1, 1, 1],
            2,
            4,
        )
        self.assertFalse(summary["certified"])
        self.assertEqual(summary["strict_weight_gap"], 0)
        self.assertFalse(summary["pairwise_union_certified"])
        self.assertTrue(summary["tight_disjoint_cover_exists"])
        self.assertFalse(summary["tight_equality_certified"])

    def test_tight_equality_without_disjoint_partition_certifies(self):
        summary = dual.summarize_exact_certificate(
            [0b0011, 0b0101, 0b1001],
            [1, 1, 1, 1],
            2,
            4,
        )
        self.assertEqual(summary["strict_weight_gap"], 0)
        self.assertTrue(summary["tight_equality_case"])
        self.assertFalse(summary["tight_disjoint_cover_exists"])
        self.assertTrue(summary["tight_equality_certified"])

    def test_pairwise_maximal_masks_compress_subsets(self):
        summary = dual.summarize_exact_certificate(
            [0b0001, 0b0011, 0b0100, 0b1100],
            [1, 1, 1, 1],
            2,
            4,
        )
        self.assertEqual(summary["inclusion_maximal_mask_count"], 2)
        self.assertEqual(
            summary["inclusion_maximal_masks"],
            [[0, 1], [2, 3]],
        )
        self.assertFalse(summary["pairwise_union_certified"])

    def test_integerize_simple_rationals(self):
        weights = dual.integerize_weights(
            [0.25, 0.5, 0.25],
            100,
        )
        self.assertEqual(weights, [1, 2, 1])


if __name__ == "__main__":
    unittest.main()
