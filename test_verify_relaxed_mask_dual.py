import unittest

import verify_relaxed_mask_dual as verifier


class RelaxedMaskDualVerifierTests(unittest.TestCase):
    def test_inclusion_maximal_masks(self):
        masks = [0b0001, 0b0011, 0b0100, 0b1100]
        self.assertEqual(
            verifier.inclusion_maximal_masks(masks),
            [0b0011, 0b1100],
        )

    def test_scalar_mask_reconstruction(self):
        rows = [
            {
                "p": 5,
                "h": 3,
                "a": 1,
                "b": 0,
                "target_modulus": 1,
                "target_residue": 0,
            },
            {
                "p": 7,
                "h": 3,
                "a": 0,
                "b": 1,
                "target_modulus": 1,
                "target_residue": 0,
            },
        ]
        masks = verifier.reconstruct_gain_masks(
            rows,
            [(1, 1), (1, 2), (2, 1)],
            {5: 0, 7: 0},
            set(),
        )
        self.assertEqual(
            masks,
            [0b010, 0b011, 0b100, 0b101],
        )


if __name__ == "__main__":
    unittest.main()
