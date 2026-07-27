import unittest

from verify_layered_exact_misses import verify_layered_misses


class VerifyLayeredExactMissesTest(unittest.TestCase):
    def test_filters_intercepted_base_hole(self):
        base = {
            "choices": [
                {
                    "p": 5,
                    "h": 2,
                    "a": 1,
                    "b": 0,
                    "target_modulus": 1,
                    "target_residue": 0,
                }
            ]
        }
        augmented = {
            "choices": [
                base["choices"][0],
                {
                    "p": 7,
                    "h": 3,
                    "a": 0,
                    "b": 1,
                    "target_modulus": 1,
                    "target_residue": 0,
                },
            ]
        }
        result = verify_layered_misses(
            base,
            augmented,
            {5: 0, 7: 1},
            [(1, 1), (1, 2)],
        )
        self.assertTrue(result["verified"])
        self.assertEqual(result["exact_augmented_misses"], [[1, 2]])
        self.assertEqual(result["intercepted_count"], 1)

    def test_rejects_false_base_hole(self):
        base = {
            "choices": [
                {
                    "p": 5,
                    "h": 2,
                    "a": 1,
                    "b": 0,
                    "target_modulus": 1,
                    "target_residue": 0,
                }
            ]
        }
        with self.assertRaisesRegex(RuntimeError, "covered by base"):
            verify_layered_misses(base, base, {5: 0}, [(0, 1)])


if __name__ == "__main__":
    unittest.main()
