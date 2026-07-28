import unittest

import quantified_phase_cover_z3 as quantified


class QuantifiedPhaseCoverZ3Tests(unittest.TestCase):
    def test_one_parity_fibre_cannot_cover(self):
        result = quantified.solve_rows(
            [
                {
                    "h": 2,
                    "a": 1,
                    "b": 0,
                    "target_modulus": 1,
                    "target_residue": 0,
                }
            ],
            10_000,
        )
        self.assertEqual(result["check"], "unsat")

    def test_two_parity_fibres_can_cover(self):
        result = quantified.solve_rows(
            [
                {
                    "h": 2,
                    "a": 1,
                    "b": 0,
                    "target_modulus": 1,
                    "target_residue": 0,
                },
                {
                    "h": 2,
                    "a": 1,
                    "b": 0,
                    "target_modulus": 1,
                    "target_residue": 0,
                },
            ],
            10_000,
        )
        self.assertEqual(result["check"], "sat")
        self.assertEqual(set(result["targets"]), {0, 1})


if __name__ == "__main__":
    unittest.main()
