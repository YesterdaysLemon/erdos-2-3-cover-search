import unittest

from materialize_axis_layered_column import materialize_column


class MaterializeAxisLayeredColumnTests(unittest.TestCase):
    def test_materializes_active_rows_and_exact_capacity(self):
        payload = {
            "schema": "axis_layered_pool_v1",
            "layer_axis": "k",
            "layer_period": 12,
            "capacity_pattern_period": 6,
            "coordinate_basis": None,
            "choices": [
                {
                    "p": 7,
                    "h": 6,
                    "a": 1,
                    "b": 2,
                    "layer_active_class": 1,
                    "target_modulus": 2,
                    "target_residue": 1,
                },
                {
                    "p": 13,
                    "h": 12,
                    "a": 1,
                    "b": 3,
                    "layer_active_class": 0,
                    "target_modulus": 3,
                    "target_residue": 0,
                },
            ],
        }
        result = materialize_column(payload, 3)
        self.assertEqual(result["coordinate"], 3)
        self.assertEqual(result["row_count"], 2)
        self.assertEqual(result["transverse_period"], 12)
        self.assertEqual(result["capacity_numerator"], 7)
        self.assertEqual(result["capacity_denominator"], 12)
        self.assertEqual(
            [(row["p"], row["h"]) for row in result["choices"]],
            [(7, 3), (13, 4)],
        )

    def test_rejects_inconsistent_target_restriction(self):
        payload = {
            "schema": "axis_layered_pool_v1",
            "layer_axis": "k",
            "layer_period": 2,
            "capacity_pattern_period": 2,
            "choices": [
                {
                    "p": 5,
                    "h": 4,
                    "a": 1,
                    "b": 2,
                    "layer_active_class": 1,
                    "target_modulus": 2,
                    "target_residue": 0,
                }
            ],
        }
        with self.assertRaisesRegex(ValueError, "target restriction"):
            materialize_column(payload, 1)


if __name__ == "__main__":
    unittest.main()
