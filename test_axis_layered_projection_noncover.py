import unittest
import json
import tempfile
from fractions import Fraction
from pathlib import Path

from certify_axis_layered_projection_noncover import (
    active_projection_rows,
    build_certificate,
    maximum_weighted_tail,
)
from verify_axis_layered_projection_noncover import replay


def toy_payload():
    return {
        "schema": "axis_layered_pool_v1",
        "layer_axis": "k",
        "layer_period": 6,
        "capacity_pattern_period": 1,
        "coordinate_basis": None,
        "choices": [
            {
                "p": 5,
                "h": 2,
                "a": 1,
                "b": 1,
                "layer_active_class": 0,
                "target_modulus": 1,
                "target_residue": 0,
            },
            {
                "p": 7,
                "h": 3,
                "a": 1,
                "b": 1,
                "layer_active_class": 0,
                "target_modulus": 1,
                "target_residue": 0,
            },
            {
                "p": 13,
                "h": 6,
                "a": 1,
                "b": 1,
                "layer_active_class": 0,
                "target_modulus": 1,
                "target_residue": 0,
            },
        ],
    }


class AxisLayeredProjectionNoncoverTests(unittest.TestCase):
    def test_projection_rows(self):
        rows = active_projection_rows(toy_payload(), 0, 6)
        self.assertEqual(
            [
                (
                    row["p"],
                    row["projected_modulus"],
                    row["conditional_denominator"],
                )
                for row in rows
            ],
            [(5, 2, 1), (7, 3, 1), (13, 6, 1)],
        )

    def test_exact_weight_bound(self):
        rows = [
            {
                "projected_modulus": 2,
                "conditional_denominator": 3,
            }
        ]
        self.assertEqual(
            maximum_weighted_tail(rows, [1, 4, 2, 3]),
            Fraction(7, 3),
        )

    def test_toy_certificate_replays(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "toy.json"
            source.write_text(json.dumps(toy_payload()))
            certificate = build_certificate(
                source,
                toy_payload(),
                coordinate=0,
                projection=6,
                base_anchor_primes=(5, 7),
                branch_anchor_prime=13,
            )
        report = replay(toy_payload(), certificate)
        self.assertTrue(report["verified"])
        self.assertEqual(report["branch_count"], 6)
        self.assertTrue(report["proved_no_declared_layered_cover"])


if __name__ == "__main__":
    unittest.main()
