import unittest

import generate_diverse_mask_repairs as diverse

try:
    import pysat  # noqa: F401
except ImportError:  # pragma: no cover - dependency-gated environment
    pysat = None


@unittest.skipIf(pysat is None, "PySAT is unavailable")
class GenerateDiverseMaskRepairsTests(unittest.TestCase):
    def test_generates_repair_separated_from_seed(self):
        rows = [
            {
                "h": 2,
                "p": prime,
                "a": 1,
                "b": 0,
                "ord2": 2,
                "ord3": 1,
                "target_modulus": 1,
                "target_residue": 0,
            }
            for prime in (5, 7)
        ]
        candidates = [
            (
                row["h"],
                row["p"],
                row["a"],
                row["b"],
                row["ord2"],
                row["ord3"],
            )
            for row in rows
        ]
        generated, results = diverse.generate_repairs(
            rows,
            candidates,
            [(0, 0), (1, 0)],
            {5: 0, 7: 0},
            set(),
            1,
            "cadical195",
            100,
            10.0,
            [[1, 0]],
            2,
            1,
        )
        self.assertEqual(len(generated), 1)
        self.assertEqual(generated[0], [0, 1])
        self.assertEqual(results[0]["status"], "INTEGER_MODEL")
        self.assertEqual(
            diverse.hamming_distance(generated[0], [1, 0]),
            2,
        )


if __name__ == "__main__":
    unittest.main()
