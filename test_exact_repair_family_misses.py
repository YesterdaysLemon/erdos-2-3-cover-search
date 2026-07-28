import unittest

from exact_repair_family_misses import build_phase_maps


class ExactRepairFamilyMissesTests(unittest.TestCase):
    def test_builds_one_phase_per_repair(self):
        base = {17: 1, 19: 2}
        repairs = [
            {"moves": [{"p": 17, "target": 3}]},
            {"moves": [{"p": 19, "target": 4}]},
        ]
        self.assertEqual(
            build_phase_maps(base, repairs),
            [{17: 3, 19: 2}, {17: 1, 19: 4}],
        )

    def test_rejects_duplicate_updates(self):
        with self.assertRaises(ValueError):
            build_phase_maps(
                {17: 1},
                [{"moves": [{"p": 17, "target": 2}, {"p": 17, "target": 3}]}],
            )


if __name__ == "__main__":
    unittest.main()
