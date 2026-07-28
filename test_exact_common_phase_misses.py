import unittest

from exact_common_phase_misses import (
    build_union_rows,
    replay_common_misses,
)


class ExactCommonPhaseMissesTest(unittest.TestCase):
    def test_union_deduplicates_shared_fibres(self):
        candidates = [
            (2, 5, 1, 0, 2, 1),
            (3, 7, 0, 1, 1, 3),
        ]
        rows = {
            5: {
                "p": 5,
                "h": 2,
                "target_modulus": 1,
                "target_residue": 0,
            },
            7: {
                "p": 7,
                "h": 3,
                "target_modulus": 1,
                "target_residue": 0,
            },
        }
        union = build_union_rows(
            candidates,
            rows,
            [{5: 0, 7: 1}, {5: 0, 7: 2}],
        )
        self.assertEqual(len(union), 3)
        self.assertEqual(
            {(int(row["p"]), int(row["c"])) for row in union},
            {(5, 0), (7, 1), (7, 2)},
        )

    def test_scalar_replay_requires_a_miss_of_every_fibre(self):
        rows = [
            {"p": 5, "h": 2, "a": 1, "b": 0, "c": 0},
            {"p": 7, "h": 3, "a": 0, "b": 1, "c": 1},
        ]
        self.assertTrue(replay_common_misses(rows, [(1, 2)], (), False))
        self.assertFalse(replay_common_misses(rows, [(0, 2)], (), False))

    def test_scalar_replay_checks_algebraic_exclusions(self):
        self.assertFalse(replay_common_misses([], [(7, 14)], (7,), False))
        self.assertFalse(replay_common_misses([], [(2, 4)], (), True))


if __name__ == "__main__":
    unittest.main()
