import unittest

from select_targeted_signature_rows import (
    best_legal_target,
    largest_prime_power_component,
    select_rows,
    target_hit_indices,
)


class TargetedSignatureRowsTest(unittest.TestCase):
    def test_best_target_respects_target_congruence(self):
        row = {
            "p": 11,
            "h": 6,
            "a": 1,
            "b": 1,
            "target_modulus": 2,
            "target_residue": 1,
        }
        points = [(0, 1), (1, 0), (2, 1), (0, 0)]
        count, target, indices = best_legal_target(row, points)
        self.assertEqual((count, target, indices), (2, 1, [0, 1]))

    def test_selector_omits_single_witness_memorizers(self):
        base = [{"p": 5, "h": 2, "a": 1, "b": 0}]
        candidates = [
            {"p": 7, "h": 3, "a": 1, "b": 0},
            {"p": 11, "h": 5, "a": 1, "b": 0},
        ]
        points = [(0, 0), (3, 1), (1, 2)]
        rows, audit = select_rows(
            base,
            candidates,
            points,
            min_hit=2,
            max_rows=0,
        )
        self.assertEqual([int(row["p"]) for row in rows], [5, 7])
        self.assertEqual(audit[0]["hit_count"], 2)
        self.assertEqual(audit[0]["point_indices"], [0, 1])

    def test_selector_prefers_smaller_modulus_on_equal_score(self):
        base = []
        candidates = [
            {"p": 13, "h": 7, "a": 0, "b": 1},
            {"p": 17, "h": 5, "a": 0, "b": 1},
        ]
        rows, audit = select_rows(
            base,
            candidates,
            [(0, 0), (1, 0)],
            min_hit=2,
            max_rows=1,
        )
        self.assertEqual([int(row["p"]) for row in rows], [17])
        self.assertEqual(audit[0]["target"], 0)

    def test_component_guard_rejects_expensive_row(self):
        candidates = [
            {"p": 19, "h": 27, "a": 0, "b": 1},
            {"p": 23, "h": 20, "a": 0, "b": 1},
        ]
        rows, audit = select_rows(
            [],
            candidates,
            [(0, 0), (1, 0)],
            min_hit=2,
            max_rows=0,
            max_component=16,
        )
        self.assertEqual([int(row["p"]) for row in rows], [23])
        self.assertEqual(audit[0]["largest_prime_power_component"], 5)
        self.assertEqual(largest_prime_power_component(884736), 32768)

    def test_validation_uses_training_selected_target(self):
        candidates = [
            {"p": 29, "h": 5, "a": 1, "b": 0},
            {"p": 31, "h": 7, "a": 1, "b": 0},
        ]
        rows, audit = select_rows(
            [],
            candidates,
            [(0, 0), (5, 1), (1, 0)],
            min_hit=2,
            max_rows=0,
            validation_points=[(5, 2), (6, 3)],
            min_validation_hit=1,
        )
        self.assertEqual([int(row["p"]) for row in rows], [29])
        self.assertEqual(audit[0]["target"], 0)
        self.assertEqual(audit[0]["validation_hit_count"], 1)
        self.assertEqual(
            target_hit_indices(candidates[0], [(5, 2), (6, 3)], 0),
            [0],
        )


if __name__ == "__main__":
    unittest.main()
