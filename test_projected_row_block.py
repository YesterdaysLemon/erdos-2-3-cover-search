import unittest

from certify_projected_row_block import (
    best_projected_target,
    normalized_shared_projection_counts,
)


class ProjectedTargetSelectionTests(unittest.TestCase):
    def test_maximizes_new_cells_before_target_number(self):
        target, new_cells = best_projected_target(
            [0b111, 0b001],
            0,
        )
        self.assertEqual((target, new_cells), (0, 3))

    def test_breaks_coverage_ties_by_smallest_target(self):
        target, new_cells = best_projected_target(
            [0b101, 0b110],
            0,
        )
        self.assertEqual((target, new_cells), (0, 2))

    def test_normalized_shared_projection_partitions_uncovered_cells(self):
        single, both = normalized_shared_projection_counts(
            0b111111,
            0b001111,
            0b110011,
        )
        self.assertEqual(single, 4)
        self.assertEqual(both, 2)


if __name__ == "__main__":
    unittest.main()
