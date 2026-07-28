import os
import sys
import unittest
from pathlib import Path

sys.path.insert(
    0,
    str(Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"),
)
import numpy as np

from enumerate_saved_mask_repairs import enumerate_owner_assignments


class EnumerateSavedMaskRepairsTests(unittest.TestCase):
    def test_enumerates_distinct_row_owner_repair(self):
        result = enumerate_owner_assignments(
            models=[[0b01, 0b10]],
            support={
                0b01: [(0, 1)],
                0b10: [(1, 1)],
            },
            assignment=[0, 0],
            targets=np.asarray([[1, 0], [0, 1]]),
            rows=[{"p": 17}, {"p": 19}],
            np=np,
        )
        self.assertEqual(result["distinct_row_owner_trials"], 1)
        self.assertEqual(result["replay_failures"], 0)
        self.assertEqual(result["exact_repair_count"], 1)
        self.assertEqual(
            result["exact_repairs"][0]["moves"],
            [{"p": 17, "target": 1}, {"p": 19, "target": 1}],
        )

    def test_rejects_same_row_owner_product(self):
        result = enumerate_owner_assignments(
            models=[[0b01, 0b10]],
            support={
                0b01: [(0, 1)],
                0b10: [(0, 2)],
            },
            assignment=[0],
            targets=np.asarray([[1], [2]]),
            rows=[{"p": 17}],
            np=np,
        )
        self.assertEqual(result["raw_owner_products"], 1)
        self.assertEqual(result["distinct_row_owner_trials"], 0)
        self.assertEqual(result["exact_repair_count"], 0)


if __name__ == "__main__":
    unittest.main()
