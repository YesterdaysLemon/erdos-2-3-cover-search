import unittest

from certify_block_plus_overlap_forest import DisjointSet
from verify_block_plus_overlap_forest import CycleChecker


class OverlapForestTests(unittest.TestCase):
    def test_generator_disjoint_set_rejects_cycle(self):
        dsu = DisjointSet([2, 3, 5, 7])
        self.assertTrue(dsu.union(2, 3))
        self.assertTrue(dsu.union(3, 5))
        self.assertFalse(dsu.union(2, 5))
        self.assertTrue(dsu.union(5, 7))

    def test_verifier_cycle_checker_is_independent(self):
        checker = CycleChecker({11, 13, 17, 19})
        self.assertTrue(checker.add_edge(11, 13))
        self.assertTrue(checker.add_edge(17, 19))
        self.assertTrue(checker.add_edge(13, 17))
        self.assertFalse(checker.add_edge(11, 19))

    def test_block_vertex_participates_in_cycle_check(self):
        checker = CycleChecker({0, 23, 29, 31})
        self.assertTrue(checker.add_edge(0, 29))
        self.assertTrue(checker.add_edge(29, 31))
        self.assertFalse(checker.add_edge(0, 31))
        self.assertTrue(checker.add_edge(0, 23))


if __name__ == "__main__":
    unittest.main()
