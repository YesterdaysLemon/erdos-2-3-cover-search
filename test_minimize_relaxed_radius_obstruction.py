import unittest

from minimize_relaxed_radius_obstruction import (
    greedy_unsat_assumption_core,
)


class MinimizeRelaxedRadiusObstructionTest(unittest.TestCase):
    def test_greedy_core_removes_redundant_assumption(self):
        try:
            from pysat.solvers import Solver
        except ModuleNotFoundError:
            self.skipTest("PySAT unavailable")
        solver = Solver(
            name="cadical195",
            bootstrap_with=[
                [-1, 4],
                [-2, -4],
                [-3, 5],
            ],
        )
        core = greedy_unsat_assumption_core(solver, [1, 2, 3])
        solver.delete()
        self.assertEqual(set(core), {1, 2})


if __name__ == "__main__":
    unittest.main()
