import math
import unittest
from fractions import Fraction

from certify_factorized_pair_extension import component_pair_union


def brute_union(shared, residual, factor, lifted):
    covered = 0
    total = shared * shared * residual * residual
    for k_shared in range(shared):
        for l_shared in range(shared):
            first = (
                factor[0] * k_shared + factor[1] * l_shared
            ) % shared == 0
            second_shared = (
                lifted[0] * k_shared + lifted[1] * l_shared
            ) % shared == 0
            for k_residual in range(residual):
                for l_residual in range(residual):
                    second = second_shared and (
                        lifted[0] * k_residual
                        + lifted[1] * l_residual
                    ) % residual == 0
                    covered += int(first or second)
    return Fraction(covered, total)


class FactorizedPairExtensionTests(unittest.TestCase):
    def test_jointly_surjective_formula_matches_brute_force(self):
        cases = (
            (5, 3, (1, 0), (0, 1)),
            (7, 4, (2, 1), (1, 3)),
            (11, 3, (3, 2), (4, 1)),
        )
        for shared, residual, factor, lifted in cases:
            determinant = (
                factor[0] * lifted[1] - factor[1] * lifted[0]
            )
            self.assertEqual(math.gcd(determinant, shared), 1)
            self.assertEqual(
                brute_union(shared, residual, factor, lifted),
                component_pair_union(shared, residual),
            )


if __name__ == "__main__":
    unittest.main()
