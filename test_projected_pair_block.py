import unittest
from fractions import Fraction

try:
    import numpy as np
except ModuleNotFoundError:
    np = None

from certify_projected_pair_block import (
    pair_union_density,
    projected_pair_density,
    residual_chain_densities,
)
from certify_projected_triple_independent_block import (
    factorized_projection_scores,
    main_chain_densities,
    main_pair_densities,
    residual_union_weights,
    residual_union_weights_with_first_projected_pair,
    split_three_events,
)
from certify_projected_endpoint_path_block import (
    endpoint_path_weights,
    endpoint_path_with_start_leaf_density,
    linear_equation_intersection_density,
)
from verify_projected_endpoint_path_block import (
    brute_endpoint_path_with_start_leaf_density,
)


class ProjectedPairBlockTests(unittest.TestCase):
    def test_affine_equation_rank_density(self):
        self.assertEqual(
            linear_equation_intersection_density(
                5,
                ((1, 0, 0), (0, 1, 0)),
            ),
            Fraction(1, 25),
        )
        self.assertEqual(
            linear_equation_intersection_density(
                5,
                ((1, 0, 0), (0, 1, 0), (1, 1, 1)),
            ),
            0,
        )

    def test_normalized_start_leaf_matches_brute_component_replay(self):
        path_rows = (
            {"a": 1, "b": 5},
            {"a": 33, "b": 37},
            {"a": 79, "b": 183},
            {"a": 95, "b": 178},
            {"a": 1, "b": 50},
        )
        leaf = {"a": 3, "b": 13}
        residuals = (7, 17, 19, 29)
        for compatible in (False, True):
            for main_active in (False, True):
                for endpoint_active in (False, True):
                    for leaf_active in (False, True):
                        arguments = {
                            "projected_enabled": main_active,
                            "endpoint_enabled": endpoint_active,
                            "start_leaf_enabled": leaf_active,
                            "compatible": compatible,
                        }
                        self.assertEqual(
                            endpoint_path_with_start_leaf_density(
                                path_rows,
                                leaf,
                                residuals,
                                **arguments,
                            ),
                            brute_endpoint_path_with_start_leaf_density(
                                path_rows,
                                leaf,
                                residuals,
                                **arguments,
                            ),
                        )

    def test_pair_component_formula(self):
        self.assertEqual(
            pair_union_density(7, 17),
            Fraction(1, 7) + Fraction(1, 119) - Fraction(1, 833),
        )

    def test_base_cell_decomposition(self):
        cells = 10
        covered = 4
        compatible = 3
        expected = (
            Fraction(covered, cells)
            + Fraction(compatible, cells)
            * (
                Fraction(1, 7)
                + Fraction(1, 119)
                - Fraction(1, 833)
            )
            + Fraction(3, cells * 119)
        )
        self.assertEqual(
            projected_pair_density(
                cells,
                covered,
                compatible,
                7,
                17,
            ),
            expected,
        )

    def test_three_row_residual_chain_inclusion_exclusion(self):
        active, inactive = residual_chain_densities(7, 17, 19)
        expected_active = (
            Fraction(1, 7)
            + Fraction(1, 7 * 17)
            + Fraction(1, 17 * 19)
            - Fraction(1, 7 * 7 * 17)
            - Fraction(1, 7 * 17 * 19)
            - Fraction(1, 7 * 17 * 17 * 19)
            + Fraction(1, 7 * 7 * 17 * 17 * 19)
        )
        expected_inactive = (
            Fraction(1, 7 * 17)
            + Fraction(1, 17 * 19)
            - Fraction(1, 7 * 17 * 17 * 19)
        )
        self.assertEqual(active, expected_active)
        self.assertEqual(inactive, expected_inactive)

    def test_four_event_path_has_product_union_density(self):
        moduli = (2, 3, 5, 7)
        covered_all = 0
        covered_without_first = 0
        total = 1
        for modulus in moduli:
            total *= modulus * modulus
        for x0 in range(2):
            for y0 in range(2):
                for x1 in range(3):
                    for y1 in range(3):
                        for x2 in range(5):
                            for y2 in range(5):
                                for x3 in range(7):
                                    for y3 in range(7):
                                        events = (
                                            x0 == 0,
                                            y0 == 0 and x1 == 0,
                                            y1 == 0 and x2 == 0,
                                            y2 == 0 and x3 == 0,
                                        )
                                        covered_all += int(any(events))
                                        covered_without_first += int(
                                            any(events[1:])
                                        )
        densities = (
            Fraction(1, 2),
            Fraction(1, 2 * 3),
            Fraction(1, 3 * 5),
            Fraction(1, 5 * 7),
        )
        expected_all_uncovered = Fraction(1)
        expected_tail_uncovered = Fraction(1)
        for density in densities:
            expected_all_uncovered *= 1 - density
        for density in densities[1:]:
            expected_tail_uncovered *= 1 - density
        self.assertEqual(
            Fraction(covered_all, total),
            1 - expected_all_uncovered,
        )
        self.assertEqual(
            Fraction(covered_without_first, total),
            1 - expected_tail_uncovered,
        )

    def test_five_event_path_with_both_endpoints(self):
        moduli = (2, 3, 5, 7)
        densities = (
            Fraction(1, 2),
            Fraction(1, 2 * 3),
            Fraction(1, 3 * 5),
            Fraction(1, 5 * 7),
            Fraction(1, 7),
        )
        weights = endpoint_path_weights(densities, Fraction(1, 11))
        expected = Fraction(1)
        for density in (*densities, Fraction(1, 11)):
            expected *= 1 - density
        self.assertEqual(weights[(True, True, True)], 1 - expected)
        expected_inactive = Fraction(1)
        for density in densities[1:4]:
            expected_inactive *= 1 - density
        self.assertEqual(
            weights[(False, False, False)],
            1 - expected_inactive,
        )

    def test_endpoint_path_with_conditional_shared_pair(self):
        densities = (
            Fraction(1, 2),
            Fraction(1, 2 * 3),
            Fraction(1, 3 * 5),
            Fraction(1, 5 * 7),
            Fraction(1, 7),
        )
        pair_active = Fraction(2, 11) - Fraction(1, 121)
        pair_inactive = Fraction(1, 11)
        weights = endpoint_path_weights(
            densities,
            Fraction(1, 13),
            pair_active,
            pair_inactive,
        )
        inactive_uncovered = Fraction(1)
        for density in densities[1:4]:
            inactive_uncovered *= 1 - density
        self.assertEqual(
            weights[(False, False, False, False)],
            1 - inactive_uncovered * (1 - pair_inactive),
        )
        active_uncovered = Fraction(1)
        for density in (*densities, Fraction(1, 13)):
            active_uncovered *= 1 - density
        self.assertEqual(
            weights[(True, True, True, True)],
            1 - active_uncovered * (1 - pair_active),
        )

    def test_endpoint_path_with_two_shared_pairs_and_extra_residual(self):
        densities = (
            Fraction(1, 2),
            Fraction(1, 2 * 3),
            Fraction(1, 3 * 5),
            Fraction(1, 5 * 7),
            Fraction(1, 7),
        )
        paired_active = Fraction(2, 11) - Fraction(1, 121)
        paired_inactive = Fraction(1, 11)
        first_pair = Fraction(2, 13) - Fraction(1, 169)
        weights = endpoint_path_weights(
            densities,
            Fraction(1, 13),
            paired_active,
            paired_inactive,
            first_pair,
            Fraction(1, 31),
        )
        inactive_uncovered = Fraction(1)
        for density in densities[1:4]:
            inactive_uncovered *= 1 - density
        inactive_uncovered *= 1 - paired_inactive
        self.assertEqual(
            weights[(False, False, False, False, False, False)],
            1 - inactive_uncovered,
        )
        active_uncovered = Fraction(1)
        for density in densities:
            active_uncovered *= 1 - density
        active_uncovered *= 1 - first_pair
        active_uncovered *= 1 - paired_active
        active_uncovered *= Fraction(30, 31)
        self.assertEqual(
            weights[(True, True, True, True, True, True)],
            1 - active_uncovered,
        )

    def test_triple_independent_residual_weights(self):
        active, inactive = main_pair_densities(2, 3)
        weights = residual_union_weights(active, inactive, (5, 7, 11))
        self.assertEqual(
            weights[(True, True, True, True)],
            1
            - (1 - active)
            * Fraction(4, 5)
            * Fraction(6, 7)
            * Fraction(10, 11),
        )
        self.assertEqual(
            weights[(False, False, False, False)],
            inactive,
        )

    def test_tail_and_shared_pair_event_weights(self):
        active_main, inactive_main = main_chain_densities(2, 3, 5)
        shared_pair_active = Fraction(2, 11) - Fraction(1, 121)
        shared_pair_inactive = Fraction(1, 11)
        weights = residual_union_weights(
            active_main,
            inactive_main,
            (7, 13, 11),
            (
                Fraction(1, 7),
                Fraction(1, 13),
                shared_pair_active,
            ),
            (
                Fraction(0),
                Fraction(0),
                shared_pair_inactive,
            ),
        )
        self.assertEqual(
            weights[(False, False, False, False)],
            1
            - (1 - inactive_main)
            * (1 - shared_pair_inactive),
        )
        self.assertEqual(
            weights[(True, True, True, True)],
            1
            - (1 - active_main)
            * Fraction(6, 7)
            * Fraction(12, 13)
            * (1 - shared_pair_active),
        )

    def test_shared_pair_weights_preserve_a_fourth_event(self):
        active_main, inactive_main = main_chain_densities(2, 3, 5)
        shared_pair_active = Fraction(2, 11) - Fraction(1, 121)
        shared_pair_inactive = Fraction(1, 11)
        weights = residual_union_weights(
            active_main,
            inactive_main,
            (7, 13, 11, 17),
            (
                Fraction(1, 7),
                Fraction(1, 13),
                shared_pair_active,
                Fraction(1, 17),
            ),
            (
                Fraction(0),
                Fraction(0),
                shared_pair_inactive,
                Fraction(0),
            ),
        )
        self.assertEqual(
            weights[(False, False, False, False, True)],
            1
            - (1 - inactive_main)
            * Fraction(10, 11)
            * Fraction(16, 17),
        )
        self.assertEqual(
            weights[(True, True, True, True, True)],
            1
            - (1 - active_main)
            * Fraction(6, 7)
            * Fraction(12, 13)
            * (1 - shared_pair_active)
            * Fraction(16, 17),
        )

    def test_two_projected_lines_share_first_residual(self):
        active_main, inactive_main = main_chain_densities(2, 3, 5)
        weights = residual_union_weights_with_first_projected_pair(
            active_main,
            inactive_main,
            (7, 11),
            (Fraction(1, 7), Fraction(2, 11) - Fraction(1, 121)),
            (Fraction(0), Fraction(1, 11)),
        )
        self.assertEqual(
            weights[(False, True, False, True)],
            1
            - (1 - inactive_main)
            * (1 - (Fraction(2, 7) - Fraction(1, 49)))
            * Fraction(10, 11),
        )

    def test_triple_independent_weights_match_small_brute_grid(self):
        active, inactive = main_pair_densities(2, 3)
        weights = residual_union_weights(active, inactive, (5, 7, 11))
        total = 2 * 2 * 3 * 5 * 7 * 11
        for state in (
            (main, first, second, third)
            for main in (False, True)
            for first in (False, True)
            for second in (False, True)
            for third in (False, True)
        ):
            covered = 0
            for shared_first in range(2):
                for shared_second in range(2):
                    for lifted_digit in range(3):
                        for first_digit in range(5):
                            for second_digit in range(7):
                                for third_digit in range(11):
                                    events = (
                                        state[0] and shared_first == 0,
                                        shared_second == 0
                                        and lifted_digit == 0,
                                        state[1] and first_digit == 0,
                                        state[2] and second_digit == 0,
                                        state[3] and third_digit == 0,
                                    )
                                    covered += int(any(events))
            self.assertEqual(
                Fraction(covered, total),
                weights[state],
            )

    def test_three_event_bitset_partition(self):
        universe = (1 << 12) - 1
        events = (
            sum(1 << index for index in range(12) if index % 2 == 0),
            sum(1 << index for index in range(12) if index % 3 == 0),
            sum(1 << index for index in range(12) if index % 5 == 0),
        )
        pieces = split_three_events(universe, events)
        self.assertEqual(sum(mask.bit_count() for mask in pieces.values()), 12)
        self.assertEqual(
            sum(pieces.values()),
            universe,
        )
        with self.assertRaises(ValueError):
            split_three_events(~universe, events)

    @unittest.skipIf(np is None, "NumPy is unavailable")
    def test_factorized_scores_match_direct_union_formula(self):
        core_denominator = 6
        extra_denominator = 5
        core_choice_weights = np.asarray(
            [[1, 3], [4, 2]],
            dtype=np.int64,
        )
        histogram = np.asarray(
            [[2, 1], [3, 4]],
            dtype=np.int64,
        )
        extra_cell_choice_weights = np.asarray(
            [[0, 2], [1, 4]],
            dtype=np.int64,
        )
        scores = factorized_projection_scores(
            core_choice_weights,
            histogram,
            extra_cell_choice_weights,
            core_denominator,
            extra_denominator,
        )
        direct = np.zeros((2, 2), dtype=np.int64)
        for core_choice in range(2):
            for extra_choice in range(2):
                for core_code in range(2):
                    for extra_code in range(2):
                        core_weight = core_choice_weights[
                            core_choice,
                            core_code,
                        ]
                        extra_weight = extra_cell_choice_weights[
                            extra_code,
                            extra_choice,
                        ]
                        combined = (
                            core_weight * extra_denominator
                            + (core_denominator - core_weight)
                            * extra_weight
                        )
                        direct[core_choice, extra_choice] += (
                            histogram[core_code, extra_code] * combined
                        )
        np.testing.assert_array_equal(scores, direct)


if __name__ == "__main__":
    unittest.main()
