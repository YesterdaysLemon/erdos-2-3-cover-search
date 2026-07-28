import unittest

import component_density_cegis


class TargetModulusTests(unittest.TestCase):
    def test_missing_target_modulus_means_unrestricted(self):
        self.assertEqual(component_density_cegis.target_modulus({"p": 5}), 1)

    def test_explicit_target_modulus_is_preserved(self):
        self.assertEqual(
            component_density_cegis.target_modulus(
                {"p": 5, "target_modulus": 7}
            ),
            7,
        )

    def test_target_restriction_splits_over_residual_component(self):
        row = {
            "h": 84,
            "p": 5,
            "target_modulus": 14,
            "target_residue": 9,
        }
        self.assertEqual(
            component_density_cegis.coarse_target_restriction(
                row,
                2,
                23,
            ),
            (2, 7),
        )
        self.assertTrue(
            component_density_cegis.coarse_target_allowed(16, (2, 7))
        )
        self.assertFalse(
            component_density_cegis.coarse_target_allowed(15, (2, 7))
        )

    def test_invalid_residual_phase_is_rejected(self):
        row = {
            "h": 84,
            "p": 5,
            "target_modulus": 14,
            "target_residue": 9,
        }
        with self.assertRaisesRegex(
            ValueError,
            "residual target restriction",
        ):
            component_density_cegis.coarse_target_restriction(
                row,
                2,
                22,
            )

    def test_algebraic_origin_reduces_the_exact_scaled_threshold(self):
        self.assertEqual(
            component_density_cegis.required_scaled_density(
                49,
                7,
                (7, 11),
            ),
            48,
        )
        self.assertEqual(
            component_density_cegis.required_scaled_density(
                343,
                7,
                (7,),
            ),
            336,
        )
        # At the first q-adic level the integer line count still rounds up
        # to q even though one affine-plane point is algebraically covered.
        self.assertEqual(
            component_density_cegis.required_scaled_density(
                7,
                7,
                (7,),
            ),
            7,
        )
        self.assertEqual(
            component_density_cegis.required_scaled_density(
                49,
                7,
                (),
            ),
            49,
        )

    def test_checker_excludes_other_algebraic_origin_cells(self):
        cells, metadata = component_density_cegis.find_low_density_cells(
            [],
            {},
            2,
            1,
            3,
            algebraic_primes=(3,),
        )
        self.assertEqual(len(cells), 1)
        k, l = cells[0]
        self.assertTrue(k % 3 or l % 3)
        self.assertEqual(metadata["components"], {3: 3})
        self.assertEqual(metadata["algebraic_primes"], [3])

    def test_residual_plane_audit_detects_density_without_cover(self):
        rows = [
            {"h": 3, "p": 7, "a": 1, "b": 0},
            {"h": 3, "p": 13, "a": 1, "b": 0},
            {"h": 3, "p": 19, "a": 1, "b": 0},
        ]
        failed = component_density_cegis.audit_residual_plane_cells(
            rows,
            {7: 0, 13: 0, 19: 1},
            3,
            [(0, 0)],
        )
        self.assertEqual(failed["failed_cells"], 1)
        self.assertEqual(failed["minimum_covered_points"], 6)
        covered = component_density_cegis.audit_residual_plane_cells(
            rows,
            {7: 0, 13: 1, 19: 2},
            3,
            [(0, 0)],
        )
        self.assertEqual(covered["failed_cells"], 0)
        self.assertEqual(covered["minimum_covered_points"], 9)

    def test_exact_density_replay_rejects_a_bad_master_phase(self):
        rows = [
            {
                "h": 6,
                "p": 5,
                "a": 1,
                "b": 0,
            },
            {
                "h": 3,
                "p": 7,
                "a": 1,
                "b": 0,
            },
        ]
        density, modulus = component_density_cegis.scaled_density_at_cell(
            rows,
            {5: 1, 7: 0},
            2,
            (1, 0),
        )
        self.assertEqual((density, modulus), (1, 2))
        with self.assertRaisesRegex(
            AssertionError,
            "failed exact residual-density replay",
        ):
            component_density_cegis.audit_density_cuts(
                rows,
                {5: 1, 7: 0},
                2,
                [(1, 0)],
            )
        replay = component_density_cegis.audit_density_cuts(
            rows,
            {5: 1, 7: 1},
            2,
            [(1, 0)],
        )
        self.assertEqual(replay["minimum_scaled_density"], 3)
        self.assertEqual(replay["violations"], 0)

    def test_all_master_engines_enforce_coarse_target_restrictions(self):
        rows = [
            {
                "h": 6,
                "p": prime,
                "a": 1,
                "b": 0,
                "target_modulus": 3,
                "target_residue": 1,
            }
            for prime in (5, 7)
        ]
        initial = {5: 1, 7: 1}
        engines = (
            lambda cells: component_density_cegis.solve_master(
                rows,
                initial,
                2,
                cells,
                {},
                "cadical195",
            ),
            lambda cells: component_density_cegis.solve_master_z3(
                rows,
                initial,
                2,
                cells,
                {},
            ),
            lambda cells: component_density_cegis.solve_master_milp(
                rows,
                initial,
                2,
                cells,
                {},
                0.0,
            ),
        )
        for solve in engines:
            with self.subTest(solve=solve):
                phases, metadata = solve([(1, 0)])
                self.assertTrue(metadata["sat"])
                self.assertTrue(
                    all(value % 3 == 1 for value in phases.values())
                )
                phases, metadata = solve([(2, 0)])
                self.assertIsNone(phases)
                self.assertFalse(metadata["sat"])

    def test_milp_can_decide_a_hard_phase_change_budget(self):
        rows = [
            {
                "h": 6,
                "p": prime,
                "a": 1,
                "b": 0,
            }
            for prime in (5, 7)
        ]
        initial = {5: 0, 7: 0}
        phases, metadata = component_density_cegis.solve_master_milp(
            rows,
            initial,
            2,
            [(1, 0)],
            {},
            0.0,
            maximum_changes=1,
        )
        self.assertIsNone(phases)
        self.assertFalse(metadata["sat"])
        phases, metadata = component_density_cegis.solve_master_milp(
            rows,
            initial,
            2,
            [(1, 0)],
            {},
            0.0,
            maximum_changes=2,
        )
        self.assertTrue(metadata["sat"])
        self.assertTrue(all(target % 3 == 1 for target in phases.values()))
        self.assertTrue(all(target % 2 == 0 for target in phases.values()))
        replay = component_density_cegis.audit_density_cuts(
            rows,
            phases,
            2,
            [(1, 0)],
        )
        self.assertEqual(replay["minimum_scaled_density"], 2)

    def test_exact_one_change_scan_finds_and_replays_a_winner(self):
        rows = [
            {
                "h": 6,
                "p": 5,
                "a": 1,
                "b": 0,
            },
            {
                "h": 3,
                "p": 7,
                "a": 1,
                "b": 0,
            },
        ]
        winners, metadata = component_density_cegis.scan_one_change(
            rows,
            {5: 0, 7: 0},
            2,
            [(1, 0)],
            {},
        )
        self.assertTrue(metadata["sat"])
        self.assertEqual(metadata["winner_count"], 1)
        phases, winner = winners[0]
        self.assertEqual(winner["prime"], 7)
        self.assertEqual(phases, {5: 0, 7: 1})
        self.assertEqual(
            winner["exact_replay"]["minimum_scaled_density"],
            2,
        )

    def test_exact_two_change_scan_finds_a_required_pair(self):
        rows = [
            {"h": 3, "p": 7, "a": 1, "b": 0},
            {"h": 3, "p": 13, "a": 1, "b": 0},
            {"h": 3, "p": 19, "a": 1, "b": 0},
        ]
        phases = {7: 0, 13: 0, 19: 0}
        winners, metadata = component_density_cegis.scan_two_changes(
            rows,
            phases,
            5,
            [(0, 0), (1, 0), (2, 0)],
            {7: 0},
        )
        self.assertTrue(metadata["sat"])
        self.assertEqual(len(winners[0][1]["changes"]), 2)
        self.assertEqual(
            {
                change["new_coarse_target"]
                for change in winners[0][1]["changes"]
            },
            {1, 2},
        )
        self.assertEqual(
            winners[0][1]["exact_replay"]["violations"],
            0,
        )

    def test_exact_two_change_scan_can_prove_no_pair(self):
        rows = [
            {"h": 4, "p": 7, "a": 1, "b": 0},
            {"h": 4, "p": 13, "a": 1, "b": 0},
            {"h": 4, "p": 19, "a": 1, "b": 0},
        ]
        phases = {7: 0, 13: 0, 19: 0}
        winners, metadata = component_density_cegis.scan_two_changes(
            rows,
            phases,
            5,
            [(0, 0), (1, 0), (2, 0), (3, 0)],
            {7: 0},
        )
        self.assertFalse(metadata["sat"])
        self.assertEqual(winners, [])


if __name__ == "__main__":
    unittest.main()
