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


if __name__ == "__main__":
    unittest.main()
