import unittest

from derive_phase_map import derive_phase, parse_update


class DerivePhaseMapTests(unittest.TestCase):
    def test_applies_updates_without_mutating_base(self):
        base = {"17": 1, "19": 2}
        result = derive_phase(base, [(19, 7)])
        self.assertEqual(result, {"17": 1, "19": 7})
        self.assertEqual(base, {"17": 1, "19": 2})

    def test_rejects_missing_and_duplicate_primes(self):
        with self.assertRaises(KeyError):
            derive_phase({"17": 1}, [(19, 2)])
        with self.assertRaises(ValueError):
            derive_phase({"17": 1}, [(17, 2), (17, 3)])

    def test_parses_update(self):
        self.assertEqual(parse_update("109=6"), (109, 6))


if __name__ == "__main__":
    unittest.main()
