#!/usr/bin/env python3
"""Regression tests for bounded homogeneous non-cover witnesses."""

from __future__ import annotations

import unittest

import certify_bounded_homogeneous_noncover as certify
import verify_bounded_homogeneous_noncover as verify


ROWS = [
    {"h": 6, "p": 7, "a": 2, "b": 1},
    {"h": 10, "p": 11, "a": 1, "b": 8},
    {"h": 30, "p": 31, "a": 24, "b": 1},
]


class BoundedHomogeneousNoncoverTests(unittest.TestCase):
    def test_signature_and_modular_membership_agree(self) -> None:
        for k in range(30):
            for l in range(30):
                signature_primes = {
                    prime
                    for _modulus, prime in certify.signature_hits(
                        ROWS, k, l
                    )
                }
                self.assertEqual(
                    signature_primes,
                    set(verify.modular_hits(ROWS, k, l)),
                )

    def test_deterministic_search_finds_a_genuine_miss(self) -> None:
        witness = certify.find_witness(
            ROWS,
            seed=203,
            bits=12,
            attempts=1_000,
        )
        self.assertIsNotNone(witness)
        assert witness is not None
        k, l, attempt = witness
        self.assertGreaterEqual(attempt, 1)
        self.assertFalse(certify.signature_hits(ROWS, k, l))
        self.assertFalse(verify.modular_hits(ROWS, k, l))

    def test_origin_is_covered_by_every_homogeneous_fibre(self) -> None:
        self.assertEqual(len(certify.signature_hits(ROWS, 0, 0)), len(ROWS))
        self.assertEqual(len(verify.modular_hits(ROWS, 0, 0)), len(ROWS))


if __name__ == "__main__":
    unittest.main()
