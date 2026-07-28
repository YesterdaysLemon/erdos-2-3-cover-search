#!/usr/bin/env python3
"""Independent verifier tests for sparse-radius obstruction certificates."""

from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from verify_sparse_radius_obstruction import verify_certificate


class SparseRadiusObstructionVerifierTests(unittest.TestCase):
    def test_one_row_phase_swap_is_an_exact_radius_one_obstruction(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pool = root / "pool.json"
            phase = root / "phase.json"
            certificate = root / "certificate.json"
            pool.write_text(
                json.dumps(
                    {
                        "choices": [
                            {
                                "h": 2,
                                "p": 5,
                                "a": 1,
                                "b": 0,
                                "ord2": 2,
                                "ord3": 1,
                                "target_modulus": 1,
                                "target_residue": 0,
                            }
                        ]
                    }
                )
                + "\n"
            )
            phase.write_text('{"5": 0}\n')

            def digest(path: Path) -> str:
                return hashlib.sha256(path.read_bytes()).hexdigest()

            certificate.write_text(
                json.dumps(
                    {
                        "complete": True,
                        "pool": str(pool),
                        "pool_sha256": digest(pool),
                        "base_phase": str(phase),
                        "base_phase_sha256": digest(phase),
                        "points": [[0, 0], [1, 0]],
                        "fixed_primes": [],
                        "max_changes": 1,
                        "discovery": {
                            "status": "UNSAT",
                            "complete_negative": True,
                            "initial_misses": 1,
                            "gain_mask_count": 1,
                            "gain_move_count": 1,
                            "minimum_relaxed_mask_count": 1,
                        },
                    }
                )
                + "\n"
            )
            result = verify_certificate(certificate)
            self.assertTrue(result["verified"])
            self.assertFalse(result["repair_exists"])


if __name__ == "__main__":
    unittest.main()
