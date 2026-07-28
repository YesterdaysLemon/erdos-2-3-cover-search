import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class RandomDerivedCoverageAuditTests(unittest.TestCase):
    def test_writes_captured_uncovered_points(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pool = root / "pool.json"
            phases = root / "phases.json"
            result = root / "result.json"
            misses = root / "misses.json"
            pool.write_text(
                '{"choices":[{"p":5,"h":2,"a":1,"b":0}],'
                '"algebraic_primes":[]}\n'
            )
            phases.write_text('{"5":0}\n')
            completed = subprocess.run(
                [
                    sys.executable,
                    "random_derived_coverage_audit.py",
                    str(pool),
                    str(phases),
                    "--draws",
                    "32",
                    "--batch",
                    "16",
                    "--seed",
                    "7",
                    "--output",
                    str(result),
                    "--uncovered-output",
                    str(misses),
                    "--uncovered-limit",
                    "3",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 1)
            import json

            result_payload = json.loads(result.read_text())
            miss_payload = json.loads(misses.read_text())
            self.assertGreater(result_payload["uncovered"], 0)
            self.assertEqual(result_payload["captured_uncovered"], 3)
            self.assertEqual(len(miss_payload["points"]), 3)
            self.assertTrue(
                all(k % 2 == 1 for k, _ell in miss_payload["points"])
            )


if __name__ == "__main__":
    unittest.main()
