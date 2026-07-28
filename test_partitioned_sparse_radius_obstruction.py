import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from compose_partitioned_sparse_radius_obstruction import compose_manifest
from verify_partitioned_sparse_radius_obstruction import verify_manifest


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PartitionedSparseRadiusObstructionTest(unittest.TestCase):
    def test_two_phase_partition_replays_both_leaves(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pool_path = root / "pool.json"
            base_path = root / "base.json"
            phase_one_path = root / "phase-one.json"
            pool_path.write_text(
                json.dumps(
                    {
                        "choices": [
                            {
                                "p": 5,
                                "h": 2,
                                "a": 1,
                                "b": 0,
                                "target_modulus": 1,
                                "target_residue": 0,
                            }
                        ]
                    }
                )
            )
            base_path.write_text(json.dumps({"5": 0}))
            phase_one_path.write_text(json.dumps({"5": 1}))
            pool_sha = file_hash(pool_path)

            def write_leaf(
                name: str,
                phase_path: Path,
                radius: int,
                point: list[int],
            ) -> Path:
                path = root / name
                path.write_text(
                    json.dumps(
                        {
                            "complete": True,
                            "pool": str(pool_path),
                            "pool_sha256": pool_sha,
                            "base_phase": str(phase_path),
                            "base_phase_sha256": file_hash(phase_path),
                            "points": [point],
                            "fixed_primes": [5],
                            "max_changes": radius,
                            "discovery": {
                                "status": "UNSAT",
                                "complete_negative": True,
                                "initial_misses": 1,
                                "gain_mask_count": 0,
                                "gain_move_count": 0,
                                "minimum_relaxed_mask_count": None,
                            },
                        }
                    )
                )
                return path

            leaf_zero = write_leaf(
                "leaf-zero.json",
                base_path,
                1,
                [1, 0],
            )
            leaf_one = write_leaf(
                "leaf-one.json",
                phase_one_path,
                0,
                [0, 0],
            )
            tree_path = root / "tree.json"
            tree_path.write_text(
                json.dumps(
                    {
                        "prime": 5,
                        "branches": {
                            "0": {"certificate": str(leaf_zero)},
                            "1": {"certificate": str(leaf_one)},
                        },
                    }
                )
            )
            manifest = compose_manifest(
                pool_path,
                base_path,
                tree_path,
                set(),
                1,
            )
            manifest_path = root / "manifest.json"
            manifest_path.write_text(json.dumps(manifest))
            result = verify_manifest(manifest_path)
            self.assertTrue(result["verified"])
            self.assertEqual(result["partition_count"], 1)
            self.assertEqual(result["leaf_count"], 2)
            self.assertEqual(result["total_leaf_point_count"], 2)


if __name__ == "__main__":
    unittest.main()
