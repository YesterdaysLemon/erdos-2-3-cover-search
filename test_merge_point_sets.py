import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import merge_point_sets


class MergePointSetsTests(unittest.TestCase):
    def test_merges_raw_and_supported_key_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw.json"
            keyed = root / "keyed.json"
            output = root / "merged.json"
            raw.write_text(json.dumps([[1, 2], [3, 4]]))
            keyed.write_text(json.dumps({"misses": [[3, 4], [5, 6]]}))
            argv = [
                "merge_point_sets.py",
                str(raw),
                str(keyed),
                "--output",
                str(output),
            ]
            with mock.patch("sys.argv", argv):
                self.assertEqual(merge_point_sets.main(), 0)
            self.assertEqual(
                json.loads(output.read_text()),
                [[1, 2], [3, 4], [5, 6]],
            )


if __name__ == "__main__":
    unittest.main()
