import itertools
import math
import random
import unittest

from certify_block_pair_anchor_star import triple_image_index_minors
from verify_block_pair_anchor_star import triple_index_explicit


def brute_image_index(rows):
    period = math.lcm(*(row["h"] for row in rows))
    image = {
        tuple(
            (row["a"] * k + row["b"] * l) % row["h"]
            for row in rows
        )
        for k in range(period)
        for l in range(period)
    }
    return math.prod(row["h"] for row in rows) // len(image)


class BlockPairAnchorStarTests(unittest.TestCase):
    def test_three_target_indices_match_brute_force(self):
        rng = random.Random(20260726)
        moduli = (4, 6, 8, 9, 10, 12)
        for _ in range(40):
            rows = tuple(
                {
                    "h": h,
                    "a": rng.randrange(h),
                    "b": rng.randrange(h),
                }
                for h in rng.sample(moduli, 3)
            )
            expected = brute_image_index(rows)
            self.assertEqual(triple_image_index_minors(rows), expected)
            self.assertEqual(triple_index_explicit(rows), expected)


if __name__ == "__main__":
    unittest.main()
