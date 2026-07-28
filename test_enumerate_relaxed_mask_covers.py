import unittest

from enumerate_relaxed_mask_covers import enumerate_exact_radius_covers


class EnumerateRelaxedMaskCoversTests(unittest.TestCase):
    def test_complete_enumeration_and_backbone(self):
        result = enumerate_exact_radius_covers(
            [0b001, 0b010, 0b100, 0b110],
            point_count=3,
            radius=2,
            solver_name="cadical195",
            model_limit=0,
            time_limit=0.0,
        )
        self.assertTrue(result["complete"])
        self.assertEqual(result["cover_count"], 1)
        self.assertEqual(result["observed_common_masks"], [0b001, 0b110])

    def test_model_limit_is_not_complete(self):
        result = enumerate_exact_radius_covers(
            [0b001, 0b010, 0b100, 0b011, 0b110],
            point_count=3,
            radius=2,
            solver_name="cadical195",
            model_limit=1,
            time_limit=0.0,
        )
        self.assertFalse(result["complete"])
        self.assertEqual(result["stop_reason"], "MODEL_LIMIT")
        self.assertEqual(result["cover_count"], 1)

    def test_models_are_saved_only_when_requested(self):
        kwargs = {
            "masks": [0b001, 0b010, 0b100, 0b110],
            "point_count": 3,
            "radius": 2,
            "solver_name": "cadical195",
            "model_limit": 0,
            "time_limit": 0.0,
        }
        without_models = enumerate_exact_radius_covers(**kwargs)
        with_models = enumerate_exact_radius_covers(
            **kwargs,
            save_models=True,
        )
        self.assertNotIn("models", without_models)
        self.assertEqual(with_models["models"], [[0b001, 0b110]])


if __name__ == "__main__":
    unittest.main()
