import json
import tempfile
import unittest
from pathlib import Path

from promote_projected_conditional_designs import (
    certificate_paths,
    collect_best_designs,
)


def fraction(numerator: int, denominator: int = 1) -> dict:
    return {"numerator": numerator, "denominator": denominator}


def payload(prime: int, intersection: tuple[int, int]) -> dict:
    numerator, denominator = intersection
    return {
        "schema": "projected_conditional_design_search_v1",
        "status": "discovery_only_requires_certificate_and_replay",
        "results": [
            {
                "outside_prime": prime,
                "best_design": {
                    "intersection": fraction(numerator, denominator),
                    "improvement": fraction(1, denominator),
                },
            }
        ],
    }


class PromoteProjectedDesignTests(unittest.TestCase):
    def test_strongest_exact_intersection_wins(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first.json"
            second = Path(directory) / "second.json"
            first.write_text(json.dumps(payload(101, (3, 10))))
            second.write_text(json.dumps(payload(101, (2, 5))))
            selected = collect_best_designs([first, second])
            self.assertEqual(
                selected[101][0]["intersection"],
                fraction(2, 5),
            )

    def test_known_invalid_discovery_is_excluded(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            invalid = (
                Path(directory)
                / "_tmp_period3139_24block_projected_design_search_"
                "paired_top13.json"
            )
            invalid.write_text(json.dumps(payload(101, (9, 10))))
            self.assertEqual(collect_best_designs([invalid]), {})

    def test_certificate_names_are_deterministic(self) -> None:
        certificate, verification = certificate_paths("frontier", 101)
        self.assertEqual(
            certificate,
            Path("frontier_conditional_fibre101_autodesign_certificate.json"),
        )
        self.assertEqual(
            verification,
            Path("frontier_conditional_fibre101_autodesign_verification.json"),
        )


if __name__ == "__main__":
    unittest.main()
