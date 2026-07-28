import json

from verify_square_quotient_period_bridge import verify


def fraction(value):
    return {
        "numerator": value[0],
        "denominator": value[1],
    }


def test_square_descent_and_forced_pair_survive_shear(tmp_path):
    rows = [
        {"p": 5, "h": 4, "a": 1, "b": 3},
        {"p": 7, "h": 6, "a": 2, "b": 1},
        {"p": 11, "h": 10, "a": 1, "b": 8},
    ]
    source_path = tmp_path / "source.json"
    source_path.write_text(json.dumps({"choices": rows}))
    period_pool_path = tmp_path / "period.json"
    period_pool_path.write_text(
        json.dumps(
            {
                "source": str(source_path),
                "period_filter": 12,
                "choices": rows[:2],
            }
        )
    )
    certificate_path = tmp_path / "certificate.json"
    certificate_path.write_text(
        json.dumps(
            {
                "pool": str(period_pool_path),
                "row_count": 2,
                "anchor_primes": [5, 7],
                "anchor_rows": rows[:2],
                "joint_image_index": 1,
                "joint_target_map_surjective": True,
                "total_reciprocal_density": fraction((5, 12)),
                "forced_pair_overlap_density": fraction((1, 24)),
                "union_density_upper_bound": fraction((3, 8)),
                "proved_no_cover": True,
            }
        )
    )
    report = verify(
        source_path,
        period_pool_path,
        certificate_path,
        period=12,
        direction=(1, 1),
        transverse=(0, 1),
    )
    assert report["verified"] is True
    assert report[
        "descending_primes_equal_divisor_period_primes"
    ] is True
    assert report["joint_image_index"] == 1
    assert report["union_upper_bound_numerator"] == 3
    assert report["union_upper_bound_denominator"] == 8
    assert report["proved_no_declared_square_quotient_cover"] is True
