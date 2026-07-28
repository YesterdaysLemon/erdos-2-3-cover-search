import json
from fractions import Fraction

from scan_reverse_group_quotients import scan


def row(prime: int, h: int) -> dict:
    return {
        "p": prime,
        "h": h,
        "a": 1,
        "b": 0,
        "ord2": h,
        "ord3": 1,
        "c": 0,
        "target_modulus": 1,
        "target_residue": 0,
    }


def test_exact_scan_ranks_without_enumerating_group_cells(tmp_path):
    payload = {
        "choices": [
            row(5, 2),
            row(7, 2),
            row(13, 3),
        ]
    }
    source = tmp_path / "pool.json"
    source.write_text(json.dumps(payload))
    result = scan(
        payload,
        source,
        [("id", (1, 0), (0, 1))],
        [1, 2, 3, 6],
        cell_limit=36,
        minimum_density=Fraction(0),
    )
    record = next(
        record
        for record in result["records"]
        if record["group"]["width"] == 6
        and record["group"]["height"] == 1
    )
    assert record["descending_row_count"] == 3
    assert record["raw_density_numerator"] == 4
    assert record["raw_density_denominator"] == 3
    assert record["forced_overlap_density_numerator"] == 1
    assert record["forced_overlap_density_denominator"] == 3
    assert record["phase_independent_union_upper_bound_numerator"] == 1
    assert record["phase_independent_union_upper_bound_denominator"] == 1
    assert record["first_moment_cover_exists"] is False
    assert result["claim"]["integer_m_found"] is False
