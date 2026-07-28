import json

import pytest

from inventory_reverse_group_lines import (
    axis_rectangle_moduli,
    canonical_projective_direction,
    inventory,
)
from verify_reverse_group_density_obstruction import verify


def row(prime: int, h: int, a: int, b: int) -> dict:
    return {
        "p": prime,
        "h": h,
        "a": a,
        "b": b,
        "ord2": h,
        "ord3": h,
        "c": 0,
        "target_modulus": 1,
        "target_residue": 0,
    }


def source(tmp_path, rows):
    path = tmp_path / "pool.json"
    payload = {"choices": rows}
    path.write_text(json.dumps(payload))
    return path, payload


def test_projective_direction_identifies_unit_multiples():
    assert canonical_projective_direction(5, 1, 2) == (
        canonical_projective_direction(5, 2, 4)
    )


def test_axis_rectangle_moduli_recovers_coprime_product():
    assert axis_rectangle_moduli(6, 3, 2) == (2, 3)
    assert axis_rectangle_moduli(6, 1, 1) is None


def test_inventory_keeps_only_rows_that_descend(tmp_path):
    path, payload = source(
        tmp_path,
        [
            row(5, 2, 1, 0),
            row(7, 3, 0, 1),
            row(11, 4, 1, 1),
        ],
    )
    result = inventory(
        payload,
        path,
        width=6,
        height=6,
        direction=(1, 0),
        transverse=(0, 1),
    )
    assert result["descending_row_count"] == 2
    assert result["excluded_not_descending_count"] == 1
    assert result["raw_density_numerator"] == 5
    assert result["raw_density_denominator"] == 6
    assert result["density_at_least_one"] is False
    assert result["expected_uncovered_cell_count_numerator"] == 12
    assert result["expected_uncovered_cell_count_denominator"] == 1
    assert result["first_moment_cover_exists"] is False
    assert result["claim"]["explicit_finite_group_cover_found"] is False


def test_first_moment_can_certify_cover_existence(tmp_path):
    path, payload = source(
        tmp_path,
        [
            row(5, 2, 1, 0),
            row(7, 2, 1, 0),
        ],
    )
    result = inventory(
        payload,
        path,
        width=2,
        height=1,
        direction=(1, 0),
        transverse=(0, 1),
    )
    assert result["expected_uncovered_cell_count_numerator"] == 1
    assert result["expected_uncovered_cell_count_denominator"] == 2
    assert result["first_moment_cover_exists"] is True
    assert result["claim"][
        "finite_group_cover_exists_by_first_moment"
    ] is True


def test_inventory_rejects_nonunimodular_basis(tmp_path):
    path, payload = source(tmp_path, [row(5, 2, 1, 0)])
    with pytest.raises(ValueError, match="unimodular"):
        inventory(
            payload,
            path,
            width=6,
            height=6,
            direction=(2, 0),
            transverse=(0, 1),
        )


def test_inventory_rejects_duplicate_primes(tmp_path):
    path, payload = source(
        tmp_path,
        [
            row(5, 2, 1, 0),
            row(5, 3, 0, 1),
        ],
    )
    with pytest.raises(ValueError, match="duplicate"):
        inventory(
            payload,
            path,
            width=6,
            height=6,
            direction=(1, 0),
            transverse=(0, 1),
        )


def test_inventory_rejects_payload_not_bound_to_source_file(tmp_path):
    path, payload = source(tmp_path, [row(5, 2, 1, 0)])
    payload["choices"][0]["p"] = 7
    with pytest.raises(ValueError, match="authenticated file"):
        inventory(
            payload,
            path,
            width=2,
            height=1,
            direction=(1, 0),
            transverse=(0, 1),
        )


def test_density_overlap_verifier_certifies_surviving_family(tmp_path):
    signatures = [
        (5, 4),
        (7, 6),
        (11, 10),
        (13, 12),
        (17, 16),
        (19, 18),
        (29, 28),
        (31, 30),
        (71, 35),
        (73, 36),
        (37, 36),
        (41, 40),
        (43, 42),
        (61, 60),
        (109, 108),
        (113, 112),
        (241, 120),
        (127, 126),
        (337, 168),
        (181, 180),
        (211, 210),
        (433, 216),
        (281, 280),
        (421, 420),
        (1009, 504),
        (631, 630),
        (2521, 1260),
        (3361, 1680),
    ]
    rows = [
        {
            **row(prime, h, 1, 0),
            "ord3": 1,
        }
        for prime, h in signatures
    ]
    path, payload = source(tmp_path, rows)
    result = inventory(
        payload,
        path,
        width=15120,
        height=2520,
        direction=(1, 0),
        transverse=(0, 1),
    )
    assert result["raw_density_numerator"] == 337
    assert result["raw_density_denominator"] == 336
    assert len(result["forced_overlap_forest"]) == 9
    assert result["forced_overlap_density_numerator"] == 11
    assert result["forced_overlap_density_denominator"] == 560
    assert (
        result["phase_independent_union_upper_bound_numerator"]
        == 413
    )
    assert (
        result["phase_independent_union_upper_bound_denominator"]
        == 420
    )
    assert result[
        "finite_group_cover_impossible_by_density_overlap"
    ] is True
    report = verify(result)
    assert report["verified"] is True
    assert report["declared_finite_group_family_noncover_proved"] is True


def test_maximum_forest_uses_jointly_surjective_noncoprime_maps(tmp_path):
    path, payload = source(
        tmp_path,
        [
            row(5, 4, 1, 0),
            row(7, 4, 1, 0),
            row(11, 4, 0, 1),
            row(13, 4, 0, 1),
        ],
    )
    result = inventory(
        payload,
        path,
        width=4,
        height=4,
        direction=(1, 0),
        transverse=(0, 1),
    )
    assert result["raw_density_numerator"] == 1
    assert result["raw_density_denominator"] == 1
    assert len(result["forced_overlap_forest"]) == 3
    assert result["forced_overlap_density_numerator"] == 3
    assert result["forced_overlap_density_denominator"] == 16
    assert (
        result["phase_independent_union_upper_bound_numerator"]
        == 13
    )
    assert (
        result["phase_independent_union_upper_bound_denominator"]
        == 16
    )
    assert result[
        "finite_group_cover_impossible_by_density_overlap"
    ] is True
    report = verify(result)
    assert report["verified"] is True
    assert report["forced_overlap_forest_edge_count"] == 3


def test_density_overlap_verifier_rejects_tampered_bound(tmp_path):
    path, payload = source(
        tmp_path,
        [
            row(5, 2, 1, 0),
            row(7, 3, 1, 0),
        ],
    )
    result = inventory(
        payload,
        path,
        width=6,
        height=1,
        direction=(1, 0),
        transverse=(0, 1),
    )
    result["phase_independent_union_upper_bound_numerator"] += 1
    with pytest.raises(AssertionError, match="bound metadata"):
        verify(result)
