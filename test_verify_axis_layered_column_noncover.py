from fractions import Fraction

import pytest

import certify_axis_layered_column_noncover as certify
import verify_axis_layered_column_noncover as verify


def source_payload() -> dict:
    def row(prime: int, modulus: int) -> dict:
        return {
            "p": prime,
            "h": modulus,
            "a": 0,
            "b": 1,
            "target_modulus": 1,
            "target_residue": 0,
            "layer_active_class": 0,
        }

    return {
        "schema": "axis_layered_pool_v1",
        "layer_axis": "k",
        "layer_period": 1,
        "coordinate_basis": None,
        "choices": [
            row(5, 2),
            row(7, 3),
            row(11, 6),
        ],
    }


def certificate_payload() -> dict:
    surviving = Fraction(1, 2)
    total = Fraction(1, 1)
    return {
        "schema": "axis_layered_column_prime_deficit_obstruction_v1",
        "layer_axis": "k",
        "coordinate_basis": None,
        "layer_period": 1,
        "column": 0,
        "active_row_count": 3,
        "total_reciprocal_density": certify.fraction_record(total),
        "obstruction_prime": 3,
        "divisible_row_count": 2,
        "divisible_rows": [
            {"p": 7, "residual_modulus": 3},
            {"p": 11, "residual_modulus": 6},
        ],
        "surviving_row_count": 1,
        "surviving_reciprocal_density": certify.fraction_record(surviving),
        "density_deficit": certify.fraction_record(1 - surviving),
        "proved_no_independent_column_cover": True,
        "proved_no_declared_layered_cover": True,
        "scope": "toy",
    }


def test_independent_replay_accepts_exact_obstruction() -> None:
    report = verify.replay(source_payload(), certificate_payload())
    assert report["verified"] is True
    assert report["prime_deficit_holds"] is True
    assert report["surviving_density_below_one"] is True


def test_independent_replay_rejects_bad_density() -> None:
    certificate = certificate_payload()
    certificate["surviving_reciprocal_density"]["numerator"] = 2
    with pytest.raises(
        AssertionError,
        match="stored decimal does not replay|surviving density mismatch",
    ):
        verify.replay(source_payload(), certificate)
