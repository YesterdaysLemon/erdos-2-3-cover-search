import copy

import pytest

import certify_axis_layered_family_noncover as certify
import verify_axis_layered_family_noncover as verify


def toy_source() -> dict:
    return {
        "schema": "axis_layered_pool_v1",
        "layer_axis": "k",
        "layer_period": 2,
        "coordinate_basis": None,
        "choices": [
            {"p": 5, "h": 2, "a": 0, "b": 1},
            {"p": 13, "h": 4, "a": 1, "b": 2},
        ],
    }


def toy_certificate() -> dict:
    source = toy_source()
    rows = certify.layer_rows(source)
    surviving, rounds = certify.prune_prime_deficits(rows)
    mdd = certify.exact_mdd_capacity_obstruction(surviving, [0, 1])
    density = certify.fraction_record(
        sum(
            (
                certify.Fraction(1, row["h"])
                for row in surviving
            ),
            certify.Fraction(),
        )
    )
    return {
        "schema": "axis_layered_family_prime_pruned_mdd_obstruction_v1",
        "layer_axis": "k",
        "coordinate_basis": None,
        "layer_period": 2,
        "initial_row_count": 2,
        "prime_deficit_pruning": rounds,
        "surviving_row_count": 2,
        "surviving_rows": surviving,
        "surviving_raw_reciprocal_density": density,
        "capacity_obstruction": mdd,
        "proved_no_active_class_placement": True,
        "proved_no_declared_layered_family_cover": True,
        "scope": "toy",
    }


def test_independent_family_replay_accepts_toy_obstruction() -> None:
    report = verify.replay(toy_source(), toy_certificate())
    assert report["verified"] is True
    assert report["terminal_root"] == 0


def test_independent_family_replay_rejects_tampered_core() -> None:
    certificate = copy.deepcopy(toy_certificate())
    certificate["capacity_obstruction"]["core_columns"] = [0]
    with pytest.raises(
        (AssertionError, RuntimeError),
        match="MDD replay mismatch|capacity infeasibility",
    ):
        verify.replay(toy_source(), certificate)
