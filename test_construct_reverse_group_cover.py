import json

import pytest

from construct_reverse_group_cover import construct, source_target
from inventory_reverse_group_lines import inventory


def row(prime: int) -> dict:
    return {
        "p": prime,
        "h": 2,
        "a": 1,
        "b": 0,
        "ord2": 2,
        "ord3": 1,
        "c": 0,
    }


def test_conditional_expectation_constructs_certified_case(tmp_path):
    pool = {"choices": [row(5), row(7)]}
    pool_path = tmp_path / "pool.json"
    pool_path.write_text(json.dumps(pool))
    line_inventory = inventory(
        pool,
        pool_path,
        width=2,
        height=1,
        direction=(1, 0),
        transverse=(0, 1),
    )
    result = construct(line_inventory, pool, pool_path)
    assert result["initial_expected_uncovered_numerator"] == 1
    assert result["initial_expected_uncovered_denominator"] == 2
    assert result["finite_group_miss_count"] == 0
    assert result["claim"]["finite_group_cover_candidate"] is True
    assert sorted(int(row["c"]) for row in result["choices"]) == [0, 1]
    assert (
        result["claim"]["independent_full_domain_verification_passed"]
        is False
    )


def test_source_target_inverts_nontrivial_projective_multiplier():
    assert source_target(
        h=5,
        source_a=2,
        source_b=4,
        canonical_a=1,
        canonical_b=2,
        canonical_target=4,
    ) == 3


def test_constructor_rejects_restricted_row_relabelled_unrestricted(
    tmp_path,
):
    restricted = row(5)
    restricted["target_modulus"] = 2
    restricted["target_residue"] = 0
    pool = {"choices": [restricted]}
    pool_path = tmp_path / "pool.json"
    pool_path.write_text(json.dumps(pool))
    line_inventory = inventory(
        pool,
        pool_path,
        width=2,
        height=1,
        direction=(1, 0),
        transverse=(0, 1),
    )
    line_inventory["line_types"][0]["target_modulus"] = 1
    line_inventory["line_types"][0]["target_residue"] = 0
    line_inventory["expected_uncovered_cell_count_numerator"] = 1
    line_inventory["expected_uncovered_cell_count_denominator"] = 1
    with pytest.raises(ValueError, match="target restriction mismatch"):
        construct(line_inventory, pool, pool_path)


def test_constructor_rejects_pool_payload_not_bound_to_file(tmp_path):
    pool = {"choices": [row(5)]}
    pool_path = tmp_path / "pool.json"
    pool_path.write_text(json.dumps(pool))
    line_inventory = inventory(
        pool,
        pool_path,
        width=2,
        height=1,
        direction=(1, 0),
        transverse=(0, 1),
    )
    tampered_pool = {"choices": [row(7)]}
    with pytest.raises(ValueError, match="authenticated file"):
        construct(line_inventory, tampered_pool, pool_path)
