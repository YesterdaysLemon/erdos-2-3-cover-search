from fractions import Fraction

import pytest

import build_axis_layered_pool as layered


def row(
    prime: int,
    h: int,
    a: int,
    b: int,
    ord2: int,
    ord3: int,
) -> dict:
    return {
        "p": prime,
        "h": h,
        "a": a,
        "b": b,
        "ord2": ord2,
        "ord3": ord3,
        "c": 0,
    }


def test_k_layer_data_and_materialized_target() -> None:
    source = row(13, 6, 2, 3, 3, 2)
    assert layered.layer_data(source, "k") == (3, 3, 2, 2)
    materialized = layered.materialize_rows([source], "k", [2])
    assert materialized[0]["target_modulus"] == 3
    assert materialized[0]["target_residue"] == 1
    assert materialized[0]["layer_active_class"] == 2


def test_l_layer_data_is_transposed() -> None:
    source = row(13, 6, 2, 3, 3, 2)
    assert layered.layer_data(source, "l") == (2, 2, 3, 3)


def test_unimodular_basis_transform() -> None:
    source = row(13, 6, 2, 3, 3, 2)
    transformed = layered.transform_basis_row(
        source,
        (3, 1),
        (-1, 0),
    )
    assert transformed["a"] == 3
    assert transformed["b"] == 4
    assert transformed["ord2"] == 2
    assert transformed["ord3"] == 3
    assert transformed["source_a"] == 2
    assert transformed["source_b"] == 3
    assert layered.layer_data(transformed, "k") == (2, 2, 3, 3)


def test_nonunimodular_basis_is_rejected() -> None:
    source = row(13, 6, 2, 3, 3, 2)
    with pytest.raises(ValueError, match="determinant"):
        layered.transform_basis_row(source, (2, 0), (0, 1))


def test_exact_capacity_replay_accepts_complete_columns() -> None:
    rows = [
        row(5, 2, 1, 0, 2, 1),
        row(7, 2, 1, 0, 2, 1),
    ]
    replay = layered.exact_capacity_replay(rows, "k", 2, [0, 1])
    assert Fraction(
        replay["minimum_numerator"],
        replay["minimum_denominator"],
    ) == 1
    assert replay["minimum_multiplicity"] == 2


def test_exact_capacity_replay_rejects_gap() -> None:
    rows = [row(5, 2, 1, 0, 2, 1)]
    with pytest.raises(RuntimeError, match="minimum=0"):
        layered.exact_capacity_replay(rows, "k", 2, [0])


def test_select_rows_uses_coordinate_order() -> None:
    payload = {
        "choices": [
            row(5, 2, 1, 0, 2, 1),
            row(7, 3, 0, 1, 1, 3),
        ]
    }
    assert [item["p"] for item in layered.select_rows(payload, "k", 2)] == [
        5,
        7,
    ]
    assert [item["p"] for item in layered.select_rows(payload, "l", 2)] == [
        5,
    ]


def test_column_relaxation_keeps_only_active_rows() -> None:
    rows = layered.materialize_rows(
        [
            row(5, 2, 1, 0, 2, 1),
            row(7, 2, 1, 0, 2, 1),
        ],
        "k",
        [0, 1],
    )
    relaxed = layered.materialize_column_relaxation(rows, "k", 0)
    assert [item["p"] for item in relaxed] == [5]
    assert relaxed[0]["h"] == 1
    assert relaxed[0]["a"] == 0
    assert relaxed[0]["b"] == 1


def test_prime_deficit_pruning_removes_rare_residual_factors() -> None:
    rows = [
        row(5, 2, 0, 1, 1, 2),
        row(7, 6, 0, 1, 1, 6),
        row(11, 3, 0, 1, 1, 3),
    ]
    surviving, rounds = layered.prune_prime_deficit_rows(rows, "k")
    assert [item["prime"] for item in rounds] == [3, 2]
    assert surviving == []


def test_unavoidable_coprime_pair_exceeds_capacity_slack() -> None:
    rows = [
        row(5, 2, 0, 1, 1, 2),
        row(7, 3, 0, 1, 1, 3),
        row(13, 6, 0, 1, 1, 6),
    ]
    violations = layered.unavoidable_pair_violations(
        rows,
        "k",
        1,
        [0, 0, 0],
    )
    assert any(item["moduli"] == [2, 3] for item in violations)
    witness = next(item for item in violations if item["moduli"] == [2, 3])
    assert Fraction(
        witness["capacity_numerator"],
        witness["capacity_denominator"],
    ) == 1
    assert Fraction(
        witness["overlap_numerator"],
        witness["overlap_denominator"],
    ) == Fraction(1, 6)


def test_non_coprime_classes_have_no_forced_pair_overlap() -> None:
    rows = [
        row(5, 2, 0, 1, 1, 2),
        row(7, 2, 0, 1, 1, 2),
    ]
    assert not layered.unavoidable_pair_violations(
        rows,
        "k",
        1,
        [0, 0],
    )
