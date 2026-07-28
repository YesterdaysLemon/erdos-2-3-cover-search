import copy

import pytest

import verify_axis_layered_pool as verifier


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


def complete_artifact() -> tuple[dict, dict]:
    source = {
        "choices": [
            row(5, 2, 1, 0, 2, 1),
            row(7, 2, 1, 0, 2, 1),
        ]
    }
    artifact = {
        "schema": "axis_layered_pool_v1",
        "layer_axis": "k",
        "coordinate_basis": None,
        "layer_period": 2,
        "row_count": 2,
        "choices": [
            {
                **source["choices"][0],
                "target_modulus": 2,
                "target_residue": 0,
                "layer_active_class": 0,
            },
            {
                **source["choices"][1],
                "target_modulus": 2,
                "target_residue": 1,
                "layer_active_class": 1,
            },
        ],
        "exact_capacity_replay": {
            "column_count": 2,
            "minimum_numerator": 1,
            "minimum_denominator": 1,
            "minimum_multiplicity": 2,
        },
        "pair_overlap_screen": {
            "complete": True,
            "pair_safe": True,
            "remaining_violations": 0,
        },
    }
    return source, artifact


def test_complete_axis_artifact_verifies() -> None:
    source, artifact = complete_artifact()
    result = verifier.verify_payload(source, artifact)
    assert result["verified"]
    assert result["unavoidable_pair_violations"] == 0


def test_target_tamper_is_rejected() -> None:
    source, artifact = complete_artifact()
    tampered = copy.deepcopy(artifact)
    tampered["choices"][1]["target_residue"] = 0
    with pytest.raises(ValueError, match="target residue"):
        verifier.verify_payload(source, tampered)


def test_sheared_signature_is_reconstructed() -> None:
    source = row(13, 6, 2, 3, 3, 2)
    transformed = verifier.transformed_signature(
        source,
        {
            "direction": [3, 1],
            "transverse": [-1, 0],
            "determinant": 1,
        },
    )
    assert transformed["a"] == 3
    assert transformed["b"] == 4
    assert transformed["ord2"] == 2
    assert transformed["ord3"] == 3


def test_prime_deficit_pruning_is_replayed() -> None:
    rows = [
        row(5, 2, 0, 1, 1, 2),
        row(7, 6, 0, 1, 1, 6),
        row(11, 3, 0, 1, 1, 3),
    ]
    surviving, rounds = verifier.replay_prime_deficit_pruning(rows, "k")
    assert [item["prime"] for item in rounds] == [3, 2]
    assert surviving == []
