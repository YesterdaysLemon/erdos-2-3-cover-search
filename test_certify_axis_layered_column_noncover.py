from fractions import Fraction

import certify_axis_layered_column_noncover as certify


def toy_payload() -> dict:
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
        "choices": [
            row(5, 2),
            row(7, 3),
            row(11, 6),
        ],
    }


def test_prime_deficit_obstruction_finds_redundant_classes() -> None:
    rows = certify.active_rows(toy_payload(), 0)
    obstructions = certify.prime_deficit_obstructions(rows)
    by_prime = {item["prime"]: item for item in obstructions}
    assert by_prime[3]["surviving_density"] == Fraction(1, 2)
    assert by_prime[3]["deficit"] == Fraction(1, 2)
    assert len(by_prime[3]["divisible_rows"]) == 2


def test_no_obstruction_when_surviving_density_is_one() -> None:
    rows = certify.active_rows(toy_payload(), 0)
    rows.append(
        {
            "p": 13,
            "active_modulus": 1,
            "active_class": 0,
            "residual_modulus": 2,
        }
    )
    assert not any(
        item["prime"] == 3
        for item in certify.prime_deficit_obstructions(rows)
    )
