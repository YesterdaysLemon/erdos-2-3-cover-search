import certify_axis_layered_family_noncover as certify


def test_exact_mdd_rejects_two_column_capacity_gap() -> None:
    rows = [
        {
            "p": 5,
            "h": 2,
            "active_modulus": 1,
            "residual_modulus": 2,
        },
        {
            "p": 13,
            "h": 4,
            "active_modulus": 2,
            "residual_modulus": 2,
        },
    ]
    result = certify.exact_mdd_capacity_obstruction(rows, [0, 1])
    assert result["capacity_period"] == 2
    assert result["terminal_root"] == 0
    assert result["proved_no_capacity_placement"] is True


def test_prime_deficit_pruning_is_iterative() -> None:
    rows = [
        {"p": 5, "h": 2, "active_modulus": 1, "residual_modulus": 2},
        {"p": 7, "h": 6, "active_modulus": 1, "residual_modulus": 6},
        {"p": 11, "h": 3, "active_modulus": 1, "residual_modulus": 3},
    ]
    surviving, rounds = certify.prune_prime_deficits(rows)
    assert rounds[0]["prime"] == 3
    assert rounds[0]["divisible_row_count"] == 2
    assert rounds[1]["prime"] == 2
    assert rounds[1]["divisible_row_count"] == 1
    assert surviving == []
