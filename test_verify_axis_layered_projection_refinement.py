from refine_axis_layered_projection import solve_anchor_projection_dual
from verify_axis_layered_projection_refinement import replay_dual_evidence


def row(prime: int, h: int) -> dict:
    return {
        "p": prime,
        "h": h,
        "a": 0,
        "b": 1,
        "ord2": 1,
        "ord3": h,
        "c": 0,
    }


def test_independent_dual_replay_accepts_expanded_obstruction():
    rows = [
        row(7, 6),
        row(11, 5),
        row(31, 30),
        row(211, 210),
    ]
    evidence = solve_anchor_projection_dual(
        rows,
        "k",
        (0, 1, 2),
        projection=30,
    )
    assert evidence is not None
    active = frozenset(evidence["infeasible_cut_active_row_indices"])
    report = replay_dual_evidence(
        rows,
        "k",
        active,
        projection=30,
        evidence=evidence,
    )
    assert report["active_row_count"] == 4
    assert report["minimum_gap_numerator"] > 0


def test_independent_dual_replay_rejects_tampered_gap():
    rows = [row(7, 6), row(11, 5), row(31, 30)]
    evidence = solve_anchor_projection_dual(
        rows,
        "k",
        (0, 1, 2),
        projection=30,
    )
    assert evidence is not None
    evidence["dual_obstruction"]["branches"][0][
        "strict_gap_numerator"
    ] += 1
    try:
        replay_dual_evidence(
            rows,
            "k",
            frozenset({0, 1, 2}),
            projection=30,
            evidence=evidence,
        )
    except AssertionError as error:
        assert "gap mismatch" in str(error)
    else:
        raise AssertionError("tampered dual gap was accepted")
