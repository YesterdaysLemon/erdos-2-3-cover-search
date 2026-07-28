import json
from fractions import Fraction

from refine_axis_layered_projection import (
    ProjectionMonotoneCache,
    active_set_groups,
    exact_projection_capacities,
    refinement_root,
    solve_anchor_projection_dual,
    solve_projection_subproblem,
)


def row(prime: int, h: int, a: int, b: int, ord2: int, ord3: int) -> dict:
    return {
        "p": prime,
        "h": h,
        "a": a,
        "b": b,
        "ord2": ord2,
        "ord3": ord3,
        "c": 0,
    }


def test_active_set_groups_identical_columns():
    rows = [
        row(5, 2, 1, 0, 2, 1),
        row(7, 2, 1, 0, 2, 1),
    ]
    groups = active_set_groups(rows, "k", 2, [0, 1])
    assert groups == {(0,): [0], (1,): [1]}


def test_exact_projection_capacities():
    rows = [
        row(5, 6, 0, 1, 1, 6),
        row(7, 10, 0, 1, 1, 10),
    ]
    capacities = exact_projection_capacities(
        rows,
        "k",
        (0, 1),
        projection=2,
        placement=[0, 1],
    )
    assert capacities == [Fraction(1, 3), Fraction(1, 5)]


def test_monotone_cache_classifies_superset_from_feasible_witness():
    rows = [
        row(5, 2, 0, 1, 1, 2),
        row(7, 2, 0, 1, 1, 2),
        row(11, 2, 0, 1, 1, 2),
    ]
    cache = ProjectionMonotoneCache()
    cache.add_feasible((0, 1), [0, 1])
    result = cache.classify((0, 1, 2), rows, "k", 2)
    assert result is not None
    assert result["status"] == "feasible"
    assert result["oracle"] == "monotone-feasible-subset"
    assert result["placement"] == [0, 1, 0]
    assert Fraction(
        result["minimum_numerator"],
        result["minimum_denominator"],
    ) >= 1


def test_monotone_cache_classifies_subset_from_infeasible_superset():
    rows = [
        row(5, 2, 0, 1, 1, 2),
        row(7, 2, 0, 1, 1, 2),
        row(11, 2, 0, 1, 1, 2),
    ]
    cache = ProjectionMonotoneCache()
    cache.add_infeasible((0, 1, 2))
    result = cache.classify((0, 2), rows, "k", 2)
    assert result == {
        "status": "infeasible",
        "oracle": "monotone-infeasible-superset",
        "infeasible_cut_active_row_indices": [0, 1, 2],
        "dominating_evidence_oracle": "unspecified-discovery",
    }


def test_monotone_cache_keeps_only_dominance_frontiers():
    cache = ProjectionMonotoneCache()
    cache.add_feasible((0, 1, 2), [0, 0, 0])
    cache.add_feasible((0, 1), [1, 0])
    cache.add_feasible((0, 1, 3), [1, 0, 0])
    assert set(cache.feasible) == {frozenset({0, 1})}

    cache.add_infeasible((0, 1))
    cache.add_infeasible((0, 1, 2))
    cache.add_infeasible((0, 2))
    assert set(cache.infeasible) == {frozenset({0, 1, 2})}
    assert cache.exact_infeasible == {}


def test_monotone_cache_round_trip():
    cache = ProjectionMonotoneCache()
    cache.add_feasible((0, 2), [1, 0])
    cache.add_infeasible((0, 1, 2))
    restored = ProjectionMonotoneCache()
    restored.import_payload(cache.export())
    assert restored.feasible == cache.feasible
    assert restored.infeasible == cache.infeasible
    assert restored.exact_infeasible == cache.exact_infeasible


def test_exact_infeasible_frontier_survives_stronger_discovery_entry():
    cache = ProjectionMonotoneCache()
    cache.add_infeasible(
        (0, 1),
        {"oracle": "exact-anchor-projection-dual"},
    )
    cache.add_infeasible(
        (0, 1, 2),
        {"oracle": "scipy-highs-discovery"},
    )
    assert set(cache.infeasible) == {frozenset({0, 1, 2})}
    assert set(cache.exact_infeasible) == {frozenset({0, 1})}
    result = cache.classify(
        (0,),
        [
            row(5, 2, 0, 1, 1, 2),
            row(7, 2, 0, 1, 1, 2),
            row(11, 2, 0, 1, 1, 2),
        ],
        "k",
        2,
        exact_infeasible_only=True,
    )
    assert result is not None
    assert result["oracle"] == "monotone-exact-infeasible-superset"


def test_anchor_projection_dual_proves_sparse_full_cell_union_infeasible():
    rows = [
        row(7, 6, 0, 1, 1, 6),
        row(11, 5, 0, 1, 1, 5),
        row(31, 30, 0, 1, 1, 30),
    ]
    result = solve_anchor_projection_dual(
        rows,
        "k",
        (0, 1, 2),
        projection=30,
    )
    assert result is not None
    assert result["status"] == "infeasible"
    assert result["oracle"] == "exact-anchor-projection-dual"
    obstruction = result["dual_obstruction"]
    assert obstruction["base_anchor_primes"] == [7, 11]
    assert obstruction["branch_anchor_prime"] == 31
    assert len(obstruction["branches"]) == 30
    assert all(
        branch["strict_gap_numerator"] > 0
        for branch in obstruction["branches"]
    )


def test_anchor_projection_dual_greedily_expands_exact_failed_set():
    rows = [
        row(7, 6, 0, 1, 1, 6),
        row(11, 5, 0, 1, 1, 5),
        row(31, 30, 0, 1, 1, 30),
        row(211, 210, 0, 1, 1, 210),
    ]
    result = solve_anchor_projection_dual(
        rows,
        "k",
        (0, 1, 2),
        projection=30,
    )
    assert result is not None
    assert result["oracle"] == "exact-anchor-projection-dual-expanded"
    assert result["infeasible_cut_active_row_indices"] == [0, 1, 2, 3]
    assert result["dual_obstruction"]["greedy_added_row_indices"] == [3]


def test_maximin_subproblem_returns_exact_feasible_witness():
    rows = [
        row(5, 2, 0, 1, 1, 2),
        row(7, 2, 0, 1, 1, 2),
    ]
    result = solve_projection_subproblem(
        rows,
        "k",
        (0, 1),
        projection=2,
        time_limit=5.0,
    )
    assert result["status"] == "feasible"
    assert result["minimum_numerator"] >= result["minimum_denominator"]


def test_maximin_subproblem_detects_simple_discovery_infeasibility():
    result = solve_projection_subproblem(
        [row(5, 2, 0, 1, 1, 2)],
        "k",
        (0,),
        projection=2,
        time_limit=5.0,
    )
    assert result["status"] == "infeasible"
    assert result["oracle"] == "scipy-highs-maximin-discovery"


def test_refinement_root_walks_hash_checked_parent_chain(tmp_path):
    root = tmp_path / "root.json"
    root.write_text(json.dumps({"schema": "axis_layered_pool_v1"}))
    from build_axis_layered_pool import sha256

    child = tmp_path / "child.json"
    child_payload = {
        "schema": "axis_layered_projected_refinement_v1",
        "refinement_source": str(root),
        "refinement_source_sha256": sha256(root),
    }
    child.write_text(json.dumps(child_payload))
    grandchild = {
        "schema": "axis_layered_projected_refinement_v1",
        "refinement_source": str(child),
        "refinement_source_sha256": sha256(child),
    }
    assert refinement_root(tmp_path / "grandchild.json", grandchild) == (
        str(root),
        sha256(root),
    )
