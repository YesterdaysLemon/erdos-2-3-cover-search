import json

from build_axis_layered_pool import sha256
from build_axis_layered_projection_opb import build_manifest_and_opb
from refine_axis_layered_projection import solve_anchor_projection_dual
from verify_axis_layered_projection_opb import replay


def row(prime: int, h: int) -> dict:
    return {
        "p": prime,
        "h": h,
        "a": 0,
        "b": 1 if h > 1 else 0,
        "ord2": 1,
        "ord3": h,
        "c": 0,
        "target_modulus": 1,
        "target_residue": 0,
        "layer_active_class": 0,
    }


def synthetic_sources(tmp_path):
    rows = [
        row(7, 6),
        row(11, 5),
        row(31, 30),
        row(211, 210),
        row(2, 1),
    ]
    evidence = solve_anchor_projection_dual(
        rows,
        "k",
        (0, 1, 2),
        projection=30,
    )
    assert evidence is not None
    active = evidence["infeasible_cut_active_row_indices"]
    base = {
        "schema": "axis_layered_pool_v1",
        "layer_axis": "k",
        "capacity_pattern_period": 1,
        "coordinate_basis": None,
        "choices": rows,
    }
    refinement = {
        **base,
        "schema": "axis_layered_projected_refinement_v1",
        "projection_refinement": {
            "projection_modulus": 30,
            "monotonicity_cache": {
                "feasible": [],
                "infeasible": [],
                "exact_infeasible": [
                    {
                        "active_row_indices": active,
                        "evidence": evidence,
                    }
                ],
            },
            "pair_benders_cuts": [
                {
                    "coordinate": 0,
                    "row_indices": [0, 1],
                }
            ],
        },
    }
    base_path = tmp_path / "base.json"
    refinement_path = tmp_path / "refinement.json"
    base_path.write_text(json.dumps(base))
    refinement_path.write_text(json.dumps(refinement))
    return base_path, refinement_path


def test_exact_opb_master_round_trip(tmp_path):
    base_path, refinement_path = synthetic_sources(tmp_path)
    opb_path = tmp_path / "master.opb"
    manifest_path = tmp_path / "manifest.json"
    manifest = build_manifest_and_opb(
        base_path,
        refinement_path,
        opb_path,
    )
    manifest_path.write_text(json.dumps(manifest))
    report = replay(json.loads(manifest_path.read_text()))
    assert report["verified"] is True
    assert report["variable_count"] == 5
    assert report["exact_infeasible_frontier_count"] == 1
    assert report["pair_cut_constraint_count"] == 1
    assert report["opb_unsat_proved"] is False


def test_exact_opb_verifier_rejects_tampered_constraint(tmp_path):
    base_path, refinement_path = synthetic_sources(tmp_path)
    opb_path = tmp_path / "master.opb"
    manifest = build_manifest_and_opb(
        base_path,
        refinement_path,
        opb_path,
    )
    text = opb_path.read_text()
    opb_path.write_text(
        text.replace(">= 1 ;", ">= 2 ;", 1),
        newline="\n",
    )
    manifest["opb_sha256"] = sha256(opb_path)
    try:
        replay(manifest)
    except AssertionError as error:
        assert "constraint mismatch" in str(error)
    else:
        raise AssertionError("tampered OPB constraint was accepted")


def test_exact_opb_verifier_rejects_tampered_extraction_source(tmp_path):
    base_path, refinement_path = synthetic_sources(tmp_path)
    opb_path = tmp_path / "master.opb"
    manifest = build_manifest_and_opb(
        base_path,
        refinement_path,
        opb_path,
    )
    refinement = json.loads(refinement_path.read_text())
    refinement["choices"][0]["p"] = 13
    refinement_path.write_text(json.dumps(refinement))
    try:
        replay(manifest)
    except AssertionError as error:
        assert "extraction source SHA-256 mismatch" in str(error)
    else:
        raise AssertionError("tampered extraction source was accepted")
