import json
import re
from itertools import product

import pytest

from build_axis_layered_pool import sha256
from build_crt_cover_qcir import build_qcir
from inventory_reverse_group_lines import inventory
from verify_crt_cover_qcir import verify


def row(prime: int, a: int = 1, b: int = 0) -> dict:
    return {
        "p": prime,
        "h": 2,
        "a": a,
        "b": b,
        "ord2": 2 if a else 1,
        "ord3": 2 if b else 1,
        "c": 0,
        "target_modulus": 1,
        "target_residue": 0,
    }


def qcir_truth(path) -> bool:
    lines = path.read_text().splitlines()

    def identifiers(line):
        body = line[line.index("(") + 1 : line.rindex(")")]
        return [int(value) for value in body.split(",")]

    existential = identifiers(lines[1])
    universal = identifiers(lines[2])
    output = identifiers(lines[3])[0]
    gates = []
    for line in lines[4:]:
        match = re.fullmatch(r"(\d+) = (and|or)\(([-,\d]+)\)", line)
        assert match is not None
        gates.append(
            (
                int(match.group(1)),
                match.group(2),
                [int(value) for value in match.group(3).split(",")],
            )
        )

    def evaluate(inputs):
        values = dict(inputs)

        def literal_value(literal):
            value = values[abs(literal)]
            return value if literal > 0 else not value

        for identifier, operation, arguments in gates:
            terms = [literal_value(argument) for argument in arguments]
            values[identifier] = (
                all(terms) if operation == "and" else any(terms)
            )
        return values[output]

    for phase_values in product((False, True), repeat=len(existential)):
        phase = dict(zip(existential, phase_values, strict=True))
        if all(
            evaluate(
                {
                    **phase,
                    **dict(zip(universal, point, strict=True)),
                }
            )
            for point in product((False, True), repeat=len(universal))
        ):
            return True
    return False


def test_qcir_uses_existential_phases_before_universal_crt_bits(tmp_path):
    pool = {"choices": [row(5), row(7)]}
    pool_path = tmp_path / "pool.json"
    pool_path.write_text(json.dumps(pool))
    artifact = inventory(
        pool,
        pool_path,
        width=2,
        height=1,
        direction=(1, 0),
        transverse=(0, 1),
    )
    inventory_path = tmp_path / "inventory.json"
    inventory_path.write_text(json.dumps(artifact))
    qcir_path = tmp_path / "cover.qcir"
    manifest = build_qcir(
        inventory_path,
        pool_path,
        qcir_path,
    )
    lines = qcir_path.read_text().splitlines()
    assert lines[0].startswith("#QCIR-G14 ")
    assert lines[1].startswith("exists(")
    assert lines[2].startswith("forall(")
    assert lines[3].startswith("output(")
    assert manifest["row_count"] == 2
    assert manifest["phase_component_count"] == 2
    assert manifest["phase_variable_count"] == 4
    assert manifest["coordinate_variable_count"] == 1
    assert manifest["claim"]["qbf_result"] == "UNRUN"
    assert manifest["claim"]["integer_m_found"] is False
    assert qcir_truth(qcir_path) is True
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    report = verify(manifest_path)
    assert report["verified"] is True
    assert report["normalized_formula_matches"] is True


def test_qcir_builder_rejects_tampered_inventory_geometry(tmp_path):
    pool = {"choices": [row(5), row(7)]}
    pool_path = tmp_path / "pool.json"
    pool_path.write_text(json.dumps(pool))
    artifact = inventory(
        pool,
        pool_path,
        width=2,
        height=1,
        direction=(1, 0),
        transverse=(0, 1),
    )
    artifact["line_types"][0]["a"] = 0
    inventory_path = tmp_path / "inventory.json"
    inventory_path.write_text(json.dumps(artifact))
    with pytest.raises(ValueError, match="line_types"):
        build_qcir(
            inventory_path,
            pool_path,
            tmp_path / "cover.qcir",
        )


def test_qcir_verifier_rejects_semantic_gate_tamper(tmp_path):
    pool = {"choices": [row(5), row(7)]}
    pool_path = tmp_path / "pool.json"
    pool_path.write_text(json.dumps(pool))
    artifact = inventory(
        pool,
        pool_path,
        width=2,
        height=1,
        direction=(1, 0),
        transverse=(0, 1),
    )
    inventory_path = tmp_path / "inventory.json"
    inventory_path.write_text(json.dumps(artifact))
    qcir_path = tmp_path / "cover.qcir"
    manifest = build_qcir(
        inventory_path,
        pool_path,
        qcir_path,
    )
    text = qcir_path.read_text()
    assert " = or(" in text
    qcir_path.write_text(text.replace(" = or(", " = and(", 1))
    manifest["qcir_sha256"] = sha256(qcir_path)
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(AssertionError, match="formula mismatch"):
        verify(manifest_path)


def test_qcir_rejects_two_crossed_parity_lines_semantically(tmp_path):
    pool = {
        "choices": [
            row(5, 1, 0),
            row(7, 0, 1),
        ]
    }
    pool_path = tmp_path / "pool.json"
    pool_path.write_text(json.dumps(pool))
    artifact = inventory(
        pool,
        pool_path,
        width=2,
        height=2,
        direction=(1, 0),
        transverse=(0, 1),
    )
    assert artifact[
        "finite_group_cover_impossible_by_density_overlap"
    ] is False
    inventory_path = tmp_path / "inventory.json"
    inventory_path.write_text(json.dumps(artifact))
    qcir_path = tmp_path / "crossed.qcir"
    build_qcir(
        inventory_path,
        pool_path,
        qcir_path,
    )
    assert qcir_truth(qcir_path) is False
