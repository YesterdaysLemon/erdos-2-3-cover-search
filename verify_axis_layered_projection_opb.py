#!/usr/bin/env python3
"""Independently reconstruct and verify an exact projected-family OPB master."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterator

from build_axis_layered_pool import sha256
from verify_axis_layered_projection_refinement import (
    replay_dual_evidence,
    row_data,
)


def offsets_and_residuals(
    rows: list[dict],
    axis: str,
) -> tuple[list[int], list[int], list[int], int]:
    offsets = []
    moduli = []
    residuals = []
    variable_count = 0
    for row in rows:
        modulus, _coefficient, residual = row_data(row, axis)
        offsets.append(variable_count)
        moduli.append(modulus)
        residuals.append(residual)
        variable_count += modulus
    return offsets, moduli, residuals, variable_count


def constraint_line(terms: list[tuple[int, int]], degree: int) -> str:
    body = " ".join(
        f"{coefficient:+d} x{variable}"
        for coefficient, variable in terms
        if coefficient
    )
    return f"{body} >= {degree} ;\n" if body else f"0 >= {degree} ;\n"


def expected_constraints(
    rows: list[dict],
    axis: str,
    period: int,
    scale: int,
    exact_sets: list[frozenset[int]],
    pair_cuts: list[dict],
) -> Iterator[str]:
    offsets, moduli, residuals, _variable_count = offsets_and_residuals(
        rows,
        axis,
    )
    for index, modulus in enumerate(moduli):
        variables = [
            (1, offsets[index] + residue + 1)
            for residue in range(modulus)
        ]
        yield constraint_line(variables, 1)
        yield constraint_line(
            [(-1, variable) for _coefficient, variable in variables],
            -1,
        )

    seen_capacity = set()
    for coordinate in range(period):
        pattern = tuple(
            (
                scale // residuals[index],
                offsets[index] + coordinate % moduli[index] + 1,
            )
            for index in range(len(rows))
        )
        if pattern in seen_capacity:
            continue
        seen_capacity.add(pattern)
        yield constraint_line(list(pattern), scale)

    all_indices = frozenset(range(len(rows)))
    for active in exact_sets:
        outside = sorted(all_indices - active)
        seen = set()
        for coordinate in range(period):
            pattern = tuple(
                offsets[index] + coordinate % moduli[index] + 1
                for index in outside
            )
            if pattern in seen:
                continue
            seen.add(pattern)
            yield constraint_line(
                [(1, variable) for variable in pattern],
                1,
            )

    for record in pair_cuts:
        coordinate = int(record["coordinate"])
        indices = [int(index) for index in record["row_indices"]]
        if len(indices) != 2:
            raise AssertionError("invalid exact pair cut")
        left, right = indices
        if (
            not 0 <= coordinate < period
            or left < 0
            or right >= len(rows)
            or left >= right
            or math.gcd(residuals[left], residuals[right]) != 1
        ):
            raise AssertionError("invalid exact pair cut")
        overlap_denominator = residuals[left] * residuals[right]
        if scale % overlap_denominator:
            raise AssertionError("invalid exact pair cut")
        overlap = scale // overlap_denominator
        terms = []
        for index in range(len(rows)):
            coefficient = scale // residuals[index]
            if index == left or index == right:
                coefficient -= overlap
            terms.append(
                (
                    coefficient,
                    offsets[index] + coordinate % moduli[index] + 1,
                )
            )
        yield constraint_line(terms, scale - overlap)


def replay(manifest: dict) -> dict:
    if (
        manifest.get("schema")
        != "axis_layered_projection_exact_opb_master_v1"
    ):
        raise AssertionError("OPB manifest schema mismatch")
    base_path = Path(manifest["base_source"])
    if sha256(base_path) != manifest["base_source_sha256"]:
        raise AssertionError("base source SHA-256 mismatch")
    base = json.loads(base_path.read_text())
    if base.get("schema") != "axis_layered_pool_v1":
        raise AssertionError("base source schema mismatch")
    rows = base["choices"]
    axis = base["layer_axis"]
    period = int(base["capacity_pattern_period"])
    projection = int(manifest["projection_modulus"])
    extraction_path = Path(manifest["extraction_source"])
    if (
        sha256(extraction_path)
        != manifest["extraction_source_sha256"]
    ):
        raise AssertionError("extraction source SHA-256 mismatch")
    extraction = json.loads(extraction_path.read_text())
    if (
        extraction.get("schema")
        != "axis_layered_projected_refinement_v1"
        or extraction.get("layer_axis") != axis
        or int(extraction.get("capacity_pattern_period", -1)) != period
        or [
            tuple(
                int(row[key])
                for key in ("p", "h", "a", "b", "ord2", "ord3")
            )
            for row in extraction.get("choices", [])
        ]
        != [
            tuple(
                int(row[key])
                for key in ("p", "h", "a", "b", "ord2", "ord3")
            )
            for row in rows
        ]
        or int(
            extraction.get("projection_refinement", {}).get(
                "projection_modulus",
                -1,
            )
        )
        != projection
    ):
        raise AssertionError("extraction source ancestry mismatch")
    if period < 1 or projection < 1:
        raise AssertionError("periods must be positive")
    if (
        axis != manifest["layer_axis"]
        or period != int(manifest["capacity_pattern_period"])
        or len(rows) != int(manifest["row_count"])
    ):
        raise AssertionError("manifest source metadata mismatch")
    _offsets, moduli, residuals, variable_count = offsets_and_residuals(
        rows,
        axis,
    )
    if any(period % modulus for modulus in moduli):
        raise AssertionError(
            "capacity pattern period does not clear an active modulus"
        )
    if variable_count != int(manifest["variable_count"]):
        raise AssertionError("manifest variable count mismatch")
    scale = math.lcm(*residuals)
    if scale != int(manifest["capacity_scale"]):
        raise AssertionError("manifest capacity scale mismatch")

    exact_sets = []
    seen_sets = set()
    for record in manifest["exact_infeasible_frontier"]:
        active = frozenset(
            int(index) for index in record["active_row_indices"]
        )
        if (
            not active
            or any(index < 0 or index >= len(rows) for index in active)
            or active in seen_sets
        ):
            raise AssertionError("invalid exact frontier set")
        replay_dual_evidence(
            rows,
            axis,
            active,
            projection,
            record["evidence"],
        )
        seen_sets.add(active)
        exact_sets.append(active)
    if len(exact_sets) != int(manifest["exact_infeasible_frontier_count"]):
        raise AssertionError("exact frontier count mismatch")
    if any(
        left < right
        for left in exact_sets
        for right in exact_sets
        if left != right
    ):
        raise AssertionError("exact infeasible frontier is not maximal")
    expected_order = sorted(
        exact_sets,
        key=lambda active: (-len(active), sorted(active)),
    )
    if exact_sets != expected_order:
        raise AssertionError("exact frontier ordering mismatch")

    pair_cuts = manifest["pair_cuts"]
    pair_keys = [
        (
            int(record["coordinate"]),
            *[int(index) for index in record["row_indices"]],
        )
        for record in pair_cuts
    ]
    if pair_keys != sorted(set(pair_keys)):
        raise AssertionError("pair cuts are not canonical")
    if len(pair_cuts) != int(manifest["pair_cut_constraint_count"]):
        raise AssertionError("pair cut count mismatch")

    opb_path = Path(manifest["opb"])
    if sha256(opb_path) != manifest["opb_sha256"]:
        raise AssertionError("OPB SHA-256 mismatch")
    expected_header = (
        f"* #variable= {variable_count} "
        f"#constraint= {int(manifest['constraint_count'])}\n"
    )
    expected_comment = (
        "* exact necessary-condition master; SAT is not a lattice cover\n"
    )
    constraint_count = 0
    exact_cut_count = 0
    capacity_count = 0
    with opb_path.open("r", newline="") as handle:
        if handle.readline() != expected_header:
            raise AssertionError("OPB header mismatch")
        if handle.readline() != expected_comment:
            raise AssertionError("OPB comment mismatch")
        one_hot_count = 2 * len(rows)
        expected = expected_constraints(
            rows,
            axis,
            period,
            scale,
            exact_sets,
            pair_cuts,
        )
        for expected_line in expected:
            actual_line = handle.readline()
            if actual_line != expected_line:
                raise AssertionError(
                    f"OPB constraint mismatch at {constraint_count + 1}"
                )
            constraint_count += 1
            if constraint_count <= one_hot_count:
                continue
            declared_capacity_count = int(
                manifest["capacity_constraint_count"]
            )
            if constraint_count <= one_hot_count + declared_capacity_count:
                capacity_count += 1
            elif constraint_count <= (
                one_hot_count
                + declared_capacity_count
                + int(manifest["exact_cut_constraint_count"])
            ):
                exact_cut_count += 1
        if handle.readline():
            raise AssertionError("OPB contains trailing constraints")

    if (
        constraint_count != int(manifest["constraint_count"])
        or int(manifest["one_hot_constraint_count"]) != 2 * len(rows)
        or capacity_count != int(manifest["capacity_constraint_count"])
        or exact_cut_count != int(manifest["exact_cut_constraint_count"])
    ):
        raise AssertionError("OPB manifest constraint counts mismatch")
    claim = manifest["claim"]
    if (
        claim.get("opb_contains_only_exact_necessary_conditions") is not True
        or claim.get("opb_unsat_proved") is not False
        or claim.get("declared_layered_family_noncover_proved") is not False
    ):
        raise AssertionError("unproved OPB claims were promoted")
    return {
        "verified": True,
        "variable_count": variable_count,
        "constraint_count": constraint_count,
        "exact_infeasible_frontier_count": len(exact_sets),
        "exact_cut_constraint_count": exact_cut_count,
        "pair_cut_constraint_count": len(pair_cuts),
        "extraction_source_sha256": manifest[
            "extraction_source_sha256"
        ],
        "opb_unsat_proved": False,
        "declared_layered_family_noncover_proved": False,
        "engine": "independent-python-exact-opb-reconstruction",
        "scope": (
            "independent exact reconstruction of the authenticated OPB "
            "and its embedded necessary-condition evidence; no UNSAT "
            "proof log was checked"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    report = replay(manifest)
    report["manifest"] = str(args.manifest)
    report["manifest_sha256"] = sha256(args.manifest)
    report["opb"] = manifest["opb"]
    report["opb_sha256"] = manifest["opb_sha256"]
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"verified={report['verified']} "
        f"variables={report['variable_count']} "
        f"constraints={report['constraint_count']} "
        f"exact_cuts={report['exact_cut_constraint_count']} "
        f"opb_unsat_proved={report['opb_unsat_proved']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
