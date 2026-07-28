#!/usr/bin/env python3
"""Build an exact pseudo-Boolean master for a projected layered family.

The input refinement may contain floating-point discovery records, but this
builder admits only:

* ordinary layer-capacity constraints, scaled to exact integers;
* exact anchor-dual infeasible supersets, replayed independently first; and
* structurally unavoidable coprime-residual pair cuts.

The resulting OPB instance asks whether any active-class placement survives
all currently certified necessary conditions.  SAT is only a survivor.
UNSAT becomes a finite-family theorem only after a proof-producing solver's
certificate and the generated manifest have both been independently checked.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterator

from build_axis_layered_pool import layer_data, sha256
from refine_axis_layered_projection import ProjectionMonotoneCache
from verify_axis_layered_projection_refinement import replay_dual_evidence


def row_identity(row: dict) -> tuple[int, ...]:
    return tuple(
        int(row[key]) for key in ("p", "h", "a", "b", "ord2", "ord3")
    )


def validate_sources(base: dict, refinement: dict) -> None:
    if base.get("schema") != "axis_layered_pool_v1":
        raise ValueError("base source is not an axis-layered pool")
    if (
        refinement.get("schema")
        != "axis_layered_projected_refinement_v1"
    ):
        raise ValueError("refinement source schema mismatch")
    if base["layer_axis"] != refinement["layer_axis"]:
        raise ValueError("layer axes differ")
    if int(base["capacity_pattern_period"]) != int(
        refinement["capacity_pattern_period"]
    ):
        raise ValueError("capacity pattern periods differ")
    base_rows = base["choices"]
    refinement_rows = refinement["choices"]
    if len(base_rows) != len(refinement_rows):
        raise ValueError("row counts differ")
    if [row_identity(row) for row in base_rows] != [
        row_identity(row) for row in refinement_rows
    ]:
        raise ValueError("refinement row identities differ from base")


def variable_layout(
    rows: list[dict],
    axis: str,
) -> tuple[list[int], int]:
    offsets = []
    variable_count = 0
    for row in rows:
        _period, modulus, _coefficient, _residual = layer_data(row, axis)
        offsets.append(variable_count)
        variable_count += modulus
    return offsets, variable_count


def literal_id(
    offsets: list[int],
    rows: list[dict],
    axis: str,
    row_index: int,
    coordinate: int,
) -> int:
    modulus = layer_data(rows[row_index], axis)[1]
    return offsets[row_index] + coordinate % modulus + 1


def format_constraint(
    terms: list[tuple[int, int]],
    degree: int,
) -> str:
    body = " ".join(
        f"{coefficient:+d} x{variable}"
        for coefficient, variable in terms
        if coefficient
    )
    return f"{body} >= {degree} ;\n" if body else f"0 >= {degree} ;\n"


def unique_capacity_patterns(
    rows: list[dict],
    axis: str,
    period: int,
    offsets: list[int],
    scale: int,
) -> Iterator[tuple[tuple[int, int], ...]]:
    data = [layer_data(row, axis) for row in rows]
    seen = set()
    for coordinate in range(period):
        pattern = tuple(
            (
                scale // residual,
                offsets[index] + coordinate % modulus + 1,
            )
            for index, (
                _coordinate_period,
                modulus,
                _coefficient,
                residual,
            ) in enumerate(data)
        )
        if pattern in seen:
            continue
        seen.add(pattern)
        yield pattern


def unique_exact_cut_patterns(
    rows: list[dict],
    axis: str,
    period: int,
    offsets: list[int],
    exact_sets: list[frozenset[int]],
) -> Iterator[tuple[int, ...]]:
    all_indices = frozenset(range(len(rows)))
    for active in exact_sets:
        outside = sorted(all_indices - active)
        seen = set()
        for coordinate in range(period):
            pattern = tuple(
                literal_id(
                    offsets,
                    rows,
                    axis,
                    row_index,
                    coordinate,
                )
                for row_index in outside
            )
            if pattern in seen:
                continue
            seen.add(pattern)
            yield pattern


def canonical_pair_cuts(
    refinement: dict,
    rows: list[dict],
    axis: str,
    period: int,
) -> list[dict]:
    records = {}
    for raw in refinement["projection_refinement"].get(
        "pair_benders_cuts",
        [],
    ):
        coordinate = int(raw["coordinate"]) % period
        left, right = sorted(int(index) for index in raw["row_indices"])
        if (
            left < 0
            or right >= len(rows)
            or left == right
            or math.gcd(
                layer_data(rows[left], axis)[3],
                layer_data(rows[right], axis)[3],
            )
            != 1
        ):
            raise ValueError("invalid persisted pair cut")
        records[(coordinate, left, right)] = {
            "coordinate": coordinate,
            "row_indices": [left, right],
        }
    return [records[key] for key in sorted(records)]


def pair_constraint(
    record: dict,
    rows: list[dict],
    axis: str,
    offsets: list[int],
    scale: int,
) -> tuple[list[tuple[int, int]], int]:
    coordinate = int(record["coordinate"])
    left, right = (int(index) for index in record["row_indices"])
    data = [layer_data(row, axis) for row in rows]
    left_residual = data[left][3]
    right_residual = data[right][3]
    overlap_denominator = left_residual * right_residual
    if scale % overlap_denominator:
        raise AssertionError("capacity scale does not clear pair overlap")
    overlap = scale // overlap_denominator
    terms = []
    for index, (
        _coordinate_period,
        modulus,
        _coefficient,
        residual,
    ) in enumerate(data):
        coefficient = scale // residual
        if index == left or index == right:
            coefficient -= overlap
        terms.append(
            (
                coefficient,
                offsets[index] + coordinate % modulus + 1,
            )
        )
    return terms, scale - overlap


def exact_frontier(
    refinement: dict,
    rows: list[dict],
    axis: str,
    projection: int,
) -> tuple[list[dict], dict]:
    cache = ProjectionMonotoneCache()
    cache.import_payload(
        refinement["projection_refinement"]["monotonicity_cache"]
    )
    strengthening = cache.strengthen_exact_infeasible(
        rows,
        axis,
        projection,
    )
    records = []
    for active, evidence in sorted(
        cache.exact_infeasible.items(),
        key=lambda item: (-len(item[0]), sorted(item[0])),
    ):
        replay_dual_evidence(rows, axis, active, projection, evidence)
        records.append(
            {
                "active_row_indices": sorted(active),
                "evidence": evidence,
            }
        )
    if not records:
        raise ValueError("refinement contains no exact infeasible frontier")
    return records, strengthening


def build_manifest_and_opb(
    base_path: Path,
    refinement_path: Path,
    opb_path: Path,
) -> dict:
    base = json.loads(base_path.read_text())
    refinement = json.loads(refinement_path.read_text())
    validate_sources(base, refinement)
    rows = base["choices"]
    axis = base["layer_axis"]
    period = int(base["capacity_pattern_period"])
    projection = int(
        refinement["projection_refinement"]["projection_modulus"]
    )
    if period < 1 or projection < 1:
        raise ValueError("periods must be positive")
    if any(period % layer_data(row, axis)[1] for row in rows):
        raise ValueError(
            "capacity pattern period does not clear an active modulus"
        )
    exact_records, strengthening = exact_frontier(
        refinement,
        rows,
        axis,
        projection,
    )
    exact_sets = [
        frozenset(int(index) for index in record["active_row_indices"])
        for record in exact_records
    ]
    offsets, variable_count = variable_layout(rows, axis)
    residuals = [layer_data(row, axis)[3] for row in rows]
    scale = math.lcm(*residuals)
    capacity_patterns = list(
        unique_capacity_patterns(rows, axis, period, offsets, scale)
    )
    exact_cut_count = sum(
        1
        for _pattern in unique_exact_cut_patterns(
            rows,
            axis,
            period,
            offsets,
            exact_sets,
        )
    )
    pair_cuts = canonical_pair_cuts(
        refinement,
        rows,
        axis,
        period,
    )
    constraint_count = (
        2 * len(rows)
        + len(capacity_patterns)
        + exact_cut_count
        + len(pair_cuts)
    )

    with opb_path.open("w", newline="\n") as handle:
        handle.write(
            f"* #variable= {variable_count} "
            f"#constraint= {constraint_count}\n"
        )
        handle.write(
            "* exact necessary-condition master; SAT is not a lattice cover\n"
        )
        for index, row in enumerate(rows):
            modulus = layer_data(row, axis)[1]
            variables = [
                (1, offsets[index] + residue + 1)
                for residue in range(modulus)
            ]
            handle.write(format_constraint(variables, 1))
            handle.write(
                format_constraint(
                    [
                        (-coefficient, variable)
                        for coefficient, variable in variables
                    ],
                    -1,
                )
            )
        for pattern in capacity_patterns:
            handle.write(format_constraint(list(pattern), scale))
        for pattern in unique_exact_cut_patterns(
            rows,
            axis,
            period,
            offsets,
            exact_sets,
        ):
            handle.write(
                format_constraint(
                    [(1, variable) for variable in pattern],
                    1,
                )
            )
        for record in pair_cuts:
            terms, degree = pair_constraint(
                record,
                rows,
                axis,
                offsets,
                scale,
            )
            handle.write(format_constraint(terms, degree))

    return {
        "schema": "axis_layered_projection_exact_opb_master_v1",
        "base_source": str(base_path),
        "base_source_sha256": sha256(base_path),
        "extraction_source": str(refinement_path),
        "extraction_source_sha256": sha256(refinement_path),
        "layer_axis": axis,
        "coordinate_basis": base.get("coordinate_basis"),
        "capacity_pattern_period": period,
        "projection_modulus": projection,
        "row_count": len(rows),
        "variable_count": variable_count,
        "constraint_count": constraint_count,
        "one_hot_constraint_count": 2 * len(rows),
        "capacity_scale": scale,
        "capacity_constraint_count": len(capacity_patterns),
        "exact_infeasible_frontier": exact_records,
        "exact_infeasible_frontier_count": len(exact_records),
        "exact_cut_constraint_count": exact_cut_count,
        "pair_cuts": pair_cuts,
        "pair_cut_constraint_count": len(pair_cuts),
        "exact_cache_strengthening": strengthening,
        "opb": str(opb_path),
        "opb_sha256": sha256(opb_path),
        "claim": {
            "opb_contains_only_exact_necessary_conditions": True,
            "opb_unsat_proved": False,
            "declared_layered_family_noncover_proved": False,
        },
        "scope": (
            "exact necessary-condition master for every active-class "
            "placement of the declared finite 81-row family; SAT would only "
            "be a projected survivor, while UNSAT still requires a checked "
            "proof log before any family-level no-cover claim"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("base", type=Path)
    parser.add_argument("refinement", type=Path)
    parser.add_argument("--opb-output", type=Path, required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    args = parser.parse_args()
    manifest = build_manifest_and_opb(
        args.base,
        args.refinement,
        args.opb_output,
    )
    args.manifest_output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"variables={manifest['variable_count']} "
        f"constraints={manifest['constraint_count']} "
        f"exact_frontier={manifest['exact_infeasible_frontier_count']} "
        f"exact_cuts={manifest['exact_cut_constraint_count']} "
        f"opb={args.opb_output} manifest={args.manifest_output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
