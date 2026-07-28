#!/usr/bin/env python3
"""Build a symbolic CRT QBF for one finite arithmetic quotient family.

The encoding is

    exists row-phase CRT digits
    forall quotient-coordinate CRT bits
        every valid coordinate is covered.

Prime-power lookup gates replace explicit quotient cells.  The output is a
QCIR instance only; SAT is a candidate phase cover, and UNSAT is not promoted
without a checked proof and independent reconstruction.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from build_axis_layered_pool import sha256
from inventory_reverse_group_lines import (
    inventory as rebuild_inventory,
    prime_power_components,
)


def prime_base(component: int) -> int:
    divisor = 2
    while component % divisor:
        divisor += 1 if divisor == 2 else 2
    value = component
    while value % divisor == 0:
        value //= divisor
    if value != 1:
        raise ValueError("component is not a prime power")
    return divisor


def valuation(value: int, prime: int) -> int:
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


class QcirBuilder:
    def __init__(self) -> None:
        self.next_identifier = 1
        self.gates: list[str] = []

    def variables(self, count: int) -> list[int]:
        identifiers = list(
            range(
                self.next_identifier,
                self.next_identifier + count,
            )
        )
        self.next_identifier += count
        return identifiers

    def gate(self, operation: str, arguments: list[int]) -> int:
        if operation not in {"and", "or"} or not arguments:
            raise ValueError("QCIR gate must be a nonempty and/or")
        identifier = self.next_identifier
        self.next_identifier += 1
        joined = ",".join(str(argument) for argument in arguments)
        self.gates.append(f"{identifier} = {operation}({joined})\n")
        return identifier


def assignment_gate(
    builder: QcirBuilder,
    bits: list[int],
    value: int,
) -> int:
    if not bits:
        raise ValueError("constant coordinate needs no assignment gate")
    literals = [
        bit if value & (1 << index) else -bit
        for index, bit in enumerate(bits)
    ]
    return builder.gate("and", literals)


def validate_inventory(
    inventory_path: Path,
    pool_path: Path,
) -> tuple[dict, dict]:
    artifact = json.loads(inventory_path.read_text())
    if artifact.get("schema") != "reverse_group_arithmetic_line_inventory_v3":
        raise ValueError("inventory schema mismatch")
    if sha256(pool_path) != artifact["source_sha256"]:
        raise ValueError("inventory source SHA-256 mismatch")
    pool = json.loads(pool_path.read_text())
    direction = tuple(
        int(value) for value in artifact["basis"]["direction"]
    )
    transverse = tuple(
        int(value) for value in artifact["basis"]["transverse"]
    )
    rebuilt = rebuild_inventory(
        pool,
        pool_path,
        int(artifact["group"]["width"]),
        int(artifact["group"]["height"]),
        direction,
        transverse,
    )
    for key in ("basis", "group", "line_types"):
        if artifact[key] != rebuilt[key]:
            raise ValueError(f"inventory {key} differs from source rebuild")
    if rebuilt[
        "finite_group_cover_impossible_by_density_overlap"
    ]:
        raise ValueError("declared family already has an exact noncover bound")
    return rebuilt, pool


def expanded_rows(artifact: dict) -> list[dict]:
    rows = []
    seen_primes = set()
    for line_type in artifact["line_types"]:
        if int(line_type["target_modulus"]) != 1:
            raise ValueError("QCIR prototype requires unrestricted targets")
        h = int(line_type["h"])
        if h <= 1:
            raise ValueError("QCIR prototype requires row moduli above one")
        a = int(line_type["a"]) % h
        b = int(line_type["b"]) % h
        for prime in line_type["primes"]:
            prime = int(prime)
            if prime in seen_primes:
                raise ValueError(f"inventory reuses prime {prime}")
            seen_primes.add(prime)
            rows.append({"p": prime, "h": h, "a": a, "b": b})
    rows.sort(key=lambda row: (row["h"], row["p"]))
    if not rows:
        raise ValueError("inventory has no descending rows")
    return rows


def build_qcir(
    inventory_path: Path,
    pool_path: Path,
    qcir_path: Path,
) -> dict:
    artifact, _pool = validate_inventory(inventory_path, pool_path)
    rows = expanded_rows(artifact)
    width = int(artifact["group"]["width"])
    height = int(artifact["group"]["height"])
    builder = QcirBuilder()

    phase_variables = []
    phase_layout = []
    row_components = []
    for row_index, row in enumerate(rows):
        components = []
        for component in prime_power_components(row["h"]):
            prime = prime_base(component)
            variables = builder.variables(component)
            phase_variables.extend(variables)
            layout = {
                "row_index": row_index,
                "p": row["p"],
                "prime": prime,
                "modulus": component,
                "variables": variables,
            }
            phase_layout.append(layout)
            components.append(layout)
        row_components.append(components)

    group_primes = sorted(
        {
            prime_base(component)
            for value in (width, height)
            for component in prime_power_components(value)
        }
    )
    coordinate_variables = []
    coordinate_layout = {}
    for prime in group_primes:
        x_modulus = prime ** valuation(width, prime)
        y_modulus = prime ** valuation(height, prime)
        x_bits = builder.variables((x_modulus - 1).bit_length())
        y_bits = builder.variables((y_modulus - 1).bit_length())
        coordinate_variables.extend(x_bits)
        coordinate_variables.extend(y_bits)
        coordinate_layout[prime] = {
            "prime": prime,
            "x_modulus": x_modulus,
            "y_modulus": y_modulus,
            "x_bits": x_bits,
            "y_bits": y_bits,
        }

    phase_constraints = []
    at_most_one_clause_count = 0
    for layout in phase_layout:
        variables = layout["variables"]
        phase_constraints.append(builder.gate("or", variables))
        for left_index, left in enumerate(variables):
            for right in variables[left_index + 1 :]:
                phase_constraints.append(
                    builder.gate("or", [-left, -right])
                )
                at_most_one_clause_count += 1
    phase_valid = builder.gate("and", phase_constraints)

    point_gates = {}
    component_valid_gates = {}
    for prime in group_primes:
        layout = coordinate_layout[prime]
        local_gates = []
        for x in range(layout["x_modulus"]):
            x_gate = (
                assignment_gate(builder, layout["x_bits"], x)
                if layout["x_bits"]
                else None
            )
            for y in range(layout["y_modulus"]):
                y_gate = (
                    assignment_gate(builder, layout["y_bits"], y)
                    if layout["y_bits"]
                    else None
                )
                arguments = [
                    gate
                    for gate in (x_gate, y_gate)
                    if gate is not None
                ]
                if len(arguments) == 1:
                    point_gate = arguments[0]
                else:
                    point_gate = builder.gate("and", arguments)
                point_gates[(prime, x, y)] = point_gate
                local_gates.append(point_gate)
        component_valid_gates[prime] = builder.gate(
            "or",
            local_gates,
        )

    bucket_cache = {}
    row_match_gates = []
    for row, components in zip(rows, row_components, strict=True):
        component_matches = []
        for phase_component in components:
            prime = phase_component["prime"]
            modulus = phase_component["modulus"]
            cache_key = (
                prime,
                modulus,
                row["a"] % modulus,
                row["b"] % modulus,
            )
            buckets = bucket_cache.get(cache_key)
            if buckets is None:
                layout = coordinate_layout[prime]
                bucket_points = [[] for _ in range(modulus)]
                for x in range(layout["x_modulus"]):
                    for y in range(layout["y_modulus"]):
                        residue = (
                            row["a"] * x + row["b"] * y
                        ) % modulus
                        bucket_points[residue].append(
                            point_gates[(prime, x, y)]
                        )
                if any(not points for points in bucket_points):
                    raise AssertionError(
                        "local affine component is not surjective"
                    )
                buckets = [
                    builder.gate("or", points)
                    for points in bucket_points
                ]
                bucket_cache[cache_key] = buckets
            target_matches = [
                builder.gate("and", [variable, buckets[residue]])
                for residue, variable in enumerate(
                    phase_component["variables"]
                )
            ]
            component_matches.append(
                builder.gate("or", target_matches)
            )
        row_match_gates.append(
            (
                component_matches[0]
                if len(component_matches) == 1
                else builder.gate("and", component_matches)
            )
        )

    coverage = builder.gate("or", row_match_gates)
    point_valid = (
        next(iter(component_valid_gates.values()))
        if len(component_valid_gates) == 1
        else builder.gate(
            "and",
            list(component_valid_gates.values()),
        )
    )
    valid_implies_covered = builder.gate(
        "or",
        [-point_valid, coverage],
    )
    output_gate = builder.gate(
        "and",
        [phase_valid, valid_implies_covered],
    )
    identifier_count = builder.next_identifier - 1
    header = [
        f"#QCIR-G14 {identifier_count}\n",
        f"exists({','.join(str(v) for v in phase_variables)})\n",
        (
            f"forall({','.join(str(v) for v in coordinate_variables)})\n"
        ),
        f"output({output_gate})\n",
    ]
    qcir_path.write_text("".join(header + builder.gates), newline="\n")
    return {
        "schema": "crt_affine_cover_qcir_manifest_v1",
        "inventory": str(inventory_path),
        "inventory_sha256": sha256(inventory_path),
        "pool": str(pool_path),
        "pool_sha256": sha256(pool_path),
        "group": artifact["group"],
        "basis": artifact["basis"],
        "row_count": len(rows),
        "rows": rows,
        "phase_component_count": len(phase_layout),
        "phase_variable_count": len(phase_variables),
        "phase_layout": phase_layout,
        "at_most_one_clause_count": at_most_one_clause_count,
        "coordinate_variable_count": len(coordinate_variables),
        "coordinate_layout": [
            coordinate_layout[prime] for prime in group_primes
        ],
        "gate_count": len(builder.gates),
        "identifier_count": identifier_count,
        "output_gate": output_gate,
        "qcir": str(qcir_path),
        "qcir_sha256": sha256(qcir_path),
        "quantifier_order": (
            "exists phase digits; forall CRT coordinate bits; "
            "deterministic lookup gates"
        ),
        "claim": {
            "encoding_generated": True,
            "qbf_result": "UNRUN",
            "finite_group_cover_proved": False,
            "finite_group_noncover_proved": False,
            "integer_m_found": False,
        },
        "scope": (
            "exact symbolic cover decision instance for the authenticated "
            "declared finite quotient family; no solver proof checked"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inventory", type=Path)
    parser.add_argument("pool", type=Path)
    parser.add_argument("--qcir-output", type=Path, required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    args = parser.parse_args()
    manifest = build_qcir(
        args.inventory,
        args.pool,
        args.qcir_output,
    )
    args.manifest_output.write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(
        f"rows={manifest['row_count']} "
        f"phase_variables={manifest['phase_variable_count']} "
        f"coordinate_variables={manifest['coordinate_variable_count']} "
        f"gates={manifest['gate_count']} "
        f"qbf_result={manifest['claim']['qbf_result']} "
        f"output={args.qcir_output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
