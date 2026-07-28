#!/usr/bin/env python3
"""Independently verify a symbolic CRT affine-cover QCIR instance.

This parser does not import the QCIR exporter.  It authenticates and rebuilds
the arithmetic inventory, parses the cleansed circuit, normalizes its Boolean
expression modulo associativity, commutativity, and idempotence, and compares
it with a separately constructed expected expression.

The verifier checks the encoding only.  It never promotes an unrun solver
result or claims an integer m.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from build_axis_layered_pool import sha256
from inventory_reverse_group_lines import inventory as rebuild_inventory


def factors(value: int) -> list[tuple[int, int]]:
    result = []
    divisor = 2
    while divisor * divisor <= value:
        component = 1
        while value % divisor == 0:
            value //= divisor
            component *= divisor
        if component > 1:
            result.append((divisor, component))
        divisor += 1 if divisor == 2 else 2
    if value > 1:
        result.append((value, value))
    return result


def exponent(value: int, prime: int) -> int:
    result = 0
    while value % prime == 0:
        value //= prime
        result += 1
    return result


class Normalizer:
    """Intern canonical Boolean DAG nodes."""

    def __init__(self) -> None:
        self.node_to_id = {}
        self.id_to_node = {}
        self.next_id = 1

    def intern(self, node: tuple) -> int:
        identifier = self.node_to_id.get(node)
        if identifier is None:
            identifier = self.next_id
            self.next_id += 1
            self.node_to_id[node] = identifier
            self.id_to_node[identifier] = node
        return identifier

    def variable(self, identifier: int) -> int:
        return self.intern(("variable", identifier))

    def negate(self, child: int) -> int:
        node = self.id_to_node[child]
        if node[0] == "not":
            return int(node[1])
        return self.intern(("not", child))

    def gate(self, operation: str, children: list[int]) -> int:
        if operation not in {"and", "or"} or not children:
            raise AssertionError("invalid nonconstant Boolean gate")
        flattened = []
        for child in children:
            node = self.id_to_node[child]
            if node[0] == operation:
                flattened.extend(int(value) for value in node[1:])
            else:
                flattened.append(child)
        canonical = sorted(set(flattened))
        if len(canonical) == 1:
            return canonical[0]
        return self.intern((operation, *canonical))


def parse_list(line: str, keyword: str) -> list[int]:
    match = re.fullmatch(
        rf"{keyword}\((\d+(?:,\d+)*)\)",
        line,
    )
    if match is None:
        raise AssertionError(f"invalid {keyword} declaration")
    return [int(value) for value in match.group(1).split(",")]


def parse_qcir(path: Path) -> dict:
    lines = path.read_text().splitlines()
    if len(lines) < 5:
        raise AssertionError("QCIR file is truncated")
    header = re.fullmatch(r"#QCIR-G14 (\d+)", lines[0])
    if header is None:
        raise AssertionError("QCIR header is not cleansed G14")
    declared_identifiers = int(header.group(1))
    existential = parse_list(lines[1], "exists")
    universal = parse_list(lines[2], "forall")
    output_values = parse_list(lines[3], "output")
    if len(output_values) != 1:
        raise AssertionError("QCIR output is not singular")
    output = output_values[0]
    quantified = set(existential) | set(universal)
    if (
        len(quantified) != len(existential) + len(universal)
        or not quantified
    ):
        raise AssertionError("QCIR prefix is duplicated or empty")

    normalizer = Normalizer()
    expressions = {
        identifier: normalizer.variable(identifier)
        for identifier in quantified
    }
    gate_identifiers = set()
    gate_pattern = re.compile(
        r"(\d+) = (and|or)\((-?\d+(?:,-?\d+)*)\)"
    )
    for line_number, line in enumerate(lines[4:], start=5):
        match = gate_pattern.fullmatch(line)
        if match is None:
            raise AssertionError(
                f"invalid QCIR gate syntax on line {line_number}"
            )
        identifier = int(match.group(1))
        if identifier in expressions:
            raise AssertionError("QCIR identifier is redefined")
        arguments = [
            int(value) for value in match.group(3).split(",")
        ]
        children = []
        for literal in arguments:
            referenced = expressions.get(abs(literal))
            if referenced is None:
                raise AssertionError(
                    "QCIR gate uses an undefined or forward identifier"
                )
            children.append(
                referenced
                if literal > 0
                else normalizer.negate(referenced)
            )
        expressions[identifier] = normalizer.gate(
            match.group(2),
            children,
        )
        gate_identifiers.add(identifier)

    all_identifiers = quantified | gate_identifiers
    if (
        all_identifiers != set(range(1, declared_identifiers + 1))
        or output not in gate_identifiers
    ):
        raise AssertionError("QCIR cleansed identifier closure mismatch")
    return {
        "declared_identifier_count": declared_identifiers,
        "existential": existential,
        "universal": universal,
        "output": output,
        "output_expression": expressions[output],
        "normalizer": normalizer,
        "gate_count": len(gate_identifiers),
    }


def authenticated_rows(manifest: dict) -> tuple[dict, list[dict]]:
    inventory_path = Path(manifest["inventory"])
    pool_path = Path(manifest["pool"])
    if sha256(inventory_path) != manifest["inventory_sha256"]:
        raise AssertionError("inventory SHA-256 mismatch")
    if sha256(pool_path) != manifest["pool_sha256"]:
        raise AssertionError("pool SHA-256 mismatch")
    artifact = json.loads(inventory_path.read_text())
    pool = json.loads(pool_path.read_text())
    if artifact.get("schema") != "reverse_group_arithmetic_line_inventory_v3":
        raise AssertionError("inventory schema mismatch")
    if artifact["source_sha256"] != manifest["pool_sha256"]:
        raise AssertionError("inventory source hash mismatch")
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
        if rebuilt[key] != artifact[key]:
            raise AssertionError(f"inventory {key} rebuild mismatch")

    rows = []
    seen_primes = set()
    for line_type in rebuilt["line_types"]:
        if int(line_type["target_modulus"]) != 1:
            raise AssertionError("QCIR contains a restricted target row")
        h = int(line_type["h"])
        if h <= 1:
            raise AssertionError("QCIR row modulus is not above one")
        for raw_prime in line_type["primes"]:
            prime = int(raw_prime)
            if prime in seen_primes:
                raise AssertionError("QCIR inventory reuses a prime")
            seen_primes.add(prime)
            rows.append(
                {
                    "p": prime,
                    "h": h,
                    "a": int(line_type["a"]) % h,
                    "b": int(line_type["b"]) % h,
                }
            )
    rows.sort(key=lambda row: (row["h"], row["p"]))
    if rows != manifest["rows"]:
        raise AssertionError("manifest row expansion mismatch")
    return rebuilt, rows


def bit_assignment(
    normalizer: Normalizer,
    bits: list[int],
    value: int,
) -> int | None:
    if not bits:
        return None
    literals = []
    for index, identifier in enumerate(bits):
        variable = normalizer.variable(identifier)
        literals.append(
            variable
            if value & (1 << index)
            else normalizer.negate(variable)
        )
    return normalizer.gate("and", literals)


def expected_expression(
    normalizer: Normalizer,
    rows: list[dict],
    group: dict,
) -> tuple[int, list[int], list[int], dict]:
    next_identifier = 1
    phase_variables = []
    row_components = []
    for row in rows:
        components = []
        for prime, modulus in factors(row["h"]):
            variables = list(
                range(next_identifier, next_identifier + modulus)
            )
            next_identifier += modulus
            phase_variables.extend(variables)
            components.append(
                {
                    "prime": prime,
                    "modulus": modulus,
                    "variables": variables,
                }
            )
        row_components.append(components)

    width = int(group["width"])
    height = int(group["height"])
    group_primes = sorted(
        {
            prime
            for value in (width, height)
            for prime, _component in factors(value)
        }
    )
    coordinate_variables = []
    coordinate_layout = {}
    for prime in group_primes:
        x_modulus = prime ** exponent(width, prime)
        y_modulus = prime ** exponent(height, prime)
        x_count = (x_modulus - 1).bit_length()
        y_count = (y_modulus - 1).bit_length()
        x_bits = list(
            range(next_identifier, next_identifier + x_count)
        )
        next_identifier += x_count
        y_bits = list(
            range(next_identifier, next_identifier + y_count)
        )
        next_identifier += y_count
        coordinate_variables.extend(x_bits)
        coordinate_variables.extend(y_bits)
        coordinate_layout[prime] = {
            "x_modulus": x_modulus,
            "y_modulus": y_modulus,
            "x_bits": x_bits,
            "y_bits": y_bits,
        }

    phase_constraints = []
    pair_count = 0
    for components in row_components:
        for component in components:
            variables = component["variables"]
            phase_constraints.append(
                normalizer.gate(
                    "or",
                    [
                        normalizer.variable(identifier)
                        for identifier in variables
                    ],
                )
            )
            for left_index, left in enumerate(variables):
                for right in variables[left_index + 1 :]:
                    phase_constraints.append(
                        normalizer.gate(
                            "or",
                            [
                                normalizer.negate(
                                    normalizer.variable(left)
                                ),
                                normalizer.negate(
                                    normalizer.variable(right)
                                ),
                            ],
                        )
                    )
                    pair_count += 1
    phase_valid = normalizer.gate("and", phase_constraints)

    point_expressions = {}
    local_valid = {}
    for prime in group_primes:
        layout = coordinate_layout[prime]
        points = []
        for x in range(layout["x_modulus"]):
            x_expression = bit_assignment(
                normalizer,
                layout["x_bits"],
                x,
            )
            for y in range(layout["y_modulus"]):
                y_expression = bit_assignment(
                    normalizer,
                    layout["y_bits"],
                    y,
                )
                arguments = [
                    expression
                    for expression in (x_expression, y_expression)
                    if expression is not None
                ]
                point = normalizer.gate("and", arguments)
                point_expressions[(prime, x, y)] = point
                points.append(point)
        local_valid[prime] = normalizer.gate("or", points)

    bucket_cache = {}
    row_matches = []
    for row, components in zip(rows, row_components, strict=True):
        component_matches = []
        for component in components:
            prime = component["prime"]
            modulus = component["modulus"]
            key = (
                prime,
                modulus,
                row["a"] % modulus,
                row["b"] % modulus,
            )
            buckets = bucket_cache.get(key)
            if buckets is None:
                layout = coordinate_layout[prime]
                bucket_points = [[] for _ in range(modulus)]
                for x in range(layout["x_modulus"]):
                    for y in range(layout["y_modulus"]):
                        residue = (
                            row["a"] * x + row["b"] * y
                        ) % modulus
                        bucket_points[residue].append(
                            point_expressions[(prime, x, y)]
                        )
                if any(not points for points in bucket_points):
                    raise AssertionError(
                        "expected local map is not surjective"
                    )
                buckets = [
                    normalizer.gate("or", points)
                    for points in bucket_points
                ]
                bucket_cache[key] = buckets
            alternatives = [
                normalizer.gate(
                    "and",
                    [
                        normalizer.variable(variable),
                        buckets[residue],
                    ],
                )
                for residue, variable in enumerate(
                    component["variables"]
                )
            ]
            component_matches.append(
                normalizer.gate("or", alternatives)
            )
        row_matches.append(
            normalizer.gate("and", component_matches)
        )

    coverage = normalizer.gate("or", row_matches)
    point_valid = normalizer.gate(
        "and",
        list(local_valid.values()),
    )
    implication = normalizer.gate(
        "or",
        [normalizer.negate(point_valid), coverage],
    )
    output = normalizer.gate(
        "and",
        [phase_valid, implication],
    )
    metadata = {
        "phase_component_count": sum(
            len(components) for components in row_components
        ),
        "phase_variable_count": len(phase_variables),
        "coordinate_variable_count": len(coordinate_variables),
        "at_most_one_clause_count": pair_count,
    }
    return output, phase_variables, coordinate_variables, metadata


def verify(manifest_path: Path) -> dict:
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema") != "crt_affine_cover_qcir_manifest_v1":
        raise AssertionError("QCIR manifest schema mismatch")
    qcir_path = Path(manifest["qcir"])
    if sha256(qcir_path) != manifest["qcir_sha256"]:
        raise AssertionError("QCIR SHA-256 mismatch")
    artifact, rows = authenticated_rows(manifest)
    parsed = parse_qcir(qcir_path)
    expected, phase_variables, coordinate_variables, metadata = (
        expected_expression(
            parsed["normalizer"],
            rows,
            artifact["group"],
        )
    )
    if parsed["existential"] != phase_variables:
        raise AssertionError("QCIR existential prefix mismatch")
    if parsed["universal"] != coordinate_variables:
        raise AssertionError("QCIR universal prefix mismatch")
    if parsed["output_expression"] != expected:
        raise AssertionError("QCIR normalized formula mismatch")
    for key, value in metadata.items():
        if int(manifest[key]) != value:
            raise AssertionError(f"QCIR manifest {key} mismatch")
    if (
        int(manifest["identifier_count"])
        != parsed["declared_identifier_count"]
        or int(manifest["gate_count"]) != parsed["gate_count"]
        or int(manifest["output_gate"]) != parsed["output"]
    ):
        raise AssertionError("QCIR identifier metadata mismatch")
    claim = manifest["claim"]
    if (
        claim.get("qbf_result") != "UNRUN"
        or claim.get("finite_group_cover_proved") is not False
        or claim.get("finite_group_noncover_proved") is not False
        or claim.get("integer_m_found") is not False
    ):
        raise AssertionError("unrun QCIR claims were promoted")
    return {
        "verified": True,
        "inventory_sha256": manifest["inventory_sha256"],
        "pool_sha256": manifest["pool_sha256"],
        "qcir_sha256": manifest["qcir_sha256"],
        "row_count": len(rows),
        **metadata,
        "identifier_count": parsed["declared_identifier_count"],
        "gate_count": parsed["gate_count"],
        "normalized_formula_matches": True,
        "qbf_result": "UNRUN",
        "finite_group_cover_proved": False,
        "finite_group_noncover_proved": False,
        "integer_m_found": False,
        "engine": "independent-normalized-Boolean-DAG-reconstruction",
        "scope": (
            "exact encoding verification for the authenticated declared "
            "finite quotient; no QBF solver result was checked"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = verify(args.manifest)
    report["manifest"] = str(args.manifest)
    report["manifest_sha256"] = sha256(args.manifest)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"verified={report['verified']} "
        f"rows={report['row_count']} "
        f"phase_variables={report['phase_variable_count']} "
        f"coordinate_variables={report['coordinate_variable_count']} "
        f"qbf_result={report['qbf_result']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
