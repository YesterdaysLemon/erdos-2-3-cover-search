#!/usr/bin/env python3
"""Audit a two-level phase-coordinate quotient with one incremental SAT model.

The low anchors are exhaustively fixed to each legal target tuple.  The repair
anchors remain free.  Exact exponent pairs are first reduced to their joint
low/repair target fingerprint after algebraic and frozen-row coverage is
discharged.  A single SAT instance is then queried under assumptions for every
low branch.

UNSAT outcomes are exact only for the supplied finite fingerprints and the
embedded frozen remainder.  They are not a global obstruction for the Erdos
problem.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

from anchor_phase_master import normalize_rows
from certify_anchor_phase_quotient import load_cells


Fingerprint = tuple[tuple[int, ...], tuple[int, ...]]


def legal_targets(row: tuple[int, int, int, int, int, int]) -> tuple[int, ...]:
    _prime, h, _a, _b, modulus, residue = row
    return tuple(
        target for target in range(h) if target % modulus == residue
    )


def point_fingerprints(
    rows: list[dict],
    base_phases: dict[int, int],
    low_anchors: tuple[int, ...],
    repair_anchors: tuple[int, ...],
    points: list[tuple[int, int]],
    algebraic_primes: tuple[int, ...] = (),
) -> tuple[set[Fingerprint], dict]:
    """Reduce exact points to fingerprints not covered by the frozen remainder."""
    if set(low_anchors) & set(repair_anchors):
        raise ValueError("low and repair anchors must be disjoint")
    free_primes = set(low_anchors) | set(repair_anchors)
    compact_rows, phases = normalize_rows(rows, base_phases)
    row_by_prime = {row[0]: row for row in compact_rows}
    missing = free_primes - row_by_prime.keys()
    if missing:
        raise ValueError(f"anchor primes are absent from the row pool: {missing}")
    frozen_rows = [row for row in compact_rows if row[0] not in free_primes]

    fingerprints: set[Fingerprint] = set()
    algebraically_discharged = 0
    frozen_discharged = 0
    started = time.monotonic()
    for k, l in points:
        if k < 0 or l < 0:
            raise ValueError("exponent pairs must be nonnegative")
        if any(
            k % prime == 0 and l % prime == 0
            for prime in algebraic_primes
        ):
            algebraically_discharged += 1
            continue
        if any(
            (a * k + b * l - phases[prime]) % h == 0
            for prime, h, a, b, _modulus, _residue in frozen_rows
        ):
            frozen_discharged += 1
            continue

        def targets(primes: tuple[int, ...]) -> tuple[int, ...]:
            return tuple(
                (row_by_prime[prime][2] * k + row_by_prime[prime][3] * l)
                % row_by_prime[prime][1]
                for prime in primes
            )

        fingerprints.add((targets(low_anchors), targets(repair_anchors)))

    return fingerprints, {
        "input_points": len(points),
        "distinct_input_points": len(set(points)),
        "algebraically_discharged_points": algebraically_discharged,
        "frozen_discharged_points": frozen_discharged,
        "eligible_fingerprints": len(fingerprints),
        "filter_seconds": time.monotonic() - started,
    }


def validate_fingerprints(
    fingerprints: set[Fingerprint],
    low_rows: tuple[tuple[int, int, int, int, int, int], ...],
    repair_rows: tuple[tuple[int, int, int, int, int, int], ...],
) -> None:
    for low, repair in fingerprints:
        if len(low) != len(low_rows) or len(repair) != len(repair_rows):
            raise ValueError("fingerprint width does not match the anchors")
        for target, row in zip(low, low_rows):
            if not 0 <= target < row[1]:
                raise ValueError("low fingerprint target is out of range")
        for target, row in zip(repair, repair_rows):
            if not 0 <= target < row[1]:
                raise ValueError("repair fingerprint target is out of range")


def binary_refinement_symmetries(
    low_rows: tuple[tuple[int, int, int, int, int, int], ...],
) -> list[dict]:
    """Find exact binary target relations that collapse low-branch cases.

    If a child linear form reduces to a binary parent's form modulo two,
    every point surviving parent phase ``s`` has child target parity
    ``1-s``.  All child phases of parity ``s`` are therefore inactive and
    mutually indistinguishable on that residual problem.
    """
    relations = []
    for parent_index, parent in enumerate(low_rows):
        parent_prime, parent_h, parent_a, parent_b, _modulus, _residue = (
            parent
        )
        if parent_h != 2 or legal_targets(parent) != (0, 1):
            continue
        for child_index, child in enumerate(low_rows):
            if parent_index == child_index:
                continue
            child_prime, child_h, child_a, child_b, _modulus, _residue = (
                child
            )
            child_options = legal_targets(child)
            if (
                child_h % 2
                or child_a % 2 != parent_a % 2
                or child_b % 2 != parent_b % 2
            ):
                continue
            quotient_classes = 0
            inactive_by_parent = {}
            for parent_target in (0, 1):
                inactive = tuple(
                    target
                    for target in child_options
                    if target % 2 == parent_target
                )
                active = tuple(
                    target
                    for target in child_options
                    if target % 2 != parent_target
                )
                quotient_classes += len(active) + bool(inactive)
                inactive_by_parent[str(parent_target)] = list(inactive)
            relations.append(
                {
                    "parent_prime": parent_prime,
                    "child_prime": child_prime,
                    "parent_form_mod_2": [
                        parent_a % 2,
                        parent_b % 2,
                    ],
                    "child_form_mod_2": [
                        child_a % 2,
                        child_b % 2,
                    ],
                    "raw_local_branches": 2 * len(child_options),
                    "quotient_local_classes": quotient_classes,
                    "inactive_child_targets_by_parent_phase": (
                        inactive_by_parent
                    ),
                }
            )
    return relations


def audit_fingerprints(
    rows: list[dict],
    base_phases: dict[int, int],
    low_anchors: tuple[int, ...],
    repair_anchors: tuple[int, ...],
    fingerprints: set[Fingerprint],
    solver_name: str = "cadical195",
) -> dict:
    """Solve every low branch by assumptions in one exact SAT instance."""
    compact_rows, _phases = normalize_rows(rows, base_phases)
    row_by_prime = {row[0]: row for row in compact_rows}
    if set(low_anchors) & set(repair_anchors):
        raise ValueError("low and repair anchors must be disjoint")
    if not low_anchors or not repair_anchors:
        raise ValueError("both anchor levels must be nonempty")
    if (set(low_anchors) | set(repair_anchors)) - row_by_prime.keys():
        raise ValueError("an anchor prime is absent from the row pool")
    low_rows = tuple(row_by_prime[prime] for prime in low_anchors)
    repair_rows = tuple(row_by_prime[prime] for prime in repair_anchors)
    validate_fingerprints(fingerprints, low_rows, repair_rows)
    refinement_symmetries = binary_refinement_symmetries(low_rows)

    dependency_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dependency_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    started = time.monotonic()
    variable_pool = IDPool()
    solver = Solver(name=solver_name)
    target_variable: dict[tuple[int, int], int] = {}
    option_sets: dict[int, tuple[int, ...]] = {}
    for row in low_rows + repair_rows:
        prime = row[0]
        options = legal_targets(row)
        option_sets[prime] = options
        variables = []
        for target in options:
            variable = variable_pool.id(("target", prime, target))
            target_variable[prime, target] = variable
            variables.append(variable)
        encoding = CardEnc.equals(
            variables,
            bound=1,
            vpool=variable_pool,
            encoding=EncType.seqcounter,
        )
        solver.append_formula(encoding.clauses)

    ordered_fingerprints = sorted(fingerprints)
    for low_targets, repair_targets in ordered_fingerprints:
        clause = []
        for prime, target in zip(low_anchors, low_targets):
            variable = target_variable.get((prime, target))
            if variable is not None:
                clause.append(variable)
        for prime, target in zip(repair_anchors, repair_targets):
            variable = target_variable.get((prime, target))
            if variable is not None:
                clause.append(variable)
        solver.add_clause(clause)

    outcomes = []
    residual_classes: Counter[str] = Counter()
    low_options = [option_sets[prime] for prime in low_anchors]
    for branch in itertools.product(*low_options):
        assumptions = [
            target_variable[prime, target]
            for prime, target in zip(low_anchors, branch)
        ]
        residual_indices = tuple(
            index
            for index, (low_targets, _repair_targets) in enumerate(
                ordered_fingerprints
            )
            if not any(
                target == point_target
                for target, point_target in zip(branch, low_targets)
            )
        )
        residual_digest = hashlib.sha256(
            ",".join(map(str, residual_indices)).encode()
        ).hexdigest()
        residual_classes[residual_digest] += 1
        sat = solver.solve(assumptions=assumptions)
        selected_targets = None
        if sat:
            model = {literal for literal in solver.get_model() if literal > 0}
            selected_targets = []
            for prime in repair_anchors:
                selected = [
                    target
                    for target in option_sets[prime]
                    if target_variable[prime, target] in model
                ]
                if len(selected) != 1:
                    raise AssertionError(
                        "SAT model did not select one repair target"
                    )
                selected_targets.append(selected[0])
        outcomes.append(
            {
                "branch": list(branch),
                "sat": sat,
                "residual_fingerprints": len(residual_indices),
                "residual_digest": residual_digest,
                "repair_targets": selected_targets,
            }
        )

    solver_clauses = solver.nof_clauses()
    solver_variables = variable_pool.top
    solver.delete()
    sat_outcomes = [outcome for outcome in outcomes if outcome["sat"]]
    top_sat = sorted(
        sat_outcomes,
        key=lambda outcome: outcome["residual_fingerprints"],
        reverse=True,
    )[:10]
    return {
        "scope_warning": (
            "Exact only for the supplied fingerprints and frozen remainder; "
            "not a global no-cover theorem."
        ),
        "low_anchors": list(low_anchors),
        "repair_anchors": list(repair_anchors),
        "binary_refinement_symmetries": refinement_symmetries,
        "summary": {
            "fingerprints": len(fingerprints),
            "low_branches": len(outcomes),
            "unsat_branches": len(outcomes) - len(sat_outcomes),
            "sat_branches": len(sat_outcomes),
            "residual_equivalence_classes": len(residual_classes),
            "largest_residual_class": max(residual_classes.values()),
            "solver_variables": solver_variables,
            "fingerprint_clauses": len(fingerprints),
            "solver_reported_clauses_after_queries": solver_clauses,
            "audit_seconds": time.monotonic() - started,
            "top_sat": top_sat,
        },
        "residual_class_sizes": sorted(
            residual_classes.values(), reverse=True
        ),
        "outcomes": outcomes,
    }


def parse_primes(value: str) -> tuple[int, ...]:
    return tuple(int(item) for item in value.split(",") if item)


def load_fingerprint_cache(path: Path | None) -> set[Fingerprint]:
    if path is None:
        return set()
    payload = json.loads(path.read_text())
    return {
        (tuple(int(value) for value in low), tuple(int(value) for value in repair))
        for low, repair in payload
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--base-phases", type=Path, required=True)
    parser.add_argument("--low-anchors", required=True)
    parser.add_argument("--repair-anchors", required=True)
    parser.add_argument("--fingerprint-input", type=Path)
    parser.add_argument("--points-file", type=Path, action="append", default=[])
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--fingerprint-output", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(args.base_phases.read_text()).items()
    }
    low_anchors = parse_primes(args.low_anchors)
    repair_anchors = parse_primes(args.repair_anchors)
    fingerprints = load_fingerprint_cache(args.fingerprint_input)
    initial_count = len(fingerprints)
    points = []
    for path in args.points_file:
        points.extend(load_cells(path))
    new_fingerprints, filter_metadata = point_fingerprints(
        rows,
        phases,
        low_anchors,
        repair_anchors,
        points,
        tuple(int(prime) for prime in payload.get("algebraic_primes", ())),
    )
    fingerprints.update(new_fingerprints)
    audit = audit_fingerprints(
        rows,
        phases,
        low_anchors,
        repair_anchors,
        fingerprints,
        args.solver,
    )
    audit["source"] = {
        "pool": str(args.pool),
        "base_phases": str(args.base_phases),
        "fingerprint_input": (
            str(args.fingerprint_input) if args.fingerprint_input else None
        ),
        "points_files": [str(path) for path in args.points_file],
    }
    audit["update"] = {
        "input_fingerprints": initial_count,
        "new_fingerprints": len(fingerprints) - initial_count,
        **filter_metadata,
    }
    args.fingerprint_output.write_text(
        json.dumps(
            [
                [list(low), list(repair)]
                for low, repair in sorted(fingerprints)
            ]
        )
        + "\n"
    )
    args.audit_output.write_text(json.dumps(audit, indent=2) + "\n")
    summary = audit["summary"]
    print(
        f"TWO_LEVEL_AUDIT fingerprints={summary['fingerprints']} "
        f"unsat={summary['unsat_branches']}/{summary['low_branches']} "
        f"classes={summary['residual_equivalence_classes']} "
        f"solve_s={summary['audit_seconds']:.6f} "
        f"output={args.audit_output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
