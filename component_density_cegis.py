#!/usr/bin/env python3
"""CEGIS for the necessary residual-density condition on one CRT component.

Fix every CRT coordinate except a selected prime-power component q^e.  Rows
compatible with the fixed coordinates restrict to affine lines of density
1/q^j on the residual plane.  Every genuine global cover must have total
residual density at least one in every such coarse cell.

The checker finds exact low-density cells with Z3.  A conservative PySAT
master retargets row phases so all accumulated cells meet the density bound.
SAT here is only a necessary-condition result, never a cover certificate.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import exact_uncovered


def target_modulus(row: dict) -> int:
    """Return the row's phase restriction, defaulting to unrestricted."""
    return int(row.get("target_modulus", 1))


def q_part(value: int, prime: int) -> int:
    result = 1
    while value % prime == 0:
        result *= prime
        value //= prime
    return result


def find_low_density_cells(
    rows: list[dict],
    phases: dict[int, int],
    prime: int,
    limit: int,
    max_component: int,
    diversity_modulus: int = 0,
) -> tuple[list[tuple[int, int]], dict]:
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import z3  # type: ignore

    maximal: dict[int, int] = {}
    row_data = []
    residual_modulus = 1
    for row in rows:
        h = int(row["h"])
        residual = q_part(h, prime)
        residual_modulus = max(residual_modulus, residual)
        other = h // residual
        factors = exact_uncovered.factor(other)
        for component_prime, exponent in factors.items():
            maximal[component_prime] = max(
                maximal.get(component_prime, 0),
                exponent,
            )
        row_data.append((other, residual, factors))
    components = {
        component_prime: component_prime**exponent
        for component_prime, exponent in maximal.items()
    }
    diversity_component = None
    if diversity_modulus:
        diversity_factors = exact_uncovered.factor(diversity_modulus)
        if len(diversity_factors) != 1:
            raise ValueError("density diversity modulus must be a prime power")
        diversity_prime = next(iter(diversity_factors))
        if (
            diversity_prime not in components
            or components[diversity_prime] % diversity_modulus
        ):
            raise ValueError(
                "density diversity modulus is absent from the coarse domain"
            )
        diversity_component = diversity_prime
    largest = max(components.values(), default=1)
    if largest > max_component:
        raise ValueError(
            f"largest coarse component {largest} exceeds guard "
            f"{max_component}"
        )

    expression_width = max(8, 2 * largest.bit_length() + 2)
    solver = z3.SolverFor("QF_BV")
    kval = {}
    lval = {}
    widths = {}
    for component_prime, modulus in components.items():
        width = max(1, (modulus - 1).bit_length())
        widths[component_prime] = width
        kval[component_prime] = z3.BitVec(
            f"k_{component_prime}", width
        )
        lval[component_prime] = z3.BitVec(
            f"l_{component_prime}", width
        )
        if modulus != 1 << width:
            solver.add(
                z3.ULT(
                    kval[component_prime],
                    z3.BitVecVal(modulus, width),
                )
            )
            solver.add(
                z3.ULT(
                    lval[component_prime],
                    z3.BitVecVal(modulus, width),
                )
            )

    def widened(value, component_prime):
        return z3.ZeroExt(
            expression_width - widths[component_prime],
            value,
        )

    activations = []
    weights = []
    for row, (other, residual, factors) in zip(rows, row_data):
        weight = residual_modulus // residual
        weights.append(weight)
        if other == 1:
            activations.append(z3.BoolVal(True))
            continue
        equalities = []
        target = phases[int(row["p"])] % other
        for component_prime, exponent in factors.items():
            modulus = component_prime**exponent
            kr = z3.URem(
                widened(kval[component_prime], component_prime),
                z3.BitVecVal(modulus, expression_width),
            )
            lr = z3.URem(
                widened(lval[component_prime], component_prime),
                z3.BitVecVal(modulus, expression_width),
            )
            expression = (
                z3.BitVecVal(
                    int(row["a"]) % modulus,
                    expression_width,
                )
                * kr
                + z3.BitVecVal(
                    int(row["b"]) % modulus,
                    expression_width,
                )
                * lr
                + z3.BitVecVal(
                    (-target) % modulus,
                    expression_width,
                )
            )
            equalities.append(
                z3.URem(
                    expression,
                    z3.BitVecVal(modulus, expression_width),
                )
                == 0
            )
        activations.append(z3.And(*equalities))

    solver.add(
        z3.PbLe(
            list(zip(activations, weights)),
            residual_modulus - 1,
        )
    )
    cells = []
    scaled_densities = []
    for _ in range(limit):
        if solver.check() != z3.sat:
            break
        model = solver.model()
        kres = []
        lres = []
        blocks = []
        for component_prime, modulus in components.items():
            kr = model.eval(
                kval[component_prime],
                model_completion=True,
            ).as_long()
            lr = model.eval(
                lval[component_prime],
                model_completion=True,
            ).as_long()
            kres.append((kr, modulus))
            lres.append((lr, modulus))
            blocks.extend(
                (
                    kval[component_prime] != kr,
                    lval[component_prime] != lr,
                )
            )
        k = exact_uncovered.crt(kres)
        l = exact_uncovered.crt(lres)
        scaled_density = 0
        for row, (other, residual, _factors) in zip(rows, row_data):
            if (
                int(row["a"]) * k
                + int(row["b"]) * l
                - phases[int(row["p"])]
            ) % other == 0:
                scaled_density += residual_modulus // residual
        if not scaled_density < residual_modulus:
            raise AssertionError("checker model does not violate density")
        cells.append((k, l))
        scaled_densities.append(scaled_density)
        if not blocks:
            break
        if diversity_component is None:
            solver.add(z3.Or(*blocks))
        else:
            width = expression_width
            modulus_value = z3.BitVecVal(diversity_modulus, width)
            solver.add(
                z3.Or(
                    z3.URem(
                        widened(
                            kval[diversity_component],
                            diversity_component,
                        ),
                        modulus_value,
                    )
                    != z3.BitVecVal(
                        k % diversity_modulus,
                        width,
                    ),
                    z3.URem(
                        widened(
                            lval[diversity_component],
                            diversity_component,
                        ),
                        modulus_value,
                    )
                    != z3.BitVecVal(
                        l % diversity_modulus,
                        width,
                    ),
                )
            )

    return cells, {
        "prime": prime,
        "residual_modulus": residual_modulus,
        "components": components,
        "period": math.prod(components.values()),
        "row_count": len(rows),
        "cell_count": len(cells),
        "scaled_densities": scaled_densities,
        "minimum_scaled_density": (
            min(scaled_densities) if scaled_densities else None
        ),
        "diversity_modulus": diversity_modulus,
        "sat": bool(cells),
        "engine": "z3-bitvector-pb",
    }


def solve_master(
    rows: list[dict],
    initial_phases: dict[int, int],
    prime: int,
    cells: list[tuple[int, int]],
    fixed_targets: dict[int, int],
    solver_name: str,
) -> tuple[dict[int, int] | None, dict]:
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    row_data = []
    residual_modulus = 1
    for row in rows:
        h = int(row["h"])
        residual = q_part(h, prime)
        residual_modulus = max(residual_modulus, residual)
        row_data.append((h // residual, residual))

    options: dict[int, set[int]] = defaultdict(set)
    required_by_cell = []
    for k, l in cells:
        required = []
        for row_index, (row, (other, _residual)) in enumerate(
            zip(rows, row_data)
        ):
            target = (
                int(row["a"]) * k + int(row["b"]) * l
            ) % other
            required.append(target)
            if other > 1 and int(row["p"]) not in fixed_targets:
                options[row_index].add(target)
        required_by_cell.append(required)

    vpool = IDPool()
    variable = {}
    variable_target = {}
    solver = Solver(name=solver_name)
    clause_count = 0
    for row_index, targets in options.items():
        literals = []
        for target in sorted(targets):
            literal = vpool.id(("coarse", row_index, target))
            variable[(row_index, target)] = literal
            variable_target[literal] = (row_index, target)
            literals.append(literal)
        if len(literals) > 1:
            encoding = CardEnc.atmost(
                literals,
                bound=1,
                vpool=vpool,
                encoding=EncType.seqcounter,
            )
            solver.append_formula(encoding.clauses)
            clause_count += len(encoding.clauses)

    cut_count = 0
    for required in required_by_cell:
        literals = []
        weights = []
        constant = 0
        for row_index, (row, (other, residual)) in enumerate(
            zip(rows, row_data)
        ):
            weight = residual_modulus // residual
            row_prime = int(row["p"])
            if other == 1:
                constant += weight
            elif row_prime in fixed_targets:
                if fixed_targets[row_prime] % other == required[row_index]:
                    constant += weight
            else:
                literals.append(
                    variable[(row_index, required[row_index])]
                )
                weights.append(weight)
        bound = residual_modulus - constant
        if bound <= 0:
            continue
        if sum(weights) < bound:
            solver.delete()
            return None, {
                "sat": False,
                "reason": "density cut lacks enough compatible row weight",
                "cells": len(cells),
            }
        if residual_modulus == prime and all(
            weight in (1, prime) for weight in weights
        ):
            # For a first-power residual component, every row either has
            # density 1 (and satisfies the cut by itself) or density 1/q.
            # Encode
            #
            #   any(full-density literal) OR at-least(bound, line literals)
            #
            # by guarding an ordinary cardinality encoding.  This is far
            # smaller than repeating each full-density literal q times.
            full_literals = [
                literal
                for literal, weight in zip(literals, weights)
                if weight == prime
            ]
            line_literals = [
                literal
                for literal, weight in zip(literals, weights)
                if weight == 1
            ]
            full_guard = None
            if full_literals:
                full_guard = vpool.id(("full_guard", cut_count))
                for literal in full_literals:
                    solver.add_clause([-literal, full_guard])
                    clause_count += 1
                solver.add_clause([-full_guard, *full_literals])
                clause_count += 1
            if len(line_literals) < bound:
                if full_guard is None:
                    solver.add_clause([])
                else:
                    solver.add_clause([full_guard])
                clause_count += 1
            else:
                encoding = CardEnc.atleast(
                    line_literals,
                    bound=bound,
                    vpool=vpool,
                    encoding=EncType.seqcounter,
                )
                for clause in encoding.clauses:
                    if full_guard is None:
                        solver.add_clause(clause)
                    else:
                        solver.add_clause([full_guard, *clause])
                clause_count += len(encoding.clauses)
        else:
            # General dependency-free fallback: repeating a literal by its
            # positive integer weight is an exact pseudo-Boolean reduction.
            expanded_literals = [
                literal
                for literal, weight in zip(literals, weights)
                for _ in range(weight)
            ]
            encoding = CardEnc.atleast(
                expanded_literals,
                bound=bound,
                vpool=vpool,
                encoding=EncType.seqcounter,
            )
            solver.append_formula(encoding.clauses)
            clause_count += len(encoding.clauses)
        cut_count += 1

    hints = []
    for (row_index, target), literal in variable.items():
        row_prime = int(rows[row_index]["p"])
        other, _residual = row_data[row_index]
        if initial_phases[row_prime] % other == target:
            hints.append(literal)
    try:
        solver.set_phases(hints)
    except NotImplementedError:
        hints = []

    started = time.monotonic()
    sat = solver.solve()
    solve_seconds = time.monotonic() - started
    result_meta = {
        "sat": sat,
        "cells": len(cells),
        "cuts": cut_count,
        "options": len(variable),
        "variables": vpool.top,
        "clauses": clause_count,
        "hints": len(hints),
        "solver": solver_name,
        "solve_seconds": solve_seconds,
    }
    if not sat:
        solver.delete()
        return None, result_meta

    model = {literal for literal in solver.get_model() if literal > 0}
    selected = {}
    for literal in model & variable_target.keys():
        row_index, target = variable_target[literal]
        if row_index in selected:
            raise AssertionError("master selected two coarse targets")
        selected[row_index] = target

    phases = dict(initial_phases)
    for row_index, target in selected.items():
        row = rows[row_index]
        row_prime = int(row["p"])
        h = int(row["h"])
        other, residual = row_data[row_index]
        if other == 1:
            continue
        if residual == 1:
            phases[row_prime] = target
        else:
            phases[row_prime] = exact_uncovered.crt(
                [
                    (target, other),
                    (phases[row_prime] % residual, residual),
                ]
            ) % h
    phases.update(fixed_targets)
    solver.delete()
    return phases, result_meta


def solve_master_z3(
    rows: list[dict],
    initial_phases: dict[int, int],
    prime: int,
    cells: list[tuple[int, int]],
    fixed_targets: dict[int, int],
    maximize_preserved: bool = True,
) -> tuple[dict[int, int] | None, dict]:
    """Solve the accumulated weighted density cuts without weight expansion.

    The PySAT fallback above is useful for small residual components, but its
    generic pseudo-Boolean reduction repeats a literal once per unit of
    weight.  For a residual component such as 2^13 this can turn one 55-row
    cut into a huge sequential counter.  Z3 can retain the same exact
    integer-weighted inequality directly.

    Each row target modulo ``h`` decomposes independently into its
    ``prime``-power part and the coprime ``other`` part.  A density cut on the
    selected residual component depends only on the latter, so one bounded
    integer variable per non-fixed row is an exact target representation.
    """

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import z3  # type: ignore

    row_data = []
    residual_modulus = 1
    for row in rows:
        h = int(row["h"])
        residual = q_part(h, prime)
        residual_modulus = max(residual_modulus, residual)
        row_data.append((h // residual, residual))

    # Preserve as much of the previous exact phase assignment as possible.
    # A merely arbitrary satisfying model tends to repair the accumulated
    # cells by moving many unrelated rows, which creates a fresh family of
    # low-density cells at the next checker call.  The Hamming objective does
    # not change feasibility or the meaning of UNSAT.
    solver = z3.SolverFor("QF_LIA")
    coarse_targets = {}
    preservation_terms = []
    for row_index, (row, (other, _residual)) in enumerate(
        zip(rows, row_data)
    ):
        row_prime = int(row["p"])
        if other == 1 or row_prime in fixed_targets:
            continue
        variable = z3.Int(f"coarse_target_{row_index}")
        coarse_targets[row_index] = variable
        solver.add(variable >= 0, variable < other)
        preservation_terms.append(
            z3.If(
                variable == initial_phases[row_prime] % other,
                1,
                0,
            )
        )

    cut_count = 0
    impossible_cell = None
    for k, l in cells:
        constant = 0
        weighted_matches = []
        available_weight = 0
        for row_index, (row, (other, residual)) in enumerate(
            zip(rows, row_data)
        ):
            weight = residual_modulus // residual
            required = (
                int(row["a"]) * k + int(row["b"]) * l
            ) % other
            row_prime = int(row["p"])
            if other == 1:
                constant += weight
            elif row_prime in fixed_targets:
                if fixed_targets[row_prime] % other == required:
                    constant += weight
            else:
                available_weight += weight
                weighted_matches.append(
                    z3.If(
                        coarse_targets[row_index] == required,
                        weight,
                        0,
                    )
                )
        bound = residual_modulus - constant
        if bound <= 0:
            continue
        if available_weight < bound:
            impossible_cell = (k, l)
            break
        solver.add(z3.Sum(*weighted_matches) >= bound)
        cut_count += 1

    if impossible_cell is not None:
        return None, {
            "sat": False,
            "reason": "density cut lacks enough compatible row weight",
            "cells": len(cells),
            "cuts": cut_count,
            "impossible_cell": list(impossible_cell),
            "engine": "z3-qf-lia-pb",
        }

    started = time.monotonic()
    status = solver.check()
    model = solver.model() if status == z3.sat else None
    maximum_preserved = None
    if status == z3.sat and preservation_terms and maximize_preserved:
        preserved_sum = z3.Sum(*preservation_terms)
        low = 0
        high = len(preservation_terms)
        best_model = model
        while low < high:
            midpoint = (low + high + 1) // 2
            solver.push()
            solver.add(preserved_sum >= midpoint)
            trial = solver.check()
            if trial == z3.sat:
                low = midpoint
                best_model = solver.model()
            elif trial == z3.unsat:
                high = midpoint - 1
            else:
                reason = solver.reason_unknown()
                solver.pop()
                raise RuntimeError(
                    f"Z3 density Hamming search returned unknown: {reason}"
                )
            solver.pop()
        maximum_preserved = low
        model = best_model
    solve_seconds = time.monotonic() - started
    result_meta = {
        "sat": status == z3.sat,
        "cells": len(cells),
        "cuts": cut_count,
        "variables": len(coarse_targets),
        "solver": "z3",
        "engine": "z3-binarymax-qf-lia-pb",
        "solve_seconds": solve_seconds,
        "maximum_preserved": maximum_preserved,
        "maximize_preserved": maximize_preserved,
    }
    if status == z3.unknown:
        result_meta["reason"] = solver.reason_unknown()
        raise RuntimeError(
            f"Z3 density master returned unknown: {solver.reason_unknown()}"
        )
    if status == z3.unsat:
        return None, result_meta

    assert model is not None
    result_meta["preserved_targets"] = sum(
        1
        for row_index, variable in coarse_targets.items()
        if model.eval(variable, model_completion=True).as_long()
        == initial_phases[int(rows[row_index]["p"])] % row_data[row_index][0]
    )
    phases = dict(initial_phases)
    for row_index, variable in coarse_targets.items():
        row = rows[row_index]
        row_prime = int(row["p"])
        h = int(row["h"])
        other, residual = row_data[row_index]
        target = model.eval(variable, model_completion=True).as_long()
        if residual == 1:
            phases[row_prime] = target
        else:
            phases[row_prime] = exact_uncovered.crt(
                [
                    (target, other),
                    (phases[row_prime] % residual, residual),
                ]
            ) % h
    phases.update(fixed_targets)
    return phases, result_meta


def solve_master_milp(
    rows: list[dict],
    initial_phases: dict[int, int],
    prime: int,
    cells: list[tuple[int, int]],
    fixed_targets: dict[int, int],
    time_limit: float,
    maximize_preserved: bool = True,
) -> tuple[dict[int, int] | None, dict]:
    """Solve exact weighted density cuts as a sparse binary MILP."""

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from scipy.optimize import Bounds, LinearConstraint, milp  # type: ignore
    from scipy.sparse import coo_matrix  # type: ignore

    row_data = []
    residual_modulus = 1
    for row in rows:
        h = int(row["h"])
        residual = q_part(h, prime)
        residual_modulus = max(residual_modulus, residual)
        row_data.append((h // residual, residual))

    mutable_rows = []
    target_options = {}
    for row_index, (row, (other, _residual)) in enumerate(
        zip(rows, row_data)
    ):
        row_prime = int(row["p"])
        if other == 1 or row_prime in fixed_targets:
            continue
        mutable_rows.append(row_index)
        target_options[row_index] = {
            initial_phases[row_prime] % other
        }

    density_cuts = 0
    impossible_cell = None
    cut_records = []
    for k, l in cells:
        constant = 0
        entries = []
        available_weight = 0
        for row_index, (row, (other, residual)) in enumerate(
            zip(rows, row_data)
        ):
            weight = residual_modulus // residual
            required = (
                int(row["a"]) * k + int(row["b"]) * l
            ) % other
            row_prime = int(row["p"])
            if other == 1:
                constant += weight
            elif row_prime in fixed_targets:
                if fixed_targets[row_prime] % other == required:
                    constant += weight
            else:
                available_weight += weight
                entries.append((row_index, required, float(weight)))
        bound = residual_modulus - constant
        if bound <= 0:
            continue
        if available_weight < bound:
            impossible_cell = (k, l)
            break
        for row_index, required, _coefficient in entries:
            target_options[row_index].add(required)
        cut_records.append((float(bound), entries))
        density_cuts += 1

    if impossible_cell is not None:
        return None, {
            "sat": False,
            "reason": "density cut lacks enough compatible row weight",
            "cells": len(cells),
            "cuts": density_cuts,
            "impossible_cell": list(impossible_cell),
            "engine": "scipy-highs-milp",
        }

    variable = {}
    variable_rows = []
    variable_count = 0
    sentinel_count = 0
    for row_index in mutable_rows:
        other, _residual = row_data[row_index]
        targets = sorted(target_options[row_index])
        target_columns = {}
        indices = []
        for target in targets:
            target_columns[target] = variable_count
            variable[(row_index, target)] = variable_count
            indices.append(variable_count)
            variable_count += 1
        sentinel_column = None
        if len(targets) < other:
            sentinel_column = variable_count
            indices.append(variable_count)
            variable_count += 1
            sentinel_count += 1
        variable_rows.append(
            {
                "row_index": row_index,
                "indices": indices,
                "target_columns": target_columns,
                "sentinel_column": sentinel_column,
                "other": other,
            }
        )

    matrix_rows = []
    matrix_cols = []
    matrix_values = []
    lower_bounds = []
    upper_bounds = []
    constraint_index = 0
    for record in variable_rows:
        for column in record["indices"]:
            matrix_rows.append(constraint_index)
            matrix_cols.append(column)
            matrix_values.append(1.0)
        lower_bounds.append(1.0)
        upper_bounds.append(1.0)
        constraint_index += 1

    for bound, entries in cut_records:
        for row_index, required, coefficient in entries:
            column = variable[(row_index, required)]
            matrix_rows.append(constraint_index)
            matrix_cols.append(column)
            matrix_values.append(coefficient)
        lower_bounds.append(bound)
        upper_bounds.append(np.inf)
        constraint_index += 1

    matrix = coo_matrix(
        (
            np.asarray(matrix_values, dtype=np.float64),
            (
                np.asarray(matrix_rows, dtype=np.int32),
                np.asarray(matrix_cols, dtype=np.int32),
            ),
        ),
        shape=(constraint_index, variable_count),
    ).tocsr()
    constraints = LinearConstraint(
        matrix,
        np.asarray(lower_bounds, dtype=np.float64),
        np.asarray(upper_bounds, dtype=np.float64),
    )
    options = {"presolve": True}
    if time_limit > 0:
        options["time_limit"] = float(time_limit)
    objective = np.zeros(variable_count, dtype=np.float64)
    if maximize_preserved:
        for record in variable_rows:
            row_index = record["row_index"]
            other = record["other"]
            initial_target = (
                initial_phases[int(rows[row_index]["p"])] % other
            )
            objective[record["target_columns"][initial_target]] = -1.0
    started = time.monotonic()
    answer = milp(
        c=objective,
        integrality=np.ones(variable_count, dtype=np.uint8),
        bounds=Bounds(
            np.zeros(variable_count, dtype=np.float64),
            np.ones(variable_count, dtype=np.float64),
        ),
        constraints=constraints,
        options=options,
    )
    solve_seconds = time.monotonic() - started
    result_meta = {
        "sat": answer.status == 0,
        "cells": len(cells),
        "cuts": density_cuts,
        "variables": variable_count,
        "sentinel_variables": sentinel_count,
        "constraints": constraint_index,
        "nonzeros": int(matrix.nnz),
        "solver": "scipy-highs",
        "engine": "scipy-highs-milp",
        "status": int(answer.status),
        "message": str(answer.message),
        "solve_seconds": solve_seconds,
        "maximize_preserved": maximize_preserved,
    }
    if answer.status == 2:
        return None, result_meta
    if answer.status != 0 or answer.x is None:
        raise RuntimeError(
            f"MILP density master did not decide feasibility: {answer.message}"
        )

    phases = dict(initial_phases)
    preserved_targets = 0
    for record in variable_rows:
        row_index = record["row_index"]
        selected_targets = [
            target
            for target, column in record["target_columns"].items()
            if answer.x[column] > 0.5
        ]
        sentinel_selected = (
            record["sentinel_column"] is not None
            and answer.x[record["sentinel_column"]] > 0.5
        )
        if len(selected_targets) + int(sentinel_selected) != 1:
            raise AssertionError(
                f"MILP selected an invalid target count for row {row_index}"
            )
        if sentinel_selected:
            represented = record["target_columns"]
            target = next(
                candidate
                for candidate in range(record["other"])
                if candidate not in represented
            )
        else:
            target = selected_targets[0]
        row = rows[row_index]
        row_prime = int(row["p"])
        h = int(row["h"])
        other, residual = row_data[row_index]
        if target == initial_phases[row_prime] % other:
            preserved_targets += 1
        if residual == 1:
            phases[row_prime] = target
        else:
            phases[row_prime] = exact_uncovered.crt(
                [
                    (target, other),
                    (phases[row_prime] % residual, residual),
                ]
            ) % h
    phases.update(fixed_targets)
    result_meta["preserved_targets"] = preserved_targets
    return phases, result_meta


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--initial-phases", type=Path, required=True)
    parser.add_argument("--fixed-targets", default="")
    parser.add_argument("--rounds", type=int, default=100)
    parser.add_argument("--batch", type=int, default=100)
    parser.add_argument("--max-component", type=int, default=1024)
    parser.add_argument(
        "--checker-diversity-modulus",
        type=int,
        default=0,
        help=(
            "after each density witness, block its whole coordinate cell "
            "modulo this coarse-domain prime power"
        ),
    )
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument(
        "--master-engine",
        choices=("pysat", "z3", "milp"),
        default="pysat",
        help=(
            "use exact Z3 weighted inequalities for large residual "
            "components instead of expanding integer weights"
        ),
    )
    parser.add_argument(
        "--milp-time-limit",
        type=float,
        default=0.0,
        help="optional positive HiGHS time limit in seconds",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="run the exact density checker once without solving a master",
    )
    parser.add_argument(
        "--direct-all-coarse-cells",
        action="store_true",
        help=(
            "solve one exact weighted master containing every cell of the "
            "coarse CRT plane; intended for modest coarse periods"
        ),
    )
    parser.add_argument(
        "--direct-max-cells",
        type=int,
        default=1_000_000,
        help="guard for --direct-all-coarse-cells",
    )
    parser.add_argument(
        "--cells-file",
        type=Path,
        help="optional JSON list of previously learned coarse cells",
    )
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--checkpoint-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    args = parser.parse_args()
    if args.prime < 2:
        raise SystemExit("--prime must be at least two")

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    if any(target_modulus(row) != 1 for row in rows):
        raise RuntimeError(
            "component-density CEGIS currently requires target_modulus=1"
        )
    phases = {
        int(row_prime): int(target)
        for row_prime, target in json.loads(
            args.initial_phases.read_text()
        ).items()
    }
    row_primes = {int(row["p"]) for row in rows}
    if row_primes - phases.keys():
        raise RuntimeError("initial phase file omits a pool prime")
    fixed_targets = {}
    for item in args.fixed_targets.split(","):
        if not item:
            continue
        row_prime, target = item.split(":", 1)
        fixed_targets[int(row_prime)] = int(target)
    if fixed_targets.keys() - row_primes:
        raise RuntimeError("a fixed prime is absent from the pool")
    phases.update(fixed_targets)

    cells: list[tuple[int, int]] = []
    if args.cells_file:
        cells_payload = json.loads(args.cells_file.read_text())
        if isinstance(cells_payload, dict):
            if "cells" not in cells_payload:
                raise RuntimeError(
                    "--cells-file checkpoint object has no cells field"
                )
            cells_payload = cells_payload["cells"]
        cells = [
            (int(k), int(l))
            for k, l in cells_payload
        ]
        if len(cells) != len(set(cells)):
            raise RuntimeError("--cells-file contains duplicate cells")
    seen = set(cells)
    history = []

    if args.direct_all_coarse_cells:
        if args.check_only:
            raise SystemExit(
                "--direct-all-coarse-cells and --check-only are exclusive"
            )
        coarse_period = math.lcm(
            *(
                int(row["h"]) // q_part(int(row["h"]), args.prime)
                for row in rows
            )
        )
        cell_count = coarse_period * coarse_period
        if cell_count > args.direct_max_cells:
            raise RuntimeError(
                f"direct coarse grid has {cell_count} cells, above guard "
                f"{args.direct_max_cells}"
            )
        cells = [
            (k, l)
            for k in range(coarse_period)
            for l in range(coarse_period)
        ]
        print(
            f"direct_all coarse_period={coarse_period} "
            f"cells={len(cells)}",
            flush=True,
        )
        if args.master_engine == "z3":
            next_phases, master = solve_master_z3(
                rows,
                phases,
                args.prime,
                cells,
                fixed_targets,
                maximize_preserved=False,
            )
        elif args.master_engine == "milp":
            next_phases, master = solve_master_milp(
                rows,
                phases,
                args.prime,
                cells,
                fixed_targets,
                args.milp_time_limit,
                maximize_preserved=False,
            )
        else:
            next_phases, master = solve_master(
                rows,
                phases,
                args.prime,
                cells,
                fixed_targets,
                args.solver,
            )
        if args.cells_file:
            args.cells_file.write_text(
                json.dumps([[k, l] for k, l in cells]) + "\n"
            )
        if next_phases is None:
            result = {
                "pool": str(args.pool),
                "prime": args.prime,
                "complete": False,
                "finite_density_core_unsat": True,
                "all_coarse_cells": True,
                "coarse_period": coarse_period,
                "cells": len(cells),
                "master": master,
            }
            args.checkpoint_output.write_text(
                json.dumps(result, indent=2) + "\n"
            )
            args.result_output.write_text(
                json.dumps(result, indent=2) + "\n"
            )
            print("UNSAT direct all-cell density master", flush=True)
            return 2

        phases = next_phases
        args.phase_output.write_text(
            json.dumps(
                {str(row_prime): target for row_prime, target in phases.items()}
            )
            + "\n"
        )
        misses, checker = find_low_density_cells(
            rows,
            phases,
            args.prime,
            1,
            args.max_component,
        )
        if misses:
            raise AssertionError(
                "direct all-cell density model failed the exact checker"
            )
        result = {
            "pool": str(args.pool),
            "prime": args.prime,
            "complete": True,
            "necessary_density_condition": True,
            "all_coarse_cells": True,
            "coarse_period": coarse_period,
            "cells": len(cells),
            "phase_file": str(args.phase_output),
            "master": master,
            "checker": checker,
        }
        args.checkpoint_output.write_text(
            json.dumps(result, indent=2) + "\n"
        )
        args.result_output.write_text(
            json.dumps(result, indent=2) + "\n"
        )
        print(
            "PASS direct all-cell residual-density condition; "
            "this is not a cover certificate",
            flush=True,
        )
        return 0

    for round_no in range(1, args.rounds + 1):
        misses, checker = find_low_density_cells(
            rows,
            phases,
            args.prime,
            args.batch,
            args.max_component,
            args.checker_diversity_modulus,
        )
        new_cells = [cell for cell in misses if cell not in seen]
        for cell in new_cells:
            seen.add(cell)
            cells.append(cell)
        print(
            f"round={round_no} violations={len(misses)} "
            f"new={len(new_cells)} total={len(cells)} "
            f"min_scaled={checker['minimum_scaled_density']}",
            flush=True,
        )
        if not misses:
            args.phase_output.write_text(
                json.dumps(
                    {str(row_prime): target for row_prime, target in phases.items()}
                )
                + "\n"
            )
            result = {
                "pool": str(args.pool),
                "prime": args.prime,
                "complete": True,
                "necessary_density_condition": True,
                "phase_file": str(args.phase_output),
                "cells": len(cells),
                "checker": checker,
                "history": history,
            }
            args.result_output.write_text(
                json.dumps(result, indent=2) + "\n"
            )
            print(
                "PASS exact residual-density condition; "
                "this is not a cover certificate",
                flush=True,
            )
            return 0
        if args.check_only:
            result = {
                "pool": str(args.pool),
                "prime": args.prime,
                "complete": False,
                "necessary_density_condition": False,
                "cells": len(cells),
                "new_cells": [[k, l] for k, l in new_cells],
                "checker": checker,
            }
            args.result_output.write_text(
                json.dumps(result, indent=2) + "\n"
            )
            print("FAIL exact residual-density condition", flush=True)
            return 1
        if not new_cells:
            raise RuntimeError("checker returned no new density cell")

        if args.master_engine == "z3":
            next_phases, master = solve_master_z3(
                rows,
                phases,
                args.prime,
                cells,
                fixed_targets,
            )
        elif args.master_engine == "milp":
            next_phases, master = solve_master_milp(
                rows,
                phases,
                args.prime,
                cells,
                fixed_targets,
                args.milp_time_limit,
            )
        else:
            next_phases, master = solve_master(
                rows,
                phases,
                args.prime,
                cells,
                fixed_targets,
                args.solver,
            )
        print(
            f"round={round_no} master_sat={master.get('sat')} "
            f"cuts={master.get('cuts', master.get('cells'))} "
            f"options={master.get('options')} "
            f"variables={master.get('variables')} "
            f"clauses={master.get('clauses')} "
            f"solve_s={master.get('solve_seconds')}",
            flush=True,
        )
        history.append(
            {
                "round": round_no,
                "checker": checker,
                "master": master,
            }
        )
        args.checkpoint_output.write_text(
            json.dumps(
                {
                    "pool": str(args.pool),
                    "prime": args.prime,
                    "cells": [[k, l] for k, l in cells],
                    "history": history,
                    "master": master,
                },
                indent=2,
            )
            + "\n"
        )
        if next_phases is None:
            result = {
                "pool": str(args.pool),
                "prime": args.prime,
                "complete": False,
                "finite_density_core_unsat": True,
                "cells": len(cells),
                "master": master,
            }
            args.result_output.write_text(
                json.dumps(result, indent=2) + "\n"
            )
            print("UNSAT finite residual-density core", flush=True)
            return 2
        phases = next_phases
        args.phase_output.write_text(
            json.dumps(
                {str(row_prime): target for row_prime, target in phases.items()}
            )
            + "\n"
        )

    args.result_output.write_text(
        json.dumps(
            {
                "pool": str(args.pool),
                "prime": args.prime,
                "complete": False,
                "round_limit": args.rounds,
                "cells": len(cells),
                "history": history,
            },
            indent=2,
        )
        + "\n"
    )
    print("INCOMPLETE round limit", flush=True)
    return 4


if __name__ == "__main__":
    raise SystemExit(main())
