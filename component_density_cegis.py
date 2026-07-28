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


def target_residue(row: dict) -> int:
    """Return the canonical allowed phase residue for a row."""
    modulus = target_modulus(row)
    if modulus < 1:
        raise ValueError("target modulus must be positive")
    return int(row.get("target_residue", 0)) % modulus


def q_part(value: int, prime: int) -> int:
    result = 1
    while value % prime == 0:
        result *= prime
        value //= prime
    return result


def coarse_target_restriction(
    row: dict,
    prime: int,
    initial_phase: int,
) -> tuple[int, int]:
    """Project a full phase restriction onto the frozen coarse coordinates.

    For the selected residual ``prime`` component, a phase modulo ``h``
    decomposes by CRT into a residual target and a coprime coarse target.
    The residual target is deliberately held at ``initial_phase`` by every
    density-master move.  This helper validates that fixed part and returns
    the remaining congruence on the coarse target.
    """
    h = int(row["h"])
    modulus = target_modulus(row)
    residue = target_residue(row)
    if h % modulus:
        raise ValueError("target modulus must divide the row modulus")
    residual = q_part(h, prime)
    residual_restriction = q_part(modulus, prime)
    coarse_modulus = modulus // residual_restriction
    if (
        residual % residual_restriction
        or (h // residual) % coarse_modulus
    ):
        raise AssertionError("target restriction does not split over CRT")
    if initial_phase % residual_restriction != residue % residual_restriction:
        raise ValueError(
            "initial phase violates the residual target restriction"
        )
    return residue % coarse_modulus, coarse_modulus


def coarse_target_allowed(
    target: int,
    restriction: tuple[int, int],
) -> bool:
    residue, modulus = restriction
    return target % modulus == residue


def required_scaled_density(
    residual_modulus: int,
    prime: int,
    algebraic_primes: tuple[int, ...] = (),
) -> int:
    """Return the union-bound weight required on one residual component.

    Ordinarily the selected fibres need total line density at least one.
    When ``prime`` is an algebraic factor of the perfect-power exponent,
    ``X^prime + 1`` already covers the origin sublattice ``prime*Z^2``.
    Its density inside the residual plane is exactly ``1 / prime^2``.
    Scaled row weights are integral with denominator ``residual_modulus``,
    so round the remaining required measure upward.
    """
    if residual_modulus < 1:
        raise ValueError("residual modulus must be positive")
    if prime not in algebraic_primes:
        return residual_modulus
    square = prime * prime
    numerator = residual_modulus * (square - 1)
    return (numerator + square - 1) // square


def scaled_density_at_cell(
    rows: list[dict],
    phases: dict[int, int],
    prime: int,
    cell: tuple[int, int],
) -> tuple[int, int]:
    """Recompute one residual-density cut using exact integer arithmetic."""
    residual_modulus = max(
        (q_part(int(row["h"]), prime) for row in rows),
        default=1,
    )
    k, l = cell
    scaled_density = 0
    for row in rows:
        h = int(row["h"])
        residual = q_part(h, prime)
        other = h // residual
        required_target = (
            int(row["a"]) * k + int(row["b"]) * l
        ) % other
        if phases[int(row["p"])] % other == required_target:
            scaled_density += residual_modulus // residual
    return scaled_density, residual_modulus


def audit_density_cuts(
    rows: list[dict],
    phases: dict[int, int],
    prime: int,
    cells: list[tuple[int, int]],
    algebraic_primes: tuple[int, ...] = (),
) -> dict:
    """Fail closed unless a proposed phase map satisfies every exact cut."""
    densities = []
    residual_modulus = 1
    for cell in cells:
        density, cell_modulus = scaled_density_at_cell(
            rows,
            phases,
            prime,
            cell,
        )
        residual_modulus = max(residual_modulus, cell_modulus)
        densities.append(density)
    if not cells:
        residual_modulus = max(
            (q_part(int(row["h"]), prime) for row in rows),
            default=1,
        )
    required_weight = required_scaled_density(
        residual_modulus,
        prime,
        algebraic_primes,
    )
    violations = [
        {
            "cell": list(cell),
            "scaled_density": density,
        }
        for cell, density in zip(cells, densities)
        if density < required_weight
    ]
    if violations:
        first = violations[0]
        raise AssertionError(
            "master phase failed exact residual-density replay at "
            f"{first['cell']}: {first['scaled_density']} < {required_weight}"
        )
    return {
        "cells": len(cells),
        "residual_modulus": residual_modulus,
        "required_scaled_density": required_weight,
        "minimum_scaled_density": min(densities) if densities else None,
        "maximum_scaled_density": max(densities) if densities else None,
        "violations": 0,
        "engine": "python-exact-integer-replay",
    }


def audit_residual_plane_cells(
    rows: list[dict],
    phases: dict[int, int],
    prime: int,
    cells: list[tuple[int, int]],
    algebraic_primes: tuple[int, ...] = (),
) -> dict:
    """Exactly test affine-line unions in selected first-level q-planes."""
    residual_modulus = max(
        (q_part(int(row["h"]), prime) for row in rows),
        default=1,
    )
    if residual_modulus != prime:
        raise ValueError(
            "residual-plane audit currently requires a first-level "
            "prime component"
        )
    full_mask = (1 << (prime * prime)) - 1
    algebraic_mask = 1 if prime in algebraic_primes else 0
    failures = []
    minimum_covered = prime * prime
    for k, l in cells:
        covered = algebraic_mask
        active_lines = 0
        active_full_rows = 0
        for row in rows:
            h = int(row["h"])
            residual = q_part(h, prime)
            other = h // residual
            coarse_target = (
                int(row["a"]) * k + int(row["b"]) * l
            ) % other
            if phases[int(row["p"])] % other != coarse_target:
                continue
            if residual == 1:
                active_full_rows += 1
                covered = full_mask
                break
            if residual != prime:
                raise ValueError("row has a higher residual prime power")
            active_lines += 1
            a = int(row["a"]) % prime
            b = int(row["b"]) % prime
            target = phases[int(row["p"])] % prime
            for x in range(prime):
                for y in range(prime):
                    if (a * x + b * y) % prime == target:
                        covered |= 1 << (x * prime + y)
        covered_points = covered.bit_count()
        minimum_covered = min(minimum_covered, covered_points)
        if covered != full_mask:
            missing = (~covered) & full_mask
            missing_index = (missing & -missing).bit_length() - 1
            failures.append(
                {
                    "cell": [k, l],
                    "covered_points": covered_points,
                    "uncovered_points": prime * prime - covered_points,
                    "first_uncovered_residual_point": [
                        missing_index // prime,
                        missing_index % prime,
                    ],
                    "active_lines": active_lines,
                    "active_full_rows": active_full_rows,
                }
            )
    return {
        "cells": len(cells),
        "prime": prime,
        "plane_points": prime * prime,
        "minimum_covered_points": minimum_covered if cells else None,
        "failed_cells": len(failures),
        "first_failure": failures[0] if failures else None,
        "engine": "python-exact-residual-plane-bitset-replay",
    }


def find_low_density_cells(
    rows: list[dict],
    phases: dict[int, int],
    prime: int,
    limit: int,
    max_component: int,
    diversity_modulus: int = 0,
    algebraic_primes: tuple[int, ...] = (),
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
    for algebraic_prime in algebraic_primes:
        if algebraic_prime != prime:
            maximal[algebraic_prime] = max(
                maximal.get(algebraic_prime, 0),
                1,
            )
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

    for algebraic_prime in algebraic_primes:
        if algebraic_prime == prime:
            continue
        width = expression_width
        modulus = z3.BitVecVal(algebraic_prime, width)
        zero = z3.BitVecVal(0, width)
        solver.add(
            z3.Or(
                z3.URem(
                    widened(kval[algebraic_prime], algebraic_prime),
                    modulus,
                )
                != zero,
                z3.URem(
                    widened(lval[algebraic_prime], algebraic_prime),
                    modulus,
                )
                != zero,
            )
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

    required_weight = required_scaled_density(
        residual_modulus,
        prime,
        algebraic_primes,
    )
    if activations:
        solver.add(
            z3.PbLe(
                list(zip(activations, weights)),
                required_weight - 1,
            )
        )
    elif required_weight <= 0:
        solver.add(z3.BoolVal(False))
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
        "required_scaled_density": required_weight,
        "components": components,
        "period": math.prod(components.values()),
        "row_count": len(rows),
        "cell_count": len(cells),
        "scaled_densities": scaled_densities,
        "minimum_scaled_density": (
            min(scaled_densities) if scaled_densities else None
        ),
        "diversity_modulus": diversity_modulus,
        "algebraic_primes": list(algebraic_primes),
        "sat": bool(cells),
        "engine": "z3-bitvector-pb",
    }


def scan_one_change(
    rows: list[dict],
    initial_phases: dict[int, int],
    prime: int,
    cells: list[tuple[int, int]],
    fixed_targets: dict[int, int],
    algebraic_primes: tuple[int, ...] = (),
    limit: int = 1,
) -> tuple[list[tuple[dict[int, int], dict]], dict]:
    """Exactly enumerate phase maps one coarse-row change away.

    The current density of every cell is computed once.  Changing one row
    subtracts its weight where its old target was active and adds the same
    weight where a candidate new target is active.  Targets absent from all
    supplied cells cannot improve any cut and need not be considered.
    """
    if limit < 0:
        raise ValueError("one-change scan limit must be nonnegative")
    started = time.monotonic()
    residual_modulus = max(
        (q_part(int(row["h"]), prime) for row in rows),
        default=1,
    )
    required_weight = required_scaled_density(
        residual_modulus,
        prime,
        algebraic_primes,
    )
    base_densities = [
        scaled_density_at_cell(rows, initial_phases, prime, cell)[0]
        for cell in cells
    ]
    row_data = []
    restrictions = []
    for row in rows:
        h = int(row["h"])
        residual = q_part(h, prime)
        row_data.append((h // residual, residual))
        restrictions.append(
            coarse_target_restriction(
                row,
                prime,
                initial_phases[int(row["p"])],
            )
        )

    scanned_rows = 0
    scanned_targets = 0
    winners = []
    for row_index, (row, (other, residual)) in enumerate(
        zip(rows, row_data)
    ):
        row_prime = int(row["p"])
        if other == 1 or row_prime in fixed_targets:
            continue
        scanned_rows += 1
        current_target = initial_phases[row_prime] % other
        weight = residual_modulus // residual
        required_targets = [
            (
                int(row["a"]) * k + int(row["b"]) * l
            ) % other
            for k, l in cells
        ]
        candidates = sorted(
            {
                target
                for target in required_targets
                if target != current_target
                and coarse_target_allowed(
                    target,
                    restrictions[row_index],
                )
            }
        )
        for target in candidates:
            scanned_targets += 1
            trial_densities = [
                density
                - (weight if required == current_target else 0)
                + (weight if required == target else 0)
                for density, required in zip(
                    base_densities,
                    required_targets,
                )
            ]
            if any(
                density < required_weight
                for density in trial_densities
            ):
                continue
            phases = dict(initial_phases)
            h = int(row["h"])
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
            exact_replay = audit_density_cuts(
                rows,
                phases,
                prime,
                cells,
                algebraic_primes,
            )
            winners.append(
                (
                    phases,
                    {
                        "row_index": row_index,
                        "prime": row_prime,
                        "h": h,
                        "old_phase": initial_phases[row_prime],
                        "new_phase": phases[row_prime],
                        "old_coarse_target": current_target,
                        "new_coarse_target": target,
                        "weight": weight,
                        "exact_replay": exact_replay,
                    },
                )
            )
            if limit and len(winners) >= limit:
                break
        if limit and len(winners) >= limit:
            break

    metadata = {
        "sat": bool(winners),
        "cells": len(cells),
        "scanned_rows": scanned_rows,
        "scanned_targets": scanned_targets,
        "winner_count": len(winners),
        "limit": limit,
        "residual_modulus": residual_modulus,
        "required_scaled_density": required_weight,
        "solve_seconds": time.monotonic() - started,
        "engine": "python-exact-one-change-scan",
    }
    return winners, metadata


def scan_two_changes(
    rows: list[dict],
    initial_phases: dict[int, int],
    prime: int,
    cells: list[tuple[int, int]],
    fixed_targets: dict[int, int],
    algebraic_primes: tuple[int, ...] = (),
    geometry_cells: list[tuple[int, int]] | None = None,
    limit: int = 1,
) -> tuple[list[tuple[dict[int, int], dict]], dict]:
    """Exactly enumerate phase maps at most two coarse-row changes away.

    Every solution must repair one selected currently failing cell.  The
    first changed row is therefore anchored to its unique useful target on
    that cell.  After applying it, a single second row must use one common
    target on every remaining failing cell.  Intersecting those target
    classes leaves only the pairs that need full exact replay.
    """
    if limit < 0:
        raise ValueError("two-change scan limit must be nonnegative")
    if not cells:
        raise ValueError("two-change scan requires at least one cell")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    started = time.monotonic()
    geometry_cells = geometry_cells or []
    row_count = len(rows)
    cell_count = len(cells)
    residual_modulus = max(
        (q_part(int(row["h"]), prime) for row in rows),
        default=1,
    )
    required_weight = required_scaled_density(
        residual_modulus,
        prime,
        algebraic_primes,
    )
    other_moduli = np.empty(row_count, dtype=np.int64)
    residual_moduli = np.empty(row_count, dtype=np.int64)
    weights = np.empty(row_count, dtype=np.int64)
    old_targets = np.empty(row_count, dtype=np.int64)
    allowed_residues = np.empty(row_count, dtype=np.int64)
    allowed_moduli = np.empty(row_count, dtype=np.int64)
    mutable = np.ones(row_count, dtype=np.bool_)
    target_matrix = np.empty((row_count, cell_count), dtype=np.int64)

    for row_index, row in enumerate(rows):
        h = int(row["h"])
        residual = q_part(h, prime)
        other = h // residual
        row_prime = int(row["p"])
        restriction = coarse_target_restriction(
            row,
            prime,
            initial_phases[row_prime],
        )
        other_moduli[row_index] = other
        residual_moduli[row_index] = residual
        weights[row_index] = residual_modulus // residual
        old_targets[row_index] = initial_phases[row_prime] % other
        allowed_residues[row_index] = restriction[0]
        allowed_moduli[row_index] = restriction[1]
        if other == 1 or row_prime in fixed_targets:
            mutable[row_index] = False
        a = int(row["a"])
        b = int(row["b"])
        if other == 1:
            target_matrix[row_index, :] = 0
        else:
            target_matrix[row_index, :] = np.fromiter(
                (
                    (
                        a * (int(k) % other)
                        + b * (int(l) % other)
                    )
                    % other
                    for k, l in cells
                ),
                dtype=np.int64,
                count=cell_count,
            )

    base_densities = np.zeros(cell_count, dtype=np.int64)
    for row_index in range(row_count):
        base_densities[
            target_matrix[row_index] == old_targets[row_index]
        ] += weights[row_index]
    failing = np.flatnonzero(base_densities < required_weight)
    if not len(failing):
        geometry_replay = (
            audit_residual_plane_cells(
                rows,
                initial_phases,
                prime,
                geometry_cells,
                algebraic_primes,
            )
            if geometry_cells
            else None
        )
        if geometry_replay and geometry_replay["failed_cells"]:
            raise ValueError(
                "two-change scan needs a failing density anchor when the "
                "initial phase fails only residual-plane geometry"
            )
        exact_replay = audit_density_cuts(
            rows,
            initial_phases,
            prime,
            cells,
            algebraic_primes,
        )
        return [
            (
                dict(initial_phases),
                {
                    "changes": [],
                    "exact_replay": exact_replay,
                    "geometry_replay": geometry_replay,
                },
            )
        ], {
            "sat": True,
            "cells": cell_count,
            "initial_violations": 0,
            "anchor_cell_index": None,
            "first_moves_scanned": 0,
            "second_candidates": 0,
            "pairs_replayed": 0,
            "winner_count": 1,
            "limit": limit,
            "residual_modulus": residual_modulus,
            "required_scaled_density": required_weight,
            "geometry_cells": len(geometry_cells),
            "solve_seconds": time.monotonic() - started,
            "engine": "numpy-exact-two-change-scan",
        }

    anchor_index = int(
        failing[np.argmin(base_densities[failing])]
    )
    first_moves_scanned = 0
    second_candidates = 0
    pairs_replayed = 0
    winners = []

    def retarget(
        phases: dict[int, int],
        row_index: int,
        target: int,
    ) -> None:
        row = rows[row_index]
        row_prime = int(row["p"])
        h = int(row["h"])
        other = int(other_moduli[row_index])
        residual = int(residual_moduli[row_index])
        if residual == 1:
            phases[row_prime] = target
        else:
            phases[row_prime] = exact_uncovered.crt(
                [
                    (target, other),
                    (phases[row_prime] % residual, residual),
                ]
            ) % h

    for first_index in np.flatnonzero(mutable):
        first_index = int(first_index)
        first_target = int(target_matrix[first_index, anchor_index])
        if first_target == int(old_targets[first_index]):
            continue
        if (
            first_target % int(allowed_moduli[first_index])
            != int(allowed_residues[first_index])
        ):
            continue
        first_moves_scanned += 1
        trial = base_densities.copy()
        first_row_targets = target_matrix[first_index]
        first_weight = int(weights[first_index])
        trial[first_row_targets == old_targets[first_index]] -= first_weight
        trial[first_row_targets == first_target] += first_weight
        deficits = np.flatnonzero(trial < required_weight)

        if not len(deficits):
            phases = dict(initial_phases)
            retarget(phases, first_index, first_target)
            phases.update(fixed_targets)
            geometry_replay = (
                audit_residual_plane_cells(
                    rows,
                    phases,
                    prime,
                    geometry_cells,
                    algebraic_primes,
                )
                if geometry_cells
                else None
            )
            if geometry_replay and geometry_replay["failed_cells"]:
                raise ValueError(
                    "geometry-aware two-change scan requires the supplied "
                    "density core to reject every one-change phase"
                )
            exact_replay = audit_density_cuts(
                rows,
                phases,
                prime,
                cells,
                algebraic_primes,
            )
            winners.append(
                (
                    phases,
                    {
                        "changes": [
                            {
                                "row_index": first_index,
                                "prime": int(rows[first_index]["p"]),
                                "old_coarse_target": int(
                                    old_targets[first_index]
                                ),
                                "new_coarse_target": first_target,
                            }
                        ],
                        "exact_replay": exact_replay,
                        "geometry_replay": geometry_replay,
                    },
                )
            )
            if limit and len(winners) >= limit:
                break
            continue

        first_deficit = int(deficits[0])
        second_targets = target_matrix[:, first_deficit]
        maximum_deficit = int(
            np.max(required_weight - trial[deficits])
        )
        candidate_mask = (
            mutable
            & (np.arange(row_count) != first_index)
            & (second_targets != old_targets)
            & (
                second_targets % allowed_moduli
                == allowed_residues
            )
            & (weights >= maximum_deficit)
        )
        candidates = np.flatnonzero(candidate_mask)
        # A useful second row must take the same target on every remaining
        # failing cell.  Intersect a few columns at a time; candidate sets
        # normally collapse after the second deficit.
        for deficit_index in deficits[1:]:
            if not len(candidates):
                break
            candidates = candidates[
                target_matrix[candidates, int(deficit_index)]
                == second_targets[candidates]
            ]
        second_candidates += len(candidates)

        for second_index in candidates:
            second_index = int(second_index)
            second_target = int(second_targets[second_index])
            second_weight = int(weights[second_index])
            second_trial = trial.copy()
            second_row_targets = target_matrix[second_index]
            second_trial[
                second_row_targets == old_targets[second_index]
            ] -= second_weight
            second_trial[
                second_row_targets == second_target
            ] += second_weight
            pairs_replayed += 1
            if np.any(second_trial < required_weight):
                continue
            phases = dict(initial_phases)
            retarget(phases, first_index, first_target)
            retarget(phases, second_index, second_target)
            phases.update(fixed_targets)
            exact_replay = audit_density_cuts(
                rows,
                phases,
                prime,
                cells,
                algebraic_primes,
            )
            geometry_replay = (
                audit_residual_plane_cells(
                    rows,
                    phases,
                    prime,
                    geometry_cells,
                    algebraic_primes,
                )
                if geometry_cells
                else None
            )
            if geometry_replay and geometry_replay["failed_cells"]:
                continue
            winners.append(
                (
                    phases,
                    {
                        "changes": [
                            {
                                "row_index": first_index,
                                "prime": int(rows[first_index]["p"]),
                                "old_coarse_target": int(
                                    old_targets[first_index]
                                ),
                                "new_coarse_target": first_target,
                            },
                            {
                                "row_index": second_index,
                                "prime": int(rows[second_index]["p"]),
                                "old_coarse_target": int(
                                    old_targets[second_index]
                                ),
                                "new_coarse_target": second_target,
                            },
                        ],
                        "exact_replay": exact_replay,
                        "geometry_replay": geometry_replay,
                    },
                )
            )
            if limit and len(winners) >= limit:
                break
        if limit and len(winners) >= limit:
            break

    metadata = {
        "sat": bool(winners),
        "cells": cell_count,
        "initial_violations": int(len(failing)),
        "anchor_cell_index": anchor_index,
        "first_moves_scanned": first_moves_scanned,
        "second_candidates": int(second_candidates),
        "pairs_replayed": pairs_replayed,
        "winner_count": len(winners),
        "limit": limit,
        "residual_modulus": residual_modulus,
        "required_scaled_density": required_weight,
        "target_matrix_bytes": int(target_matrix.nbytes),
        "geometry_cells": len(geometry_cells),
        "solve_seconds": time.monotonic() - started,
        "engine": "numpy-exact-two-change-scan",
    }
    return winners, metadata


def solve_master(
    rows: list[dict],
    initial_phases: dict[int, int],
    prime: int,
    cells: list[tuple[int, int]],
    fixed_targets: dict[int, int],
    solver_name: str,
    algebraic_primes: tuple[int, ...] = (),
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
    coarse_restrictions = [
        coarse_target_restriction(
            row,
            prime,
            initial_phases[int(row["p"])],
        )
        for row in rows
    ]
    required_weight = required_scaled_density(
        residual_modulus,
        prime,
        algebraic_primes,
    )

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
            if (
                other > 1
                and int(row["p"]) not in fixed_targets
                and coarse_target_allowed(
                    target,
                    coarse_restrictions[row_index],
                )
            ):
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
                if not coarse_target_allowed(
                    required[row_index],
                    coarse_restrictions[row_index],
                ):
                    continue
                literals.append(
                    variable[(row_index, required[row_index])]
                )
                weights.append(weight)
        bound = required_weight - constant
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
        "required_scaled_density": required_weight,
        "algebraic_primes": list(algebraic_primes),
        "restricted_rows": sum(
            modulus > 1 for _residue, modulus in coarse_restrictions
        ),
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
    algebraic_primes: tuple[int, ...] = (),
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
    coarse_restrictions = [
        coarse_target_restriction(
            row,
            prime,
            initial_phases[int(row["p"])],
        )
        for row in rows
    ]
    required_weight = required_scaled_density(
        residual_modulus,
        prime,
        algebraic_primes,
    )

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
        allowed_residue, allowed_modulus = coarse_restrictions[row_index]
        if allowed_modulus > 1:
            solver.add(variable % allowed_modulus == allowed_residue)
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
                if not coarse_target_allowed(
                    required,
                    coarse_restrictions[row_index],
                ):
                    continue
                available_weight += weight
                weighted_matches.append(
                    z3.If(
                        coarse_targets[row_index] == required,
                        weight,
                        0,
                    )
                )
        bound = required_weight - constant
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
        "required_scaled_density": required_weight,
        "algebraic_primes": list(algebraic_primes),
        "restricted_rows": sum(
            modulus > 1 for _residue, modulus in coarse_restrictions
        ),
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
    algebraic_primes: tuple[int, ...] = (),
    maximum_changes: int = 0,
) -> tuple[dict[int, int] | None, dict]:
    """Solve exact weighted density cuts as a sparse binary MILP."""

    if maximum_changes < 0:
        raise ValueError("maximum changes must be nonnegative")

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
    coarse_restrictions = [
        coarse_target_restriction(
            row,
            prime,
            initial_phases[int(row["p"])],
        )
        for row in rows
    ]
    required_weight = required_scaled_density(
        residual_modulus,
        prime,
        algebraic_primes,
    )

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
                if not coarse_target_allowed(
                    required,
                    coarse_restrictions[row_index],
                ):
                    continue
                available_weight += weight
                entries.append((row_index, required, float(weight)))
        bound = required_weight - constant
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
        allowed_residue, allowed_modulus = coarse_restrictions[row_index]
        allowed_target_count = other // allowed_modulus
        target_columns = {}
        indices = []
        for target in targets:
            target_columns[target] = variable_count
            variable[(row_index, target)] = variable_count
            indices.append(variable_count)
            variable_count += 1
        sentinel_column = None
        if len(targets) < allowed_target_count:
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
                "allowed_residue": allowed_residue,
                "allowed_modulus": allowed_modulus,
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

    if maximum_changes:
        for record in variable_rows:
            row_index = record["row_index"]
            other = record["other"]
            initial_target = (
                initial_phases[int(rows[row_index]["p"])] % other
            )
            matrix_rows.append(constraint_index)
            matrix_cols.append(
                record["target_columns"][initial_target]
            )
            matrix_values.append(1.0)
        lower_bounds.append(
            float(max(0, len(variable_rows) - maximum_changes))
        )
        upper_bounds.append(np.inf)
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
    if maximize_preserved and not maximum_changes:
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
        "maximum_changes": maximum_changes,
        "required_scaled_density": required_weight,
        "algebraic_primes": list(algebraic_primes),
        "restricted_rows": sum(
            modulus > 1 for _residue, modulus in coarse_restrictions
        ),
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
            allowed_residue = record["allowed_residue"]
            allowed_modulus = record["allowed_modulus"]
            target = next(
                candidate
                for candidate in range(
                    allowed_residue,
                    record["other"],
                    allowed_modulus,
                )
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
        "--maximum-changes",
        type=int,
        default=0,
        help=(
            "MILP-only hard bound on phase changes from --initial-phases; "
            "zero disables the bound and retains the preservation objective"
        ),
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="run the exact density checker once without solving a master",
    )
    parser.add_argument(
        "--master-first",
        action="store_true",
        help=(
            "solve the supplied --cells-file cuts before asking the checker "
            "for additional violations"
        ),
    )
    parser.add_argument(
        "--one-change-scan",
        action="store_true",
        help=(
            "exactly scan every legal one-row retarget against the supplied "
            "cells, without invoking a general master"
        ),
    )
    parser.add_argument(
        "--two-change-scan",
        action="store_true",
        help=(
            "exactly scan every legal phase map at most two row retargets "
            "from the initial phase against the supplied cells"
        ),
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
        action="append",
        help=(
            "optional JSON list or checker object containing learned coarse "
            "cells; repeat to union several sources"
        ),
    )
    parser.add_argument(
        "--geometry-cells-file",
        type=Path,
        action="append",
        help=(
            "exact-hole file whose coarse cells must also have complete "
            "first-level residual affine-plane line unions"
        ),
    )
    parser.add_argument("--phase-output", type=Path, required=True)
    parser.add_argument("--checkpoint-output", type=Path, required=True)
    parser.add_argument("--result-output", type=Path, required=True)
    args = parser.parse_args()
    if args.prime < 2:
        raise SystemExit("--prime must be at least two")
    if args.maximum_changes < 0:
        parser.error("--maximum-changes must be nonnegative")
    if args.maximum_changes and args.master_engine != "milp":
        parser.error(
            "--maximum-changes currently requires --master-engine milp"
        )
    if args.one_change_scan and args.two_change_scan:
        parser.error("choose at most one exact change-scan mode")
    if (args.one_change_scan or args.two_change_scan) and (
        args.check_only
        or args.master_first
        or args.direct_all_coarse_cells
    ):
        parser.error(
            "exact change scans are incompatible with "
            "checker/direct/master-first modes"
        )
    if args.geometry_cells_file and not args.two_change_scan:
        parser.error(
            "--geometry-cells-file currently requires --two-change-scan"
        )

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    algebraic_primes = tuple(
        int(value) for value in payload.get("algebraic_primes", ())
    )
    if len(set(algebraic_primes)) != len(algebraic_primes):
        raise RuntimeError("pool repeats an algebraic prime")
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
    for row in rows:
        row_prime = int(row["p"])
        h = int(row["h"])
        modulus = target_modulus(row)
        residue = target_residue(row)
        if h % modulus:
            raise RuntimeError(
                f"target modulus does not divide h for prime {row_prime}"
            )
        if phases[row_prime] % modulus != residue:
            raise RuntimeError(
                f"phase for prime {row_prime} violates its target restriction"
            )
        # Validate that the restriction decomposes across the selected
        # residual component before invoking either checker or master.
        coarse_target_restriction(
            row,
            args.prime,
            phases[row_prime],
        )

    cells: list[tuple[int, int]] = []
    loaded_cells = set()
    for cells_file in args.cells_file or ():
        cells_payload = json.loads(cells_file.read_text())
        if isinstance(cells_payload, dict):
            cell_key = next(
                (
                    key
                    for key in (
                        "cells",
                        "exact_misses",
                        "misses",
                        "new_cells",
                    )
                    if key in cells_payload
                ),
                None,
            )
            if cell_key is None:
                raise RuntimeError(
                    "--cells-file checkpoint object has no supported "
                    "cell-list field"
                )
            cells_payload = cells_payload[cell_key]
        for raw_k, raw_l in cells_payload:
            cell = (int(raw_k), int(raw_l))
            if not all(
                algebraic_prime == args.prime
                or cell[0] % algebraic_prime
                or cell[1] % algebraic_prime
                for algebraic_prime in algebraic_primes
            ):
                continue
            if cell not in loaded_cells:
                loaded_cells.add(cell)
                cells.append(cell)
    if args.master_first and not cells:
        raise RuntimeError("--master-first requires nonempty supplied cells")
    if (args.one_change_scan or args.two_change_scan) and not cells:
        raise RuntimeError("exact change scan requires supplied cells")
    seen = set(cells)
    history = []
    geometry_cells = []
    seen_geometry_cells = set()
    for geometry_file in args.geometry_cells_file or ():
        geometry_payload = json.loads(geometry_file.read_text())
        if isinstance(geometry_payload, dict):
            geometry_key = next(
                (
                    key
                    for key in (
                        "cells",
                        "exact_misses",
                        "misses",
                        "new_cells",
                    )
                    if key in geometry_payload
                ),
                None,
            )
            if geometry_key is None:
                raise RuntimeError(
                    "--geometry-cells-file has no supported cell-list field"
                )
            geometry_payload = geometry_payload[geometry_key]
        for raw_k, raw_l in geometry_payload:
            cell = (int(raw_k), int(raw_l))
            if cell not in seen_geometry_cells:
                seen_geometry_cells.add(cell)
                geometry_cells.append(cell)

    if args.one_change_scan:
        winners, scan = scan_one_change(
            rows,
            phases,
            args.prime,
            cells,
            fixed_targets,
            algebraic_primes,
        )
        checkpoint = {
            "pool": str(args.pool),
            "prime": args.prime,
            "cells": [[k, l] for k, l in cells],
            "scan": scan,
            "winner": winners[0][1] if winners else None,
        }
        args.checkpoint_output.write_text(
            json.dumps(checkpoint, indent=2) + "\n"
        )
        if not winners:
            args.phase_output.write_text(
                json.dumps(
                    {
                        str(row_prime): target
                        for row_prime, target in phases.items()
                    }
                )
                + "\n"
            )
            result = {
                "pool": str(args.pool),
                "prime": args.prime,
                "complete": False,
                "one_change_feasible": False,
                "scan": scan,
            }
            args.result_output.write_text(
                json.dumps(result, indent=2) + "\n"
            )
            print(
                f"NO_ONE_CHANGE cells={len(cells)} "
                f"targets={scan['scanned_targets']} "
                f"solve_s={scan['solve_seconds']}",
                flush=True,
            )
            return 2
        phases, winner = winners[0]
        args.phase_output.write_text(
            json.dumps(
                {
                    str(row_prime): target
                    for row_prime, target in phases.items()
                }
            )
            + "\n"
        )
        result = {
            "pool": str(args.pool),
            "prime": args.prime,
            "complete": False,
            "one_change_feasible": True,
            "scan": scan,
            "winner": winner,
        }
        args.result_output.write_text(
            json.dumps(result, indent=2) + "\n"
        )
        print(
            f"ONE_CHANGE_SAT cells={len(cells)} "
            f"prime={winner['prime']} h={winner['h']} "
            f"targets={scan['scanned_targets']} "
            f"solve_s={scan['solve_seconds']}",
            flush=True,
        )
        return 0

    if args.two_change_scan:
        winners, scan = scan_two_changes(
            rows,
            phases,
            args.prime,
            cells,
            fixed_targets,
            algebraic_primes,
            geometry_cells,
        )
        winner = winners[0][1] if winners else None
        args.checkpoint_output.write_text(
            json.dumps(
                {
                    "pool": str(args.pool),
                    "prime": args.prime,
                    "cells": [[k, l] for k, l in cells],
                    "scan": scan,
                    "winner": winner,
                },
                indent=2,
            )
            + "\n"
        )
        if winners:
            phases = winners[0][0]
        args.phase_output.write_text(
            json.dumps(
                {
                    str(row_prime): target
                    for row_prime, target in phases.items()
                }
            )
            + "\n"
        )
        result = {
            "pool": str(args.pool),
            "prime": args.prime,
            "complete": False,
            "two_change_feasible": bool(winners),
            "scan": scan,
            "winner": winner,
        }
        args.result_output.write_text(
            json.dumps(result, indent=2) + "\n"
        )
        label = "TWO_CHANGE_SAT" if winners else "NO_TWO_CHANGE"
        print(
            f"{label} cells={len(cells)} "
            f"first={scan['first_moves_scanned']} "
            f"pairs={scan['pairs_replayed']} "
            f"solve_s={scan['solve_seconds']}",
            flush=True,
        )
        return 0 if winners else 2

    def invoke_master(
        current_phases: dict[int, int],
        current_cells: list[tuple[int, int]],
        maximize_preserved: bool = True,
    ) -> tuple[dict[int, int] | None, dict]:
        if args.master_engine == "z3":
            next_phases, master = solve_master_z3(
                rows,
                current_phases,
                args.prime,
                current_cells,
                fixed_targets,
                maximize_preserved=maximize_preserved,
                algebraic_primes=algebraic_primes,
            )
        elif args.master_engine == "milp":
            next_phases, master = solve_master_milp(
                rows,
                current_phases,
                args.prime,
                current_cells,
                fixed_targets,
                args.milp_time_limit,
                maximize_preserved=maximize_preserved,
                algebraic_primes=algebraic_primes,
                maximum_changes=args.maximum_changes,
            )
        else:
            next_phases, master = solve_master(
                rows,
                current_phases,
                args.prime,
                current_cells,
                fixed_targets,
                args.solver,
                algebraic_primes,
            )
        if next_phases is not None:
            master["exact_replay"] = audit_density_cuts(
                rows,
                next_phases,
                args.prime,
                current_cells,
                algebraic_primes,
            )
        return next_phases, master

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
            if all(
                algebraic_prime == args.prime
                or k % algebraic_prime
                or l % algebraic_prime
                for algebraic_prime in algebraic_primes
            )
        ]
        print(
            f"direct_all coarse_period={coarse_period} "
            f"cells={len(cells)}",
            flush=True,
        )
        next_phases, master = invoke_master(
            phases,
            cells,
            maximize_preserved=False,
        )
        if args.cells_file and len(args.cells_file) == 1:
            args.cells_file[0].write_text(
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
            algebraic_primes=algebraic_primes,
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

    if args.master_first:
        next_phases, master = invoke_master(phases, cells)
        print(
            f"seed_master master_sat={master.get('sat')} "
            f"cuts={master.get('cuts', master.get('cells'))} "
            f"options={master.get('options')} "
            f"variables={master.get('variables')} "
            f"clauses={master.get('clauses')} "
            f"solve_s={master.get('solve_seconds')}",
            flush=True,
        )
        history.append(
            {
                "round": 0,
                "source": "supplied_cells",
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
            print("UNSAT supplied residual-density core", flush=True)
            return 2
        phases = next_phases
        args.phase_output.write_text(
            json.dumps(
                {
                    str(row_prime): target
                    for row_prime, target in phases.items()
                }
            )
            + "\n"
        )

    for round_no in range(1, args.rounds + 1):
        misses, checker = find_low_density_cells(
            rows,
            phases,
            args.prime,
            args.batch,
            args.max_component,
            args.checker_diversity_modulus,
            algebraic_primes,
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

        next_phases, master = invoke_master(phases, cells)
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
