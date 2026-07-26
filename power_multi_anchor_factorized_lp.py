#!/usr/bin/env python3
"""CRT-factorized exact dual for power constructions and 2--3 anchors."""

from __future__ import annotations

import argparse
import functools
import itertools
import json
import math
import os
import random
import sys
from fractions import Fraction
from pathlib import Path

import exact_uncovered
from power_anchor_capacity_lp import power_target_congruence


def valuation(value: int, prime: int) -> int:
    exponent = 0
    while value % prime == 0:
        exponent += 1
        value //= prime
    return exponent


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument(
        "--period",
        type=int,
        default=0,
        help="if nonzero, retain only rows whose subgroup order divides it",
    )
    parser.add_argument(
        "--max-component",
        type=int,
        default=0,
        help="if nonzero, retain rows whose prime-power components are at most this",
    )
    parser.add_argument("--anchor-primes", default="5,7,13")
    parser.add_argument(
        "--anchor-targets",
        help=(
            "comma-separated explicit power-compatible anchor targets; "
            "required when the valid target space has multiple translation "
            "orbits"
        ),
    )
    parser.add_argument("--weight-denominator", type=int, default=100000)
    parser.add_argument(
        "--float-only",
        action="store_true",
        help="report the HiGHS primal ratio without constructing a rational dual",
    )
    parser.add_argument(
        "--linprog-method",
        choices=("highs", "highs-ds", "highs-ipm"),
        default="highs",
        help="SciPy HiGHS algorithm used for the anchor-capacity LP",
    )
    parser.add_argument(
        "--fractional-output",
        type=Path,
        help=(
            "write the primal LP mixture over each row's distinct local "
            "target patterns"
        ),
    )
    parser.add_argument(
        "--evaluate-phase-file",
        type=Path,
        help=(
            "skip optimization and report exact conditional-density "
            "statistics for this saved phase assignment"
        ),
    )
    parser.add_argument(
        "--profile-output",
        type=Path,
        help="write per-cell conditional densities in evaluation mode",
    )
    parser.add_argument(
        "--count-required-only",
        action="store_true",
        help=(
            "report anchor/residual/required cell counts after algebraic "
            "identities, then exit before phase evaluation or LP construction"
        ),
    )
    parser.add_argument(
        "--cell-profile-input",
        type=Path,
        help=(
            "restrict the exact cell calculation to cells below a threshold "
            "in a previously exported compatible profile"
        ),
    )
    parser.add_argument(
        "--cell-profile-below",
        type=float,
        default=1.2,
        help="coverage cutoff used with --cell-profile-input",
    )
    parser.add_argument(
        "--optimize-phase-output",
        type=Path,
        help=(
            "in evaluation mode, coordinate-optimize the exact conditional "
            "capacity profile and write the resulting full phase assignment"
        ),
    )
    parser.add_argument(
        "--capacity-sweeps",
        type=int,
        default=0,
        help="number of exact conditional-capacity coordinate sweeps",
    )
    parser.add_argument(
        "--capacity-threshold",
        type=float,
        default=1.0,
        help="coverage threshold used by the squared-deficit objective",
    )
    parser.add_argument("--capacity-seed", type=int, default=203)
    args = parser.parse_args()
    if args.power < 1:
        raise SystemExit("--power must be positive")
    if args.capacity_sweeps < 0:
        raise SystemExit("--capacity-sweeps must be nonnegative")
    if args.capacity_threshold <= 0:
        raise SystemExit("--capacity-threshold must be positive")
    if args.cell_profile_below <= 0:
        raise SystemExit("--cell-profile-below must be positive")
    if args.optimize_phase_output and not args.evaluate_phase_file:
        raise SystemExit("--optimize-phase-output requires --evaluate-phase-file")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore
    from scipy.optimize import linprog  # type: ignore
    from scipy.sparse import coo_matrix  # type: ignore

    candidates = json.loads(args.pool.read_text())["choices"]
    if args.period:
        candidates = [
            row for row in candidates if args.period % int(row["h"]) == 0
        ]
    if args.max_component:
        candidates = [
            row
            for row in candidates
            if max(
                (
                    prime**exponent
                    for prime, exponent in exact_uncovered.factor(
                        int(row["h"])
                    ).items()
                ),
                default=1,
            )
            <= args.max_component
        ]
    requested = tuple(int(value) for value in args.anchor_primes.split(",") if value)
    if not 2 <= len(requested) <= 6:
        raise SystemExit("two to six anchor primes are required")
    by_prime = {int(row["p"]): (index, row) for index, row in enumerate(candidates)}
    if any(prime not in by_prime for prime in requested):
        raise RuntimeError("an anchor prime is missing from the pool")
    anchors = [by_prime[prime][1] for prime in requested]
    anchor_indices = {by_prime[prime][0] for prime in requested}
    anchor_period = math.lcm(*(int(row["h"]) for row in anchors))
    if anchor_period > 10000:
        raise RuntimeError("anchor period exceeds the explicit-cell guard")

    algebraic_primes = tuple(
        prime for prime in exact_uncovered.factor(args.power) if prime % 2
    )
    algebraic_set = set(algebraic_primes)
    sophie_germain = args.power % 4 == 0
    period = math.lcm(
        anchor_period,
        *algebraic_primes,
        *(4,) if sophie_germain else (),
    )
    components = {
        prime: prime**exponent
        for prime, exponent in exact_uncovered.factor(period).items()
    }

    valid_anchor_targets = []
    for row in anchors:
        residue, modulus = power_target_congruence(
            int(row["h"]), int(row["p"]), args.power
        )
        valid_anchor_targets.append(
            tuple(
                target
                for target in range(int(row["h"]))
                if target % modulus == residue
            )
        )
    valid_phase_space = set(itertools.product(*valid_anchor_targets))
    translation_image = {
        tuple(
            (
                int(row["a"]) * args.power * shift_k
                + int(row["b"]) * args.power * shift_l
            )
            % int(row["h"])
            for row in anchors
        )
        for shift_k in range(anchor_period)
        for shift_l in range(anchor_period)
    }
    if args.anchor_targets:
        normalized_targets = tuple(
            int(value) for value in args.anchor_targets.split(",") if value
        )
        if len(normalized_targets) != len(anchors):
            raise SystemExit("--anchor-targets must match --anchor-primes")
        if normalized_targets not in valid_phase_space:
            raise SystemExit("--anchor-targets are not power-compatible")
    else:
        normalized_targets = min(valid_phase_space)
    normalized_orbit = {
        tuple(
            (target + shift) % int(row["h"])
            for target, shift, row in zip(normalized_targets, image, anchors)
        )
        for image in translation_image
    }
    if not args.anchor_targets and normalized_orbit != valid_phase_space:
        raise RuntimeError(
            "anchor targets are not one translation orbit; specify "
            "--anchor-targets for an orbit representative"
        )

    anchor_cells = tuple(
        sorted(
            {
                tuple(
                    (int(row["a"]) * k + int(row["b"]) * l) % int(row["h"])
                    for row in anchors
                )
                for k in range(anchor_period)
                for l in range(anchor_period)
            }
        )
    )
    residual_cells = tuple(
        cell
        for cell in anchor_cells
        if all(
            value != normalized_targets[index]
            for index, value in enumerate(cell)
        )
    )

    @functools.lru_cache(maxsize=None)
    def local_count(
        prime: int,
        component: int,
        row_modulus: int,
        row_a: int,
        row_b: int,
        row_target: int,
        anchor_conditions: tuple[tuple[int, int, int, int], ...],
    ) -> int:
        count = 0
        for k in range(component):
            for l in range(component):
                if row_modulus and (
                    row_a * k + row_b * l - row_target
                ) % row_modulus:
                    continue
                if any(
                    (a * k + b * l - target) % modulus
                    for modulus, a, b, target in anchor_conditions
                ):
                    continue
                if prime in algebraic_set and k % prime == 0 and l % prime == 0:
                    continue
                if prime == 2 and sophie_germain and k % 4 == 2 and l % 4 == 0:
                    continue
                count += 1
        return count

    def count_for_cell(
        cell: tuple[int, ...], shared: int, a: int, b: int, target: int
    ) -> int:
        count = 1
        for prime, component in components.items():
            row_exponent = valuation(shared, prime)
            row_modulus = prime**row_exponent if row_exponent else 0
            conditions = []
            for row, value in zip(anchors, cell):
                exponent = valuation(int(row["h"]), prime)
                if not exponent:
                    continue
                modulus = prime**exponent
                conditions.append(
                    (
                        modulus,
                        int(row["a"]) % modulus,
                        int(row["b"]) % modulus,
                        value % modulus,
                    )
                )
            count *= local_count(
                prime,
                component,
                row_modulus,
                a % row_modulus if row_modulus else 0,
                b % row_modulus if row_modulus else 0,
                target % row_modulus if row_modulus else 0,
                tuple(conditions),
            )
        return count

    raw_demand = tuple(count_for_cell(cell, 1, 0, 0, 0) for cell in residual_cells)
    required_cells = tuple(
        cell for cell, demand in zip(residual_cells, raw_demand) if demand
    )
    demand_vector = tuple(demand for demand in raw_demand if demand)
    if not required_cells:
        raise RuntimeError("anchors and algebraic identities already cover the plane")
    if args.count_required_only:
        print(
            f"anchor_cells={len(anchor_cells)} residual={len(residual_cells)} "
            f"required={len(required_cells)} algebraic_zero="
            f"{len(residual_cells) - len(required_cells)} period={period}"
        )
        return 0
    if args.cell_profile_input:
        prior_profile = json.loads(args.cell_profile_input.read_text())
        expected_profile = {
            "power": args.power,
            "max_component": args.max_component,
            "anchor_primes": list(requested),
            "anchor_targets": list(normalized_targets),
        }
        for key, expected in expected_profile.items():
            if prior_profile.get(key) != expected:
                raise RuntimeError(
                    f"cell profile {key} does not match this calculation"
                )
        selected_cells = {
            tuple(int(value) for value in record["cell"])
            for record in prior_profile["cells"]
            if float(record["coverage"]) < args.cell_profile_below
        }
        filtered = [
            (cell, demand)
            for cell, demand in zip(required_cells, demand_vector)
            if cell in selected_cells
        ]
        if not filtered:
            raise RuntimeError("cell profile restriction selected no cells")
        required_cells = tuple(cell for cell, _demand in filtered)
        demand_vector = tuple(demand for _cell, demand in filtered)
        print(
            f"cell_profile_restriction selected={len(required_cells)} "
            f"below={args.cell_profile_below}"
        )

    def vector(shared: int, a: int, b: int, target: int) -> tuple[int, ...]:
        return tuple(
            count_for_cell(cell, shared, a, b, target) for cell in required_cells
        )

    if args.evaluate_phase_file:
        phases = {
            int(prime): int(target)
            for prime, target in json.loads(
                args.evaluate_phase_file.read_text()
            ).items()
        }
        for row, target in zip(anchors, normalized_targets):
            prime = int(row["p"])
            if phases.get(prime) != target:
                raise RuntimeError(
                    f"phase target for anchor {prime} does not match "
                    f"--anchor-targets"
                )
        coverage = np.zeros(len(required_cells), dtype=np.float64)
        compatible = 0
        evaluation_patterns: dict[
            tuple[int, int, int, int], "np.ndarray"
        ] = {}
        phase_rows = []
        for index, row in enumerate(candidates):
            if index in anchor_indices:
                continue
            h = int(row["h"])
            p = int(row["p"])
            a = int(row["a"])
            b = int(row["b"])
            try:
                residue, modulus = power_target_congruence(
                    h, p, args.power
                )
            except RuntimeError:
                continue
            if p not in phases:
                raise RuntimeError(f"phase assignment is missing prime {p}")
            target = phases[p] % h
            if target % modulus != residue:
                raise RuntimeError(f"phase for prime {p} is power-incompatible")
            shared = math.gcd(h, period)
            if not args.capacity_sweeps and not args.optimize_phase_output:
                key = (
                    shared,
                    a % shared,
                    b % shared,
                    target % shared,
                )
                if key not in evaluation_patterns:
                    evaluation_patterns[key] = (
                        np.asarray(vector(*key), dtype=np.float64)
                        / np.asarray(demand_vector, dtype=np.float64)
                    )
                coverage += evaluation_patterns[key] * (shared / h)
                compatible += 1
                continue
            target_period = math.lcm(shared, modulus)
            local_targets = sorted(
                {
                    value % shared
                    for value in range(residue, target_period, modulus)
                }
            )
            options = []
            option_by_counts = {}
            chosen = None
            for local_target in local_targets:
                key = (
                    shared,
                    a % shared,
                    b % shared,
                    local_target,
                )
                counts = vector(*key)
                if counts in option_by_counts:
                    if local_target == target % shared:
                        chosen = option_by_counts[counts]
                        _old_target, old_key = options[chosen]
                        options[chosen] = (local_target, old_key)
                    continue
                option_by_counts[counts] = len(options)
                if key not in evaluation_patterns:
                    evaluation_patterns[key] = (
                        np.asarray(counts, dtype=np.float64)
                        / np.asarray(demand_vector, dtype=np.float64)
                    )
                if local_target == target % shared:
                    chosen = len(options)
                options.append((local_target, key))
            if chosen is None:
                raise RuntimeError(f"current phase has no local option for prime {p}")
            scale = shared / h
            coverage += evaluation_patterns[options[chosen][1]] * scale
            phase_rows.append(
                {
                    "p": p,
                    "h": h,
                    "residue": residue,
                    "modulus": modulus,
                    "shared": shared,
                    "options": options,
                    "chosen": chosen,
                    "scale": scale,
                }
            )
            compatible += 1
        if args.capacity_sweeps:
            rng = random.Random(args.capacity_seed)

            def deficit_score(values: "np.ndarray") -> float:
                deficits = np.maximum(args.capacity_threshold - values, 0.0)
                return float(deficits @ deficits)

            mutable = [
                row_index
                for row_index, phase_row in enumerate(phase_rows)
                if len(phase_row["options"]) > 1
            ]
            for sweep in range(1, args.capacity_sweeps + 1):
                rng.shuffle(mutable)
                moves = 0
                score = deficit_score(coverage)
                for row_index in mutable:
                    phase_row = phase_rows[row_index]
                    old_index = int(phase_row["chosen"])
                    old_key = phase_row["options"][old_index][1]
                    scale = float(phase_row["scale"])
                    base = coverage - evaluation_patterns[old_key] * scale
                    best_index = old_index
                    best_score = score
                    for option_index, (_target, key) in enumerate(
                        phase_row["options"]
                    ):
                        if option_index == old_index:
                            continue
                        trial = base + evaluation_patterns[key] * scale
                        trial_score = deficit_score(trial)
                        if trial_score < best_score - 1e-15:
                            best_score = trial_score
                            best_index = option_index
                    if best_index != old_index:
                        best_key = phase_row["options"][best_index][1]
                        coverage = base + evaluation_patterns[best_key] * scale
                        phase_row["chosen"] = best_index
                        score = best_score
                        moves += 1
                print(
                    f"capacity_sweep={sweep} moves={moves} "
                    f"score={score:.15g} min={float(coverage.min()):.15f} "
                    f"below1={int(np.count_nonzero(coverage < 1.0 - 1e-12))}"
                )
                if not moves:
                    break
        if args.optimize_phase_output:
            from round_fractional_phases import combine_congruences

            optimized = dict(phases)
            for phase_row in phase_rows:
                p = int(phase_row["p"])
                h = int(phase_row["h"])
                shared = int(phase_row["shared"])
                residue = int(phase_row["residue"])
                modulus = int(phase_row["modulus"])
                current = phases[p] % h
                local_target = int(
                    phase_row["options"][int(phase_row["chosen"])][0]
                )
                hidden_coprime = 1
                visible_primes = set(exact_uncovered.factor(shared))
                for prime, exponent in exact_uncovered.factor(h).items():
                    if prime not in visible_primes:
                        hidden_coprime *= prime**exponent
                combined, combined_modulus = combine_congruences(
                    local_target, shared, current % hidden_coprime, hidden_coprime
                )
                combined, combined_modulus = combine_congruences(
                    combined, combined_modulus, residue, modulus
                )
                representatives = range(combined, h, combined_modulus)
                optimized[p] = min(
                    representatives,
                    key=lambda value: min(
                        (value - current) % h, (current - value) % h
                    ),
                )
            args.optimize_phase_output.write_text(
                json.dumps(
                    {str(prime): target for prime, target in optimized.items()},
                    indent=2,
                    sort_keys=True,
                )
                + "\n"
            )
        quantiles = np.quantile(
            coverage, [0.0, 0.01, 0.05, 0.10, 0.50, 0.90, 1.0]
        )
        if args.profile_output:
            records = [
                {
                    "cell": list(cell),
                    "coverage": float(value),
                    "demand": int(demand),
                }
                for cell, value, demand in zip(
                    required_cells, coverage, demand_vector
                )
            ]
            records.sort(key=lambda item: item["coverage"])
            args.profile_output.write_text(
                json.dumps(
                    {
                        "pool": str(args.pool),
                        "power": args.power,
                        "max_component": args.max_component,
                        "anchor_primes": list(requested),
                        "anchor_targets": list(normalized_targets),
                        "anchor_period": anchor_period,
                        "algebraic_primes": list(algebraic_primes),
                        "cells": records,
                    },
                    indent=2,
                )
                + "\n"
            )
        print(
            f"pool={args.pool} power={args.power} candidates={len(candidates)} "
            f"compatible_nonanchors={compatible} "
            f"anchors={[(row['p'], row['h']) for row in anchors]} "
            f"normalized_targets={normalized_targets} "
            f"required={len(required_cells)} period={period}"
        )
        print(
            "density_profile "
            f"min={quantiles[0]:.15f} p01={quantiles[1]:.15f} "
            f"p05={quantiles[2]:.15f} p10={quantiles[3]:.15f} "
            f"median={quantiles[4]:.15f} p90={quantiles[5]:.15f} "
            f"max={quantiles[6]:.15f} mean={float(coverage.mean()):.15f} "
            f"below1={int(np.count_nonzero(coverage < 1.0 - 1e-12))}"
        )
        return 0

    pattern_cache: dict[tuple[int, int, int, int], tuple[int, ...]] = {}
    rows = []
    incompatible = 0
    for index, row in enumerate(candidates):
        if index in anchor_indices:
            continue
        h = int(row["h"])
        p = int(row["p"])
        a = int(row["a"])
        b = int(row["b"])
        shared = math.gcd(h, period)
        try:
            residue, modulus = power_target_congruence(h, p, args.power)
        except RuntimeError:
            incompatible += 1
            continue
        target_period = math.lcm(shared, modulus)
        targets = sorted(
            {
                target % shared
                for target in range(target_period)
                if target % modulus == residue
            }
        )
        options = []
        seen = set()
        for target in targets:
            key = (shared, a % shared, b % shared, target)
            counts = pattern_cache.get(key)
            if counts is None:
                counts = vector(shared, a, b, target)
                pattern_cache[key] = counts
            if counts in seen:
                continue
            seen.add(counts)
            options.append((shared, h, target, counts))
        if not options:
            raise RuntimeError(f"no valid target for prime {p}")
        rows.append((p, options))

    columns = [
        (row_index, option_index)
        for row_index, (_p, options) in enumerate(rows)
        for option_index in range(len(options))
    ]
    column_by_option = {
        option: column for column, option in enumerate(columns)
    }
    eq = coo_matrix(
        (
            np.ones(len(columns)),
            ([row_index for row_index, _ in columns], range(len(columns))),
        ),
        shape=(len(rows), len(columns) + 1),
    ).tocsr()
    ub_rows = []
    ub_columns = []
    ub_values = []
    for column, (row_index, option_index) in enumerate(columns):
        shared, h, _target, counts = rows[row_index][1][option_index]
        for cell_index, count in enumerate(counts):
            if count:
                ub_rows.append(cell_index)
                ub_columns.append(column)
                ub_values.append(-count * shared / (h * demand_vector[cell_index]))
    for cell_index in range(len(required_cells)):
        ub_rows.append(cell_index)
        ub_columns.append(len(columns))
        ub_values.append(1.0)
    ub = coo_matrix(
        (ub_values, (ub_rows, ub_columns)),
        shape=(len(required_cells), len(columns) + 1),
    ).tocsr()
    objective = np.zeros(len(columns) + 1)
    objective[-1] = -1.0
    result = linprog(
        objective,
        A_ub=ub,
        b_ub=np.zeros(len(required_cells)),
        A_eq=eq,
        b_eq=np.ones(len(rows)),
        bounds=[(0.0, 1.0)] * len(columns) + [(None, None)],
        method=args.linprog_method,
    )
    if not result.success:
        raise RuntimeError(result.message)
    if args.fractional_output:
        primal_rows = []
        for row_index, (prime, options) in enumerate(rows):
            mixture = []
            for option_index, (shared, h, target, _counts) in enumerate(options):
                column = column_by_option[(row_index, option_index)]
                weight = float(result.x[column])
                if weight > 1e-12:
                    mixture.append(
                        {
                            "shared": shared,
                            "h": h,
                            "target_mod_shared": target,
                            "weight": weight,
                        }
                    )
            primal_rows.append({"p": prime, "options": mixture})
        args.fractional_output.write_text(
            json.dumps(
                {
                    "pool": str(args.pool),
                    "power": args.power,
                    "anchor_primes": requested,
                    "normalized_targets": normalized_targets,
                    "lp_ratio": float(result.x[-1]),
                    "rows": primal_rows,
                },
                indent=2,
            )
            + "\n"
        )
    if args.float_only:
        print(
            f"pool={args.pool} power={args.power} candidates={len(candidates)} "
            f"compatible={len(rows) + len(anchors)} incompatible={incompatible} "
            f"anchors={[(row['p'], row['h']) for row in anchors]} "
            f"normalized_targets={normalized_targets} cells={len(anchor_cells)} "
            f"required={len(required_cells)} algebraic_primes={algebraic_primes} "
            f"sophie_germain={sophie_germain} period={period} "
            f"local_patterns={local_count.cache_info().currsize} "
            f"patterns={len(pattern_cache)} columns={len(columns)}"
        )
        print(f"lp_ratio={result.x[-1]:.15f} result=FLOAT_ONLY")
        return 0

    weights = [
        Fraction(-float(value)).limit_denominator(args.weight_denominator)
        for value in result.ineqlin.marginals
    ]
    weight_sum = sum(weights, Fraction())
    weights = [weight / weight_sum for weight in weights]
    exact_upper = Fraction()
    for _p, options in rows:
        exact_upper += max(
            sum(
                (
                    weights[cell_index]
                    * Fraction(count * shared, h * demand_vector[cell_index])
                    for cell_index, count in enumerate(counts)
                ),
                Fraction(),
            )
            for shared, h, _target, counts in options
        )
    impossible = exact_upper < 1
    print(
        f"pool={args.pool} power={args.power} candidates={len(candidates)} "
        f"compatible={len(rows) + len(anchors)} incompatible={incompatible} "
        f"anchors={[(row['p'], row['h']) for row in anchors]} "
        f"normalized_targets={normalized_targets} cells={len(anchor_cells)} "
        f"required={len(required_cells)} algebraic_primes={algebraic_primes} "
        f"sophie_germain={sophie_germain} period={period} "
        f"local_patterns={local_count.cache_info().currsize} "
        f"patterns={len(pattern_cache)} columns={len(columns)}"
    )
    print(
        f"lp_ratio={result.x[-1]:.15f} "
        f"exact_ratio_upper={float(exact_upper):.15f} "
        f"bits=({exact_upper.numerator.bit_length()},"
        f"{exact_upper.denominator.bit_length()})"
    )
    print(f"result={'IMPOSSIBLE' if impossible else 'INCONCLUSIVE'}")
    return 2 if impossible else 0


if __name__ == "__main__":
    raise SystemExit(main())
