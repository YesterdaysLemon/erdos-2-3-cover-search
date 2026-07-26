#!/usr/bin/env python3
"""Exact lazy master for a product of prime-power affine components.

Rows with the same normalized component directions are interchangeable, so
the master chooses bounded sets of affine CRT flats.  It learns cover clauses
only for concrete counterexamples returned by the exact CRT checker.  A SAT
result is emitted solely after both the PySAT checker and an independent Z3
bit-vector checker prove that no uncovered point exists.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import exact_uncovered
import exact_uncovered_z3_bv


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument(
        "--components",
        required=True,
        help="comma-separated pairwise-coprime prime-power moduli",
    )
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--checker-solver", default="cadical195")
    parser.add_argument("--rounds", type=int, default=10000)
    parser.add_argument("--batch", type=int, default=500)
    parser.add_argument(
        "--diversity-coordinate-moduli",
        default="",
        help="coordinate moduli used to diversify each checker batch",
    )
    parser.add_argument("--diversity-quota", type=int, default=1)
    parser.add_argument(
        "--no-translation-symmetry",
        action="store_true",
        help=(
            "disable the proof-safe zero-phase anchors for two independent "
            "simple direction groups in each component"
        ),
    )
    parser.add_argument(
        "--resume-checkpoint",
        type=Path,
        help="restore learned points and use saved phases as solver hints",
    )
    parser.add_argument("--checkpoint-output", type=Path)
    parser.add_argument("--phase-output", type=Path)
    parser.add_argument("--result-output", type=Path)
    args = parser.parse_args()

    components = tuple(
        int(value) for value in args.components.split(",") if value
    )
    if not components or len(set(components)) != len(components):
        raise SystemExit("--components must be a nonempty unique list")
    if any(value < 2 for value in components):
        raise SystemExit("component moduli must be at least two")
    for index, left in enumerate(components):
        factors = exact_uncovered.factor(left)
        if len(factors) != 1:
            raise SystemExit(f"component {left} is not a prime power")
        if any(math.gcd(left, right) != 1 for right in components[index + 1 :]):
            raise SystemExit("component moduli must be pairwise coprime")
    period = math.prod(components)
    diversity_moduli = tuple(
        int(value)
        for value in args.diversity_coordinate_moduli.split(",")
        if value
    )
    if args.rounds < 1 or args.batch < 1 or args.diversity_quota < 1:
        raise SystemExit("rounds, batch, and diversity quota must be positive")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    payload = json.loads(args.pool.read_text())
    rows = [
        row
        for row in payload["choices"]
        if int(row["h"]) > 1 and period % int(row["h"]) == 0
    ]
    groups: dict[tuple, list[tuple[dict, tuple]]] = defaultdict(list)
    for row in rows:
        h = int(row["h"])
        if int(row["target_modulus"]) != 1:
            raise RuntimeError("symmetry reduction requires target_modulus=1")
        component_data = []
        for modulus in components:
            if h % modulus:
                continue
            a = int(row["a"]) % modulus
            b = int(row["b"]) % modulus
            if math.gcd(a, modulus) == 1:
                scale = pow(a, -1, modulus)
            elif math.gcd(b, modulus) == 1:
                scale = pow(b, -1, modulus)
            else:
                raise RuntimeError(
                    f"row {row['p']} has no unit normal coefficient "
                    f"modulo component {modulus}"
                )
            direction = (
                a * scale % modulus,
                b * scale % modulus,
            )
            component_data.append((modulus, direction, scale))
        component_product = math.prod(
            modulus for modulus, _direction, _scale in component_data
        )
        if component_product != h:
            raise RuntimeError(
                f"row modulus {h} is not a product of whole components"
            )
        signature = tuple(
            (modulus, direction)
            for modulus, direction, _scale in component_data
        )
        groups[signature].append((row, tuple(component_data)))

    vpool = IDPool()
    solver = Solver(name=args.solver)
    flat_var: dict[tuple[tuple, tuple[int, ...]], int] = {}
    variable_flat: dict[int, tuple[tuple, tuple[int, ...]]] = {}
    group_variables: dict[tuple, list[int]] = {}
    clause_count = 0
    for signature, members in groups.items():
        target_space = itertools.product(
            *(range(modulus) for modulus, _direction in signature)
        )
        variables = []
        for targets in target_space:
            variable = vpool.id(("flat", signature, targets))
            flat_var[(signature, targets)] = variable
            variable_flat[variable] = (signature, targets)
            variables.append(variable)
        group_variables[signature] = variables
        # We may assume every interchangeable row uses a distinct phase until
        # the phase space is exhausted.  If two members duplicate a flat while
        # an unused flat exists, moving one duplicate to the unused flat keeps
        # the old flat selected and only adds covered points.  Thus a cover
        # exists iff one exists with exactly min(members, flats) selections in
        # each group.
        selected_count = min(len(members), len(variables))
        if selected_count == len(variables):
            for variable in variables:
                solver.add_clause([variable])
                clause_count += 1
        else:
            encoding = CardEnc.equals(
                variables,
                bound=selected_count,
                vpool=vpool,
                encoding=EncType.seqcounter,
            )
            solver.append_formula(encoding.clauses)
            clause_count += len(encoding.clauses)

    symmetry_fixed_groups = []
    if not args.no_translation_symmetry:
        # Translation acts independently on each pairwise-coprime CRT
        # component.  Every simple direction group is nonempty.  For two
        # directions with unit determinant, translate one selected flat from
        # each group to target zero.  Repeating this independently for every
        # component is WLOG and removes the two-dimensional translation orbit.
        for modulus in components:
            simple_signatures = sorted(
                signature
                for signature in groups
                if len(signature) == 1 and signature[0][0] == modulus
            )
            anchor_pair = None
            for left_index, left in enumerate(simple_signatures):
                left_direction = left[0][1]
                for right in simple_signatures[left_index + 1 :]:
                    right_direction = right[0][1]
                    determinant = (
                        left_direction[0] * right_direction[1]
                        - left_direction[1] * right_direction[0]
                    ) % modulus
                    if math.gcd(determinant, modulus) == 1:
                        anchor_pair = (left, right)
                        break
                if anchor_pair is not None:
                    break
            if anchor_pair is None:
                continue
            for signature in anchor_pair:
                zero_targets = tuple(0 for _entry in signature)
                solver.add_clause([flat_var[(signature, zero_targets)]])
                clause_count += 1
                symmetry_fixed_groups.append(signature)

    def cover_clause(k: int, l: int) -> list[int]:
        clause = []
        for signature in groups:
            targets = tuple(
                (direction[0] * k + direction[1] * l) % modulus
                for modulus, direction in signature
            )
            clause.append(flat_var[(signature, targets)])
        return clause

    learned_points: set[tuple[int, int]] = set()

    resume_payload = None
    if args.resume_checkpoint and args.resume_checkpoint.exists():
        resume_payload = json.loads(args.resume_checkpoint.read_text())
        if (
            resume_payload.get("pool") != str(args.pool)
            or tuple(resume_payload.get("components", ())) != components
        ):
            raise RuntimeError("resume checkpoint does not match this problem")
        for raw_k, raw_l in resume_payload.get("learned_point_values", ()):
            point = (int(raw_k) % period, int(raw_l) % period)
            if point in learned_points:
                continue
            solver.add_clause(cover_clause(*point))
            learned_points.add(point)
            clause_count += 1
        saved_phases = {
            int(prime): int(target)
            for prime, target in resume_payload.get("phases", {}).items()
        }
        hints = []
        for signature, members in groups.items():
            for row, component_data in members:
                prime = int(row["p"])
                if prime not in saved_phases:
                    continue
                target = saved_phases[prime]
                canonical_targets = tuple(
                    scale * target % modulus
                    for modulus, _direction, scale in component_data
                )
                hints.append(flat_var[(signature, canonical_targets)])
        try:
            solver.set_phases(hints)
        except NotImplementedError:
            pass

    def decode(model: list[int]) -> tuple[dict[str, int], list[dict]]:
        positive = {literal for literal in model if literal > 0}
        selected: dict[tuple, list[tuple[int, ...]]] = defaultdict(list)
        for variable in positive & variable_flat.keys():
            signature, targets = variable_flat[variable]
            selected[signature].append(targets)
        phases: dict[str, int] = {}
        checked_rows = []
        for signature, members in groups.items():
            targets_list = selected[signature]
            expected = min(len(members), len(group_variables[signature]))
            if len(targets_list) != expected:
                raise AssertionError(
                    "model does not saturate a direction group"
                )
            for index, (row, component_data) in enumerate(members):
                canonical_targets = targets_list[
                    min(index, len(targets_list) - 1)
                ]
                congruences = []
                for canonical_target, (
                    modulus,
                    _direction,
                    scale,
                ) in zip(canonical_targets, component_data):
                    congruences.append(
                        (
                            canonical_target
                            * pow(scale, -1, modulus)
                            % modulus,
                            modulus,
                        )
                    )
                target = exact_uncovered.crt(congruences)
                phases[str(int(row["p"]))] = target
                checked_row = dict(row)
                checked_row["c"] = target
                checked_rows.append(checked_row)
        return phases, checked_rows

    started = time.monotonic()
    print(
        f"components={components} period={period} rows={len(rows)} "
        f"groups={len(groups)} variables={vpool.top} "
        f"capacity_clauses={clause_count} "
        f"symmetry_anchors={len(symmetry_fixed_groups)}",
        flush=True,
    )
    final_phases = None
    final_rows = None
    last_meta = None
    for round_no in range(1, args.rounds + 1):
        solve_started = time.monotonic()
        master_sat = solver.solve()
        solve_seconds = time.monotonic() - solve_started
        if not master_sat:
            result = {
                "pool": str(args.pool),
                "components": list(components),
                "period": period,
                "row_count": len(rows),
                "group_count": len(groups),
                "symmetry_fixed_groups": symmetry_fixed_groups,
                "master_sat": False,
                "cover_sat": False,
                "round": round_no,
                "learned_points": len(learned_points),
                "clause_count": clause_count,
                "elapsed_seconds": time.monotonic() - started,
            }
            if args.result_output:
                args.result_output.write_text(
                    json.dumps(result, indent=2) + "\n"
                )
            print(
                f"result=UNSAT round={round_no} "
                f"learned={len(learned_points)}",
                flush=True,
            )
            solver.delete()
            return 2

        phases, checked_rows = decode(solver.get_model())
        misses, last_meta = exact_uncovered.find_uncovered(
            checked_rows,
            max_component=max(components),
            limit=args.batch,
            solver_name=args.checker_solver,
            diversity_coordinate_moduli=diversity_moduli,
            diversity_quota=args.diversity_quota,
        )
        if not misses:
            z3_misses, z3_meta = exact_uncovered_z3_bv.find_uncovered(
                checked_rows,
                max_component=max(components),
                limit=1,
            )
            if z3_misses:
                raise AssertionError(
                    "independent Z3 checker found a missed point"
                )
            final_phases = phases
            final_rows = checked_rows
            if args.phase_output:
                args.phase_output.write_text(json.dumps(phases) + "\n")
            result = {
                "pool": str(args.pool),
                "components": list(components),
                "period": period,
                "row_count": len(rows),
                "group_count": len(groups),
                "symmetry_fixed_groups": symmetry_fixed_groups,
                "master_sat": True,
                "cover_sat": True,
                "round": round_no,
                "learned_points": len(learned_points),
                "clause_count": clause_count,
                "primary_verification": last_meta,
                "independent_verification": z3_meta,
                "elapsed_seconds": time.monotonic() - started,
            }
            if args.result_output:
                args.result_output.write_text(
                    json.dumps(result, indent=2) + "\n"
                )
            print(
                f"result=SAT round={round_no} "
                f"learned={len(learned_points)}",
                flush=True,
            )
            solver.delete()
            return 0

        added = 0
        for k, l in misses:
            point = (k % period, l % period)
            if point in learned_points:
                continue
            clause = cover_clause(*point)
            solver.add_clause(clause)
            learned_points.add(point)
            clause_count += 1
            added += 1
        if not added:
            raise AssertionError("checker returned no new master lesson")
        if args.checkpoint_output:
            args.checkpoint_output.write_text(
                json.dumps(
                    {
                        "pool": str(args.pool),
                        "components": list(components),
                        "period": period,
                        "round": round_no,
                        "learned_points": len(learned_points),
                        "learned_point_values": sorted(learned_points),
                        "clause_count": clause_count,
                        "last_master_solve_seconds": solve_seconds,
                        "last_checker": last_meta,
                        "phases": phases,
                    },
                    indent=2,
                )
                + "\n"
            )
        print(
            f"round={round_no} master_s={solve_seconds:.3f} "
            f"misses={len(misses)} added={added} "
            f"learned={len(learned_points)}",
            flush=True,
        )

    result = {
        "pool": str(args.pool),
        "components": list(components),
        "period": period,
        "row_count": len(rows),
        "group_count": len(groups),
        "symmetry_fixed_groups": symmetry_fixed_groups,
        "master_sat": True,
        "cover_sat": None,
        "round": args.rounds,
        "learned_points": len(learned_points),
        "clause_count": clause_count,
        "last_checker": last_meta,
        "elapsed_seconds": time.monotonic() - started,
    }
    if args.result_output:
        args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"result=INCOMPLETE rounds={args.rounds} "
        f"learned={len(learned_points)}",
        flush=True,
    )
    solver.delete()
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
