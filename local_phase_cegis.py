#!/usr/bin/env python3
"""Compact multivalued CEGIS master for affine-line covers.

Each prime owns one integer target c modulo h.  Unlike the binary SAT master,
this keeps sampled target values in a dense uint32 matrix and repairs uncovered
sample points by min-conflicts moves.  It is heuristic only on the synthesis
side: a construction is emitted solely after exact_uncovered proves UNSAT.
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import os
import random
import sys
import time
from pathlib import Path

import exact_greedy
import exact_uncovered

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def build_targets(points, candidates, np):
    started = time.monotonic()
    # Repair is column-oriented: each move inspects all retained points for
    # one candidate.  Fortran order keeps those candidate columns contiguous.
    targets = np.empty(
        (len(points), len(candidates)), dtype=np.uint32, order="F"
    )
    uint64_flags = [
        all(0 <= coordinate < (1 << 64) for coordinate in point)
        for point in points
    ]
    uint64_count = sum(uint64_flags)
    uint64_prefix = (
        uint64_count == len(points)
        or uint64_flags == [True] * uint64_count + [False] * (
            len(points) - uint64_count
        )
    )
    if uint64_count:
        uint64_points = (
            points[:uint64_count]
            if uint64_prefix
            else [point for point, flag in zip(points, uint64_flags) if flag]
        )
        ks = np.fromiter(
            (point[0] for point in uint64_points),
            dtype=np.uint64,
            count=uint64_count,
        )
        ls = np.fromiter(
            (point[1] for point in uint64_points),
            dtype=np.uint64,
            count=uint64_count,
        )
        uint64_selector = (
            slice(0, uint64_count)
            if uint64_prefix
            else np.flatnonzero(uint64_flags)
        )
    big_points = (
        points[uint64_count:]
        if uint64_prefix
        else [point for point, flag in zip(points, uint64_flags) if not flag]
    )
    big_selector = (
        slice(uint64_count, len(points))
        if uint64_prefix
        else np.flatnonzero(np.logical_not(uint64_flags))
    )
    for column, (h, _p, a, b, _ord2, _ord3) in enumerate(candidates):
        if uint64_count:
            targets[uint64_selector, column] = (
                (a * (ks % h) + b * (ls % h)) % h
            ).astype(np.uint32)
        if big_points:
            targets[big_selector, column] = np.fromiter(
                ((a * (k % h) + b * (l % h)) % h for k, l in big_points),
                dtype=np.uint32,
                count=len(big_points),
            )
    return targets, time.monotonic() - started


def merge_congruences(residue_a, modulus_a, residue_b, modulus_b):
    """Return the intersection of two congruences, or None if it is empty."""
    common = math.gcd(modulus_a, modulus_b)
    if (residue_b - residue_a) % common:
        return None
    left = modulus_a // common
    right = modulus_b // common
    if right == 1:
        multiplier = 0
    else:
        multiplier = (
            ((residue_b - residue_a) // common) * pow(left, -1, right)
        ) % right
    modulus = math.lcm(modulus_a, modulus_b)
    return (residue_a + modulus_a * multiplier) % modulus, modulus


def repair(
    targets,
    assignment,
    mutable,
    rng,
    np,
    max_steps,
    sample_size,
    valid_moduli=None,
    valid_residues=None,
    candidate_moduli=None,
    required_coverage=1,
):
    cover = np.count_nonzero(targets == assignment, axis=1).astype(np.int32)
    start_uncovered = int(np.count_nonzero(cover < required_coverage))
    best_uncovered = start_uncovered
    best_deficit = int(
        np.maximum(required_coverage - cover, 0).sum()
    )
    best_assignment = assignment.copy()
    stagnant = 0
    for step in range(max_steps + 1):
        uncovered = np.flatnonzero(cover < required_coverage)
        if not len(uncovered):
            return True, step, start_uncovered, 0
        if step == max_steps:
            break
        point = int(uncovered[rng.randrange(len(uncovered))])
        size = min(sample_size, len(mutable))
        selected = rng.sample(mutable, size)
        selected_array = np.asarray(selected, dtype=np.int32)
        new_values = targets[point, selected_array]
        changed = new_values != assignment[selected_array]
        if valid_moduli is not None:
            changed &= (
                new_values % valid_moduli[selected_array]
                == valid_residues[selected_array]
            )
        if not np.any(changed):
            continue
        selected_array = selected_array[changed]
        new_values = new_values[changed]
        slab = targets[:, selected_array]
        old_matches = slab == assignment[selected_array]
        new_matches = slab == new_values
        losses = np.count_nonzero(
            old_matches & (cover[:, None] <= required_coverage),
            axis=0,
        )
        gains = np.count_nonzero(
            new_matches & (cover[:, None] < required_coverage),
            axis=0,
        )
        scores = gains - losses
        best_score = int(scores.max())
        choices = np.flatnonzero(scores == best_score)
        choice = int(choices[rng.randrange(len(choices))])
        candidate = int(selected_array[choice])
        old_value = int(assignment[candidate])
        new_value = int(new_values[choice])
        old_mask = targets[:, candidate] == old_value
        new_mask = targets[:, candidate] == new_value
        cover -= old_mask.astype(np.int32)
        cover += new_mask.astype(np.int32)
        assignment[candidate] = new_value

        current = int(np.count_nonzero(cover < required_coverage))
        current_deficit = int(
            np.maximum(required_coverage - cover, 0).sum()
        )
        if current_deficit < best_deficit:
            best_deficit = current_deficit
            best_uncovered = current
            best_assignment[:] = assignment
            stagnant = 0
        else:
            stagnant += 1
        if stagnant >= 2500:
            # A small exact-domain-neutral perturbation breaks repair cycles.
            for candidate in rng.sample(mutable, min(25, len(mutable))):
                old_value = int(assignment[candidate])
                if valid_moduli is None:
                    new_value = rng.randrange(int(candidate_moduli[candidate]))
                else:
                    modulus = int(valid_moduli[candidate])
                    residue = int(valid_residues[candidate])
                    choices = (int(candidate_moduli[candidate]) - 1 - residue) // modulus + 1
                    new_value = residue + modulus * rng.randrange(choices)
                if new_value == old_value:
                    continue
                old_mask = targets[:, candidate] == old_value
                new_mask = targets[:, candidate] == new_value
                cover -= old_mask.astype(np.int32)
                cover += new_mask.astype(np.int32)
                assignment[candidate] = new_value
            stagnant = 0
    assignment[:] = best_assignment
    return False, max_steps, start_uncovered, best_uncovered


def repair_coordinate(
    targets,
    assignment,
    mutable,
    rng,
    np,
    max_steps,
    sample_size,
    valid_moduli,
    valid_residues,
    candidate_moduli,
    required_coverage=1,
):
    """Coordinate descent using the best valid target on all missed points."""
    cover = np.count_nonzero(targets == assignment, axis=1).astype(np.int32)
    start_uncovered = int(np.count_nonzero(cover < required_coverage))
    best_uncovered = start_uncovered
    best_deficit = int(
        np.maximum(required_coverage - cover, 0).sum()
    )
    best_assignment = assignment.copy()
    for step in range(max_steps + 1):
        uncovered_mask = cover < required_coverage
        if not np.any(uncovered_mask):
            return True, step, start_uncovered, 0
        if step == max_steps:
            break
        selected = rng.sample(mutable, min(sample_size, len(mutable)))
        best_move = None
        best_score = None
        for candidate in selected:
            values = targets[:, candidate]
            h = int(candidate_moduli[candidate])
            modulus = int(valid_moduli[candidate])
            residue = int(valid_residues[candidate])
            observed = values[uncovered_mask]
            valid_observed = observed[
                observed % modulus == residue
            ]
            if not len(valid_observed):
                continue
            # Large subgroup orders made the old minlength=h bincount allocate
            # mostly empty arrays for every coordinate.  Count only targets
            # that actually occur on the current missed points unless the
            # domain is already compact.
            if h <= max(256, 4 * len(valid_observed)):
                counts = np.bincount(valid_observed, minlength=h)
                valid_counts = counts[residue::modulus]
                target_index = int(np.argmax(valid_counts))
                new_value = residue + modulus * target_index
                gains = int(valid_counts[target_index])
            else:
                observed_targets, observed_counts = np.unique(
                    valid_observed, return_counts=True
                )
                target_index = int(np.argmax(observed_counts))
                new_value = int(observed_targets[target_index])
                gains = int(observed_counts[target_index])
            old_value = int(assignment[candidate])
            if new_value == old_value:
                continue
            losses = int(
                np.count_nonzero(
                    (cover <= required_coverage) & (values == old_value)
                )
            )
            score = gains - losses
            if best_score is None or score > best_score:
                best_score = score
                best_move = (candidate, old_value, new_value, values)
        if best_move is None:
            continue
        candidate, old_value, new_value, values = best_move
        cover -= (values == old_value).astype(np.int32)
        cover += (values == new_value).astype(np.int32)
        assignment[candidate] = new_value
        current = int(np.count_nonzero(cover < required_coverage))
        current_deficit = int(
            np.maximum(required_coverage - cover, 0).sum()
        )
        if current_deficit < best_deficit:
            best_deficit = current_deficit
            best_uncovered = current
            best_assignment[:] = assignment
    assignment[:] = best_assignment
    return False, max_steps, start_uncovered, best_uncovered


def repair_component(
    targets,
    assignment,
    mutable,
    rng,
    np,
    max_steps,
    sample_size,
    valid_moduli,
    valid_residues,
    candidate_moduli,
    required_coverage=1,
):
    """Greedily retarget one prime-power CRT digit, then fall back.

    A component move preserves every coprime digit of a phase.  This gives a
    much smoother neighborhood than replacing the entire phase by one point's
    target, especially when a candidate modulus has several components.
    """
    cover = np.count_nonzero(targets == assignment, axis=1).astype(np.int32)
    start_uncovered = int(np.count_nonzero(cover < required_coverage))
    if not start_uncovered:
        return True, 0, 0, 0
    components = {
        candidate: tuple(
            prime**exponent
            for prime, exponent in exact_uncovered.factor(
                int(candidate_moduli[candidate])
            ).items()
        )
        for candidate in mutable
    }
    moves = 0
    while moves < max_steps:
        selected = rng.sample(
            mutable, min(sample_size, len(mutable))
        )
        rng.shuffle(selected)
        sweep_moves = 0
        for candidate in selected:
            h = int(candidate_moduli[candidate])
            values = targets[:, candidate]
            powers = list(components[candidate])
            rng.shuffle(powers)
            for component in powers:
                if moves >= max_steps:
                    break
                other = h // component
                old_value = int(assignment[candidate])
                old_hits = values == old_value
                base_cover = cover - old_hits
                needed = base_cover < required_coverage
                other_residue = old_value % other if other > 1 else 0
                compatible = needed
                if other > 1:
                    compatible = compatible & (
                        values % other == other_residue
                    )
                component_constraint = math.gcd(
                    component, int(valid_moduli[candidate])
                )
                legal_residue = (
                    int(valid_residues[candidate])
                    % component_constraint
                )
                component_values = values % component
                if component_constraint > 1:
                    compatible = compatible & (
                        component_values % component_constraint
                        == legal_residue
                    )
                observed = component_values[compatible]
                if not len(observed):
                    continue
                choices, counts = np.unique(
                    observed, return_counts=True
                )
                best_score = int(counts.max())
                old_component = old_value % component
                old_score = int(
                    np.count_nonzero(
                        compatible
                        & (component_values == old_component)
                    )
                )
                if best_score <= old_score:
                    continue
                best_components = choices[counts == best_score]
                new_component = int(
                    best_components[
                        rng.randrange(len(best_components))
                    ]
                )
                if other == 1:
                    new_value = new_component
                elif component == 1:
                    new_value = other_residue
                else:
                    multiplier = (
                        (new_component - other_residue)
                        * pow(other, -1, component)
                    ) % component
                    new_value = (
                        other_residue + other * multiplier
                    ) % h
                if (
                    new_value % int(valid_moduli[candidate])
                    != int(valid_residues[candidate])
                ):
                    raise AssertionError(
                        "component move violates target restriction"
                    )
                new_hits = values == new_value
                cover = base_cover + new_hits
                assignment[candidate] = new_value
                moves += 1
                sweep_moves += 1
                if not np.any(cover < required_coverage):
                    return True, moves, start_uncovered, 0
        if not sweep_moves:
            break

    # Whole-phase min-conflicts remains useful for escaping a component-local
    # optimum.  Component moves were strictly deficit-improving, so the
    # current assignment is a safe warm start for this fallback.
    remaining = max(0, max_steps - moves)
    solved, fallback_steps, _before, after = repair(
        targets,
        assignment,
        mutable,
        rng,
        np,
        remaining,
        sample_size,
        valid_moduli,
        valid_residues,
        candidate_moduli,
        required_coverage,
    )
    return solved, moves + fallback_steps, start_uncovered, after


def random_misses(
    candidates,
    assignment,
    coordinate_bits,
    rng,
    checks,
    batch,
    algebraic_primes=(),
    sophie_germain=False,
):
    misses = []
    tested = []
    for _ in range(checks):
        # A huge exact LCM is expensive to construct for large pools.  Fresh
        # fixed-width integers are ample for heuristic discovery; the exact
        # checker remains the only acceptance test.
        k = rng.getrandbits(coordinate_bits)
        l = rng.getrandbits(coordinate_bits)
        if any(k % prime == 0 and l % prime == 0 for prime in algebraic_primes):
            continue
        if sophie_germain and k % 4 == 2 and l % 4 == 0:
            continue
        tested.append((k, l))
        if all(
            (a * (k % h) + b * (l % h) - int(assignment[index])) % h
            for index, (h, _p, a, b, _ord2, _ord3) in enumerate(candidates)
        ):
            misses.append((k, l))
            if len(misses) >= batch:
                break
    return misses, tested


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_pool", type=Path)
    parser.add_argument("--derived-pool", action="store_true")
    parser.add_argument(
        "--derived-targets",
        action="store_true",
        help=(
            "read target_residue and target_modulus from each derived-pool "
            "row instead of recomputing a perfect-power condition"
        ),
    )
    parser.add_argument(
        "--period",
        type=int,
        default=0,
        help="if nonzero, retain only candidates whose h divides this period",
    )
    parser.add_argument(
        "--min-parent-common",
        type=int,
        default=0,
        help=(
            "for a conditioned derived pool, retain only rows whose recorded "
            "parent_common is at least this value"
        ),
    )
    parser.add_argument("--candidate-count", type=int)
    parser.add_argument(
        "--fixed-prefix",
        type=int,
        default=0,
        help=(
            "keep this many leading compatible candidate phases fixed at "
            "their saved values"
        ),
    )
    parser.add_argument(
        "--fixed-primes",
        default="",
        help="comma-separated primes to keep fixed at their saved phases",
    )
    parser.add_argument(
        "--drop-fixed-covered-points",
        action="store_true",
        help=(
            "omit training points already covered by a fixed candidate; "
            "useful for residual high-band repair"
        ),
    )
    parser.add_argument("--normalize-primes", default="")
    parser.add_argument(
        "--diversity-primes",
        default="",
        help=(
            "comma-separated small fibres whose exact target fingerprint is "
            "blocked between checker witnesses, without fixing their phases"
        ),
    )
    parser.add_argument(
        "--diversity-quota",
        type=int,
        default=1,
        help=(
            "retain this many exact CRT representatives per diversity "
            "fingerprint before blocking the whole fingerprint"
        ),
    )
    parser.add_argument(
        "--diversity-coordinate-moduli",
        default="",
        help=(
            "comma-separated prime-power moduli whose (k,l) residue cells "
            "are diversified by the exact checker"
        ),
    )
    parser.add_argument(
        "--diversity-coordinate-schedule",
        default="",
        help=(
            "semicolon-separated coordinate-modulus stages, for example "
            "'7;7,11;7,13;none'; reuses the retained target matrix"
        ),
    )
    parser.add_argument(
        "--diversity-stage-rounds",
        type=int,
        default=0,
        help="number of exact rounds per scheduled diversity stage",
    )
    parser.add_argument(
        "--diversity-schedule-cycle",
        action="store_true",
        help="repeat the coordinate-diversity schedule after its last stage",
    )
    parser.add_argument(
        "--freeze-target-coordinate-moduli",
        default="",
        help=(
            "comma-separated pairwise-compatible coordinate moduli; preserve "
            "each saved target modulo gcd(h,lcm(moduli)) while repairing the "
            "remaining CRT components"
        ),
    )
    parser.add_argument(
        "--checker-coordinate-cell",
        default="",
        help=(
            "comma-separated modulus:k:l restrictions for a focused exact "
            "cell search, for example 16:0:1,27:8:0; random checks must be 0"
        ),
    )
    parser.add_argument("--points-file", type=Path, required=True)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument("--initial-phase-file", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--checker-checkpoint-file",
        type=Path,
        help=(
            "overwrite a JSON record of the latest random and exact checker "
            "witnesses; with --rounds 1 the phase file is the tested phase"
        ),
    )
    parser.add_argument("--rounds", type=int, default=1000)
    parser.add_argument(
        "--seed",
        type=int,
        help="override the deterministic repair/random-audit seed",
    )
    parser.add_argument(
        "--randomize-mutable-initial",
        action="store_true",
        help=(
            "replace saved mutable phases by seed-deterministic random valid "
            "targets before the first repair; fixed phases are unchanged"
        ),
    )
    parser.add_argument("--max-steps", type=int, default=50000)
    parser.add_argument("--sample-size", type=int, default=256)
    parser.add_argument(
        "--repair-mode",
        choices=("point", "coordinate", "component"),
        default="point",
        help=(
            "point-wise min-conflicts, global best-target coordinate descent, "
            "or prime-power CRT-component descent with point fallback"
        ),
    )
    parser.add_argument(
        "--required-coverage",
        type=int,
        default=1,
        help="minimum number of selected fibres covering every retained point",
    )
    parser.add_argument(
        "--repair-only",
        action="store_true",
        help=(
            "write a repaired finite-sample phase and exit before random or "
            "exact checking; this is heuristic preprocessing, never a cover "
            "certificate"
        ),
    )
    parser.add_argument("--random-checks", type=int, default=100000)
    parser.add_argument("--random-coordinate-bits", type=int, default=128)
    parser.add_argument("--random-batch", type=int, default=500)
    parser.add_argument(
        "--retain-random-tested",
        action="store_true",
        help="retain covered as well as missed random tests to prevent phase overfitting",
    )
    parser.add_argument("--exact-batch", type=int, default=200)
    parser.add_argument(
        "--exact-when-random-misses-at-most",
        type=int,
        default=0,
        help=(
            "invoke the exact checker when the random audit has at most this "
            "many misses; exact witnesses and random misses are both retained"
        ),
    )
    parser.add_argument("--max-component", type=int, default=300000)
    parser.add_argument("--checker-solver", default="cadical195")
    parser.add_argument(
        "--retain-target-matrix-during-exact",
        action="store_true",
        help=(
            "keep the synthesis matrix while the exact checker runs and "
            "append only new lesson columns afterward; trades memory for "
            "faster repeated exact rounds"
        ),
    )
    parser.add_argument("--power", type=int, default=1)
    parser.add_argument(
        "--algebraic-primes",
        help=(
            "comma-separated exact origin sublattices for the checker; with "
            "--derived-targets, defaults to candidate-pool metadata"
        ),
    )
    parser.add_argument(
        "--through-origin-mod",
        type=int,
        default=1,
        help=(
            "restrict every fibre to target c=0 modulo gcd(h,Q); for Q=2 "
            "all even-order lines pass through the origin parity cell"
        ),
    )
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    source_payload = json.loads(args.candidate_pool.read_text())
    if args.derived_targets and not args.derived_pool:
        raise SystemExit("--derived-targets requires --derived-pool")
    derived_target_by_prime = {}
    if args.derived_targets:
        for row in source_payload["choices"]:
            prime = int(row["p"])
            if prime in derived_target_by_prime:
                raise RuntimeError(f"duplicate derived target prime {prime}")
            modulus = int(row["target_modulus"])
            residue = int(row["target_residue"]) % modulus
            if int(row["h"]) % modulus:
                raise RuntimeError(
                    f"derived target modulus does not divide h for prime {prime}"
                )
            derived_target_by_prime[prime] = (residue, modulus)

    candidates = exact_greedy.load_candidates(args.candidate_pool, args.derived_pool)
    if args.candidate_count is not None:
        candidates = candidates[: args.candidate_count]
    if args.period:
        candidates = [item for item in candidates if args.period % item[0] == 0]
    if args.min_parent_common:
        if not args.derived_pool:
            raise SystemExit("--min-parent-common requires --derived-pool")
        parent_common_by_prime = {
            int(row["p"]): int(row.get("parent_common", 1))
            for row in source_payload["choices"]
        }
        candidates = [
            item
            for item in candidates
            if parent_common_by_prime.get(item[1], 1)
            >= args.min_parent_common
        ]
    if args.power < 1:
        raise SystemExit("--power must be positive")
    if args.through_origin_mod < 1:
        raise SystemExit("--through-origin-mod must be positive")
    if args.random_coordinate_bits < 32:
        raise SystemExit("--random-coordinate-bits must be at least 32")
    if args.diversity_quota < 1:
        raise SystemExit("--diversity-quota must be positive")
    if args.required_coverage < 1:
        raise SystemExit("--required-coverage must be positive")
    if args.exact_when_random_misses_at_most < 0:
        raise SystemExit("--exact-when-random-misses-at-most must be nonnegative")
    raw_count = len(candidates)
    if args.derived_targets:
        candidates = [
            item for item in candidates if item[1] in derived_target_by_prime
        ]
    else:
        candidates = [
            item
            for item in candidates
            if ((item[1] - 1) // 2)
            % math.gcd(
                (item[1] - 1) // item[0],
                math.gcd(args.power, item[1] - 1),
            )
            == 0
        ]
    if not 0 <= args.fixed_prefix <= len(candidates):
        raise SystemExit("--fixed-prefix is outside the compatible candidate pool")
    by_prime = {p: index for index, (_h, p, _a, _b, _o2, _o3) in enumerate(candidates)}
    normalizers = tuple(
        int(value) for value in args.normalize_primes.split(",") if value
    )
    diversity_primes = tuple(
        int(value) for value in args.diversity_primes.split(",") if value
    )
    diversity_coordinate_moduli = tuple(
        int(value)
        for value in args.diversity_coordinate_moduli.split(",")
        if value
    )
    diversity_coordinate_schedule = []
    if args.diversity_coordinate_schedule:
        if args.diversity_stage_rounds < 1:
            raise SystemExit(
                "--diversity-stage-rounds must be positive with a schedule"
            )
        for raw_stage in args.diversity_coordinate_schedule.split(";"):
            raw_stage = raw_stage.strip()
            if raw_stage.lower() in ("", "-", "none"):
                diversity_coordinate_schedule.append(())
                continue
            stage = tuple(
                int(value) for value in raw_stage.split(",") if value
            )
            if not stage or any(value < 1 for value in stage):
                raise SystemExit("invalid diversity coordinate schedule stage")
            diversity_coordinate_schedule.append(stage)
    elif args.diversity_stage_rounds:
        raise SystemExit(
            "--diversity-stage-rounds requires a coordinate schedule"
        )
    if args.diversity_schedule_cycle and not diversity_coordinate_schedule:
        raise SystemExit(
            "--diversity-schedule-cycle requires a coordinate schedule"
        )
    freeze_target_coordinate_moduli = tuple(
        int(value)
        for value in args.freeze_target_coordinate_moduli.split(",")
        if value
    )
    if any(value < 1 for value in freeze_target_coordinate_moduli):
        raise SystemExit("--freeze-target-coordinate-moduli must be positive")
    freeze_target_period = math.lcm(
        1, *freeze_target_coordinate_moduli
    )
    checker_coordinate_cell = []
    for raw_spec in args.checker_coordinate_cell.split(","):
        if not raw_spec:
            continue
        parts = raw_spec.split(":")
        if len(parts) != 3:
            raise SystemExit(
                "--checker-coordinate-cell entries must be modulus:k:l"
            )
        modulus, kr, lr = map(int, parts)
        if modulus < 1:
            raise SystemExit("checker coordinate moduli must be positive")
        checker_coordinate_cell.append(
            (modulus, kr % modulus, lr % modulus)
        )
    checker_coordinate_cell = tuple(checker_coordinate_cell)
    if checker_coordinate_cell and args.random_checks:
        raise SystemExit(
            "--random-checks must be 0 with --checker-coordinate-cell"
        )
    fixed_primes = tuple(
        int(value) for value in args.fixed_primes.split(",") if value
    )
    if args.algebraic_primes is not None:
        algebraic_primes = tuple(
            int(value)
            for value in args.algebraic_primes.split(",")
            if value
        )
    elif args.derived_targets:
        algebraic_primes = tuple(
            int(value)
            for value in source_payload.get("algebraic_primes", ())
        )
    else:
        algebraic_primes = tuple(
            prime for prime in exact_uncovered.factor(args.power) if prime % 2
        )
    sophie_germain = (
        bool(source_payload.get("sophie_germain", False))
        if args.derived_targets
        else args.power % 4 == 0
    )
    if not all(prime in by_prime for prime in normalizers):
        raise SystemExit("normalization prime missing from selected candidates")
    if not all(prime in by_prime for prime in diversity_primes):
        raise SystemExit("diversity prime missing from selected candidates")
    if not all(prime in by_prime for prime in fixed_primes):
        raise SystemExit("fixed prime missing from selected candidates")
    zero_fixed = {by_prime[prime] for prime in normalizers}
    saved_fixed = {by_prime[prime] for prime in fixed_primes}
    saved_fixed.update(range(args.fixed_prefix))
    fixed = zero_fixed | saved_fixed
    mutable = [index for index in range(len(candidates)) if index not in fixed]

    phase_source = (
        args.phase_file if args.phase_file.exists() else args.initial_phase_file
    )
    saved_phases = (
        json.loads(phase_source.read_text())
        if phase_source and phase_source.exists()
        else {}
    )
    if freeze_target_period > 1:
        missing_saved = [
            p for _h, p, _a, _b, _ord2, _ord3 in candidates
            if str(p) not in saved_phases
        ]
        if missing_saved:
            raise RuntimeError(
                "component freezing requires a saved phase for every "
                f"candidate; missing {missing_saved[:5]}"
            )
    valid_moduli = []
    valid_residues = []
    initial_values = []
    for h, p, _a, _b, _ord2, _ord3 in candidates:
        if args.derived_targets:
            residue, modulus = derived_target_by_prime[p]
        else:
            d = math.gcd(args.power, p - 1)
            step = (p - 1) // h
            divisor = math.gcd(step, d)
            half = (p - 1) // 2
            if half % divisor:
                raise RuntimeError(f"prime {p} has no power-compatible target")
            modulus = d // divisor
            residue = (
                0
                if modulus == 1
                else (half // divisor)
                * pow(step // divisor, -1, modulus)
                % modulus
            )
        origin_modulus = math.gcd(h, args.through_origin_mod)
        merged = merge_congruences(residue, modulus, 0, origin_modulus)
        if merged is None:
            raise RuntimeError(
                f"prime {p} has incompatible power and origin target constraints"
            )
        residue, modulus = merged
        frozen_modulus = math.gcd(h, freeze_target_period)
        if frozen_modulus > 1:
            frozen_residue = int(saved_phases[str(p)]) % frozen_modulus
            merged = merge_congruences(
                residue,
                modulus,
                frozen_residue,
                frozen_modulus,
            )
            if merged is None:
                raise RuntimeError(
                    f"saved frozen components for {p} are incompatible"
                )
            residue, modulus = merged
        if by_prime[p] in saved_fixed and str(p) not in saved_phases:
            raise RuntimeError(f"fixed prime {p} has no saved phase")
        value = int(saved_phases.get(str(p), residue)) % h
        if value % modulus != residue:
            if by_prime[p] in saved_fixed:
                raise RuntimeError(f"saved fixed phase for {p} is incompatible")
            value = residue
        valid_moduli.append(modulus)
        valid_residues.append(residue)
        initial_values.append(value)
    assignment = np.asarray(initial_values, dtype=np.uint32)
    valid_moduli = np.asarray(valid_moduli, dtype=np.uint32)
    valid_residues = np.asarray(valid_residues, dtype=np.uint32)
    candidate_moduli = np.asarray([item[0] for item in candidates], dtype=np.uint32)
    for index in zero_fixed:
        if valid_residues[index] != 0:
            raise RuntimeError("normalization target zero is not power-compatible")
        assignment[index] = 0

    points = []
    seen = set()

    def fixed_covered(point):
        k, l = point
        return any(
            (
                candidates[index][2] * (k % candidates[index][0])
                + candidates[index][3] * (l % candidates[index][0])
                - int(assignment[index])
            )
            % candidates[index][0]
            == 0
            for index in fixed
        )

    if args.points_file.exists():
        for raw_k, raw_l in json.loads(args.points_file.read_text()):
            point = (int(raw_k), int(raw_l))
            if any(
                point[0] % modulus != kr or point[1] % modulus != lr
                for modulus, kr, lr in checker_coordinate_cell
            ):
                continue
            if any(
                point[0] % prime == 0 and point[1] % prime == 0
                for prime in algebraic_primes
            ):
                continue
            if sophie_germain and point[0] % 4 == 2 and point[1] % 4 == 0:
                continue
            if args.drop_fixed_covered_points and fixed_covered(point):
                continue
            if point not in seen:
                seen.add(point)
                points.append(point)
    rng = random.Random(
        args.seed if args.seed is not None else 203_5000 + len(points)
    )
    if args.randomize_mutable_initial:
        for index in mutable:
            modulus = int(valid_moduli[index])
            residue = int(valid_residues[index])
            choices = (
                int(candidate_moduli[index]) - 1 - residue
            ) // modulus + 1
            assignment[index] = residue + modulus * rng.randrange(choices)
    print(
        f"candidates={len(candidates)}/{raw_count} points={len(points)} "
        f"density={sum(1 / item[0] for item in candidates):.12f} "
        f"randomized_initial={args.randomize_mutable_initial} "
        f"random_coordinate_bits={args.random_coordinate_bits}",
        flush=True,
    )

    targets, build_seconds = build_targets(points, candidates, np)
    for round_no in range(1, args.rounds + 1):
        repair_started = time.monotonic()
        repair_function = {
            "point": repair,
            "coordinate": repair_coordinate,
            "component": repair_component,
        }[args.repair_mode]
        solved, steps, before, after = repair_function(
            targets,
            assignment,
            mutable,
            rng,
            np,
            args.max_steps,
            args.sample_size,
            valid_moduli,
            valid_residues,
            candidate_moduli,
            args.required_coverage,
        )
        repair_seconds = time.monotonic() - repair_started
        print(
            f"round={round_no} points={len(points)} matrix_s={build_seconds:.3f} "
            f"repair={solved} steps={steps} misses={before}->{after} "
            f"repair_s={repair_seconds:.3f}",
            flush=True,
        )
        if not solved:
            args.phase_file.write_text(
                json.dumps(
                    {str(item[1]): int(assignment[index]) for index, item in enumerate(candidates)}
                )
                + "\n"
            )
            return 3
        args.phase_file.write_text(
            json.dumps(
                {str(item[1]): int(assignment[index]) for index, item in enumerate(candidates)}
            )
            + "\n"
        )
        if args.repair_only:
            print(
                f"round={round_no} repair_only=True checker=skipped",
                flush=True,
            )
            return 0

        random_miss_points, tested = random_misses(
            candidates,
            assignment,
            args.random_coordinate_bits,
            rng,
            args.random_checks,
            args.random_batch,
            algebraic_primes,
            sophie_germain,
        )
        misses = random_miss_points
        exact_misses = []
        engine = "random"
        meta = {"sat": True, "engine": engine}
        round_diversity_coordinate_moduli = diversity_coordinate_moduli
        if diversity_coordinate_schedule:
            stage_index = (
                (round_no - 1) // args.diversity_stage_rounds
            )
            if args.diversity_schedule_cycle:
                stage_index %= len(diversity_coordinate_schedule)
            else:
                stage_index = min(
                    stage_index,
                    len(diversity_coordinate_schedule) - 1,
                )
            round_diversity_coordinate_moduli = (
                diversity_coordinate_schedule[stage_index]
            )
        if len(random_miss_points) <= args.exact_when_random_misses_at_most:
            rows = [
                exact_greedy.as_row(item, int(assignment[index]))
                for index, item in enumerate(candidates)
            ]
            if not args.retain_target_matrix_during_exact:
                del targets
                gc.collect()
            exact_misses, exact_meta = exact_uncovered.find_uncovered(
                rows,
                max_component=args.max_component,
                limit=args.exact_batch,
                solver_name=args.checker_solver,
                diversity_primes=(
                    diversity_primes
                    if diversity_primes
                    or round_diversity_coordinate_moduli
                    else normalizers
                ),
                diversity_coordinate_moduli=(
                    round_diversity_coordinate_moduli
                ),
                diversity_quota=args.diversity_quota,
                algebraic_primes=algebraic_primes,
                sophie_germain=sophie_germain,
                fixed_coordinate_residues=checker_coordinate_cell,
            )
            if not exact_misses and random_miss_points:
                raise AssertionError(
                    "exact checker proved coverage despite random misses"
                )
            misses = list(random_miss_points) + list(exact_misses)
            engine = "exact" if not random_miss_points else "random+exact"
            meta = {
                "sat": bool(exact_misses),
                "engine": engine,
                "random_misses": len(random_miss_points),
                "exact": exact_meta,
            }
            if not exact_misses:
                payload = {
                    "candidate_pool": str(args.candidate_pool),
                    "candidate_count": len(candidates),
                    "choices": rows,
                    "checker": meta,
                    "synthesis": "local_phase_cegis",
                    "power": args.power,
                    "derived_targets": args.derived_targets,
                    "algebraic_primes": algebraic_primes,
                    "sophie_germain": sophie_germain,
                    "checker_coordinate_cell": checker_coordinate_cell,
                }
                args.output.write_text(json.dumps(payload, indent=2) + "\n")
                result_label = (
                    "CELL_COVER" if checker_coordinate_cell else "COVER"
                )
                print(
                    f"{result_label} rows={len(rows)} output={args.output}",
                    flush=True,
                )
                return 0
        if args.checker_checkpoint_file:
            coordinate_fingerprints = []
            if round_diversity_coordinate_moduli:
                coordinate_fingerprints = sorted(
                    {
                        tuple(
                            coordinate
                            for modulus in round_diversity_coordinate_moduli
                            for coordinate in (
                                int(k) % modulus,
                                int(l) % modulus,
                            )
                        )
                        for k, l in exact_misses
                    }
                )
            args.checker_checkpoint_file.write_text(
                json.dumps(
                    {
                        "round": round_no,
                        "engine": engine,
                        "phase_file": str(args.phase_file),
                        "random_misses": random_miss_points,
                        "exact_misses": exact_misses,
                        "diversity_coordinate_moduli": list(
                            round_diversity_coordinate_moduli
                        ),
                        "exact_coordinate_fingerprints": (
                            coordinate_fingerprints
                        ),
                        "checker": meta,
                    },
                    indent=2,
                )
                + "\n"
            )
        lesson_points = (
            tested if engine == "random" and args.retain_random_tested else misses
        )
        new_points = []
        for point in lesson_points:
            point = (int(point[0]), int(point[1]))
            if args.drop_fixed_covered_points and fixed_covered(point):
                continue
            if point not in seen:
                seen.add(point)
                points.append(point)
                new_points.append(point)
        args.points_file.write_text(json.dumps(points) + "\n")
        print(
            f"round={round_no} checker={engine} "
            f"misses={len(misses)}/"
            f"{len(tested) if engine == 'random' else 'exact'} "
            f"random={len(random_miss_points)} "
            f"exact={len(exact_misses) if engine != 'random' else 0} "
            f"diversity={','.join(map(str, round_diversity_coordinate_moduli)) or '-'} "
            f"new={len(new_points)} total={len(points)}",
            flush=True,
        )
        if not new_points:
            raise RuntimeError("checker returned no new point")
        if engine != "random":
            if args.retain_target_matrix_during_exact:
                new_targets, build_seconds = build_targets(
                    new_points, candidates, np
                )
                targets = np.asfortranarray(
                    np.vstack((targets, new_targets))
                )
            else:
                targets, build_seconds = build_targets(points, candidates, np)
        else:
            new_targets, build_seconds = build_targets(new_points, candidates, np)
            targets = np.vstack((targets, new_targets))
    return 4


if __name__ == "__main__":
    raise SystemExit(main())
