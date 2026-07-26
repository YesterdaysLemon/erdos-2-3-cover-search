#!/usr/bin/env python3
"""Coordinate descent for maximizing affine-line union on a fixed sample."""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
import time
from pathlib import Path

import exact_greedy
import exact_uncovered
from local_phase_cegis import build_targets, merge_congruences
from power_anchor_capacity_lp import power_target_congruence

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def audit(
    candidates,
    assignment,
    algebraic_primes,
    sophie_germain,
    seed,
    count,
):
    rng = random.Random(seed)
    missed = 0
    tested = 0
    while tested < count:
        k = rng.getrandbits(64)
        l = rng.getrandbits(64)
        if any(k % prime == 0 and l % prime == 0 for prime in algebraic_primes):
            continue
        if sophie_germain and k % 4 == 2 and l % 4 == 0:
            continue
        tested += 1
        if all(
            (a * (k % h) + b * (l % h) - int(assignment[index])) % h
            for index, (h, _p, a, b, _o2, _o3) in enumerate(candidates)
        ):
            missed += 1
    return missed


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
    parser.add_argument("--candidate-count", type=int)
    parser.add_argument("--period", type=int, default=0)
    parser.add_argument(
        "--max-component",
        type=int,
        default=0,
        help="if nonzero, retain rows whose prime-power components are at most this",
    )
    parser.add_argument("--power", type=int, default=1)
    parser.add_argument(
        "--normalize-primes",
        default="",
        help="comma-separated jointly normalized primes fixed at target zero",
    )
    parser.add_argument(
        "--fixed-primes",
        default="",
        help="comma-separated primes whose saved phases must not move",
    )
    parser.add_argument("--freeze-first", type=int, default=0)
    parser.add_argument("--phase-file", type=Path, required=True)
    parser.add_argument("--initial-phase-file", type=Path)
    parser.add_argument("--overlay-phase-file", type=Path, action="append", default=[])
    parser.add_argument(
        "--points-file",
        type=Path,
        help="ordered JSON point list; when present, do not generate a sample",
    )
    parser.add_argument(
        "--fill-random-sample",
        action="store_true",
        help=(
            "after loading --points-file, append admissible random points "
            "until --sample-size is reached"
        ),
    )
    parser.add_argument("--sample-size", type=int, default=20000)
    parser.add_argument("--audit-size", type=int, default=5000)
    parser.add_argument("--sweeps", type=int, default=100)
    parser.add_argument(
        "--stream-targets",
        action="store_true",
        help=(
            "recompute candidate columns during each sweep instead of "
            "retaining the full point-by-candidate matrix"
        ),
    )
    parser.add_argument("--zero-move-rate", type=float, default=0.0)
    parser.add_argument(
        "--mutable-limit",
        type=int,
        default=0,
        help=(
            "if positive, examine only this many randomly ordered mutable "
            "fibres per sweep; useful for bounded soft-margin passes"
        ),
    )
    parser.add_argument(
        "--coverage-rho",
        type=float,
        default=0.0,
        help=(
            "optimize sum(rho**coverage); zero recovers the exact hole-count "
            "objective, while a small positive value rewards coverage margin"
        ),
    )
    parser.add_argument(
        "--preserve-hole-count",
        action="store_true",
        help=(
            "with a soft objective, reject moves that increase the number "
            "of uncovered retained points on either sample split"
        ),
    )
    parser.add_argument(
        "--validation-fraction",
        type=float,
        default=0.0,
        help="require each move to be non-worsening on both fixed sample splits",
    )
    parser.add_argument(
        "--split-index",
        type=int,
        help="explicit first validation index, overriding --validation-fraction",
    )
    parser.add_argument("--seed", type=int, default=2032026)
    parser.add_argument(
        "--algebraic-primes",
        help=(
            "comma-separated exact origin sublattices; with "
            "--derived-targets, defaults to candidate-pool metadata"
        ),
    )
    parser.add_argument(
        "--through-origin-mod",
        type=int,
        default=1,
        help=(
            "restrict every fibre to target c=0 modulo gcd(h,Q), matching "
            "local_phase_cegis"
        ),
    )
    args = parser.parse_args()
    if not 0.0 <= args.coverage_rho < 1.0:
        raise SystemExit("--coverage-rho must lie in [0,1)")
    if args.mutable_limit < 0:
        raise SystemExit("--mutable-limit must be nonnegative")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import numpy as np  # type: ignore

    source_payload = json.loads(args.candidate_pool.read_text())
    if args.derived_targets and not args.derived_pool:
        raise SystemExit("--derived-targets requires --derived-pool")
    if args.fill_random_sample and not args.points_file:
        raise SystemExit("--fill-random-sample requires --points-file")
    if args.through_origin_mod < 1:
        raise SystemExit("--through-origin-mod must be positive")
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
    if args.period:
        candidates = [row for row in candidates if args.period % row[0] == 0]
    if args.max_component:
        candidates = [
            row
            for row in candidates
            if max(
                (
                    prime**exponent
                    for prime, exponent in exact_uncovered.factor(row[0]).items()
                ),
                default=1,
            )
            <= args.max_component
        ]
    if args.derived_targets:
        candidates = [
            row for row in candidates if row[1] in derived_target_by_prime
        ]
    else:
        candidates = [
            row
            for row in candidates
            if (
                ((row[1] - 1) // 2)
                % math.gcd(
                    (row[1] - 1) // row[0],
                    math.gcd(args.power, row[1] - 1),
                )
                == 0
            )
        ]
    if args.candidate_count is not None:
        candidates = candidates[: args.candidate_count]
    by_prime = {p: index for index, (_h, p, *_rest) in enumerate(candidates)}
    normalizers = tuple(
        int(value) for value in args.normalize_primes.split(",") if value
    )
    normalizer_indices = {by_prime[prime] for prime in normalizers}
    fixed = set(normalizer_indices)
    fixed_primes = tuple(
        int(value) for value in args.fixed_primes.split(",") if value
    )
    missing_fixed = [prime for prime in fixed_primes if prime not in by_prime]
    if missing_fixed:
        raise SystemExit(f"fixed primes are absent: {missing_fixed}")
    fixed.update(by_prime[prime] for prime in fixed_primes)
    phase_source = (
        args.phase_file
        if args.phase_file.exists()
        else args.initial_phase_file
    )
    saved = json.loads(phase_source.read_text()) if phase_source and phase_source.exists() else {}
    for overlay in args.overlay_phase_file:
        saved.update(json.loads(overlay.read_text()))
    valid_residues = []
    valid_moduli = []
    initial_targets = []
    for h, p, *_rest in candidates:
        if args.derived_targets:
            residue, modulus = derived_target_by_prime[p]
        else:
            residue, modulus = power_target_congruence(h, p, args.power)
        origin_modulus = math.gcd(h, args.through_origin_mod)
        merged = merge_congruences(residue, modulus, 0, origin_modulus)
        if merged is None:
            raise RuntimeError(
                f"prime {p} has incompatible target and origin constraints"
            )
        residue, modulus = merged
        valid_residues.append(residue)
        valid_moduli.append(modulus)
        target = int(saved.get(str(p), residue)) % h
        if target % modulus != residue:
            target = residue
        initial_targets.append(target)
    valid_residues = np.asarray(valid_residues, dtype=np.uint32)
    valid_moduli = np.asarray(valid_moduli, dtype=np.uint32)
    assignment = np.asarray(initial_targets, dtype=np.uint32)
    for index in normalizer_indices:
        assignment[index] = 0
    fixed.update(range(min(args.freeze_first, len(candidates))))

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
    rng = random.Random(args.seed)
    points = []
    seen_points = set()
    if args.points_file:
        for raw_k, raw_l in json.loads(args.points_file.read_text()):
            point = (int(raw_k), int(raw_l))
            if point in seen_points:
                continue
            if any(
                point[0] % prime == 0 and point[1] % prime == 0
                for prime in algebraic_primes
            ):
                continue
            if sophie_germain and point[0] % 4 == 2 and point[1] % 4 == 0:
                continue
            seen_points.add(point)
            points.append(point)
    if not args.points_file or args.fill_random_sample:
        while len(points) < args.sample_size:
            point = (rng.getrandbits(64), rng.getrandbits(64))
            if any(
                point[0] % prime == 0 and point[1] % prime == 0
                for prime in algebraic_primes
            ):
                continue
            if sophie_germain and point[0] % 4 == 2 and point[1] % 4 == 0:
                continue
            if point in seen_points:
                continue
            seen_points.add(point)
            points.append(point)
    if args.stream_targets:
        build_started = time.monotonic()
        uint64_flags = [
            all(0 <= coordinate < (1 << 64) for coordinate in point)
            for point in points
        ]
        uint64_indices = np.flatnonzero(uint64_flags)
        big_indices = np.flatnonzero(np.logical_not(uint64_flags))
        ks = np.fromiter(
            (points[index][0] for index in uint64_indices),
            dtype=np.uint64,
            count=len(uint64_indices),
        )
        ls = np.fromiter(
            (points[index][1] for index in uint64_indices),
            dtype=np.uint64,
            count=len(uint64_indices),
        )

        def stream_column(item):
            h, _p, a, b, _ord2, _ord3 = item
            column = np.empty(len(points), dtype=np.uint32)
            if len(uint64_indices):
                column[uint64_indices] = (
                    (a * (ks % h) + b * (ls % h)) % h
                ).astype(np.uint32)
            if len(big_indices):
                column[big_indices] = np.fromiter(
                    (
                        (
                            a * (points[index][0] % h)
                            + b * (points[index][1] % h)
                        )
                        % h
                        for index in big_indices
                    ),
                    dtype=np.uint32,
                    count=len(big_indices),
                )
            return column

        cover = np.zeros(len(points), dtype=np.int32)
        for index, (h, _p, a, b, _ord2, _ord3) in enumerate(candidates):
            values = stream_column(candidates[index])
            cover += (values == assignment[index]).astype(np.int32)
        targets = None
        build_seconds = time.monotonic() - build_started
    else:
        targets, build_seconds = build_targets(points, candidates, np)
        cover = np.count_nonzero(targets == assignment, axis=1).astype(np.int32)
    initial_histogram = np.bincount(cover)
    initial_singletons = (
        int(initial_histogram[1]) if len(initial_histogram) > 1 else 0
    )
    print(
        f"candidates={len(candidates)} sample={len(points)} "
        f"density={sum(1 / item[0] for item in candidates):.12f} "
        f"matrix_s={build_seconds:.3f} "
        f"initial_missed={int(np.count_nonzero(cover == 0))} "
        f"initial_singletons={initial_singletons}",
        flush=True,
    )

    mutable = [index for index in range(len(candidates)) if index not in fixed]
    split = (
        args.split_index
        if args.split_index is not None
        else int(len(points) * (1.0 - args.validation_fraction))
    )
    if not 0 <= split <= len(points):
        raise SystemExit("--split-index is outside the point list")
    use_split = 0 < split < len(points)
    for sweep in range(1, args.sweeps + 1):
        rng.shuffle(mutable)
        sweep_mutable = (
            mutable
            if not args.mutable_limit
            else mutable[: args.mutable_limit]
        )
        moves = 0
        zero_moves = 0
        total_gain = 0.0
        started = time.monotonic()
        for index in sweep_mutable:
            if targets is None:
                column = stream_column(candidates[index])
            else:
                column = targets[:, index]
            current = int(assignment[index])
            old_mask = column == current
            uncovered = cover == 0
            if args.coverage_rho == 0.0 and not np.any(uncovered):
                break
            h = candidates[index][0]
            if use_split:
                if args.coverage_rho == 0.0:
                    train_uncovered = uncovered[:split]
                    valid_uncovered = uncovered[split:]
                    train_loss = float(
                        np.count_nonzero(
                            old_mask[:split] & (cover[:split] == 1)
                        )
                    )
                    valid_loss = float(
                        np.count_nonzero(
                            old_mask[split:] & (cover[split:] == 1)
                        )
                    )
                    train_gains = np.bincount(
                        column[:split][train_uncovered], minlength=h
                    ).astype(np.float64)
                    valid_gains = np.bincount(
                        column[split:][valid_uncovered], minlength=h
                    ).astype(np.float64)
                else:
                    rho = args.coverage_rho
                    upward = (1.0 - rho) * np.power(
                        rho, cover, dtype=np.float64
                    )
                    old_indices = np.flatnonzero(old_mask)
                    downward = (1.0 - rho) * np.power(
                        rho,
                        cover[old_indices] - 1,
                        dtype=np.float64,
                    )
                    train_old = old_indices < split
                    train_loss = float(np.sum(downward[train_old]))
                    valid_loss = float(np.sum(downward[~train_old]))
                    train_gains = np.bincount(
                        column[:split],
                        weights=upward[:split],
                        minlength=h,
                    )
                    valid_gains = np.bincount(
                        column[split:],
                        weights=upward[split:],
                        minlength=h,
                    )
                train_delta = train_gains - train_loss
                valid_delta = valid_gains - valid_loss
                scores = train_delta + valid_delta
                invalid_moves = (
                    (train_delta < -1e-12) | (valid_delta < -1e-12)
                )
                if args.preserve_hole_count and args.coverage_rho != 0.0:
                    hard_train_loss = float(
                        np.count_nonzero(
                            old_mask[:split] & (cover[:split] == 1)
                        )
                    )
                    hard_valid_loss = float(
                        np.count_nonzero(
                            old_mask[split:] & (cover[split:] == 1)
                        )
                    )
                    hard_train_gains = np.bincount(
                        column[:split][uncovered[:split]], minlength=h
                    )
                    hard_valid_gains = np.bincount(
                        column[split:][uncovered[split:]], minlength=h
                    )
                    invalid_moves |= (
                        hard_train_gains < hard_train_loss
                    ) | (hard_valid_gains < hard_valid_loss)
                scores[invalid_moves] = -float("inf")
                scores[current] = -float("inf")
                valid_targets = np.arange(
                    int(valid_residues[index]),
                    h,
                    int(valid_moduli[index]),
                    dtype=np.int64,
                )
                valid_scores = scores[valid_targets]
                best_gain = float(valid_scores.max())
                best_targets = valid_targets[
                    np.flatnonzero(valid_scores == best_gain)
                ]
                delta = best_gain
            else:
                if args.coverage_rho == 0.0:
                    loss = float(
                        np.count_nonzero(old_mask & (cover == 1))
                    )
                    gains = np.bincount(
                        column[uncovered], minlength=h
                    ).astype(np.float64)
                else:
                    rho = args.coverage_rho
                    upward = (1.0 - rho) * np.power(
                        rho, cover, dtype=np.float64
                    )
                    loss = float(
                        np.sum(
                            (1.0 - rho)
                            * np.power(
                                rho,
                                cover[old_mask] - 1,
                                dtype=np.float64,
                            )
                        )
                    )
                    gains = np.bincount(
                        column, weights=upward, minlength=h
                    )
                scores = gains.astype(np.float64) - loss
                if args.preserve_hole_count and args.coverage_rho != 0.0:
                    hard_loss = float(
                        np.count_nonzero(old_mask & (cover == 1))
                    )
                    hard_gains = np.bincount(
                        column[uncovered], minlength=h
                    )
                    scores[hard_gains < hard_loss] = -float("inf")
                scores[current] = -float("inf")
                valid_targets = np.arange(
                    int(valid_residues[index]),
                    h,
                    int(valid_moduli[index]),
                    dtype=np.int64,
                )
                valid_scores = scores[valid_targets]
                delta = float(valid_scores.max())
                best_targets = valid_targets[
                    np.flatnonzero(valid_scores == delta)
                ]
            target = int(best_targets[rng.randrange(len(best_targets))])
            accept_zero = (
                abs(delta) <= 1e-12
                and target != current
                and rng.random() < args.zero_move_rate
            )
            if delta <= 1e-12 and not accept_zero:
                continue
            new_mask = column == target
            cover -= old_mask.astype(np.int32)
            cover += new_mask.astype(np.int32)
            assignment[index] = target
            moves += 1
            total_gain += delta
            zero_moves += int(accept_zero)

        args.phase_file.write_text(
            json.dumps(
                {
                    str(item[1]): int(assignment[index])
                    for index, item in enumerate(candidates)
                }
            )
            + "\n"
        )
        missed = int(np.count_nonzero(cover == 0))
        final_histogram = np.bincount(cover)
        singletons = (
            int(final_histogram[1]) if len(final_histogram) > 1 else 0
        )
        audit_missed = audit(
            candidates,
            assignment,
            algebraic_primes,
            sophie_germain,
            args.seed + 10_000 + sweep,
            args.audit_size,
        )
        print(
            f"sweep={sweep} missed={missed}/{len(points)} moves={moves} "
            f"zero={zero_moves} scanned={len(sweep_mutable)} "
            f"singletons={singletons} gain={total_gain:.6f} "
            f"s={time.monotonic()-started:.3f} "
            f"audit={audit_missed}/{args.audit_size}",
            flush=True,
        )
        if missed == 0 and args.coverage_rho == 0.0:
            return 0
        if moves == 0:
            return 2
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
