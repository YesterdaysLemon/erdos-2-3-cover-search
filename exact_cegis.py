#!/usr/bin/env python3
"""Exact CEGIS synthesis of a periodic affine-line cover for Erdős #203.

The master SAT problem chooses exactly one fibre for every eligible prime.
The checker in exact_uncovered.py either proves that assignment covers the
entire CRT torus or returns genuine missed residue classes, which are added to
the master.  SAT with checker-UNSAT is a construction; master-UNSAT rules out
the declared finite prime pool.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
import time
from pathlib import Path

import exact_uncovered
import exact_uncovered_z3
import cegis_cover
import local_cover
import search_cover
from power_anchor_capacity_lp import power_target_congruence

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument("--prime-limit", type=int, default=5_000_000)
    parser.add_argument(
        "--candidate-pool",
        type=Path,
        help="reuse signatures from a JSON file containing a choices array",
    )
    parser.add_argument(
        "--derived-pool",
        action="store_true",
        help="use h,a,b from a rigorously transformed affine-line pool",
    )
    parser.add_argument(
        "--derived-targets",
        action="store_true",
        help=(
            "read target_residue and target_modulus from each derived-pool "
            "row instead of recomputing a perfect-power condition"
        ),
    )
    parser.add_argument(
        "--normalize-primes",
        help="comma-separated jointly-surjective fibres to fix at target zero",
    )
    parser.add_argument(
        "--diversity-primes",
        default="",
        help=(
            "comma-separated checker fingerprints to block between exact "
            "witnesses without fixing their targets"
        ),
    )
    parser.add_argument(
        "--diversity-coordinate-moduli",
        default="",
        help=(
            "comma-separated prime-power coordinate moduli whose (k,l) "
            "residue cells are diversified between exact witnesses"
        ),
    )
    parser.add_argument("--diversity-quota", type=int, default=1)
    parser.add_argument("--candidate-count", type=int)
    parser.add_argument("--initial-witness", type=Path)
    parser.add_argument("--fix-initial-witness", action="store_true")
    parser.add_argument(
        "--complete-period-primes",
        action="store_true",
        help="factor gcd(2^L-1,3^L-1) and use every eligible prime",
    )
    parser.add_argument("--seed-points", type=int, default=100)
    parser.add_argument("--batch", type=int, default=100)
    parser.add_argument("--rounds", type=int, default=1000)
    parser.add_argument("--max-component", type=int, default=64)
    parser.add_argument("--checker", choices=("pysat", "z3"), default="pysat")
    parser.add_argument("--master-solver", default="cadical195")
    parser.add_argument(
        "--target-encoding", choices=("binary", "onehot"), default="binary"
    )
    parser.add_argument(
        "--random-checks",
        type=int,
        default=0,
        help="try this many genuine random points before invoking the exact checker",
    )
    parser.add_argument("--random-batch", type=int, default=1)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--progress-file",
        type=Path,
        help="overwrite this JSON checkpoint after every completed round",
    )
    parser.add_argument("--points-file", type=Path)
    parser.add_argument(
        "--phase-file",
        type=Path,
        help="persist the last satisfying per-prime target assignment as restart phases",
    )
    parser.add_argument(
        "--initial-phase-file",
        type=Path,
        help="read restart phases here when --phase-file does not yet exist",
    )
    parser.add_argument("--checkpoint-every", type=int, default=100)
    parser.add_argument(
        "--power",
        type=int,
        default=1,
        help="require m=M^power and use the exact X^q+1 algebraic cover",
    )
    parser.add_argument(
        "--algebraic-primes",
        help=(
            "comma-separated exact origin sublattices for the checker; with "
            "--derived-targets, defaults to candidate-pool metadata"
        ),
    )
    args = parser.parse_args()

    if args.power < 1:
        raise SystemExit("--power must be a positive integer")
    if args.derived_targets and (not args.derived_pool or not args.candidate_pool):
        raise SystemExit(
            "--derived-targets requires --derived-pool and --candidate-pool"
        )
    algebraic_primes = (
        tuple(
            int(value)
            for value in args.algebraic_primes.split(",")
            if value
        )
        if args.algebraic_primes is not None
        else tuple(
            prime for prime in exact_uncovered.factor(args.power) if prime % 2
        )
    )
    sophie_germain = args.power % 4 == 0 and not args.derived_targets

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    common = None
    factorization = None
    initial_targets = None
    incompatible_period_candidates = 0
    derived_target_by_prime = {}
    if args.candidate_pool:
        source = json.loads(args.candidate_pool.read_text())
        if args.derived_targets:
            if args.algebraic_primes is None:
                algebraic_primes = tuple(
                    int(value)
                    for value in source.get("algebraic_primes", ())
                )
            sophie_germain = bool(source.get("sophie_germain", False))
            for row in source["choices"]:
                p = int(row["p"])
                if p in derived_target_by_prime:
                    raise RuntimeError(f"duplicate derived target prime {p}")
                modulus = int(row["target_modulus"])
                residue = int(row["target_residue"]) % modulus
                if int(row["h"]) % modulus:
                    raise RuntimeError(
                        f"derived target modulus does not divide h for prime {p}"
                    )
                derived_target_by_prime[p] = (residue, modulus)
        candidates = []
        initial_targets = {}
        for row in source["choices"]:
            p = int(row["p"])
            if args.derived_pool:
                h = int(row["h"])
                a = int(row["a"]) % h
                b = int(row["b"]) % h
                if math.gcd(a, b, h) != 1:
                    raise RuntimeError(f"non-surjective derived row for prime {p}")
                item = (
                    h,
                    a,
                    b,
                    int(row.get("ord2", h)),
                    int(row.get("ord3", h)),
                )
            else:
                ord2 = search_cover.multiplicative_order(2, p)
                ord3 = search_cover.multiplicative_order(3, p)
                item = cegis_cover.signature(p, ord2, ord3)
            if args.period and args.period % item[0]:
                incompatible_period_candidates += 1
                continue
            candidates.append((item[0], p, item[1], item[2], item[3], item[4]))
            old_signature = (
                int(row["h"]),
                int(row["a"]),
                int(row["b"]),
                int(row["ord2"]),
                int(row["ord3"]),
            )
            initial_targets[p] = int(row["c"]) if old_signature == item else 0
        candidates.sort()
        if args.candidate_count is not None:
            candidates = candidates[: args.candidate_count]
    elif args.complete_period_primes:
        candidates, common, factorization = local_cover.get_complete_period_candidates(
            args.period
        )
    else:
        candidates = local_cover.get_candidates(args.period, args.prime_limit, None)
    power_targets = {}
    compatible_candidates = []
    incompatible_power_candidates = 0
    for item in candidates:
        h, p = item[0], item[1]
        if args.derived_targets:
            if p not in derived_target_by_prime:
                incompatible_power_candidates += 1
                continue
            residue, modulus = derived_target_by_prime[p]
        else:
            try:
                residue, modulus = power_target_congruence(h, p, args.power)
            except RuntimeError:
                incompatible_power_candidates += 1
                continue
        power_targets[p] = (residue, modulus)
        compatible_candidates.append(item)
    candidates = compatible_candidates
    if not candidates:
        raise RuntimeError("no power-compatible candidates")
    if args.period == 0:
        args.period = math.lcm(*(item[0] for item in candidates))
    by_prime_candidate = {item[1]: item for item in candidates}
    diversity_primes = tuple(
        int(value)
        for value in args.diversity_primes.split(",")
        if value
    )
    diversity_coordinate_moduli = tuple(
        int(value)
        for value in args.diversity_coordinate_moduli.split(",")
        if value
    )
    if args.diversity_quota < 1:
        raise SystemExit("--diversity-quota must be positive")
    if not all(p in by_prime_candidate for p in diversity_primes):
        raise RuntimeError("diversity prime is missing from the pool")
    if args.normalize_primes:
        normalize_primes = tuple(
            int(value) for value in args.normalize_primes.split(",") if value
        )
    elif args.derived_pool:
        normalize_primes = ()
    else:
        normalize_primes = (5, 7)
    if normalize_primes:
        if not all(p in by_prime_candidate for p in normalize_primes):
            raise RuntimeError("normalization prime is missing from the pool")
        if any(power_targets[p][0] != 0 for p in normalize_primes):
            raise RuntimeError(
                "a normalization prime cannot have power-compatible target zero"
            )
        normalizers = [
            (
                by_prime_candidate[p][0],
                by_prime_candidate[p][2],
                by_prime_candidate[p][3],
            )
            for p in normalize_primes
        ]
        normalizer_period = math.lcm(*(row[0] for row in normalizers))
        image_size = len(
            {
                tuple((a * k + b * l) % h for h, a, b in normalizers)
                for k in range(normalizer_period)
                for l in range(normalizer_period)
            }
        )
        if image_size != math.prod(row[0] for row in normalizers):
            raise RuntimeError("declared normalization maps are not jointly surjective")
    fixed_targets = {}
    if args.initial_witness:
        witness = json.loads(args.initial_witness.read_text())
        for row in witness["rows"]:
            fixed_targets[int(row["p"])] = int(row["line"][2])
        if initial_targets is None:
            initial_targets = {}
        initial_targets.update(
            {p: c for p, c in fixed_targets.items() if any(item[1] == p for item in candidates)}
        )
    phase_source = (
        args.phase_file
        if args.phase_file and args.phase_file.exists()
        else args.initial_phase_file
    )
    if phase_source and phase_source.exists():
        saved_phases = json.loads(phase_source.read_text())
        if initial_targets is None:
            initial_targets = {}
        initial_targets.update(
            {
                int(p): int(c)
                for p, c in saved_phases.items()
                if int(p) in by_prime_candidate
            }
        )
    period_label = (
        str(args.period)
        if args.period.bit_length() <= 1024
        else f"<{args.period.bit_length()} bits>"
    )
    print(
        f"period={period_label} candidates={len(candidates)} "
        f"density={sum(1/item[0] for item in candidates):.12f} "
        f"complete={args.complete_period_primes} "
        f"outside_period={incompatible_period_candidates}"
    )
    pool = IDPool()
    solver = Solver(name=args.master_solver)
    target_bits = {}
    target_counts = {}
    target_values = {}
    for h, p, a, b, ord2, ord3 in candidates:
        residue, modulus = power_targets[p]
        target_count = (h - 1 - residue) // modulus + 1
        target_counts[p] = target_count
        if args.target_encoding == "binary":
            width = max(1, (target_count - 1).bit_length())
            target_bits[p] = [
                pool.id(f"target_index_{p}_bit_{bit}") for bit in range(width)
            ]
        else:
            values = [
                pool.id(f"target_index_{p}_{index}")
                for index in range(target_count)
            ]
            target_values[p] = values
            solver.append_formula(
                CardEnc.equals(
                    values, 1, vpool=pool, encoding=EncType.seqcounter
                ).clauses
            )
    # Both encodings store the index j in c = residue + modulus*j, so every
    # represented value is power-compatible by construction.
    # A jointly-surjective tuple of line maps can have all its targets
    # normalized to zero by translating the exponent lattice.  For the
    # original power mode the default p=5,p=7 normalization uses translations
    # by power*(s,t), which preserve qZ^2 and D-th-power residue classes.
    for anchor_prime in normalize_primes:
        if anchor_prime in target_values:
            solver.add_clause([target_values[anchor_prime][0]])
        elif anchor_prime in target_bits:
            for bit in target_bits[anchor_prime]:
                solver.add_clause([-bit])
    if args.fix_initial_witness:
        for p, target in fixed_targets.items():
            if p in by_prime_candidate:
                h = by_prime_candidate[p][0]
                residue, modulus = power_targets[p]
                if not 0 <= target < h or target % modulus != residue:
                    raise RuntimeError(
                        f"fixed target {target} for prime {p} is not power-compatible"
                    )
            if p in target_values:
                residue, modulus = power_targets[p]
                solver.add_clause(
                    [target_values[p][(target - residue) // modulus]]
                )
            elif p in target_bits:
                residue, modulus = power_targets[p]
                target = (target - residue) // modulus
                for bit_no, bit in enumerate(target_bits[p]):
                    solver.add_clause([bit if (target >> bit_no) & 1 else -bit])
    if initial_targets is not None:
        by_prime = {item[1]: item for item in candidates}
        shift = (0, 0)
        if not args.derived_pool and all(
            p in by_prime and p in initial_targets for p in (5, 7)
        ):
            shift = next(
                (args.power * k, args.power * l)
                for k in range(12)
                for l in range(12)
                if all(
                    (
                        by_prime[p][2] * args.power * k
                        + by_prime[p][3] * args.power * l
                        - initial_targets[p]
                    )
                    % by_prime[p][0]
                    == 0
                    for p in (5, 7)
                )
            )
        phases = []
        for h, p, a, b, ord2, ord3 in candidates:
            suggested = initial_targets.get(p, power_targets[p][0])
            target = (suggested - a * shift[0] - b * shift[1]) % h
            residue, modulus = power_targets[p]
            if target % modulus != residue:
                target = residue
            if args.target_encoding == "binary":
                target = (target - residue) // modulus
                for bit_no, bit in enumerate(target_bits[p]):
                    phases.append(bit if (target >> bit_no) & 1 else -bit)
            else:
                target_index = (target - residue) // modulus
                phases.extend(
                    variable if index == target_index else -variable
                    for index, variable in enumerate(target_values[p])
                )
        try:
            solver.set_phases(phases)
        except NotImplementedError:
            # Some portfolio backends (notably Kissat through PySAT) do not
            # expose phase hints.  Hints affect search order, not the formula.
            pass

    points: set[tuple[int, int]] = set()
    equality_indicators: dict[tuple[int, int], int] = {}

    def indicator_for(p: int, target: int) -> int | None:
        residue, modulus = power_targets[p]
        if target % modulus != residue:
            return None
        target_index = (target - residue) // modulus
        key = (p, target_index)
        if key in equality_indicators:
            return equality_indicators[key]
        indicator = pool.id(f"target_index_{p}_equals_{target_index}")
        equality_indicators[key] = indicator
        matching_literals = []
        for bit_no, bit in enumerate(target_bits[p]):
            literal = bit if (target_index >> bit_no) & 1 else -bit
            matching_literals.append(literal)
            solver.add_clause([-indicator, literal])
        # Complete the equality gadget in both directions.  The former
        # one-way encoding was satisfiability-correct because point clauses
        # could choose a matching indicator, but it left the solver to guess
        # those indicators even after every target bit was fixed.
        solver.add_clause(
            [indicator, *(-literal for literal in matching_literals)]
        )
        return indicator

    def add_point(k: int, l: int) -> None:
        k %= args.period
        l %= args.period
        if any(k % prime == 0 and l % prime == 0 for prime in algebraic_primes):
            return
        if sophie_germain and k % 4 == 2 and l % 4 == 0:
            return
        if (k, l) in points:
            return
        points.add((k, l))
        covering = []
        for h, p, a, b, ord2, ord3 in candidates:
            target = (a * k + b * l) % h
            if args.target_encoding == "onehot":
                residue, modulus = power_targets[p]
                if target % modulus == residue:
                    covering.append(
                        target_values[p][(target - residue) // modulus]
                    )
                continue
            # Reuse the same equality gadget whenever several sampled points
            # ask a prime for the same target.  This is exact and saves large
            # amounts of memory once the counterexample set grows.
            indicator = indicator_for(p, target)
            if indicator is not None:
                covering.append(indicator)
        if not covering:
            raise RuntimeError(f"no compatible fibre can cover point {(k, l)}")
        solver.add_clause(covering)

    rng = random.Random(203)
    if args.points_file and args.points_file.exists():
        for k, l in json.loads(args.points_file.read_text()):
            add_point(int(k), int(l))
    for _ in range(args.seed_points):
        add_point(rng.randrange(args.period), rng.randrange(args.period))
    if args.points_file:
        args.points_file.write_text(json.dumps(sorted(points)) + "\n")

    for round_no in range(1, args.rounds + 1):
        master_started = time.monotonic()
        if not solver.solve():
            if args.points_file:
                args.points_file.write_text(json.dumps(sorted(points)) + "\n")
            payload = {
                "period": args.period,
                "candidate_pool": (
                    str(args.candidate_pool) if args.candidate_pool else None
                ),
                "derived_pool": args.derived_pool,
                "derived_targets": args.derived_targets,
                "candidate_count": len(candidates),
                "master_points": len(points),
                "normalize_primes": normalize_primes,
                "fixed_targets": (
                    {
                        str(prime): target
                        for prime, target in sorted(fixed_targets.items())
                        if prime in by_prime_candidate
                    }
                    if args.fix_initial_witness
                    else {}
                ),
                "target_encoding": args.target_encoding,
                "master_solver": args.master_solver,
                "master_unsat": True,
                "scope": (
                    "the declared finite candidate pool subject to the "
                    "recorded target restrictions"
                ),
            }
            if args.output:
                args.output.write_text(json.dumps(payload, indent=2) + "\n")
            if args.progress_file:
                args.progress_file.write_text(
                    json.dumps(payload, indent=2) + "\n"
                )
            print(f"MASTER_UNSAT round={round_no} points={len(points)}")
            return 2
        master_seconds = time.monotonic() - master_started
        model = {literal for literal in solver.get_model() if literal > 0}
        rows = []
        for h, p, a, b, ord2, ord3 in candidates:
            if args.target_encoding == "onehot":
                target_index = next(
                    index
                    for index, variable in enumerate(target_values[p])
                    if variable in model
                )
                residue, modulus = power_targets[p]
                c = residue + modulus * target_index
            else:
                raw = sum(
                    (1 << bit_no) if bit in model else 0
                    for bit_no, bit in enumerate(target_bits[p])
                )
                residue, modulus = power_targets[p]
                # An out-of-range target-index code means this prime is unused
                # by the current master assignment. Fill it with the first
                # valid fibre for checking; checker-UNSAT is still a real cover.
                c = (
                    residue + modulus * raw
                    if raw < target_counts[p]
                    else residue
                )
            if args.power != 1:
                d = math.gcd(args.power, p - 1)
                exponent = (p - 1) // 2 - ((p - 1) // h) * c
                assert exponent % d == 0, (p, h, c, d)
            rows.append(
                {
                    "h": h,
                    "p": p,
                    "a": a,
                    "b": b,
                    "ord2": ord2,
                    "ord3": ord3,
                    "c": c,
                }
            )
        if args.phase_file:
            args.phase_file.write_text(
                json.dumps({str(row["p"]): int(row["c"]) for row in rows}) + "\n"
            )
        checker_started = time.monotonic()
        missed = []
        for _ in range(args.random_checks):
            point = (rng.randrange(args.period), rng.randrange(args.period))
            if any(
                point[0] % prime == 0 and point[1] % prime == 0
                for prime in algebraic_primes
            ):
                continue
            if (
                sophie_germain
                and point[0] % 4 == 2
                and point[1] % 4 == 0
            ):
                continue
            if all(
                (int(row["a"]) * point[0] + int(row["b"]) * point[1] - int(row["c"]))
                % int(row["h"])
                for row in rows
            ):
                missed.append(point)
                if len(missed) >= args.random_batch:
                    break
        if missed:
            meta = {"engine": "random", "sat": True, "rows": len(rows)}
        else:
            checker = (
                exact_uncovered_z3.find_uncovered
                if args.checker == "z3"
                else exact_uncovered.find_uncovered
            )
            checker_kwargs = {
                "max_component": args.max_component,
                "limit": args.batch,
                "algebraic_primes": algebraic_primes,
                "sophie_germain": sophie_germain,
            }
            if checker is exact_uncovered.find_uncovered:
                checker_kwargs["diversity_primes"] = (
                    diversity_primes if diversity_primes else normalize_primes
                )
                checker_kwargs["diversity_coordinate_moduli"] = (
                    diversity_coordinate_moduli
                )
                checker_kwargs["diversity_quota"] = args.diversity_quota
            missed, meta = checker(rows, **checker_kwargs)
        checker_seconds = time.monotonic() - checker_started
        print(
            f"round={round_no} master_points={len(points)} "
            f"checker={'SAT' if missed else 'UNSAT'} new={len(missed)} "
            f"master_s={master_seconds:.3f} checker_s={checker_seconds:.3f}"
        )
        if args.progress_file:
            args.progress_file.write_text(
                json.dumps(
                    {
                        "round": round_no,
                        "master_points": len(points),
                        "missed": missed,
                        "master_seconds": master_seconds,
                        "checker_seconds": checker_seconds,
                    },
                    indent=2,
                )
                + "\n"
            )
        if not missed:
            payload = {
                "period": args.period,
                "prime_limit": args.prime_limit,
                "candidate_pool": str(args.candidate_pool) if args.candidate_pool else None,
                "derived_pool": args.derived_pool,
                "complete_period_primes": args.complete_period_primes,
                "common_gcd": str(common) if common is not None else None,
                "factorization": factorization,
                "choices": rows,
                "power": args.power,
                "derived_targets": args.derived_targets,
                "algebraic_primes": algebraic_primes,
                "sophie_germain": sophie_germain,
                "incompatible_power_candidates": incompatible_power_candidates,
                "diversity_coordinate_moduli": diversity_coordinate_moduli,
                "checker": meta,
            }
            if args.output:
                args.output.write_text(json.dumps(payload, indent=2) + "\n")
            else:
                print(json.dumps(payload, indent=2))
            return 0
        before = len(points)
        for k, l in missed:
            add_point(k, l)
        if (
            args.points_file
            and args.checkpoint_every > 0
            and round_no % args.checkpoint_every == 0
        ):
            args.points_file.write_text(json.dumps(sorted(points)) + "\n")
        if len(points) == before:
            raise RuntimeError("checker returned no new counterexample")
    print(f"ROUND_LIMIT points={len(points)}")
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
