#!/usr/bin/env python3
"""Incremental exact greedy affine-line cover search.

Unlike exact_greedy.py, this keeps one uncovered-point SAT solver alive while
new line constraints and prime-power residue components are added.  Every SAT
model is a genuine missed exponent pair; solver UNSAT is independently rebuilt
and checked before any cover is written.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
from pathlib import Path

import exact_greedy
import exact_uncovered

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


class IncrementalChecker:
    def __init__(self, max_component: int, solver_name: str = "cadical195"):
        dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
        sys.path.insert(0, str(dep_path))
        from pysat.card import CardEnc, EncType  # type: ignore
        from pysat.formula import IDPool  # type: ignore
        from pysat.solvers import Solver  # type: ignore

        self.CardEnc = CardEnc
        self.EncType = EncType
        self.pool = IDPool()
        self.solver = Solver(name=solver_name)
        self.max_component = max_component
        self.components: dict[int, int] = {}
        self.klevels: dict[tuple[int, int], list[int]] = {}
        self.llevels: dict[tuple[int, int], list[int]] = {}
        self.rows: list[dict] = []

    def _exactly_one(self, variables: list[int]) -> None:
        self.solver.append_formula(
            self.CardEnc.equals(
                variables,
                1,
                vpool=self.pool,
                encoding=self.EncType.seqcounter,
            ).clauses
        )

    def _new_level(self, prime: int, modulus: int) -> tuple[list[int], list[int]]:
        ks = [self.pool.id(f"k_{prime}_{modulus}_{r}") for r in range(modulus)]
        ls = [self.pool.id(f"l_{prime}_{modulus}_{r}") for r in range(modulus)]
        self._exactly_one(ks)
        self._exactly_one(ls)
        self.klevels[prime, modulus] = ks
        self.llevels[prime, modulus] = ls
        return ks, ls

    def _link_refinement(
        self,
        coarse: list[int],
        coarse_modulus: int,
        fine: list[int],
        fine_modulus: int,
    ) -> None:
        if fine_modulus % coarse_modulus:
            raise AssertionError((coarse_modulus, fine_modulus))
        for residue, variable in enumerate(coarse):
            refinements = [
                fine[value]
                for value in range(residue, fine_modulus, coarse_modulus)
            ]
            self.solver.add_clause([-variable, *refinements])
        for value, variable in enumerate(fine):
            self.solver.add_clause([-variable, coarse[value % coarse_modulus]])

    def ensure_component(self, prime: int, modulus: int) -> None:
        if modulus > self.max_component:
            raise ValueError(
                f"prime-power component {modulus} exceeds guard {self.max_component}"
            )
        old_modulus = self.components.get(prime)
        if old_modulus is None:
            self._new_level(prime, modulus)
            self.components[prime] = modulus
            return
        if modulus <= old_modulus:
            if old_modulus % modulus:
                raise AssertionError((old_modulus, modulus))
            return
        old_k = self.klevels[prime, old_modulus]
        old_l = self.llevels[prime, old_modulus]
        new_k, new_l = self._new_level(prime, modulus)
        self._link_refinement(old_k, old_modulus, new_k, modulus)
        self._link_refinement(old_l, old_modulus, new_l, modulus)
        self.components[prime] = modulus

    def level(self, prime: int, modulus: int) -> tuple[list[int], list[int]]:
        key = (prime, modulus)
        if key in self.klevels:
            return self.klevels[key], self.llevels[key]
        fine_modulus = self.components[prime]
        fine_k = self.klevels[prime, fine_modulus]
        fine_l = self.llevels[prime, fine_modulus]
        coarse_k, coarse_l = self._new_level(prime, modulus)
        self._link_refinement(coarse_k, modulus, fine_k, fine_modulus)
        self._link_refinement(coarse_l, modulus, fine_l, fine_modulus)
        return coarse_k, coarse_l

    def add_row(self, row: dict) -> None:
        factors = exact_uncovered.factor(int(row["h"]))
        for prime, exponent in factors.items():
            self.ensure_component(prime, prime**exponent)
        index = len(self.rows)
        a, b, c = int(row["a"]), int(row["b"]), int(row["c"])
        failures = []
        for prime, exponent in factors.items():
            modulus = prime**exponent
            ks, ls = self.level(prime, modulus)
            failure = self.pool.id(f"failure_{index}_{prime}")
            failures.append(failure)
            matches = []
            if math.gcd(b, modulus) == 1:
                inverse = pow(b, -1, modulus)
                pairs = [
                    (kr, ((c - a * kr) * inverse) % modulus)
                    for kr in range(modulus)
                ]
            elif math.gcd(a, modulus) == 1:
                inverse = pow(a, -1, modulus)
                pairs = [
                    (((c - b * lr) * inverse) % modulus, lr)
                    for lr in range(modulus)
                ]
            else:
                raise AssertionError("derived line is not surjective")
            for pair_no, (kr, lr) in enumerate(pairs):
                match = self.pool.id(f"match_{index}_{prime}_{pair_no}")
                matches.append(match)
                self.solver.add_clause([-match, ks[kr]])
                self.solver.add_clause([-match, ls[lr]])
                self.solver.add_clause([-ks[kr], -ls[lr], match])
                self.solver.add_clause([-match, -failure])
            self.solver.add_clause([failure, *matches])
        self.solver.add_clause(failures)
        self.rows.append(row)

    def solve(self) -> tuple[int, int] | None:
        if not self.solver.solve():
            return None
        model = {literal for literal in self.solver.get_model() if literal > 0}
        kres = []
        lres = []
        for prime, modulus in self.components.items():
            ks = self.klevels[prime, modulus]
            ls = self.llevels[prime, modulus]
            kr = next(residue for residue, variable in enumerate(ks) if variable in model)
            lr = next(residue for residue, variable in enumerate(ls) if variable in model)
            kres.append((kr, modulus))
            lres.append((lr, modulus))
        k = exact_uncovered.crt(kres)
        l = exact_uncovered.crt(lres)
        if not all(
            (int(row["a"]) * k + int(row["b"]) * l - int(row["c"]))
            % int(row["h"])
            for row in self.rows
        ):
            raise AssertionError("incremental checker returned a covered point")
        return k, l

    def initialize_phases(self) -> None:
        seed = sum(
            (index + 1) * (int(row["c"]) + 1) * int(row["p"])
            for index, row in enumerate(self.rows)
        )
        rng = random.Random(seed)
        phases = []
        for prime, modulus in self.components.items():
            preferred_k = rng.randrange(modulus)
            preferred_l = rng.randrange(modulus)
            ks = self.klevels[prime, modulus]
            ls = self.llevels[prime, modulus]
            for residue in range(modulus):
                phases.append(ks[residue] if residue == preferred_k else -ks[residue])
                phases.append(ls[residue] if residue == preferred_l else -ls[residue])
        self.solver.set_phases(phases)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate_pool", type=Path)
    parser.add_argument("--derived-pool", action="store_true")
    parser.add_argument("--normalize-primes", required=True)
    parser.add_argument("--max-component", type=int, default=300000)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--checkpoint-every", type=int, default=25)
    parser.add_argument("--rounds", type=int, default=50000)
    parser.add_argument("--add-per-round", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.add_per_round < 1:
        raise SystemExit("--add-per-round must be positive")

    candidates = exact_greedy.load_candidates(
        args.candidate_pool, args.derived_pool
    )
    by_prime = {item[1]: item for item in candidates}
    normalizers = tuple(int(value) for value in args.normalize_primes.split(","))
    anchors = [(by_prime[p][0], by_prime[p][2], by_prime[p][3]) for p in normalizers]
    period = math.lcm(*(row[0] for row in anchors))
    image = {
        tuple((a * k + b * l) % h for h, a, b in anchors)
        for k in range(period)
        for l in range(period)
    }
    if len(image) != math.prod(row[0] for row in anchors):
        raise RuntimeError("normalizer maps are not jointly surjective")

    selected = []
    used = set()
    start_round = 0
    if args.checkpoint and args.checkpoint.exists():
        saved = json.loads(args.checkpoint.read_text())
        start_round = int(saved.get("round", 0))
        for old in saved["choices"]:
            p = int(old["p"])
            selected.append(exact_greedy.as_row(by_prime[p], int(old["c"])))
            used.add(p)
    else:
        for p in normalizers:
            selected.append(exact_greedy.as_row(by_prime[p], 0))
            used.add(p)

    checker = IncrementalChecker(args.max_component, args.solver)
    for row in selected:
        checker.add_row(row)
    checker.initialize_phases()
    print(
        f"loaded={len(selected)} candidates={len(candidates)} "
        f"components={len(checker.components)}",
        flush=True,
    )

    ordered_unused = (item for item in candidates if item[1] not in used)
    for round_no in range(start_round + 1, args.rounds + 1):
        witness = checker.solve()
        if witness is None:
            missed, meta = exact_uncovered.find_uncovered(
                selected,
                max_component=args.max_component,
                limit=1,
                solver_name=args.solver,
            )
            if missed:
                raise RuntimeError(
                    f"incremental UNSAT disagrees with rebuilt checker: {missed[0]}"
                )
            payload = {
                "candidate_pool": str(args.candidate_pool),
                "derived_pool": args.derived_pool,
                "choices": selected,
                "checker": meta,
                "incremental_checker": "UNSAT",
            }
            args.output.write_text(json.dumps(payload, indent=2) + "\n")
            print(f"COVER rows={len(selected)} output={args.output}", flush=True)
            return 0
        k, l = witness
        added = []
        for _ in range(args.add_per_round):
            try:
                item = next(ordered_unused)
            except StopIteration:
                if not added:
                    print("EXHAUSTED candidates without cover", flush=True)
                    return 2
                break
            h, p, a, b, ord2, ord3 = item
            row = exact_greedy.as_row(item, (a * k + b * l) % h)
            selected.append(row)
            used.add(p)
            checker.add_row(row)
            added.append((p, h))
        added_summary = (
            f"count={len(added)} first={added[0]} last={added[-1]}"
            if added
            else "count=0"
        )
        print(
            f"round={round_no} selected={len(selected)} added={added_summary} "
            f"components={len(checker.components)}",
            flush=True,
        )
        if (
            args.checkpoint
            and args.checkpoint_every > 0
            and round_no % args.checkpoint_every == 0
        ):
            args.checkpoint.write_text(
                json.dumps({"round": round_no, "choices": selected}, indent=2)
                + "\n"
            )
    print(f"ROUND_LIMIT selected={len(selected)}", flush=True)
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
