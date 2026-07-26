#!/usr/bin/env python3
"""Prune fibres that cannot occur essentially in a finite affine cover.

At the largest q-adic exponent among the surviving moduli, freeze every
lower digit and every other CRT component.  Each fibre at that exponent then
cuts out one affine line in the top-digit plane F_q^2.  If even freely choosing
the offsets of all available directions cannot cover that plane, those fibres
are redundant in every cover and may all be removed.  Repeating the argument
gives a rigorously equivalent finite candidate core.

The optional SAT screen is deliberately a relaxation: it forgets all lower
digit compatibility between offsets.  UNSAT is therefore a safe pruning
certificate; SAT only means that this local obstruction is inconclusive.

For a declared odd algebraic prime q, the q-th-power identity covers qZ^2.
At a first q-adic top digit this is exactly the origin of F_q^2, so a
punctured-plane SAT relaxation replaces the full-plane test.  At higher
q-adic exponents, lower-digit states outside qZ^2 still require a full
top-plane cover.
"""

from __future__ import annotations

import argparse
import collections
import itertools
import json
import math
import os
import sys
from pathlib import Path

import exact_uncovered


def smallest_prime_factors(limit: int) -> list[int]:
    """Return an SPF table for every integer through limit."""
    spf = [0] * (limit + 1)
    if limit >= 1:
        spf[1] = 1
    for prime in range(2, limit + 1):
        if spf[prime]:
            continue
        spf[prime] = prime
        if prime * prime > limit:
            continue
        for multiple in range(prime * prime, limit + 1, prime):
            if not spf[multiple]:
                spf[multiple] = prime
    return spf


def factor_from_spf(value: int, spf: list[int]) -> dict[int, int]:
    out: dict[int, int] = {}
    while value > 1:
        prime = spf[value]
        if not prime:
            raise AssertionError(f"missing SPF entry for {value}")
        out[prime] = out.get(prime, 0) + 1
        value //= prime
    return out


def valuation(value: int, prime: int) -> int:
    exponent = 0
    while value % prime == 0:
        exponent += 1
        value //= prime
    return exponent


def direction(row: dict, prime: int) -> tuple[int, int]:
    a = int(row["a"]) % prime
    b = int(row["b"]) % prime
    if a:
        return 1, b * pow(a, -1, prime) % prime
    if not b:
        raise RuntimeError(f"non-surjective q-component in row {row}")
    return 0, 1


def tangent_capacity_impossible(
    prime: int, counts: collections.Counter
) -> bool:
    """Apply Blokhuis and Blokhuis--Brouwer with direction capacities."""
    total = sum(counts.values())
    if prime % 2 == 0:
        return False
    if max(counts.values()) >= prime or len(counts) >= prime + 1:
        return False
    # Dualize the affine line cover and add the point P dual to the line at
    # infinity.  Because at least one direction is absent, P is essential.
    # Reduce to a minimal blocking subset B containing P.  It is nontrivial:
    # a projective line through P would require q points in one available
    # direction, while a projective line away from P would require all q+1
    # directions.  Thus |B| >= 3(q+1)/2.  Moreover P lies on at least
    # 2q+1-|B| tangents,
    # hence at most |B|-q direction rays through P contain the other |B|-1
    # points.  Those rays retain the declared direction capacities.
    minimum_blocking_size = 3 * (prime + 1) // 2
    capacities = sorted(counts.values(), reverse=True)
    prefix = [0]
    for capacity in capacities:
        prefix.append(prefix[-1] + capacity)
    for blocking_size in range(minimum_blocking_size, total + 2):
        occupied_directions = min(blocking_size - prime, len(capacities))
        if prefix[occupied_directions] >= blocking_size - 1:
            return False
    return True


def cheap_impossible(prime: int, counts: collections.Counter) -> bool:
    """Finite-plane line-cover obstructions before invoking SAT."""
    total = sum(counts.values())
    if total < prime:
        return True
    # At least q rows in one direction may choose all q distinct parallel
    # lines and cover F_q^2 by themselves.  The overlap inequality below
    # assumes a proper parallel class of size s < q; applying it to the raw
    # row multiplicity s >= q would count forced duplicate offsets as if
    # they created additional overlap and can yield a false obstruction.
    if max(counts.values()) >= prime:
        return False
    if tangent_capacity_impossible(prime, counts):
        return True
    # s parallel lines cover sq points.  Every nonparallel line meets all s,
    # hence adds at most q-s new points.
    return any(
        total * prime - size * (total - size) < prime * prime
        for size in counts.values()
    )


def normalized_plane_cover_enumeration(
    prime: int, counts: collections.Counter, max_states: int
):
    """Exactly test a manageable direction-capacity plane cover.

    We may use every available row and choose distinct offsets within a
    direction: duplicates never help, while adding a new parallel line never
    hurts.  Translation of F_q^2 independently normalizes one chosen offset
    in any two distinct directions to zero.  The remaining finite search is
    exhaustive.  Return None when its raw state count exceeds the limit.
    """
    if max(counts.values()) >= prime or len(counts) >= prime + 1:
        return True
    directions = sorted(counts, key=lambda item: (counts[item], item))
    normalized = set(directions[:2])
    option_offsets = {}
    states = 1
    for item in directions:
        capacity = counts[item]
        if item in normalized:
            option_count = math.comb(prime - 1, capacity - 1)
            if states * option_count > max_states:
                return None
            choices = [
                (0,) + rest
                for rest in itertools.combinations(range(1, prime), capacity - 1)
            ]
        else:
            option_count = math.comb(prime, capacity)
            if states * option_count > max_states:
                return None
            choices = list(itertools.combinations(range(prime), capacity))
        option_offsets[item] = choices
        states *= option_count

    full = (1 << (prime * prime)) - 1
    option_masks = {}
    for kind, slope in directions:
        masks = []
        for offsets in option_offsets[(kind, slope)]:
            allowed = set(offsets)
            mask = 0
            for x in range(prime):
                for y in range(prime):
                    value = y if kind == 0 else (x + slope * y) % prime
                    if value in allowed:
                        mask |= 1 << (x * prime + y)
            masks.append(mask)
        option_masks[(kind, slope)] = masks
    # Small option sets first improves both discovery and pruning.
    directions.sort(key=lambda item: len(option_masks[item]))

    def search(position: int, covered: int) -> bool:
        if covered == full:
            return True
        if position == len(directions):
            return False
        uncovered = full ^ covered
        remaining_upper = 0
        for item in directions[position:]:
            remaining_upper += max(
                (mask & uncovered).bit_count() for mask in option_masks[item]
            )
        if remaining_upper < uncovered.bit_count():
            return False
        item = directions[position]
        return any(
            search(position + 1, covered | mask) for mask in option_masks[item]
        )

    return search(0, 0)


def relaxed_plane_cover_sat(
    prime: int,
    counts: collections.Counter,
    solver_name: str,
    normalize_translation: bool = True,
    excluded_points: frozenset[tuple[int, int]] = frozenset(),
) -> tuple[bool, dict]:
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import CNF, IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    pool = IDPool()
    variables = {
        (kind, slope, offset): pool.id((kind, slope, offset))
        for kind, slope in counts
        for offset in range(prime)
    }
    cnf = CNF()
    for (kind, slope), capacity in counts.items():
        literals = [variables[(kind, slope, offset)] for offset in range(prime)]
        if capacity < prime:
            # Selecting additional parallel lines cannot destroy a cover, so
            # every feasible at-most-capacity solution extends to one using
            # the full capacity.  Equality propagates far better in SAT.
            encoded = CardEnc.equals(
                literals,
                bound=capacity,
                vpool=pool,
                encoding=EncType.seqcounter,
            )
            cnf.extend(encoded.clauses)
    normalized_directions = []
    if excluded_points:
        normalize_translation = False
    if normalize_translation and len(counts) >= 2:
        # Every direction has positive capacity.  Pick one selected line from
        # each of two distinct directions and translate their intersection to
        # the origin.  Translation acts transitively on such ordered pairs,
        # so requiring offset zero in both directions loses no covers and
        # removes a q^2 symmetry orbit from the SAT search.
        normalized_directions = sorted(
            counts, key=lambda item: (counts[item], item)
        )[:2]
        for kind, slope in normalized_directions:
            cnf.append([variables[(kind, slope, 0)]])
    directions = list(counts)
    for x in range(prime):
        for y in range(prime):
            if (x, y) in excluded_points:
                continue
            clause = []
            for kind, slope in directions:
                offset = y if kind == 0 else (x + slope * y) % prime
                clause.append(variables[(kind, slope, offset)])
            cnf.append(clause)
    with Solver(name=solver_name, bootstrap_with=cnf.clauses) as solver:
        sat = solver.solve()
        stats = solver.accum_stats()
    return sat, {
        "variables": pool.top,
        "clauses": len(cnf.clauses),
        "solver": solver_name,
        "normalized_directions": normalized_directions,
        "excluded_points": sorted(excluded_points),
        "stats": stats,
    }


def brute_essential_rows(
    prime: int,
    row_directions: list[tuple[int, int]],
    excluded_points: frozenset[tuple[int, int]] = frozenset(),
):
    """Return row positions that can be essential in some top-plane cover."""
    essential = set()
    points = [
        (x, y)
        for x in range(prime)
        for y in range(prime)
        if (x, y) not in excluded_points
    ]
    line_points = {}
    for row_index, (kind, slope) in enumerate(row_directions):
        for offset in range(prime):
            line_points[(row_index, offset)] = {
                point_index
                for point_index, (x, y) in enumerate(points)
                if (y if kind == 0 else (x + slope * y) % prime) == offset
            }
    for offsets in itertools.product(range(prime), repeat=len(row_directions)):
        covers = [0] * len(points)
        for row_index, offset in enumerate(offsets):
            for point_index in line_points[(row_index, offset)]:
                covers[point_index] += 1
        if 0 in covers:
            continue
        for row_index, offset in enumerate(offsets):
            if any(
                covers[point_index] == 1
                for point_index in line_points[(row_index, offset)]
            ):
                essential.add(row_index)
        if len(essential) == len(row_directions):
            break
    return essential


def sat_essential_rows(
    prime: int,
    row_directions: list[tuple[int, int]],
    solver_name: str,
    excluded_points: frozenset[tuple[int, int]] = frozenset(),
):
    """Find rows that can uniquely cover a point in some plane cover."""
    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    from pysat.card import CardEnc, EncType  # type: ignore
    from pysat.formula import CNF, IDPool  # type: ignore
    from pysat.solvers import Solver  # type: ignore

    pool = IDPool()
    variable = {
        (row, offset): pool.id((row, offset))
        for row in range(len(row_directions))
        for offset in range(prime)
    }
    cnf = CNF()
    for row in range(len(row_directions)):
        cnf.extend(
            CardEnc.equals(
                [variable[(row, offset)] for offset in range(prime)],
                bound=1,
                vpool=pool,
                encoding=EncType.seqcounter,
            ).clauses
        )
    points = [
        (x, y)
        for x in range(prime)
        for y in range(prime)
        if (x, y) not in excluded_points
    ]
    point_offsets = []
    for x, y in points:
        offsets = [
            y if kind == 0 else (x + slope * y) % prime
            for kind, slope in row_directions
        ]
        point_offsets.append(offsets)
        cnf.append(
            [variable[(row, offset)] for row, offset in enumerate(offsets)]
        )
    essential = set()
    with Solver(name=solver_name, bootstrap_with=cnf.clauses) as solver:
        for distinguished in range(len(row_directions)):
            for offsets in point_offsets:
                assumptions = [variable[(distinguished, offsets[distinguished])]]
                assumptions.extend(
                    -variable[(row, offsets[row])]
                    for row in range(len(row_directions))
                    if row != distinguished
                )
                if solver.solve(assumptions=assumptions):
                    essential.add(distinguished)
                    break
    return essential


def derived_lower_digit_potentially_essential(
    rows: list[dict],
    indices: list[int],
    prime: int,
    exponent: int,
    solver_name: str,
    algebraic: bool,
):
    """Relaxed essential-row screen using a fixed derived target mod q.

    For a maximal q^e component with e >= 2, freeze the lowest q-adic
    digits of (k,l).  A row whose target is fixed modulo q is available only
    on its declared line in that lower-digit plane.  If the available rows
    cannot cover the free top-digit plane even after choosing every top
    offset independently, none of them can be essential in that state.

    Rows whose target fixes the top digit itself are deliberately unsupported:
    treating their top offset as free would still be a relaxation, but marking
    all rows in a SAT state as potentially essential is only implemented for
    the common square-free-power case v_q(target_modulus) < e.
    """
    if exponent < 2:
        return set(indices), {"skipped": "exponent_below_two"}
    target_valuations = {
        index: valuation(int(rows[index]["target_modulus"]), prime)
        for index in indices
    }
    if not any(value for value in target_valuations.values()):
        return set(indices), {"skipped": "no_fixed_lower_digit"}
    if any(value >= exponent for value in target_valuations.values()):
        return set(indices), {"skipped": "fixed_top_digit_present"}

    possible_states = []
    impossible_states = []
    potentially_essential = set()
    for x in range(prime):
        for y in range(prime):
            if algebraic and x == 0 and y == 0:
                continue
            available = []
            for index in indices:
                target_valuation = target_valuations[index]
                if not target_valuation:
                    available.append(index)
                    continue
                row = rows[index]
                if (
                    int(row["a"]) * x
                    + int(row["b"]) * y
                    - int(row["target_residue"])
                ) % prime == 0:
                    available.append(index)
            counts = collections.Counter(
                direction(rows[index], prime) for index in available
            )
            if not counts or cheap_impossible(prime, counts):
                possible = False
                sat_meta = None
            elif max(counts.values()) >= prime or len(counts) >= prime + 1:
                possible = True
                sat_meta = None
            else:
                possible, sat_meta = relaxed_plane_cover_sat(
                    prime,
                    counts,
                    solver_name,
                )
            state_record = {
                "state": [x, y],
                "available_rows": len(available),
                "directions": len(counts),
                "max_parallel": max(counts.values(), default=0),
                "sat": sat_meta,
            }
            if possible:
                potentially_essential.update(available)
                possible_states.append(state_record)
            else:
                impossible_states.append(state_record)
    return potentially_essential, {
        "required_states": prime * prime - int(algebraic),
        "possible_state_count": len(possible_states),
        "impossible_state_count": len(impossible_states),
        "possible_states": possible_states,
        "impossible_states": impossible_states,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sat-max-prime", type=int, default=0)
    parser.add_argument(
        "--skip-primes",
        default="",
        help="comma-separated components not to prune (for algebraic domains)",
    )
    parser.add_argument(
        "--algebraic-primes",
        default="",
        help=(
            "comma-separated odd q for which qZ^2 is covered by the "
            "X^q+1 identity"
        ),
    )
    parser.add_argument(
        "--algebraic-sat-max-prime",
        type=int,
        default=100,
        help="largest algebraic q on which to run the punctured-plane SAT test",
    )
    parser.add_argument(
        "--algebraic-essential-max-prime",
        type=int,
        default=19,
        help="largest algebraic q on which to SAT-prune inessential rows",
    )
    parser.add_argument(
        "--brute-essential-states",
        type=int,
        default=1_000_000,
        help="exactly prune top-plane rows that can never be essential",
    )
    parser.add_argument("--sat-essential-max-prime", type=int, default=3)
    parser.add_argument(
        "--derived-lower-digit-pruning",
        action="store_true",
        help=(
            "for derived pools, use target_residue/target_modulus on the "
            "lowest q-adic digit to prune maximal q^e rows that can never "
            "be essential"
        ),
    )
    parser.add_argument("--enumeration-max-states", type=int, default=10_000_000)
    parser.add_argument("--solver", default="cadical195")
    args = parser.parse_args()
    skipped_primes = {
        int(value) for value in args.skip_primes.split(",") if value
    }
    algebraic_primes = {
        int(value) for value in args.algebraic_primes.split(",") if value
    }
    if any(prime % 2 == 0 for prime in algebraic_primes):
        raise SystemExit("--algebraic-primes currently supports odd primes only")

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    max_modulus = max((int(row["h"]) for row in rows), default=1)
    # An SPF table is excellent for dense pools of modest moduli, but its
    # memory is proportional to the largest h rather than to the row count.
    # A handful of newly discovered high-order fibres can otherwise allocate
    # hundreds of megabytes merely to factor a few pure-power products.
    if max_modulus <= 5_000_000:
        spf = smallest_prime_factors(max_modulus)
        factors = [
            factor_from_spf(int(row["h"]), spf)
            for row in rows
        ]
    else:
        factors = [
            exact_uncovered.factor(int(row["h"]))
            for row in rows
        ]

    alive = set(range(len(rows)))
    audit = []
    round_no = 0
    while True:
        round_no += 1
        remove = set()
        inessential_findings = []
        derived_inessential_findings = []
        groups: dict[int, dict[int, list[int]]] = collections.defaultdict(
            lambda: collections.defaultdict(list)
        )
        for index in alive:
            for prime, exponent in factors[index].items():
                groups[prime][exponent].append(index)
        for prime in sorted(groups, reverse=True):
            if prime in skipped_primes:
                continue
            exponent = max(groups[prime])
            indices = groups[prime][exponent]
            counts = collections.Counter(direction(rows[index], prime) for index in indices)
            reason = None
            sat_meta = None
            punctured = prime in algebraic_primes and exponent == 1
            if punctured:
                # For a q-th-power construction the algebraic identity covers
                # exactly the origin of the first top-digit plane.  If even
                # freely chosen offsets cannot cover the other q^2-1 points,
                # these maximal fibres cannot be essential.
                if sum(counts.values()) < prime:
                    reason = "punctured_plane_line_count_bound"
                elif (
                    prime <= args.algebraic_sat_max_prime
                    and max(counts.values()) < prime
                    and len(counts) < prime + 1
                ):
                    sat, sat_meta = relaxed_plane_cover_sat(
                        prime,
                        counts,
                        args.solver,
                        normalize_translation=False,
                        excluded_points=frozenset({(0, 0)}),
                    )
                    if not sat:
                        reason = "punctured_plane_relaxation_unsat"
            else:
                if tangent_capacity_impossible(prime, counts):
                    reason = "blokhuis_brouwer_tangent_capacity_bound"
                elif cheap_impossible(prime, counts):
                    reason = "parallel_class_union_bound"
                else:
                    enumerated = normalized_plane_cover_enumeration(
                        prime, counts, args.enumeration_max_states
                    )
                    if enumerated is False:
                        reason = "normalized_plane_cover_exhaustion"
            if not punctured and reason is None and (
                args.sat_max_prime
                and prime <= args.sat_max_prime
                and max(counts.values()) < prime
            ):
                sat, sat_meta = relaxed_plane_cover_sat(prime, counts, args.solver)
                if not sat:
                    reason = "relaxed_plane_cover_unsat"
            if reason:
                remove.update(indices)
                audit.append(
                    {
                        "round": round_no,
                        "prime": prime,
                        "exponent": exponent,
                        "rows": len(indices),
                        "directions": len(counts),
                        "max_parallel": max(counts.values()),
                        "reason": reason,
                        "sat": sat_meta,
                    }
                )
            elif args.derived_lower_digit_pruning and exponent >= 2:
                if not all(
                    "target_residue" in rows[index]
                    and "target_modulus" in rows[index]
                    for index in indices
                ):
                    raise RuntimeError(
                        "--derived-lower-digit-pruning requires derived target fields"
                    )
                potentially_essential, derived_meta = (
                    derived_lower_digit_potentially_essential(
                        rows,
                        indices,
                        prime,
                        exponent,
                        args.solver,
                        prime in algebraic_primes,
                    )
                )
                inessential = [
                    index
                    for index in indices
                    if index not in potentially_essential
                ]
                if (
                    inessential
                    and not potentially_essential
                    and derived_meta.get("possible_state_count") == 0
                ):
                    # This group can be peeled in bulk.  If a cover used one
                    # of these rows essentially, freezing its lower state
                    # would produce one of the impossible top-plane covers
                    # just certified above.  Hence one row can be deleted
                    # without losing the cover; the state obstruction is
                    # monotone under deletion, so repeating removes the
                    # entire group.
                    remove.update(indices)
                    audit.append(
                        {
                            "round": round_no,
                            "prime": prime,
                            "exponent": exponent,
                            "rows": len(indices),
                            "directions": len(counts),
                            "max_parallel": max(counts.values()),
                            "reason": (
                                "derived_lower_digit_all_states_impossible"
                            ),
                            "derived_lower_digit": derived_meta,
                            "sat": None,
                        }
                    )
                elif inessential:
                    derived_inessential_findings.append(
                        (
                            prime,
                            exponent,
                            indices,
                            counts,
                            inessential,
                            derived_meta,
                        )
                    )
            elif prime ** len(indices) <= args.brute_essential_states:
                row_directions = [direction(rows[index], prime) for index in indices]
                essential = brute_essential_rows(
                    prime,
                    row_directions,
                    frozenset({(0, 0)}) if punctured else frozenset(),
                )
                inessential = [
                    index for position, index in enumerate(indices) if position not in essential
                ]
                if inessential:
                    inessential_findings.append(
                        (prime, exponent, indices, counts, inessential)
                    )
            elif (
                (punctured and prime <= args.algebraic_essential_max_prime)
                or (not punctured and prime <= args.sat_essential_max_prime)
            ):
                row_directions = [direction(rows[index], prime) for index in indices]
                essential = sat_essential_rows(
                    prime,
                    row_directions,
                    args.solver,
                    frozenset({(0, 0)}) if punctured else frozenset(),
                )
                inessential = [
                    index for position, index in enumerate(indices) if position not in essential
                ]
                if inessential:
                    inessential_findings.append(
                        (prime, exponent, indices, counts, inessential)
                    )
        # Inessentiality of one row can rely on another currently available
        # row, so remove at most one such row before recomputing.  Whole-group
        # impossibility is monotone and may still be applied in parallel.
        if not remove and derived_inessential_findings:
            (
                prime,
                exponent,
                indices,
                counts,
                inessential,
                derived_meta,
            ) = derived_inessential_findings[0]
            remove.add(inessential[0])
            audit.append(
                {
                    "round": round_no,
                    "prime": prime,
                    "exponent": exponent,
                    "rows": len(indices),
                    "directions": len(counts),
                    "max_parallel": max(counts.values()),
                    "reason": (
                        "never_essential_under_derived_lower_digit_states"
                    ),
                    "removed_rows": 1,
                    "removed_prime": int(rows[inessential[0]]["p"]),
                    "derived_lower_digit": derived_meta,
                    "sat": None,
                }
            )
        elif not remove and inessential_findings:
            prime, exponent, indices, counts, inessential = inessential_findings[0]
            remove.add(inessential[0])
            audit.append(
                {
                    "round": round_no,
                    "prime": prime,
                    "exponent": exponent,
                    "rows": len(indices),
                    "directions": len(counts),
                    "max_parallel": max(counts.values()),
                    "reason": "never_essential_in_any_top_plane_cover",
                    "removed_rows": 1,
                    "removed_prime": int(rows[inessential[0]]["p"]),
                    "sat": None,
                }
            )
        if not remove:
            break
        alive.difference_update(remove)
        print(
            f"round={round_no} removed={len(remove)} alive={len(alive)} "
            f"density={sum(1 / int(rows[i]['h']) for i in alive):.12f}",
            flush=True,
        )

    kept = [row for index, row in enumerate(rows) if index in alive]
    result = dict(payload)
    result["source"] = str(args.pool)
    result["choices"] = kept
    result["component_core"] = {
        "input_rows": len(rows),
        "output_rows": len(kept),
        "density": sum(1 / int(row["h"]) for row in kept),
        "skipped_primes": sorted(skipped_primes),
        "algebraic_primes": sorted(algebraic_primes),
        "audit": audit,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"stable rows={len(kept)} density={result['component_core']['density']:.12f} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
