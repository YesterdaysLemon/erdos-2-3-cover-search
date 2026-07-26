#!/usr/bin/env python3
"""Discover pure-(2,3)-order fibres by scanning small group multipliers.

The exhaustive order scanner factors gcd(2^h-1, 3^h-1), which becomes
expensive when h has hundreds of millions of bits.  This complementary
construction scan tests primes of the form

    p = multiplier * h + 1

for selected h = 2^a 3^b.  If both 2 and 3 have joint order exactly h
modulo p, the resulting affine fibre is emitted.  The multiplier range is
finite and explicit; finding no row is not an exclusion theorem.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

def pure23_values(start: int, end: int) -> list[int]:
    values = set()
    power2 = 1
    while power2 <= end:
        power3 = 1
        while power2 * power3 <= end:
            value = power2 * power3
            if value >= start:
                values.add(value)
            power3 *= 3
        power2 *= 2
    return sorted(values)


def order_from_factorization(
    base: int,
    prime: int,
    group_order: int,
    prime_factors: tuple[int, ...],
) -> int:
    order = group_order
    for divisor in prime_factors:
        while order % divisor == 0 and pow(base, order // divisor, prime) == 1:
            order //= divisor
    return order


def primitive_root_from_factorization(
    prime: int,
    group_order: int,
    prime_factors: tuple[int, ...],
) -> int:
    for generator in range(2, prime):
        if all(
            pow(generator, group_order // divisor, prime) != 1
            for divisor in prime_factors
        ):
            return generator
    raise AssertionError(f"no primitive root found modulo {prime}")


def discrete_log_smooth23(
    generator: int,
    target: int,
    order: int,
    modulus: int,
) -> int:
    """Return log_generator(target) for an order of the form 2^a*3^b.

    This is Pohlig--Hellman digit lifting specialized to the only two prime
    factors that occur in this scanner.  It avoids the O(sqrt(order)) memory
    of a generic baby-step/giant-step table.
    """

    remainder = order
    prime_powers = []
    for prime in (2, 3):
        exponent = 0
        prime_power = 1
        while remainder % prime == 0:
            remainder //= prime
            exponent += 1
            prime_power *= prime
        if exponent:
            prime_powers.append((prime, exponent, prime_power))
    if remainder != 1:
        raise ValueError("discrete_log_smooth23 requires a pure-(2,3) order")
    if order == 1:
        if target % modulus != 1:
            raise ValueError("target is outside the trivial subgroup")
        return 0

    inverse_generator = pow(generator, -1, modulus)
    residues = []
    for prime, exponent, prime_power in prime_powers:
        digit_generator = pow(generator, order // prime, modulus)
        digit_values = {
            pow(digit_generator, digit, modulus): digit
            for digit in range(prime)
        }
        residue = 0
        place = 1
        for _ in range(exponent):
            adjusted = (
                target
                * pow(inverse_generator, residue, modulus)
            ) % modulus
            digit_target = pow(
                adjusted,
                order // (place * prime),
                modulus,
            )
            if digit_target not in digit_values:
                raise ValueError("target is outside the generated subgroup")
            residue += digit_values[digit_target] * place
            place *= prime
        if (
            pow(
                target * pow(inverse_generator, residue, modulus) % modulus,
                order // prime_power,
                modulus,
            )
            != 1
        ):
            raise AssertionError("Pohlig--Hellman prime-power lift failed")
        residues.append((residue, prime_power))

    value, value_modulus = residues[0]
    for residue, prime_power in residues[1:]:
        step = (
            (residue - value)
            * pow(value_modulus, -1, prime_power)
        ) % prime_power
        value += value_modulus * step
        value_modulus *= prime_power
        value %= value_modulus
    if value_modulus != order or pow(generator, value, modulus) != target:
        raise AssertionError("Pohlig--Hellman CRT reconstruction failed")
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start-h", type=int, required=True)
    parser.add_argument("--end-h", type=int, required=True)
    parser.add_argument("--start-multiplier", type=int, default=1)
    parser.add_argument("--end-multiplier", type=int, required=True)
    parser.add_argument(
        "--h-values",
        help=(
            "optional comma-separated pure-(2,3) values; otherwise scan "
            "every such value in the declared h interval"
        ),
    )
    parser.add_argument("--known-pool", type=Path)
    parser.add_argument("--progress-every", type=int, default=100_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not (
        1 <= args.start_h <= args.end_h
        and 1 <= args.start_multiplier <= args.end_multiplier
    ):
        raise SystemExit("invalid h or multiplier interval")

    if args.h_values:
        h_values = sorted(
            {
                int(value)
                for value in args.h_values.split(",")
                if value
            }
        )
        if not h_values:
            raise SystemExit("--h-values is empty")
        for value in h_values:
            remainder = value
            while remainder % 2 == 0:
                remainder //= 2
            while remainder % 3 == 0:
                remainder //= 3
            if (
                not args.start_h <= value <= args.end_h
                or remainder != 1
            ):
                raise SystemExit(
                    f"--h-values entry {value} is outside the interval "
                    "or is not pure-(2,3)"
                )
    else:
        h_values = pure23_values(args.start_h, args.end_h)
    if not h_values:
        raise SystemExit("the h interval contains no pure-(2,3) value")

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import gmpy2  # type: ignore
    import sympy  # type: ignore

    known_primes = set()
    if args.known_pool:
        known = json.loads(args.known_pool.read_text())
        if int(known.get("unresolved", 0)):
            raise RuntimeError("known pool has unresolved cofactors")
        known_primes = {int(row["p"]) for row in known["choices"]}

    rows_by_prime = {}
    primality_tests = 0
    prime_candidates = 0
    joint_subgroup_candidates = 0
    next_progress = args.progress_every
    for h_index, h in enumerate(h_values, 1):
        for multiplier in range(
            args.start_multiplier,
            args.end_multiplier + 1,
        ):
            p = multiplier * h + 1
            if p in known_primes or p in rows_by_prime or p in (2, 3):
                continue
            # Every odd prime p has even p-1.  This eliminates all even p
            # before invoking the stronger probable-prime test.
            if p % 2 == 0 and p != 2:
                continue
            primality_tests += 1
            if not gmpy2.is_prime(p):
                continue
            prime_candidates += 1
            if pow(2, h, p) != 1 or pow(3, h, p) != 1:
                continue
            joint_subgroup_candidates += 1

            multiplier_factors = sympy.factorint(multiplier)
            factor_primes = tuple(
                sorted({2, 3, *(int(q) for q in multiplier_factors)})
            )
            group_order = p - 1
            ord2 = order_from_factorization(
                2, p, group_order, factor_primes
            )
            ord3 = order_from_factorization(
                3, p, group_order, factor_primes
            )
            exact_h = math.lcm(ord2, ord3)
            if exact_h != h:
                continue

            primitive_root = primitive_root_from_factorization(
                p, group_order, factor_primes
            )
            subgroup_generator = pow(
                primitive_root,
                group_order // h,
                p,
            )
            log2 = discrete_log_smooth23(
                subgroup_generator, 2, h, p
            )
            log3 = discrete_log_smooth23(
                subgroup_generator, 3, h, p
            )
            rows_by_prime[p] = {
                "h": h,
                "p": p,
                "a": log2,
                "b": log3,
                "ord2": ord2,
                "ord3": ord3,
                "c": 0,
                "multiplier": multiplier,
            }
        if args.progress_every > 0 and (
            h_index == len(h_values)
            or primality_tests >= next_progress
        ):
            print(
                f"h={h} h_index={h_index}/{len(h_values)} "
                f"tests={primality_tests} primes={prime_candidates} "
                f"joint={joint_subgroup_candidates} "
                f"rows={len(rows_by_prime)}",
                flush=True,
            )
            while next_progress <= primality_tests:
                next_progress += args.progress_every

    rows = sorted(
        rows_by_prime.values(),
        key=lambda row: (int(row["h"]), int(row["p"])),
    )
    payload = {
        "scan_kind": "pure23_group_multiplier_scan",
        "complete_order_interval": False,
        "complete_multiplier_box": True,
        "start_h": args.start_h,
        "end_h": args.end_h,
        "h_values": h_values,
        "start_multiplier": args.start_multiplier,
        "end_multiplier": args.end_multiplier,
        "primality_tests": primality_tests,
        "prime_candidates": prime_candidates,
        "joint_subgroup_candidates": joint_subgroup_candidates,
        "known_pool": (
            str(args.known_pool) if args.known_pool else None
        ),
        "unresolved": 0,
        "choices": rows,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"wrote={args.output} h_values={len(h_values)} "
        f"tests={primality_tests} rows={len(rows)} "
        f"density={sum(1 / int(row['h']) for row in rows):.12f}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
