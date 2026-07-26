#!/usr/bin/env python3
"""Round a power-anchor LP's local pattern mixture to valid full phases.

The factorized LP only distinguishes a target modulo the part of a fibre
visible in its anchor/algebraic period.  This script samples that local option
and then samples uniformly among the compatible full targets modulo h.  Its
JSON output can seed ``local_phase_cegis.py``.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

from power_anchor_capacity_lp import power_target_congruence


def combine_congruences(
    residue_a: int, modulus_a: int, residue_b: int, modulus_b: int
) -> tuple[int, int]:
    """Return the least nonnegative solution and lcm modulus."""
    common = math.gcd(modulus_a, modulus_b)
    delta = residue_b - residue_a
    if delta % common:
        raise RuntimeError("inconsistent local and power target congruences")
    reduced_b = modulus_b // common
    if reduced_b == 1:
        multiplier = 0
    else:
        multiplier = (
            (delta // common)
            * pow(modulus_a // common, -1, reduced_b)
        ) % reduced_b
    combined_modulus = math.lcm(modulus_a, modulus_b)
    return (
        (residue_a + modulus_a * multiplier) % combined_modulus,
        combined_modulus,
    )


def weighted_option(options: list[dict], rng: random.Random) -> dict:
    total = sum(float(option["weight"]) for option in options)
    if total <= 0:
        raise RuntimeError("LP row has no positive option weight")
    threshold = rng.random() * total
    running = 0.0
    for option in options:
        running += float(option["weight"])
        if running >= threshold:
            return option
    return options[-1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("fractional_solution", type=Path)
    parser.add_argument("--power", type=int, required=True)
    parser.add_argument("--seed", type=int, default=203)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = json.loads(args.pool.read_text())
    fractional = json.loads(args.fractional_solution.read_text())
    if int(fractional["power"]) != args.power:
        raise RuntimeError("fractional solution power does not match --power")
    if Path(fractional["pool"]).name != args.pool.name:
        raise RuntimeError("fractional solution pool does not match input pool")

    pool_by_prime = {
        int(row["p"]): row for row in source["choices"]
    }
    fractional_by_prime = {
        int(row["p"]): row for row in fractional["rows"]
    }
    anchors = tuple(int(p) for p in fractional["anchor_primes"])
    anchor_targets = tuple(
        int(target) for target in fractional["normalized_targets"]
    )
    if len(anchors) != len(anchor_targets):
        raise RuntimeError("anchor metadata is inconsistent")

    rng = random.Random(args.seed)
    phases: dict[str, int] = {}
    for prime, target in zip(anchors, anchor_targets):
        row = pool_by_prime.get(prime)
        if row is None:
            raise RuntimeError(f"anchor prime {prime} is missing")
        phases[str(prime)] = target

    compatible = 0
    for prime, row in pool_by_prime.items():
        h = int(row["h"])
        try:
            power_residue, power_modulus = power_target_congruence(
                h, prime, args.power
            )
        except RuntimeError:
            continue
        compatible += 1
        if prime in anchors:
            target = phases[str(prime)]
            if target % power_modulus != power_residue:
                raise RuntimeError(f"anchor {prime} is power-incompatible")
            continue
        local = fractional_by_prime.get(prime)
        if local is None:
            raise RuntimeError(f"compatible prime {prime} is absent from LP output")
        option = weighted_option(local["options"], rng)
        shared = int(option["shared"])
        local_target = int(option["target_mod_shared"])
        base, step = combine_congruences(
            power_residue,
            power_modulus,
            local_target,
            shared,
        )
        if h % step:
            raise RuntimeError(f"combined target modulus does not divide h for {prime}")
        count = h // step
        target = base + step * rng.randrange(count)
        if not 0 <= target < h:
            raise AssertionError("rounded target is out of range")
        phases[str(prime)] = target

    if len(phases) != compatible:
        raise RuntimeError(
            f"rounded {len(phases)} phases but expected {compatible}"
        )
    args.output.write_text(json.dumps(phases) + "\n")
    print(
        f"power={args.power} compatible={compatible} anchors={anchors} "
        f"seed={args.seed} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
