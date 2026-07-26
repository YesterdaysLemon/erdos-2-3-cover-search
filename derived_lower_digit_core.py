#!/usr/bin/env python3
"""Fast certified peel using fixed lower q-adic target digits.

This is the target-aware companion to component_core.py.  It considers only
the new derived-lower-digit obstruction, so it can be interleaved with the
ordinary component cascade without re-running every expensive plane
enumeration after each target-aware deletion.
"""

from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path

import component_core
import exact_uncovered


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--algebraic-primes", default="")
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    algebraic_primes = {
        int(value) for value in args.algebraic_primes.split(",") if value
    }
    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    if not all(
        "target_residue" in row and "target_modulus" in row for row in rows
    ):
        raise RuntimeError("derived target fields are required")
    factors = [
        {
            prime: exponent
            for prime, exponent in exact_uncovered.factor(int(row["h"])).items()
        }
        for row in rows
    ]

    alive = set(range(len(rows)))
    audit = []
    round_no = 0
    while True:
        round_no += 1
        groups: dict[int, dict[int, list[int]]] = collections.defaultdict(
            lambda: collections.defaultdict(list)
        )
        for index in alive:
            for prime, exponent in factors[index].items():
                groups[prime][exponent].append(index)
        remove = set()
        round_records = []
        for prime in sorted(groups, reverse=True):
            exponent = max(groups[prime])
            if exponent < 2:
                continue
            indices = groups[prime][exponent]
            potentially_essential, meta = (
                component_core.derived_lower_digit_potentially_essential(
                    rows,
                    indices,
                    prime,
                    exponent,
                    args.solver,
                    prime in algebraic_primes,
                )
            )
            if (
                not potentially_essential
                and meta.get("possible_state_count") == 0
            ):
                remove.update(indices)
                round_records.append(
                    {
                        "round": round_no,
                        "prime": prime,
                        "exponent": exponent,
                        "rows": len(indices),
                        "reason": (
                            "derived_lower_digit_all_states_impossible"
                        ),
                        "derived_lower_digit": meta,
                    }
                )
        if not remove:
            break
        alive.difference_update(remove)
        audit.extend(round_records)
        print(
            f"round={round_no} removed={len(remove)} alive={len(alive)} "
            f"groups={[(item['prime'], item['exponent'], item['rows']) for item in round_records]}",
            flush=True,
        )

    kept = [row for index, row in enumerate(rows) if index in alive]
    result = dict(payload)
    result["source"] = str(args.pool)
    result["choices"] = kept
    result["derived_lower_digit_core"] = {
        "input_rows": len(rows),
        "output_rows": len(kept),
        "algebraic_primes": sorted(algebraic_primes),
        "audit": audit,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"stable rows={len(kept)} output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
