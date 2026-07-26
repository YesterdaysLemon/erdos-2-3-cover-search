#!/usr/bin/env python3
"""Compose component reductions into one locked-branch obstruction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("grand_pool", type=Path)
    parser.add_argument("independent_pool", type=Path)
    parser.add_argument("coarse_reduction", type=Path)
    parser.add_argument("coarse_obstruction", type=Path)
    parser.add_argument("--component-prime", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    grand = json.loads(args.grand_pool.read_text())
    independent = json.loads(args.independent_pool.read_text())
    reduction = json.loads(args.coarse_reduction.read_text())
    obstruction = json.loads(args.coarse_obstruction.read_text())
    component = args.component_prime

    independent_rows = [
        row
        for row in grand["choices"]
        if int(row["h"]) % component != 0
    ]
    incident_rows = [
        row
        for row in grand["choices"]
        if int(row["h"]) % component == 0
    ]
    independent_primes = sorted(
        int(row["p"]) for row in independent_rows
    )
    materialized_primes = sorted(
        int(row["p"]) for row in independent["choices"]
    )
    if independent_primes != materialized_primes:
        raise RuntimeError("independent pool does not match grand pool")
    for row in incident_rows:
        if math.gcd(
            int(row["a"]),
            int(row["b"]),
            component,
        ) != 1:
            raise RuntimeError(
                f"incident row p={row['p']} is degenerate"
            )
    if len(incident_rows) >= component:
        raise RuntimeError(
            "incident component lines can meet the density threshold"
        )
    if Path(reduction["pool"]).name != args.independent_pool.name:
        raise RuntimeError("coarse reduction uses another independent pool")
    if not reduction.get("proved_equivalence"):
        raise RuntimeError("coarse reduction is not proved")
    if (
        Path(obstruction["pool"]).name
        != Path(reduction["coarse_pool"]).name
    ):
        raise RuntimeError("coarse obstruction uses another pool")
    if not obstruction.get("proved_no_cover"):
        raise RuntimeError("coarse obstruction is not proved")

    fixed_inactive = [
        int(prime)
        for prime in grand.get("fixed_inactive_primes", ())
    ]
    result = {
        "grand_pool": str(args.grand_pool),
        "independent_pool": str(args.independent_pool),
        "coarse_reduction": str(args.coarse_reduction),
        "coarse_obstruction": str(args.coarse_obstruction),
        "component_prime": component,
        "grand_row_count": len(grand["choices"]),
        "component_independent_rows": len(independent_rows),
        "component_incident_rows": len(incident_rows),
        "component_plane_line_threshold": component,
        "component_hole_lower_bound": (
            component * component
            - len(incident_rows) * component
        ),
        "inherited_fixed_inactive_primes": fixed_inactive,
        "inherited_fixed_inactive_count": len(fixed_inactive),
        "proved_no_cover": True,
        "scope": (
            "all phase assignments of rows present in the grand conditioned "
            "pool; parent rows listed as fixed inactive remain omitted"
        ),
        "proof_chain": [
            (
                f"A point missed by the {len(independent_rows)} "
                f"component-{component}-independent rows leaves at most "
                f"{len(incident_rows)} affine lines on F_{component}^2, "
                f"so it has an uncovered refinement."
            ),
            (
                "The independently certified coarse reduction makes a cover "
                "of the independent pool equivalent to a cover by its "
                "squarefree coarse rows."
            ),
            (
                "The independently certified three-component obstruction "
                "proves that the squarefree coarse rows have no phase cover."
            ),
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"grand={len(grand['choices'])} "
        f"independent={len(independent_rows)} "
        f"incident_{component}={len(incident_rows)} "
        f"fixed_inactive={len(fixed_inactive)}",
        flush=True,
    )
    print("PROVED locked grand-branch obstruction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
