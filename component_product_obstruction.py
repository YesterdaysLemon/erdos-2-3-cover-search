#!/usr/bin/env python3
"""Certify a product-capacity obstruction to one component-density condition.

The supported geometry has single-component affine fibres modulo x, y, z,
one x*z cross-fibre family, one y*z cross-fibre family, and residual-q
fibres.  Fixed z-fibres plus conservative per-line coverage bounds give a
universal lower bound on the remaining z holes.  A product/capacity argument
then proves that the y*z cross-fibres cannot cover all remaining pairs.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path


def parse_targets(text: str) -> dict[int, int]:
    result = {}
    for item in text.split(","):
        if not item:
            continue
        row_prime, target = item.split(":", 1)
        result[int(row_prime)] = int(target)
    return result


def is_on_line(
    row: dict,
    target: int,
    k: int,
    l: int,
    modulus: int,
) -> bool:
    return (
        int(row["a"]) * k + int(row["b"]) * l - target
    ) % modulus == 0


def require_nondegenerate(rows: list[dict], modulus: int) -> None:
    for row in rows:
        if (
            int(row["a"]) % modulus == 0
            and int(row["b"]) % modulus == 0
        ):
            raise RuntimeError(
                f"row p={row['p']} is degenerate modulo {modulus}"
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--density-prime", type=int, required=True)
    parser.add_argument("--x-prime", type=int, required=True)
    parser.add_argument("--y-prime", type=int, required=True)
    parser.add_argument("--z-prime", type=int, required=True)
    parser.add_argument("--fixed-targets", default="")
    parser.add_argument("--result-output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.pool.read_text())
    rows = payload["choices"]
    if any(int(row["target_modulus"]) != 1 for row in rows):
        raise RuntimeError("certificate requires target_modulus=1")
    q = args.density_prime
    x = args.x_prime
    y = args.y_prime
    z = args.z_prime
    if len({q, x, y, z}) != 4:
        raise RuntimeError("component primes must be distinct")
    allowed_moduli = {
        q,
        x,
        y,
        z,
        q * z,
        x * z,
        y * z,
    }
    actual_moduli = {int(row["h"]) for row in rows}
    if not actual_moduli <= allowed_moduli:
        raise RuntimeError(
            f"unsupported residual moduli: "
            f"{sorted(actual_moduli - allowed_moduli)}"
        )

    by_modulus = {
        modulus: [
            row for row in rows if int(row["h"]) == modulus
        ]
        for modulus in allowed_moduli
    }
    require_nondegenerate(by_modulus[q], q)
    require_nondegenerate(by_modulus[q * z], q)
    require_nondegenerate(by_modulus[x], x)
    require_nondegenerate(by_modulus[y], y)
    require_nondegenerate(by_modulus[z], z)
    require_nondegenerate(by_modulus[x * z], x)
    require_nondegenerate(by_modulus[x * z], z)
    require_nondegenerate(by_modulus[y * z], y)
    require_nondegenerate(by_modulus[y * z], z)
    fixed_targets = parse_targets(args.fixed_targets)
    row_by_prime = {int(row["p"]): row for row in rows}
    if fixed_targets.keys() - row_by_prime.keys():
        raise RuntimeError("a fixed prime is absent from the pool")
    fixed_z_rows = [
        row
        for row in by_modulus[z]
        if int(row["p"]) in fixed_targets
    ]
    if any(
        int(row_by_prime[row_prime]["h"]) != z
        for row_prime in fixed_targets
    ):
        raise RuntimeError(
            "this certificate permits fixed targets only on z-fibres"
        )
    mutable_z_rows = [
        row
        for row in by_modulus[z]
        if int(row["p"]) not in fixed_targets
    ]

    # With no q-free fibre active, q and q*z rows contribute at most one
    # residual line apiece.  Their total scaled density must stay below q.
    residual_line_upper_bound = (
        len(by_modulus[q]) + len(by_modulus[q * z])
    )
    if residual_line_upper_bound >= q:
        raise RuntimeError(
            "residual q-lines alone can meet the density threshold"
        )

    # A union of r affine lines in F_p^2 contains at most r*p points.
    # The bounds are valid even when directions repeat or lines overlap.
    x_holes_lower_bound = max(
        0, x * x - len(by_modulus[x]) * x
    )
    y_holes_lower_bound = max(
        0, y * y - len(by_modulus[y]) * y
    )

    fixed_z_holes = [
        (k, l)
        for k in range(z)
        for l in range(z)
        if not any(
            is_on_line(
                row,
                fixed_targets[int(row["p"])],
                k,
                l,
                z,
            )
            for row in fixed_z_rows
        )
    ]
    mutable_z_maxima = []
    for row in mutable_z_rows:
        target_counts = Counter(
            (
                int(row["a"]) * k + int(row["b"]) * l
            )
            % z
            for k, l in fixed_z_holes
        )
        mutable_z_maxima.append(
            {
                "p": int(row["p"]),
                "maximum_new_coverage": max(
                    target_counts.values(), default=0
                ),
            }
        )
    z_holes_lower_bound = max(
        0,
        len(fixed_z_holes)
        - sum(
            item["maximum_new_coverage"]
            for item in mutable_z_maxima
        ),
    )

    # Every x*z row projects to at most one x-line.  If the universal x-hole
    # lower bound exceeds their total projected capacity, one x-hole avoids
    # every x*z row, whatever their phases.
    xz_projection_capacity = len(by_modulus[x * z]) * x
    x_escape = x_holes_lower_bound > xz_projection_capacity

    # At that escaped x, every pair in U_y x U_z must be covered by a y*z
    # row.  Each such row covers at most y*z pairs.
    required_yz_pairs_lower_bound = (
        y_holes_lower_bound * z_holes_lower_bound
    )
    yz_cross_capacity = len(by_modulus[y * z]) * y * z
    capacity_contradiction = (
        required_yz_pairs_lower_bound > yz_cross_capacity
    )
    proved = x_escape and capacity_contradiction

    result = {
        "pool": str(args.pool),
        "density_prime": q,
        "components": {"x": x, "y": y, "z": z},
        "row_counts": {
            str(modulus): len(by_modulus[modulus])
            for modulus in sorted(allowed_moduli)
        },
        "fixed_targets": {
            str(row_prime): target
            for row_prime, target in sorted(fixed_targets.items())
        },
        "fixed_z_rows": len(fixed_z_rows),
        "mutable_z_rows": len(mutable_z_rows),
        "residual_line_scaled_density_upper_bound": (
            residual_line_upper_bound
        ),
        "residual_density_threshold": q,
        "x_holes_lower_bound": x_holes_lower_bound,
        "xz_projection_capacity": xz_projection_capacity,
        "x_escape": x_escape,
        "y_holes_lower_bound": y_holes_lower_bound,
        "fixed_z_holes": len(fixed_z_holes),
        "mutable_z_maxima": mutable_z_maxima,
        "z_holes_lower_bound": z_holes_lower_bound,
        "required_yz_pairs_lower_bound": (
            required_yz_pairs_lower_bound
        ),
        "yz_cross_capacity": yz_cross_capacity,
        "capacity_contradiction": capacity_contradiction,
        "proved_component_density_obstruction": proved,
        "scope": (
            "all phase assignments preserving the supplied fixed z-fibres"
        ),
    }
    args.result_output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"fixed_z_holes={len(fixed_z_holes)} "
        f"mutable_z_max_sum="
        f"{sum(item['maximum_new_coverage'] for item in mutable_z_maxima)} "
        f"z_holes_lb={z_holes_lower_bound}",
        flush=True,
    )
    print(
        f"x_holes_lb={x_holes_lower_bound} "
        f"xz_capacity={xz_projection_capacity} "
        f"yz_required_lb={required_yz_pairs_lower_bound} "
        f"yz_capacity={yz_cross_capacity}",
        flush=True,
    )
    print(
        "PROVED locked component-density obstruction"
        if proved
        else "NO obstruction from these bounds",
        flush=True,
    )
    return 0 if proved else 1


if __name__ == "__main__":
    raise SystemExit(main())
