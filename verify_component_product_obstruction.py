#!/usr/bin/env python3
"""Independent replay of a component-product obstruction certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    pool_path = args.certificate.parent / certificate["pool"]
    rows = json.loads(pool_path.read_text())["choices"]
    q = int(certificate["density_prime"])
    x = int(certificate["components"]["x"])
    y = int(certificate["components"]["y"])
    z = int(certificate["components"]["z"])
    fixed = {
        int(row_prime): int(target)
        for row_prime, target in certificate["fixed_targets"].items()
    }

    grouped: dict[int, list[dict]] = {}
    for row in rows:
        grouped.setdefault(int(row["h"]), []).append(row)
    expected_moduli = {q, x, y, z, q * z, x * z, y * z}
    if set(grouped) != expected_moduli:
        raise AssertionError("pool modulus support changed")
    for modulus in expected_moduli:
        if len(grouped[modulus]) != int(
            certificate["row_counts"][str(modulus)]
        ):
            raise AssertionError("row count changed")

    line_components = [
        (q, q),
        (q * z, q),
        (x, x),
        (y, y),
        (z, z),
        (x * z, x),
        (x * z, z),
        (y * z, y),
        (y * z, z),
    ]
    for row_modulus, component in line_components:
        for row in grouped[row_modulus]:
            if (
                int(row["a"]) % component == 0
                and int(row["b"]) % component == 0
            ):
                raise AssertionError("degenerate projected fibre")

    z_rows = grouped[z]
    fixed_rows = [row for row in z_rows if int(row["p"]) in fixed]
    mutable_rows = [row for row in z_rows if int(row["p"]) not in fixed]
    if set(fixed) != {int(row["p"]) for row in fixed_rows}:
        raise AssertionError("certificate fixes a non-z row")

    holes = []
    for k in range(z):
        for l in range(z):
            covered = False
            for row in fixed_rows:
                value = (
                    int(row["a"]) * k + int(row["b"]) * l
                ) % z
                if value == fixed[int(row["p"])] % z:
                    covered = True
                    break
            if not covered:
                holes.append((k, l))

    maxima = {}
    for row in mutable_rows:
        best = 0
        for target in range(z):
            count = sum(
                (
                    int(row["a"]) * k
                    + int(row["b"]) * l
                )
                % z
                == target
                for k, l in holes
            )
            best = max(best, count)
        maxima[int(row["p"])] = best

    residual_upper = len(grouped[q]) + len(grouped[q * z])
    x_holes = max(0, x * x - len(grouped[x]) * x)
    y_holes = max(0, y * y - len(grouped[y]) * y)
    z_holes = max(0, len(holes) - sum(maxima.values()))
    xz_capacity = len(grouped[x * z]) * x
    yz_required = y_holes * z_holes
    yz_capacity = len(grouped[y * z]) * y * z
    checks = {
        "residual_lines_below_threshold": residual_upper < q,
        "x_escape": x_holes > xz_capacity,
        "yz_capacity_contradiction": yz_required > yz_capacity,
    }
    recomputed = {
        "residual_line_scaled_density_upper_bound": residual_upper,
        "x_holes_lower_bound": x_holes,
        "y_holes_lower_bound": y_holes,
        "fixed_z_holes": len(holes),
        "mutable_z_maxima": {
            str(row_prime): value
            for row_prime, value in sorted(maxima.items())
        },
        "z_holes_lower_bound": z_holes,
        "xz_projection_capacity": xz_capacity,
        "required_yz_pairs_lower_bound": yz_required,
        "yz_cross_capacity": yz_capacity,
    }
    expected_maxima = {
        str(item["p"]): int(item["maximum_new_coverage"])
        for item in certificate["mutable_z_maxima"]
    }
    comparisons = {
        "residual_upper": (
            residual_upper
            == int(
                certificate[
                    "residual_line_scaled_density_upper_bound"
                ]
            )
        ),
        "x_holes": (
            x_holes == int(certificate["x_holes_lower_bound"])
        ),
        "y_holes": (
            y_holes == int(certificate["y_holes_lower_bound"])
        ),
        "fixed_z_holes": (
            len(holes) == int(certificate["fixed_z_holes"])
        ),
        "mutable_z_maxima": (
            recomputed["mutable_z_maxima"] == expected_maxima
        ),
        "z_holes": (
            z_holes == int(certificate["z_holes_lower_bound"])
        ),
        "xz_capacity": (
            xz_capacity
            == int(certificate["xz_projection_capacity"])
        ),
        "yz_required": (
            yz_required
            == int(
                certificate["required_yz_pairs_lower_bound"]
            )
        ),
        "yz_capacity": (
            yz_capacity == int(certificate["yz_cross_capacity"])
        ),
    }
    passed = (
        all(checks.values())
        and all(comparisons.values())
        and bool(
            certificate["proved_component_density_obstruction"]
        )
    )
    output = {
        "certificate": str(args.certificate),
        "passed": passed,
        "checks": checks,
        "comparisons": comparisons,
        "recomputed": recomputed,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print("PASS" if passed else "FAIL", flush=True)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
