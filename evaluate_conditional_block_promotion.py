#!/usr/bin/env python3
"""Evaluate promoting a verified conditional fibre into an anchor block."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

from certify_block_plus_overlap_forest import (
    read_fraction,
    recorded_anchor_rows,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("base_certificate", type=Path)
    parser.add_argument("conditional_certificate", type=Path)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument("--block-output", type=Path, required=True)
    parser.add_argument("--period-output", type=Path, required=True)
    args = parser.parse_args()

    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).with_name(
                "certify_conditional_block_extension.py"
            )),
            str(args.pool),
            str(args.base_certificate),
            str(args.conditional_certificate),
            "--output",
            str(args.block_output),
        ],
        check=True,
    )

    source = json.loads(args.pool.read_text())
    block = json.loads(args.block_output.read_text())
    rows = source["choices"]
    selected_primes = {
        int(row["p"])
        for row in rows
        if args.period % int(row["h"]) == 0
    }
    anchor_primes = {
        int(row["p"]) for row in recorded_anchor_rows(block)
    }

    best: dict[int, tuple[Fraction, Path]] = {}
    for certificate_path in Path(".").glob(
        "order_pool_1050000_conditional_fibre*_certificate.json"
    ):
        try:
            certificate = json.loads(certificate_path.read_text())
            outside_prime = int(certificate["outside_prime"])
            conditional_anchors = {
                int(prime) for prime in certificate["anchor_primes"]
            }
            if (
                outside_prime in anchor_primes
                or outside_prime not in selected_primes
                or not conditional_anchors <= anchor_primes
            ):
                continue
            verification_path = Path(
                str(certificate_path).replace(
                    "_certificate.json",
                    "_verification.json",
                )
            )
            if not verification_path.exists():
                continue
            verification = json.loads(verification_path.read_text())
            if (
                verification.get("verified") is not True
                or Path(verification.get("certificate", "")).name
                != certificate_path.name
            ):
                continue
            lower = read_fraction(
                certificate["forced_intersection_density"]
            )
            if outside_prime not in best or lower > best[outside_prime][0]:
                best[outside_prime] = (lower, certificate_path)
        except (KeyError, TypeError, ValueError, json.JSONDecodeError):
            continue

    command = [
        sys.executable,
        str(Path(__file__).with_name(
            "certify_ranked_period_conditional_star.py"
        )),
        str(args.pool),
        str(args.block_output),
        "--period",
        str(args.period),
    ]
    for _, certificate_path in sorted(best.values(), key=lambda item: str(
        item[1]
    )):
        command.extend([
            "--conditional-certificate",
            str(certificate_path),
        ])
    command.extend(["--output", str(args.period_output)])
    result = subprocess.run(command, check=False)
    if not args.period_output.exists():
        return result.returncode
    period_certificate = json.loads(args.period_output.read_text())
    upper = read_fraction(
        period_certificate["union_density_upper_bound"]
    )
    print(
        f"promoted={block['extra_row']['p']} "
        f"anchors={len(anchor_primes)} conditionals={len(best)} "
        f"upper={upper} proved_no_cover={upper < 1}",
        flush=True,
    )
    return 0 if upper < 1 else 1


if __name__ == "__main__":
    raise SystemExit(main())
