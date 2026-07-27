#!/usr/bin/env python3
"""Screen recorded small-anchor templates for stronger conditional edges.

This is a discovery tool.  A reported improvement must still be emitted by
the dedicated certificate generator and replayed by its independent verifier
before it is used in a proof.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
from fractions import Fraction
from pathlib import Path

try:
    import numpy as np
except ImportError:  # pragma: no cover - explicit runtime guard
    np = None

from certify_block_plus_overlap_forest import recorded_anchor_rows
from certify_conditional_fibre_overlap import kernel_image
from certify_forced_pair_overlap import joint_index


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def conditional_intersection(
    outside: dict,
    anchors: list[dict],
    normalizer_index: int,
) -> tuple[Fraction, int, int]:
    image, _generators = kernel_image(outside, anchors)
    moduli = tuple(int(row["h"]) for row in anchors)
    shape = moduli[:normalizer_index] + moduli[normalizer_index + 1:]
    counts = np.zeros(shape, dtype=np.int64)
    for values in image:
        if values[normalizer_index] == 0:
            continue
        index = values[:normalizer_index] + values[normalizer_index + 1:]
        counts[index] += 1
    uncovered = counts
    for axis in range(uncovered.ndim):
        uncovered = (
            np.expand_dims(uncovered.sum(axis=axis), axis) - uncovered
        )
    minimum_covered = len(image) - int(uncovered.max())
    return (
        Fraction(
            minimum_covered,
            int(outside["h"]) * len(image),
        ),
        len(image),
        math.prod(shape),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("block_certificate", type=Path)
    parser.add_argument("--outside-primes", required=True)
    parser.add_argument("--period-certificate", type=Path)
    parser.add_argument(
        "--extra-anchor-templates",
        default="",
        help=(
            "semicolon-separated anchor lists with @normalizer, for example "
            "'5,11,13,19,23,601@5;7,11,13,17,23,599@7'"
        ),
    )
    parser.add_argument(
        "--extra-only",
        action="store_true",
        help="screen only command-line templates, not recorded templates",
    )
    parser.add_argument(
        "--max-target-combinations",
        type=int,
        default=5_000_000,
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if np is None:
        raise RuntimeError("NumPy is required for template screening")

    source = json.loads(args.pool.read_text())
    block = json.loads(args.block_certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    block_primes = {
        int(row["p"]) for row in recorded_anchor_rows(block)
    }
    outside_primes = [
        int(value) for value in args.outside_primes.split(",") if value
    ]
    missing = set(outside_primes) - by_prime.keys()
    if missing:
        raise RuntimeError(f"outside primes absent from pool: {sorted(missing)}")

    baseline_by_prime: dict[int, Fraction] = {}
    if args.period_certificate:
        period = json.loads(args.period_certificate.read_text())
        baseline_by_prime = {
            int(record["outside_prime"]): read_fraction(
                record["used_intersection_lower_bound"]
            )
            for record in period["outside_row_lowers"]
        }

    templates: dict[tuple[tuple[int, ...], int], str] = {}
    if not args.extra_only:
        for raw_path in glob.glob("*conditional*fibre*certificate.json"):
            path = Path(raw_path)
            try:
                certificate = json.loads(path.read_text())
                if certificate.get("schema") != "conditional_fibre_overlap_v1":
                    continue
                anchor_primes = tuple(
                    int(prime) for prime in certificate["anchor_primes"]
                )
                normalizer = int(certificate["normalizer_anchor_prime"])
                if set(anchor_primes) <= block_primes:
                    templates[(anchor_primes, normalizer)] = str(path)
            except (KeyError, TypeError, ValueError, json.JSONDecodeError):
                continue
    for raw_template in args.extra_anchor_templates.split(";"):
        if not raw_template:
            continue
        try:
            raw_anchors, raw_normalizer = raw_template.split("@", 1)
            anchor_primes = tuple(
                int(value) for value in raw_anchors.split(",") if value
            )
            normalizer = int(raw_normalizer)
        except ValueError as error:
            raise SystemExit(
                f"invalid extra anchor template: {raw_template}"
            ) from error
        if (
            len(anchor_primes) < 2
            or len(set(anchor_primes)) != len(anchor_primes)
            or normalizer not in anchor_primes
            or not set(anchor_primes) <= block_primes
        ):
            raise SystemExit(
                f"invalid extra anchor template: {raw_template}"
            )
        templates[(anchor_primes, normalizer)] = "command_line"

    results = []
    for outside_prime in outside_primes:
        outside = by_prime[outside_prime]
        baseline = baseline_by_prime.get(outside_prime, Fraction(0))
        best = baseline
        for (anchor_primes, normalizer), source_path in templates.items():
            anchors = [by_prime[prime] for prime in anchor_primes]
            normalizer_index = anchor_primes.index(normalizer)
            if joint_index(outside, anchors[normalizer_index]) != 1:
                continue
            target_combinations = math.prod(
                int(row["h"])
                for index, row in enumerate(anchors)
                if index != normalizer_index
            )
            if target_combinations > args.max_target_combinations:
                continue
            value, image_size, checked_targets = conditional_intersection(
                outside,
                anchors,
                normalizer_index,
            )
            if value <= best:
                continue
            best = value
            results.append(
                {
                    "outside_prime": outside_prime,
                    "anchor_primes": list(anchor_primes),
                    "normalizer_anchor_prime": normalizer,
                    "template_source": source_path,
                    "kernel_image_size": image_size,
                    "target_combinations": checked_targets,
                    "baseline": fraction_payload(baseline),
                    "intersection": fraction_payload(value),
                    "improvement": fraction_payload(value - baseline),
                }
            )
            print(
                f"outside={outside_prime} anchors={anchor_primes} "
                f"intersection={value} improvement={value - baseline}",
                flush=True,
            )

    best_by_prime = {}
    for record in results:
        prime = int(record["outside_prime"])
        value = read_fraction(record["intersection"])
        if (
            prime not in best_by_prime
            or value > read_fraction(best_by_prime[prime]["intersection"])
        ):
            best_by_prime[prime] = record
    output = {
        "pool": str(args.pool),
        "block_certificate": str(args.block_certificate),
        "period_certificate": (
            str(args.period_certificate)
            if args.period_certificate
            else None
        ),
        "outside_primes": outside_primes,
        "template_count": len(templates),
        "improved_outside_count": len(best_by_prime),
        "best_improvements": [
            best_by_prime[prime] for prime in sorted(best_by_prime)
        ],
        "status": "discovery_only_requires_certificate_and_replay",
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"templates={len(templates)} "
        f"improved={len(best_by_prime)}/{len(outside_primes)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
