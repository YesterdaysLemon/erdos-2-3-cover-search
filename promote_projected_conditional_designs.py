#!/usr/bin/env python3
"""Promote discovery-only projected designs to replayed proof artifacts.

The discovery search may visit the same outside fibre in several batches.
This script selects the strongest exact design for each outside prime,
regenerates it with the dedicated certificate program, and then invokes the
independently written verifier.  A resumable manifest is rewritten after each
successful replay.
"""

from __future__ import annotations

import argparse
import glob
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path


INVALID_DISCOVERY_BASENAMES = {
    "_tmp_period3139_24block_projected_design_search_paired_top13.json",
}


def read_fraction(payload: dict) -> Fraction:
    return Fraction(int(payload["numerator"]), int(payload["denominator"]))


def fraction_payload(value: Fraction) -> dict:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def collect_best_designs(
    paths: list[Path],
) -> dict[int, tuple[dict, Path]]:
    best: dict[int, tuple[dict, Path]] = {}
    for path in sorted(paths):
        if path.name in INVALID_DISCOVERY_BASENAMES:
            print(f"excluding known-invalid discovery file: {path}", flush=True)
            continue
        payload = json.loads(path.read_text())
        if payload.get("schema") != "projected_conditional_design_search_v1":
            raise RuntimeError(f"unexpected discovery schema in {path}")
        if not str(payload.get("status", "")).startswith("discovery_"):
            raise RuntimeError(f"incomplete discovery status in {path}")
        for result in payload.get("results", []):
            design = result.get("best_design")
            if design is None:
                continue
            outside_prime = int(result["outside_prime"])
            if read_fraction(design["improvement"]) <= 0:
                raise RuntimeError(
                    f"nonpositive winning improvement for p={outside_prime}"
                )
            current = best.get(outside_prime)
            if (
                current is None
                or read_fraction(design["intersection"])
                > read_fraction(current[0]["intersection"])
            ):
                best[outside_prime] = design, path
    return best


def certificate_paths(prefix: str, outside_prime: int) -> tuple[Path, Path]:
    stem = f"{prefix}_conditional_fibre{outside_prime}_autodesign"
    return Path(f"{stem}_certificate.json"), Path(
        f"{stem}_verification.json"
    )


def generator_command(
    pool: Path,
    design: dict,
    outside_prime: int,
    certificate: Path,
) -> list[str]:
    projected = ",".join(
        (
            f"{int(record['prime'])}:"
            f"{int(record['projection_modulus'])}:"
            f"{int(record['residual_modulus'])}"
        )
        for record in design["projected_anchors"]
    )
    base_primes = ",".join(
        str(int(prime)) for prime in design["base_anchor_primes"]
    )
    command = [
        sys.executable,
        "certify_projected_conditional_fibre_overlap.py",
        str(pool),
        "--outside-prime",
        str(outside_prime),
        "--normalizer-anchor-prime",
        str(int(design["normalizer_anchor_prime"])),
        f"--base-anchor-primes={base_primes}",
        "--projected-anchors",
        projected,
        "--base-period",
        str(int(design["base_period"])),
        "--max-tensor-cells",
        str(int(design["tensor_cells"])),
        "--output",
        str(certificate),
    ]
    paired_prime = design.get("paired_projected_prime")
    if paired_prime is not None:
        command.extend(
            [
                "--paired-projected-prime",
                str(int(paired_prime)),
                "--paired-shared-prime",
                str(int(design["paired_shared_prime"])),
                "--paired-projection-modulus",
                str(int(design["paired_projection_modulus"])),
                "--paired-residual-modulus",
                str(int(design["paired_residual_modulus"])),
            ]
        )
    return command


def manifest_payload(
    *,
    pool: Path,
    prefix: str,
    discovery_paths: list[Path],
    selected_count: int,
    records: list[dict],
    status: str,
) -> dict:
    total_improvement = sum(
        (
            read_fraction(record["improvement"])
            for record in records
        ),
        Fraction(0),
    )
    return {
        "schema": "promoted_projected_conditional_designs_v1",
        "pool": str(pool),
        "artifact_prefix": prefix,
        "discovery_files": [str(path) for path in discovery_paths],
        "selected_outside_count": selected_count,
        "completed_outside_count": len(records),
        "records": records,
        "total_improvement": fraction_payload(total_improvement),
        "status": status,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("--discovery-glob", required=True)
    parser.add_argument("--artifact-prefix", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--outside-primes",
        help="optional comma-separated subset after strongest-design selection",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="reuse already replayed records from a matching output manifest",
    )
    args = parser.parse_args()

    discovery_paths = [
        Path(path) for path in glob.glob(args.discovery_glob)
    ]
    if not discovery_paths:
        raise RuntimeError("discovery glob matched no files")
    selected = collect_best_designs(discovery_paths)
    if args.outside_primes:
        requested = {
            int(item) for item in args.outside_primes.split(",") if item
        }
        missing = requested - selected.keys()
        if missing:
            raise RuntimeError(
                f"requested outside primes lack a design: {sorted(missing)}"
            )
        selected = {
            prime: selected[prime] for prime in requested
        }

    existing_records: dict[int, dict] = {}
    if args.resume and args.output.exists():
        existing = json.loads(args.output.read_text())
        if (
            existing.get("schema")
            != "promoted_projected_conditional_designs_v1"
            or existing.get("pool") != str(args.pool)
            or existing.get("artifact_prefix") != args.artifact_prefix
            or int(existing["selected_outside_count"]) != len(selected)
        ):
            raise RuntimeError("resume manifest metadata does not match")
        existing_records = {
            int(record["outside_prime"]): record
            for record in existing.get("records", [])
        }

    records = []
    for index, outside_prime in enumerate(
        sorted(selected),
        start=1,
    ):
        design, discovery_path = selected[outside_prime]
        certificate, verification = certificate_paths(
            args.artifact_prefix,
            outside_prime,
        )
        expected_intersection = read_fraction(design["intersection"])
        existing = existing_records.get(outside_prime)
        if existing is not None:
            certificate_payload = json.loads(certificate.read_text())
            verification_payload = json.loads(verification.read_text())
            if (
                existing.get("certificate") != str(certificate)
                or existing.get("verification") != str(verification)
                or read_fraction(
                    certificate_payload["forced_intersection_density"]
                )
                != expected_intersection
                or not verification_payload.get("verified")
                or verification_payload.get("certificate")
                != str(certificate)
                or read_fraction(
                    verification_payload["forced_intersection_density"]
                )
                != expected_intersection
            ):
                raise RuntimeError(
                    f"p={outside_prime} resume artifact does not match"
                )
            records.append(existing)
            print(
                f"promoted={index}/{len(selected)} "
                f"outside={outside_prime} resumed=True",
                flush=True,
            )
            continue
        subprocess.run(
            generator_command(
                args.pool,
                design,
                outside_prime,
                certificate,
            ),
            check=True,
        )
        certificate_payload = json.loads(certificate.read_text())
        actual_intersection = read_fraction(
            certificate_payload["forced_intersection_density"]
        )
        if actual_intersection != expected_intersection:
            raise RuntimeError(
                f"p={outside_prime} regenerated a different intersection"
            )
        subprocess.run(
            [
                sys.executable,
                "verify_projected_conditional_fibre_overlap.py",
                str(args.pool),
                str(certificate),
                "--output",
                str(verification),
            ],
            check=True,
        )
        verification_payload = json.loads(verification.read_text())
        if (
            not verification_payload.get("verified")
            or read_fraction(
                verification_payload["forced_intersection_density"]
            )
            != expected_intersection
        ):
            raise RuntimeError(
                f"p={outside_prime} failed independent replay"
            )
        record = {
            "outside_prime": outside_prime,
            "discovery_file": str(discovery_path),
            "certificate": str(certificate),
            "verification": str(verification),
            "baseline": design["baseline"],
            "forced_intersection_density": fraction_payload(
                expected_intersection
            ),
            "improvement": design["improvement"],
            "target_combinations": int(
                verification_payload["target_combinations"]
            ),
            "verified": True,
        }
        records.append(record)
        args.output.write_text(
            json.dumps(
                manifest_payload(
                    pool=args.pool,
                    prefix=args.artifact_prefix,
                    discovery_paths=discovery_paths,
                    selected_count=len(selected),
                    records=records,
                    status="promotion_in_progress",
                ),
                indent=2,
            )
            + "\n"
        )
        print(
            f"promoted={index}/{len(selected)} outside={outside_prime} "
            f"improvement={read_fraction(design['improvement'])}",
            flush=True,
        )

    args.output.write_text(
        json.dumps(
            manifest_payload(
                pool=args.pool,
                prefix=args.artifact_prefix,
                discovery_paths=discovery_paths,
                selected_count=len(selected),
                records=records,
                status="all_selected_designs_independently_replayed",
            ),
            indent=2,
        )
        + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
