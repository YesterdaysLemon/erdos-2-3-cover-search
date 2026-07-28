#!/usr/bin/env python3
"""Certify a period block-star bound with stronger conditional fibre edges."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from certify_block_plus_overlap_forest import (
    read_fraction,
    recorded_anchor_rows,
)
from scan_all_ranked_pairanchor_star import (
    best_anchor_lower,
    fraction_payload,
    joint_index,
)
from search_mitm_pairwise_conditional_bounds import (
    exact_pairwise_subset_lower,
)


def select_star_anchors(
    outside: dict,
    anchors: list[dict],
    limit: int | None,
) -> list[dict]:
    """Choose a proof-safe compatible anchor subset for one star edge."""
    if limit is None:
        return anchors
    compatible = [
        anchor for anchor in anchors
        if joint_index(outside, anchor) == 1
    ]
    compatible.sort(
        key=lambda anchor: (int(anchor["h"]), int(anchor["p"]))
    )
    return compatible[:limit]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pool", type=Path)
    parser.add_argument("block_certificate", type=Path)
    parser.add_argument("--period", type=int, required=True)
    parser.add_argument(
        "--conditional-certificate",
        action="append",
        type=Path,
        default=[],
    )
    parser.add_argument(
        "--conditional-manifest",
        action="append",
        type=Path,
        default=[],
        help=(
            "promotion manifest whose replayed certificate records are "
            "added to the conditional star"
        ),
    )
    parser.add_argument("--fourway-triples", action="store_true")
    parser.add_argument("--fourterm-quads", action="store_true")
    parser.add_argument(
        "--star-anchor-limit",
        type=int,
        help=(
            "use at most this many smallest-modulus compatible anchors "
            "for each outside star edge"
        ),
    )
    parser.add_argument(
        "--star-witness-file",
        type=Path,
        help=(
            "MITM discovery output whose recorded anchor subsets are "
            "replayed exactly instead of optimizing over every subset"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.star_anchor_limit is not None and args.star_anchor_limit < 1:
        raise SystemExit("--star-anchor-limit must be positive")

    source = json.loads(args.pool.read_text())
    block = json.loads(args.block_certificate.read_text())
    rows = source["choices"]
    by_prime = {int(row["p"]): row for row in rows}
    recorded = recorded_anchor_rows(block)
    anchors = []
    for row in recorded:
        prime = int(row["p"])
        if prime not in by_prime or any(
            int(by_prime[prime][key]) != int(row[key])
            for key in ("h", "a", "b")
        ):
            raise RuntimeError("block anchor differs from the pool")
        anchors.append(by_prime[prime])
    anchor_primes = {int(row["p"]) for row in anchors}
    selected = [
        row for row in rows if args.period % int(row["h"]) == 0
    ]
    selected_primes = {int(row["p"]) for row in selected}
    if not anchor_primes <= selected_primes:
        raise RuntimeError("block does not lie in the requested family")

    conditional_paths = list(args.conditional_certificate)
    for manifest_path in args.conditional_manifest:
        manifest = json.loads(manifest_path.read_text())
        if (
            manifest.get("schema")
            != "promoted_projected_conditional_designs_v1"
            or manifest.get("status")
            != "all_selected_designs_independently_replayed"
            or int(manifest["completed_outside_count"])
            != int(manifest["selected_outside_count"])
        ):
            raise RuntimeError(
                f"conditional manifest is incomplete: {manifest_path}"
            )
        conditional_paths.extend(
            Path(record["certificate"])
            for record in manifest["records"]
        )

    conditionals = {}
    for path in conditional_paths:
        cert = json.loads(path.read_text())
        outside_prime = int(cert["outside_prime"])
        if outside_prime in conditionals:
            raise RuntimeError("repeated conditional outside prime")
        if (
            outside_prime in anchor_primes
            or outside_prime not in selected_primes
        ):
            raise RuntimeError(
                "conditional outside prime is not an outside family row"
            )
        conditional_anchor_primes = {
            int(prime) for prime in cert["anchor_primes"]
        }
        if not conditional_anchor_primes <= anchor_primes:
            raise RuntimeError(
                "conditional certificate uses a non-block anchor"
            )
        for row in cert["anchor_rows"]:
            prime = int(row["p"])
            if prime not in by_prime or any(
                int(by_prime[prime][key]) != int(row[key])
                for key in ("h", "a", "b")
            ):
                raise RuntimeError(
                    "conditional anchor differs from the pool"
                )
        outside_record = cert["outside_row"]
        if any(
            int(by_prime[outside_prime][key])
            != int(outside_record[key])
            for key in ("p", "h", "a", "b")
        ):
            raise RuntimeError(
                "conditional outside row differs from the pool"
            )
        conditionals[outside_prime] = (
            path,
            read_fraction(cert["forced_intersection_density"]),
        )

    star_witnesses = {}
    if args.star_witness_file is not None:
        witness_payload = json.loads(args.star_witness_file.read_text())
        witness_schema = witness_payload.get("schema")
        if witness_schema not in {
            "mitm_pairwise_conditional_search_v1",
            "pairwise_star_witnesses_v1",
        }:
            raise RuntimeError("unexpected star witness schema")
        witness_records = witness_payload[
            "records"
            if witness_schema == "pairwise_star_witnesses_v1"
            else "results"
        ]
        for record in witness_records:
            outside_prime = int(record["outside_prime"])
            if outside_prime in star_witnesses:
                raise RuntimeError("repeated outside prime in star witnesses")
            star_witnesses[outside_prime] = tuple(
                int(prime)
                for prime in record["selected_anchor_primes"]
            )
        expected_outside = selected_primes - anchor_primes
        if set(star_witnesses) != expected_outside:
            missing = expected_outside - star_witnesses.keys()
            extra = star_witnesses.keys() - expected_outside
            raise RuntimeError(
                "star witnesses do not match outside family rows: "
                f"missing={sorted(missing)} extra={sorted(extra)}"
            )

    total_density = sum(
        (Fraction(1, int(row["h"])) for row in selected),
        Fraction(0),
    )
    row_lowers = []
    star_loss = Fraction(0)
    for row in selected:
        prime = int(row["p"])
        if prime in anchor_primes:
            continue
        witness_primes = None
        if star_witnesses:
            witness_primes = star_witnesses[prime]
            if (
                len(set(witness_primes)) != len(witness_primes)
                or not set(witness_primes) <= anchor_primes
            ):
                raise RuntimeError(
                    f"invalid star witness anchors for p={prime}"
                )
            edge_anchors = tuple(
                by_prime[anchor_prime]
                for anchor_prime in witness_primes
            )
            if any(
                joint_index(row, anchor) != 1
                for anchor in edge_anchors
            ):
                raise RuntimeError(
                    f"incompatible star witness anchor for p={prime}"
                )
            baseline = exact_pairwise_subset_lower(
                row,
                edge_anchors,
            )
        else:
            edge_anchors = select_star_anchors(
                row,
                anchors,
                args.star_anchor_limit,
            )
            baseline = best_anchor_lower(
                row,
                edge_anchors,
                fourway_triples=args.fourway_triples,
                fourterm_quads=args.fourterm_quads,
            )
        path = None
        enhanced = baseline
        if prime in conditionals:
            path, conditional_lower = conditionals[prime]
            enhanced = max(baseline, conditional_lower)
        star_loss += enhanced
        row_lowers.append(
            {
                "outside_prime": prime,
                "baseline_intersection_lower_bound":
                fraction_payload(baseline),
                "baseline_anchor_primes": (
                    list(witness_primes)
                    if witness_primes is not None
                    else None
                ),
                "conditional_certificate": (
                    str(path) if path is not None else None
                ),
                "used_intersection_lower_bound":
                fraction_payload(enhanced),
            }
        )
    block_loss = read_fraction(block["forced_overlap_loss"])
    upper = total_density - block_loss - star_loss
    strengthened_note = ""
    if args.fourway_triples:
        strengthened_note += (
            " Surjective outside-plus-three-anchor systems restore their "
            "exact positive third-order term."
        )
    if args.fourterm_quads:
        strengthened_note += (
            " Four-anchor subsystems also use a fourth-order Bonferroni "
            "bound with independently checkable target-map indices."
        )
    result = {
        "schema": "ranked_period_conditional_star_v1",
        "pool": str(args.pool),
        "period": args.period,
        "row_count": len(selected),
        "block_certificate": str(args.block_certificate),
        "block_anchor_primes": sorted(anchor_primes),
        "block_overlap_loss": fraction_payload(block_loss),
        "fourway_triples": args.fourway_triples,
        "fourterm_quads": args.fourterm_quads,
        "star_anchor_limit": args.star_anchor_limit,
        "star_witness_file": (
            str(args.star_witness_file)
            if args.star_witness_file is not None
            else None
        ),
        "conditional_certificates": [
            str(path) for path in conditional_paths
        ],
        "conditional_manifests": [
            str(path) for path in args.conditional_manifest
        ],
        "outside_row_lowers": row_lowers,
        "star_overlap_sum": fraction_payload(star_loss),
        "total_density": fraction_payload(total_density),
        "union_density_upper_bound": fraction_payload(upper),
        "proved_no_cover": upper < 1,
        "argument": (
            "The exact anchor block is one node. Each outside fibre is a "
            "leaf of a star and contributes a certified lower bound for its "
            "intersection with the block union. Conditional certificates may "
            "replace the generic Bonferroni edge by an exact stronger edge. "
            "Subtracting every star edge is pointwise valid."
            + (
                " Each generic star edge records an explicit anchor subset "
                "whose exact pairwise Bonferroni value is sufficient; no "
                "optimality claim about the discovery search is used."
                if star_witnesses
                else ""
            )
            + strengthened_note
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"period={args.period} rows={len(selected)} anchors={len(anchors)} "
        f"conditional_edges={len(conditionals)} block_loss={block_loss} "
        f"star_loss={star_loss} upper={upper} "
        f"proved_no_cover={upper < 1}",
        flush=True,
    )
    return 0 if upper < 1 else 1


if __name__ == "__main__":
    raise SystemExit(main())
