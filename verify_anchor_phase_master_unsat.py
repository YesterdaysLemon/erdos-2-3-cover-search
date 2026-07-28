#!/usr/bin/env python3
"""Independent Z3 replay of an anchor-phase finite UNSAT certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path


EXPECTED_COLUMNS = [
    "p",
    "h",
    "a",
    "b",
    "base_phase",
    "target_modulus",
    "target_residue",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_certificate(certificate: dict) -> dict:
    if certificate.get("schema_version") != 1:
        raise ValueError("unsupported certificate schema")
    if certificate.get("row_columns") != EXPECTED_COLUMNS:
        raise ValueError("unexpected compact-row schema")

    rows = []
    row_by_prime = {}
    for raw_row in certificate["rows"]:
        if not isinstance(raw_row, list) or len(raw_row) != len(EXPECTED_COLUMNS):
            raise ValueError("malformed compact row")
        prime, h, a, b, phase, modulus, residue = map(int, raw_row)
        if prime in row_by_prime:
            raise ValueError(f"duplicate row prime {prime}")
        if h < 1 or modulus < 1 or h % modulus:
            raise ValueError(f"invalid row restriction for p={prime}")
        a %= h
        b %= h
        phase %= h
        residue %= modulus
        if math.gcd(math.gcd(a, b), h) != 1:
            raise ValueError(f"nonsurjective affine row for p={prime}")
        if phase % modulus != residue:
            raise ValueError(f"forbidden base phase for p={prime}")
        row = (prime, h, a, b, phase, modulus, residue)
        rows.append(row)
        row_by_prime[prime] = row

    anchors = certificate["anchors"]
    anchor_primes = tuple(int(anchor["p"]) for anchor in anchors)
    if not anchor_primes or len(set(anchor_primes)) != len(anchor_primes):
        raise ValueError("anchor primes must be nonempty and distinct")
    if set(anchor_primes) - row_by_prime.keys():
        raise ValueError("an anchor row is absent")
    legal_spaces = {}
    for anchor, prime in zip(anchors, anchor_primes):
        _p, h, a, b, _phase, modulus, residue = row_by_prime[prime]
        legal = tuple(
            target for target in range(h) if target % modulus == residue
        )
        expected = {
            "p": prime,
            "h": h,
            "a": a,
            "b": b,
            "target_modulus": modulus,
            "target_residue": residue,
            "legal_target_count": len(legal),
        }
        if anchor != expected:
            raise ValueError(f"anchor metadata mismatch for p={prime}")
        legal_spaces[prime] = legal

    raw_points = certificate["points"]
    points = []
    for raw_point in raw_points:
        if not isinstance(raw_point, list) or len(raw_point) != 2:
            raise ValueError("malformed core point")
        point = tuple(map(int, raw_point))
        if point[0] < 0 or point[1] < 0:
            raise ValueError("core points must be nonnegative")
        points.append(point)
    if len(set(points)) != len(points):
        raise ValueError("core point list contains duplicates")

    algebraic_primes = tuple(
        int(prime) for prime in certificate.get("algebraic_primes", ())
    )
    if len(set(algebraic_primes)) != len(algebraic_primes):
        raise ValueError("duplicate algebraic prime")
    anchor_set = set(anchor_primes)
    frozen_rows = [row for row in rows if row[0] not in anchor_set]
    for k, l in points:
        if any(
            k % prime == 0 and l % prime == 0
            for prime in algebraic_primes
        ):
            raise ValueError("core point is algebraically covered")
        for prime, h, a, b, phase, _modulus, _residue in frozen_rows:
            if (a * k + b * l - phase) % h == 0:
                raise ValueError(
                    f"core point is covered by frozen row p={prime}"
                )

    dependency_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dependency_path))
    import z3  # type: ignore

    solver = z3.Solver()
    phase_variables = {
        prime: z3.Int(f"phase_{prime}") for prime in anchor_primes
    }
    for prime in anchor_primes:
        solver.add(
            z3.Or(
                [
                    phase_variables[prime] == target
                    for target in legal_spaces[prime]
                ]
            )
        )
    for k, l in points:
        covering_targets = []
        for prime in anchor_primes:
            _p, h, a, b, _phase, modulus, residue = row_by_prime[prime]
            target = (a * k + b * l) % h
            if target % modulus == residue:
                covering_targets.append(
                    phase_variables[prime] == target
                )
        if not covering_targets:
            raise ValueError("core point has no legal anchor cover")
        solver.add(z3.Or(covering_targets))
    solve_started = time.monotonic()
    verdict = solver.check()
    solve_seconds = time.monotonic() - solve_started
    if verdict != z3.unsat:
        raise ValueError(f"independent anchor master returned {verdict}")

    target_space_sizes = [
        len(legal_spaces[prime]) for prime in anchor_primes
    ]
    summary = certificate["summary"]
    expected_summary = {
        "rows": len(rows),
        "anchor_rows": len(anchor_primes),
        "frozen_rows": len(frozen_rows),
        "core_points": len(points),
        "legal_target_space_sizes": target_space_sizes,
        "joint_legal_target_assignments": math.prod(target_space_sizes),
        "pysat_unsat": True,
    }
    if summary != expected_summary:
        raise ValueError("certificate summary mismatch")
    master = certificate["master"]
    if (
        master.get("sat") is not False
        or int(master.get("eligible_points", -1)) != len(points)
        or master.get("empty_clause_point") is not None
    ):
        raise ValueError("embedded PySAT result is inconsistent")

    return {
        "verified": True,
        "rows": len(rows),
        "anchor_rows": len(anchor_primes),
        "frozen_rows": len(frozen_rows),
        "core_points": len(points),
        "joint_legal_target_assignments": math.prod(target_space_sizes),
        "z3_verdict": str(verdict),
        "z3_assertions": len(solver.assertions()),
        "z3_solve_seconds": solve_seconds,
        "engine": "independent-z3-integer-anchor-master",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    result = verify_certificate(certificate)
    result["certificate"] = args.certificate.name
    result["certificate_sha256"] = sha256_file(args.certificate)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"VERIFIED_UNSAT anchors={result['anchor_rows']} "
        f"points={result['core_points']} "
        f"assignments={result['joint_legal_target_assignments']} "
        f"output={args.output}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
