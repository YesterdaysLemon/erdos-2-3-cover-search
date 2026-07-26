#!/usr/bin/env python3
"""Independently replay primary-SAT low-cell row cores with Z3."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import exact_uncovered


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()

    dep_path = Path(os.environ.get("TEMP", ".")) / "erdos203-pydeps"
    sys.path.insert(0, str(dep_path))
    import z3  # type: ignore

    certificate = json.loads(args.certificate.read_text())
    source = json.loads(Path(certificate["candidate_pool"]).read_text())
    phases = {
        int(prime): int(target)
        for prime, target in json.loads(
            Path(certificate["phase_file"]).read_text()
        ).items()
    }
    all_rows = source["choices"]
    algebraic_primes = tuple(
        int(prime) for prime in source.get("algebraic_primes", ())
    )
    records = certificate["checker"]["core_coordinate_cells"]
    for record_no, record in enumerate(records, 1):
        rows = []
        maximal = {2: 4, 3: 3}
        for index in record["row_indices"]:
            raw = all_rows[int(index)]
            row = {
                key: int(raw[key]) for key in ("h", "p", "a", "b")
            }
            row["c"] = phases[row["p"]] % row["h"]
            row["factors"] = exact_uncovered.factor(row["h"])
            rows.append(row)
            for prime, exponent in row["factors"].items():
                maximal[prime] = max(maximal.get(prime, 0), exponent)
        for prime in algebraic_primes:
            maximal[prime] = max(maximal.get(prime, 0), 1)
        components = {
            prime: prime**exponent for prime, exponent in maximal.items()
        }
        solver = z3.Solver()
        kval = {prime: z3.Int(f"k_{record_no}_{prime}") for prime in components}
        lval = {prime: z3.Int(f"l_{record_no}_{prime}") for prime in components}
        for prime, modulus in components.items():
            solver.add(kval[prime] >= 0, kval[prime] < modulus)
            solver.add(lval[prime] >= 0, lval[prime] < modulus)
        for row in rows:
            failures = []
            for prime, exponent in row["factors"].items():
                modulus = prime**exponent
                failures.append(
                    (
                        row["a"] * kval[prime]
                        + row["b"] * lval[prime]
                        - row["c"]
                    )
                    % modulus
                    != 0
                )
            solver.add(z3.Or(*failures))
        for prime in algebraic_primes:
            solver.add(
                z3.Or(
                    kval[prime] % prime != 0,
                    lval[prime] % prime != 0,
                )
            )
        if source.get("sophie_germain", False):
            solver.add(z3.Or(kval[2] % 4 != 2, lval[2] % 4 != 0))
        for modulus, kr, lr in record["cell"]:
            factors = exact_uncovered.factor(int(modulus))
            if len(factors) != 1:
                raise RuntimeError(f"non-prime-power cell modulus {modulus}")
            prime = next(iter(factors))
            solver.add(kval[prime] % int(modulus) == int(kr))
            solver.add(lval[prime] % int(modulus) == int(lr))
        status = solver.check()
        if status != z3.unsat:
            raise RuntimeError(
                f"cell core {record_no} failed independent replay: {status}"
            )
        print(
            f"cell={record_no}/{len(records)} rows={len(rows)} UNSAT",
            flush=True,
        )
    print(
        f"PASS cells={len(records)} engine=z3 certificate={args.certificate}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
