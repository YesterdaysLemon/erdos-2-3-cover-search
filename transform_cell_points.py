#!/usr/bin/env python3
"""Map original lattice points into a conditioned cell's canonical coordinates."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("conditioned_pool", type=Path)
    parser.add_argument("points", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.conditioned_pool.read_text())
    cell = payload["cell_condition"]
    period = int(cell["period"])
    k0 = int(cell["k_residue"])
    l0 = int(cell["l_residue"])
    shift_x = int(cell["coordinate_shift"]["x"])
    shift_y = int(cell["coordinate_shift"]["y"])
    algebraic_primes = tuple(int(value) for value in payload["algebraic_primes"])

    transformed = []
    seen = set()
    inspected = 0
    in_cell = 0
    for path in args.points:
        for raw_k, raw_l in json.loads(path.read_text()):
            inspected += 1
            k, l = int(raw_k), int(raw_l)
            if k % period != k0 or l % period != l0:
                continue
            in_cell += 1
            x = (k - k0) // period - shift_x
            y = (l - l0) // period - shift_y
            if any(x % prime == 0 and y % prime == 0 for prime in algebraic_primes):
                raise AssertionError("canonical point lies in an algebraic cell")
            point = (x, y)
            if point not in seen:
                seen.add(point)
                transformed.append(point)

    args.output.write_text(json.dumps(transformed) + "\n")
    print(
        f"inspected={inspected} in_cell={in_cell} "
        f"deduplicated={len(transformed)} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
