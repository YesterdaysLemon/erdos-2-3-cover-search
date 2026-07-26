#!/usr/bin/env python3
"""Overlay one or more per-prime phase maps without changing other phases."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("base", type=Path)
    parser.add_argument("overlays", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    merged = {
        str(prime): int(target)
        for prime, target in json.loads(args.base.read_text()).items()
    }
    for overlay_path in args.overlays:
        overlay = json.loads(overlay_path.read_text())
        unknown = sorted(set(overlay) - set(merged), key=int)
        if unknown:
            raise RuntimeError(
                f"{overlay_path} has primes absent from base: {unknown[:5]}"
            )
        merged.update(
            {str(prime): int(target) for prime, target in overlay.items()}
        )
    args.output.write_text(json.dumps(merged) + "\n")
    print(
        f"base={len(json.loads(args.base.read_text()))} "
        f"overlays={len(args.overlays)} output={len(merged)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
