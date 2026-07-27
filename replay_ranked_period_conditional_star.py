#!/usr/bin/env python3
"""Replay an assembled conditional-star certificate and its named reports."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


def verification_path(certificate: str) -> str:
    suffix = "_certificate.json"
    if not certificate.endswith(suffix):
        raise RuntimeError(
            f"conditional certificate has an unexpected name: {certificate}"
        )
    return certificate.removesuffix(suffix) + "_verification.json"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    pool = str(certificate["pool"])
    block = str(certificate["block_certificate"])
    block_verification = verification_path(block)
    conditional_verifications = [
        verification_path(str(path))
        for path in certificate["conditional_certificates"]
    ]

    command = [
        sys.executable,
        str(Path(__file__).with_name(
            "verify_ranked_period_conditional_star.py"
        )),
        pool,
        str(args.certificate),
        block_verification,
    ]
    for verification in conditional_verifications:
        command.extend(["--conditional-verification", verification])
    command.extend(["--output", str(args.output)])
    return subprocess.run(command, check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
