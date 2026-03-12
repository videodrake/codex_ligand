#!/usr/bin/env python3
"""Shared pre-qsub guard for production PBS submission."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def should_enforce_precheck(mode: str, skip_precheck_guard: str) -> bool:
    """Return True when production should require a precheck success marker."""
    return mode != "status" and skip_precheck_guard != "1"


def validate_precheck_status(
    workdir: Path,
    mode: str,
    skip_precheck_guard: str,
) -> tuple[bool, str]:
    """Validate the precheck marker and return (ok, message)."""
    status_file = workdir / "output" / "pre_qsub_status" / "last_pass.json"

    if not should_enforce_precheck(mode, skip_precheck_guard):
        return True, f"Precheck guard bypassed (mode={mode}, skip={skip_precheck_guard})."

    if status_file.exists():
        return True, f"Found pre-qsub success marker: {status_file}"

    lines = [
        "[ERROR] Missing pre-qsub success marker:",
        f"        {status_file}",
        "",
        "Run the validation lane first:",
        "  qsub config/run_pre_qsub_checks.pbs",
        "",
        "Best practice is to chain them with PBS dependency:",
        "  PRECHECK_JOB=$(qsub config/run_pre_qsub_checks.pbs)",
        "  qsub -W depend=afterok:${PRECHECK_JOB} config/run_production.pbs",
        "",
        "If you intentionally want to bypass this guard, set:",
        "  qsub -v SKIP_PRECHECK_GUARD=1 config/run_production.pbs",
    ]
    return False, "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate pre-qsub success marker.")
    parser.add_argument("--workdir", required=True, help="PBS working directory")
    parser.add_argument("--mode", default="auto", help="Production mode")
    parser.add_argument(
        "--skip-precheck-guard",
        default="0",
        help="Use 1 to bypass the precheck marker requirement",
    )
    args = parser.parse_args()

    ok, message = validate_precheck_status(
        Path(args.workdir),
        mode=args.mode,
        skip_precheck_guard=args.skip_precheck_guard,
    )
    print(message)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
