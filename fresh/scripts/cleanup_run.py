#!/usr/bin/env python3
"""M1 Phase 1: thin wrapper around `python -m egfr_myo1d.cli cleanup`.

Delegates all cleanup logic to ``egfr_myo1d.core.cleanup.run_cleanup`` via the
CLI subcommand. Sets PYTHONPATH so the package is importable from repo root
without requiring an editable install.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Cleanup wrapper. Forwards to `python -m egfr_myo1d.cli cleanup`. "
            "test mode deletes intermediates by default; production mode is "
            "dry-run by default."
        )
    )
    parser.add_argument(
        "--run-id",
        required=True,
        help="Existing run identifier under fresh/runs/.",
    )
    parser.add_argument(
        "--mode",
        required=True,
        choices=["test", "production"],
        help="test deletes intermediates; production defaults to dry-run.",
    )
    parser.add_argument(
        "--dry-run",
        choices=["true", "false"],
        default=None,
        help=(
            "Override default. Default is 'false' for --mode test, 'true' for "
            "--mode production."
        ),
    )
    parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
        help="Profile label recorded in cleanup_report.json.",
    )
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    repo_root = Path(__file__).resolve().parents[2]
    src_dir = repo_root / "fresh" / "src"

    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = (
        str(src_dir) + os.pathsep + existing if existing else str(src_dir)
    )

    cli_args = [
        sys.executable,
        "-m",
        "egfr_myo1d.cli",
        "cleanup",
        "--run-id",
        args.run_id,
        "--mode",
        args.mode,
        "--profile",
        args.profile,
    ]
    if args.dry_run is not None:
        cli_args.extend(["--dry-run", args.dry_run])

    proc = subprocess.run(cli_args, cwd=str(repo_root), env=env)
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
