#!/usr/bin/env python3
"""M1 Phase 6: thin wrapper around `python -m egfr_myo1d.cli prepare-pbs`.

Forwards arguments to the CLI subcommand. Sets PYTHONPATH so the package
is importable from repo root without an editable install.

Does NOT call qsub. Submission is a user-side HPC step.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Generate a concrete PBS job file under "
            "fresh/runs/<run_id>/scripts/. Forwards to "
            "`python -m egfr_myo1d.cli prepare-pbs`."
        )
    )
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--job-name", required=True)
    parser.add_argument(
        "--mode",
        required=True,
        choices=["smoke_env", "smoke_input", "mini", "scaling", "production"],
    )
    parser.add_argument(
        "--node",
        default=None,
        choices=["node04", "node05", "node06"],
    )
    parser.add_argument("--ppn", type=int, default=None)
    parser.add_argument("--walltime", default=None)
    parser.add_argument("--output-path", default=None)
    parser.add_argument("--input-root", default="fresh/data/raw")
    parser.add_argument(
        "--profile",
        default="codex_dev",
        choices=["codex_dev", "hpc_strict"],
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
        "prepare-pbs",
        "--run-id",
        args.run_id,
        "--job-name",
        args.job_name,
        "--mode",
        args.mode,
        "--profile",
        args.profile,
        "--input-root",
        args.input_root,
    ]
    if args.node is not None:
        cli_args.extend(["--node", args.node])
    if args.ppn is not None:
        cli_args.extend(["--ppn", str(args.ppn)])
    if args.walltime is not None:
        cli_args.extend(["--walltime", args.walltime])
    if args.output_path is not None:
        cli_args.extend(["--output-path", args.output_path])

    proc = subprocess.run(cli_args, cwd=str(repo_root), env=env)
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
