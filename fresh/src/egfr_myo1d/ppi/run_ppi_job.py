"""Dry-run entry point for M2.2 PyRosetta PPI jobs.

This module validates job-manifest wiring and writes a dry-run status record.
It refuses real execution in M2.2.
"""

import argparse
import os
import subprocess
import sys


def _candidate_python_commands():
    commands = []
    env_python = os.environ.get("EGFR_MYO1D_PYTHON")
    if env_python:
        commands.append([env_python])

    path_value = os.environ.get("PATH", "")
    for directory in path_value.split(os.pathsep):
        pytest_exe = os.path.join(directory, "pytest.exe")
        if os.path.exists(pytest_exe):
            candidate = os.path.abspath(os.path.join(directory, os.pardir, "python.exe"))
            if os.path.exists(candidate):
                commands.append([candidate])

    commands.extend([["py", "-3"], ["python3"], ["python"]])
    return commands


def _candidate_supports_task(command):
    if command[0] == sys.executable:
        return False
    try:
        proc = subprocess.Popen(
            command
            + [
                "-c",
                "import sys; raise SystemExit(0 if sys.version_info[:2] >= (3, 9) else 1)",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        proc.communicate()
        return proc.returncode == 0
    except OSError:
        return False


def _bootstrap_python3_if_needed():
    if sys.version_info[:2] >= (3, 9):
        return
    for command in _candidate_python_commands():
        if _candidate_supports_task(command):
            raise SystemExit(
                subprocess.call(command + ["-m", "egfr_myo1d.ppi.run_ppi_job"] + sys.argv[1:])
            )
    sys.stderr.write("Python 3.9+ is required for M2.2 PyRosetta harness commands.\n")
    raise SystemExit(1)


_bootstrap_python3_if_needed()

from pathlib import Path

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import now_iso, write_json
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.ppi.pyrosetta_adapter import load_job_by_name


def _bool_text(value):
    text = value.lower().strip()
    if text == "true":
        return True
    if text == "false":
        return False
    raise argparse.ArgumentTypeError("expected true or false")


def _repo_or_abs(ctx, text):
    path = Path(text)
    if path.is_absolute():
        return path.resolve()
    return (ctx.repo_root / path).resolve()


def dry_run_job(ctx, job_name):
    job = load_job_by_name(ctx, job_name)
    output_dir = ctx.require_within_run_dir(_repo_or_abs(ctx, job["output_dir"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    status_path = output_dir / "dry_run_status.json"
    payload = {
        "run_id": ctx.run_id,
        "job_name": job_name,
        "timestamp": now_iso(),
        "status": "DRY_RUN_PASS",
        "execution_allowed": False,
        "pyrosetta_imported": False,
        "docking_executed": False,
        "relaxation_executed": False,
        "validated_inputs": {
            "input_spec_json": job["input_spec_json"],
            "receptor_pdb": job["receptor_pdb"],
            "partner_pdb": job["partner_pdb"],
            "mapping_csv": job["mapping_csv"],
        },
        "notes": "M2.2 dry-run only; no scientific engine executed.",
    }
    write_json(status_path, payload, ctx)
    append_job_status(
        ctx,
        job_name,
        "DRY_RUN_PASS",
        stdout=job.get("stdout_path"),
        stderr=job.get("stderr_path"),
        details={
            "dry_run_status": ctx.relative_to_repo(status_path),
            "execution_allowed": False,
        },
    )
    append_phase_status(
        ctx,
        "m2.2-pyrosetta-job-dry-run",
        "PASS",
        "PyRosetta adapter job dry-run passed for {0}".format(job_name),
        {"job_name": job_name, "dry_run_status": ctx.relative_to_repo(status_path)},
    )
    return status_path


def build_parser():
    parser = argparse.ArgumentParser(prog="python -m egfr_myo1d.ppi.run_ppi_job")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--job-name", required=True)
    parser.add_argument(
        "--dry-run",
        type=_bool_text,
        default=True,
        help="M2.2 supports only true. Real execution is deferred to a later explicitly scoped phase.",
    )
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.dry_run:
        sys.stderr.write(
            "ERROR: M2.2 PyRosetta adapter refuses real execution; pass --dry-run true.\n"
        )
        return 2
    ctx = RunContext.for_existing(args.run_id)
    status_path = dry_run_job(ctx, args.job_name)
    print("dry_run_status={0}".format(status_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
