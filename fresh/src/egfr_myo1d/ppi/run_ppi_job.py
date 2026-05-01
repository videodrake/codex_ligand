"""Entry point for M2.2 PyRosetta PPI jobs.

This module validates job-manifest wiring and writes a dry-run status record.
Real execution is allowed only when explicitly requested with
``--dry-run false --allow-real true``. The scientific chain convention for real
jobs is EGFR dimer chains A/B docked against MYO1D chain C (Rosetta partners
``AB_C``).
"""

import argparse
import csv
import importlib
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
from egfr_myo1d.core.run_context import RunContext, ensure_within
from egfr_myo1d.ppi.pyrosetta_adapter import load_job_by_name


REAL_SCORE_FIELDS = [
    "run_id",
    "job_name",
    "state_id",
    "seed",
    "model_index",
    "pose_id",
    "score",
    "pose_path",
    "status",
]


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


def _repo_path(ctx, path):
    return ctx.relative_to_repo(Path(path))


def _renumber_atom_line(line, serial):
    text = line.rstrip("\n")
    if len(text) < 80:
        text = text.ljust(80)
    return "{0}{1:5d}{2}".format(text[:6], serial, text[11:])


def _set_chain_id(line, chain_id):
    text = line.rstrip("\n")
    if len(text) < 22:
        text = text.ljust(22)
    return "{0}{1}{2}".format(text[:21], chain_id, text[22:])


def _copy_atoms_with_chain(source_path, target_handle, serial_start, chain_id=None):
    serial = serial_start
    atom_count = 0
    chains = set()
    with Path(source_path).open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            output = _renumber_atom_line(line, serial)
            if chain_id is not None:
                output = _set_chain_id(output, chain_id)
            chains.add((output[21:22].strip() or "_"))
            target_handle.write(output.rstrip("\n") + "\n")
            serial += 1
            atom_count += 1
    return serial, atom_count, sorted(chains)


def build_ab_c_input(ctx, job, output_dir):
    """Combine EGFR A/B receptor and MYO1D partner into a run-local AB_C PDB."""

    receptor_path = ensure_within(_repo_or_abs(ctx, job["receptor_pdb"]), ctx.run_dir)
    partner_path = ensure_within(_repo_or_abs(ctx, job["partner_pdb"]), ctx.run_dir)
    input_dir = ctx.require_within_run_dir(output_dir / "inputs")
    input_dir.mkdir(parents=True, exist_ok=True)
    target = ctx.require_within_run_dir(
        input_dir / "{0}_AB_C_input.pdb".format(job["job_name"])
    )
    serial = 1
    with target.open("w", encoding="utf-8") as handle:
        serial, receptor_atom_count, receptor_chains = _copy_atoms_with_chain(
            receptor_path, handle, serial
        )
        handle.write("TER\n")
        serial, partner_atom_count, partner_chains = _copy_atoms_with_chain(
            partner_path, handle, serial, chain_id="C"
        )
        handle.write("TER\nEND\n")
    return target, {
        "receptor_pdb": _repo_path(ctx, receptor_path),
        "partner_pdb": _repo_path(ctx, partner_path),
        "combined_input_pdb": _repo_path(ctx, target),
        "receptor_atom_count": receptor_atom_count,
        "partner_atom_count": partner_atom_count,
        "receptor_chains": receptor_chains,
        "partner_input_chains": partner_chains,
        "partner_output_chain": "C",
    }


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


def _write_score_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REAL_SCORE_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _import_pyrosetta():
    return importlib.import_module("pyrosetta")


def real_run_job(
    ctx,
    job_name,
    *,
    allow_real=False,
    max_models=None,
    dock_partners="AB_C",
):
    """Execute one explicit PyRosetta AB_C docking job and dump pose PDBs."""

    if not allow_real:
        raise PermissionError("real execution requires --allow-real true")
    job = load_job_by_name(ctx, job_name)
    output_dir = ctx.require_within_run_dir(_repo_or_abs(ctx, job["output_dir"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    poses_dir = ctx.require_within_run_dir(output_dir / "poses")
    poses_dir.mkdir(parents=True, exist_ok=True)
    combined_input, input_summary = build_ab_c_input(ctx, job, output_dir)

    seed = int(job.get("seed", 0))
    planned_models = int(job.get("models_per_seed", 1))
    model_count = min(planned_models, int(max_models)) if max_models else planned_models
    if model_count <= 0:
        raise ValueError("model count must be positive")

    pyrosetta = _import_pyrosetta()
    pyrosetta.init("-mute all -constant_seed -jran {0}".format(seed))
    rosetta = pyrosetta.rosetta
    movable_jumps = rosetta.utility.vector1_int()
    movable_jumps.append(1)
    scorefxn = pyrosetta.create_score_function("ref2015")
    base_pose = pyrosetta.pose_from_pdb(str(combined_input))

    rows = []
    failures = []
    for model_index in range(1, model_count + 1):
        pose_id = "{0}_model_{1:05d}".format(job_name, model_index)
        pose_path = ctx.require_within_run_dir(poses_dir / "{0}.pdb".format(pose_id))
        try:
            pose = base_pose.clone()
            rosetta.protocols.docking.setup_foldtree(pose, dock_partners, movable_jumps)
            perturb = rosetta.protocols.rigid.RigidBodyPerturbMover(1, 360.0, 100.0)
            perturb.apply(pose)
            slide = rosetta.protocols.docking.DockingSlideIntoContact(1)
            slide.apply(pose)
            docking = rosetta.protocols.docking.DockMCMProtocol()
            docking.set_scorefxn(scorefxn)
            docking.apply(pose)
            score = float(scorefxn(pose))
            pose.dump_pdb(str(pose_path))
            status = "PASS"
        except Exception as exc:  # pragma: no cover - exercised on HPC with PyRosetta
            score = ""
            status = "FAIL"
            failures.append("model_{0:05d}:{1}".format(model_index, exc))
        rows.append(
            {
                "run_id": ctx.run_id,
                "job_name": job_name,
                "state_id": str(job.get("state_id", "")),
                "seed": str(seed),
                "model_index": model_index,
                "pose_id": pose_id,
                "score": "{0:.6f}".format(score) if isinstance(score, float) else "",
                "pose_path": _repo_path(ctx, pose_path) if pose_path.is_file() else "",
                "status": status,
            }
        )

    score_csv = ctx.require_within_run_dir(output_dir / "pose_scores.csv")
    _write_score_csv(score_csv, rows)
    status = "REAL_RUN_PASS" if not failures else "REAL_RUN_FAIL"
    status_path = output_dir / "real_run_status.json"
    payload = {
        "run_id": ctx.run_id,
        "job_name": job_name,
        "timestamp": now_iso(),
        "status": status,
        "execution_allowed": True,
        "pyrosetta_imported": True,
        "docking_executed": not failures,
        "relaxation_executed": False,
        "dock_partners": dock_partners,
        "planned_models": planned_models,
        "executed_models": model_count,
        "successful_models": sum(1 for row in rows if row["status"] == "PASS"),
        "failed_models": len(failures),
        "input_summary": input_summary,
        "outputs": {
            "poses_dir": _repo_path(ctx, poses_dir),
            "score_csv": _repo_path(ctx, score_csv),
        },
        "failures": failures,
    }
    write_json(status_path, payload, ctx)
    append_job_status(
        ctx,
        job_name,
        status,
        stdout=job.get("stdout_path"),
        stderr=job.get("stderr_path"),
        details={
            "real_run_status": _repo_path(ctx, status_path),
            "pose_score_csv": _repo_path(ctx, score_csv),
            "docking_executed": not failures,
        },
    )
    append_phase_status(
        ctx,
        "m2.2-pyrosetta-job-real-run",
        "PASS" if not failures else "FAIL",
        "PyRosetta AB_C real job {0} for {1}".format(status, job_name),
        {
            "job_name": job_name,
            "real_run_status": _repo_path(ctx, status_path),
            "dock_partners": dock_partners,
        },
    )
    if failures:
        raise RuntimeError("; ".join(failures[:3]))
    return status_path


def build_parser():
    parser = argparse.ArgumentParser(prog="python -m egfr_myo1d.ppi.run_ppi_job")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--job-name", required=True)
    parser.add_argument(
        "--dry-run",
        type=_bool_text,
        default=True,
        help="Use true for wiring checks. Use false only with --allow-real true.",
    )
    parser.add_argument(
        "--allow-real",
        type=_bool_text,
        default=False,
        help="Required with --dry-run false to execute PyRosetta AB_C docking.",
    )
    parser.add_argument(
        "--max-models",
        type=int,
        default=None,
        help="Optional cap on models for this one job. Default uses the manifest models_per_seed.",
    )
    parser.add_argument(
        "--dock-partners",
        default="AB_C",
        help="Rosetta docking partners string. Default docks EGFR dimer A/B against MYO1D C.",
    )
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.dry_run and not args.allow_real:
        sys.stderr.write(
            "ERROR: M2.2 PyRosetta adapter refuses real execution unless --allow-real true is supplied.\n"
        )
        return 2
    ctx = RunContext.for_existing(args.run_id)
    if args.dry_run:
        status_path = dry_run_job(ctx, args.job_name)
        print("dry_run_status={0}".format(status_path))
    else:
        status_path = real_run_job(
            ctx,
            args.job_name,
            allow_real=args.allow_real,
            max_models=args.max_models,
            dock_partners=args.dock_partners,
        )
        print("real_run_status={0}".format(status_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
