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
from collections import Counter, defaultdict


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
STANDARD_AA3 = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}
ROSETTA_RESNAME_NORMALIZE = {
    "HSD": "HIS",
    "HSE": "HIS",
    "HSP": "HIS",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "CYX": "CYS",
}
ROSETTA_ATOM_NAME_NORMALIZE = {
    ("ILE", "CD"): "CD1",
}
ROSETTA_ALTLOC_KEEP = {"", "A"}
BACKBONE_ATOMS = {"N", "CA", "C", "O"}


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


def _set_resname(line, resname):
    text = line.rstrip("\n")
    if len(text) < 20:
        text = text.ljust(20)
    return "{0}{1:>3}{2}".format(text[:17], resname, text[20:])


def _set_atom_name(line, atom_name):
    text = line.rstrip("\n")
    if len(text) < 16:
        text = text.ljust(16)
    return "{0}{1:>4}{2}".format(text[:12], atom_name, text[16:])


def _clear_altloc(line):
    text = line.rstrip("\n")
    if len(text) < 17:
        text = text.ljust(17)
    return "{0} {1}".format(text[:16], text[17:])


def _get_atom_name(line):
    return line[12:16].strip()


def _get_resname(line):
    return line[17:20].strip()


def _get_altloc(line):
    return line[16:17].strip()


def _get_occupancy(line):
    value = line[54:60].strip()
    if not value:
        return 1.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def _get_element(line):
    return line[76:78].strip().upper()


def _is_hydrogen_atom(line):
    element = _get_element(line)
    if element:
        return element == "H"
    atom_name = _get_atom_name(line).lstrip("0123456789").upper()
    return atom_name.startswith("H")


def _atom_identity(line):
    return (
        line[21:22].strip() or "_",
        line[22:26].strip(),
        line[26:27].strip(),
        line[12:16].strip(),
    )


def _residue_identity(line):
    return (
        line[21:22].strip() or "_",
        line[22:26].strip(),
        line[26:27].strip(),
    )


def _new_sanitize_stats():
    return {
        "atom_records_in": 0,
        "atoms_out": 0,
        "hetatm_removed": 0,
        "hydrogens_removed": 0,
        "zero_occupancy_removed": 0,
        "alternate_location_removed": 0,
        "alternate_locations_collapsed": 0,
        "nonstandard_residue_removed": 0,
        "duplicate_atoms_removed": 0,
        "incomplete_residues_removed": 0,
        "incomplete_residue_atoms_removed": 0,
        "residue_renames": Counter(),
        "atom_renames": Counter(),
    }


def _json_sanitize_stats(stats):
    return {
        key: value
        for key, value in stats.items()
        if key not in {"residue_renames", "atom_renames"}
    } | {
        "residue_renames": dict(sorted(stats["residue_renames"].items())),
        "atom_renames": dict(sorted(stats["atom_renames"].items())),
    }


def _normalize_rosetta_atom_line(line, stats, chain_id=None):
    if line.startswith("HETATM"):
        stats["hetatm_removed"] += 1
        return None
    if not line.startswith("ATOM  "):
        return None

    stats["atom_records_in"] += 1
    text = line.rstrip("\n")
    if len(text) < 80:
        text = text.ljust(80)

    altloc = _get_altloc(text)
    if altloc and altloc not in ROSETTA_ALTLOC_KEEP:
        stats["alternate_location_removed"] += 1
        return None
    if altloc:
        stats["alternate_locations_collapsed"] += 1
        text = _clear_altloc(text)

    if _get_occupancy(text) <= 0.0:
        stats["zero_occupancy_removed"] += 1
        return None

    if _is_hydrogen_atom(text):
        stats["hydrogens_removed"] += 1
        return None

    resname = _get_resname(text)
    normalized_resname = ROSETTA_RESNAME_NORMALIZE.get(resname, resname)
    if normalized_resname != resname:
        stats["residue_renames"]["{0}->{1}".format(resname, normalized_resname)] += 1
        text = _set_resname(text, normalized_resname)
        resname = normalized_resname
    if resname not in STANDARD_AA3:
        stats["nonstandard_residue_removed"] += 1
        return None

    atom_name = _get_atom_name(text)
    normalized_atom_name = ROSETTA_ATOM_NAME_NORMALIZE.get((resname, atom_name), atom_name)
    if normalized_atom_name != atom_name:
        stats["atom_renames"]["{0}:{1}->{2}".format(resname, atom_name, normalized_atom_name)] += 1
        text = _set_atom_name(text, normalized_atom_name)

    if chain_id is not None:
        text = _set_chain_id(text, chain_id)
    return text


def _strip_backbone_incomplete_residues(lines, stats):
    residue_atoms = defaultdict(set)
    for line in lines:
        residue_atoms[_residue_identity(line)].add(_get_atom_name(line))

    incomplete = {
        residue
        for residue, atoms in residue_atoms.items()
        if not BACKBONE_ATOMS.issubset(atoms)
    }
    if not incomplete:
        return lines

    cleaned = []
    for line in lines:
        if _residue_identity(line) in incomplete:
            stats["incomplete_residue_atoms_removed"] += 1
            continue
        cleaned.append(line)
    stats["incomplete_residues_removed"] += len(incomplete)
    return cleaned


def _sanitized_rosetta_atom_lines(source_path, chain_id=None):
    stats = _new_sanitize_stats()
    lines = []
    seen_atom_identities = set()
    with Path(source_path).open("r", encoding="utf-8") as handle:
        for line in handle:
            output = _normalize_rosetta_atom_line(line, stats, chain_id=chain_id)
            if output is None:
                continue
            atom_identity = _atom_identity(output)
            if atom_identity in seen_atom_identities:
                stats["duplicate_atoms_removed"] += 1
                continue
            seen_atom_identities.add(atom_identity)
            lines.append(output)
    lines = _strip_backbone_incomplete_residues(lines, stats)
    stats["atoms_out"] = len(lines)
    return lines, stats


def _copy_atoms_with_chain(source_path, target_handle, serial_start, chain_id=None):
    serial = serial_start
    atom_count = 0
    chains = set()
    current_chain = None
    lines, stats = _sanitized_rosetta_atom_lines(source_path, chain_id=chain_id)
    for line in lines:
        output_chain = line[21:22].strip() or "_"
        if chain_id is None and current_chain is not None and output_chain != current_chain:
            target_handle.write("TER\n")
        current_chain = output_chain
        output = _renumber_atom_line(line, serial)
        chains.add(output_chain)
        target_handle.write(output.rstrip("\n") + "\n")
        serial += 1
        atom_count += 1
    return serial, atom_count, sorted(chains), _json_sanitize_stats(stats)


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
        serial, receptor_atom_count, receptor_chains, receptor_sanitize = _copy_atoms_with_chain(
            receptor_path, handle, serial
        )
        handle.write("TER\n")
        serial, partner_atom_count, partner_chains, partner_sanitize = _copy_atoms_with_chain(
            partner_path, handle, serial, chain_id="C"
        )
        handle.write("TER\nEND\n")
    if receptor_atom_count == 0 or partner_atom_count == 0:
        raise ValueError(
            "Rosetta input sanitation removed all atoms for receptor or partner: "
            "receptor_atoms={0}, partner_atoms={1}".format(
                receptor_atom_count, partner_atom_count
            )
        )
    return target, {
        "receptor_pdb": _repo_path(ctx, receptor_path),
        "partner_pdb": _repo_path(ctx, partner_path),
        "combined_input_pdb": _repo_path(ctx, target),
        "receptor_atom_count": receptor_atom_count,
        "partner_atom_count": partner_atom_count,
        "receptor_chains": receptor_chains,
        "partner_input_chains": partner_chains,
        "partner_output_chain": "C",
        "receptor_rosetta_sanitize": receptor_sanitize,
        "partner_rosetta_sanitize": partner_sanitize,
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


def _real_status_payload(
    *,
    ctx,
    job,
    job_name,
    status,
    dock_partners,
    planned_models,
    model_count,
    input_summary,
    poses_dir,
    score_csv=None,
    failures=None,
):
    failures = failures or []
    return {
        "run_id": ctx.run_id,
        "job_name": job_name,
        "timestamp": now_iso(),
        "status": status,
        "execution_allowed": True,
        "pyrosetta_imported": True,
        "docking_executed": status == "REAL_RUN_PASS",
        "relaxation_executed": False,
        "dock_partners": dock_partners,
        "planned_models": planned_models,
        "executed_models": model_count,
        "successful_models": 0,
        "failed_models": len(failures),
        "input_summary": input_summary,
        "outputs": {
            "poses_dir": _repo_path(ctx, poses_dir),
            "score_csv": _repo_path(ctx, score_csv) if score_csv else "",
        },
        "failures": failures,
        "notes": "Real PyRosetta execution uses sanitized protein ATOM-only AB_C input.",
        "validated_inputs": {
            "receptor_pdb": job.get("receptor_pdb", ""),
            "partner_pdb": job.get("partner_pdb", ""),
            "mapping_csv": job.get("mapping_csv", ""),
        },
    }


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

    status_path = output_dir / "real_run_status.json"
    pyrosetta = _import_pyrosetta()
    pyrosetta.init("-mute all -constant_seed -jran {0}".format(seed))
    rosetta = pyrosetta.rosetta
    movable_jumps = rosetta.utility.vector1_int()
    movable_jumps.append(1)
    scorefxn = pyrosetta.create_score_function("ref2015")
    try:
        base_pose = pyrosetta.pose_from_pdb(str(combined_input))
    except Exception as exc:
        failure = "pose_import_failed:{0}".format(exc)
        payload = _real_status_payload(
            ctx=ctx,
            job=job,
            job_name=job_name,
            status="REAL_RUN_INPUT_IMPORT_FAIL",
            dock_partners=dock_partners,
            planned_models=planned_models,
            model_count=0,
            input_summary=input_summary,
            poses_dir=poses_dir,
            failures=[failure],
        )
        payload["docking_executed"] = False
        write_json(status_path, payload, ctx)
        append_job_status(
            ctx,
            job_name,
            "REAL_RUN_INPUT_IMPORT_FAIL",
            stdout=job.get("stdout_path"),
            stderr=job.get("stderr_path"),
            details={
                "real_run_status": _repo_path(ctx, status_path),
                "combined_input_pdb": input_summary["combined_input_pdb"],
                "docking_executed": False,
            },
        )
        raise

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
    payload = _real_status_payload(
        ctx=ctx,
        job=job,
        job_name=job_name,
        status=status,
        dock_partners=dock_partners,
        planned_models=planned_models,
        model_count=model_count,
        input_summary=input_summary,
        poses_dir=poses_dir,
        score_csv=score_csv,
        failures=failures,
    )
    payload["docking_executed"] = not failures
    payload["successful_models"] = sum(1 for row in rows if row["status"] == "PASS")
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
