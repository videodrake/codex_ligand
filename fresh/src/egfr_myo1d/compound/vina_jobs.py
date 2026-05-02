"""M3-T5 PBS generator and compound docking job manifest.

This module plans focused_pocket_first Vina jobs and renders PBS wrappers. It
does not submit qsub jobs, execute Vina, parse poses, score compounds, or make
candidate claims.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import shlex
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from egfr_myo1d.compound.confidentiality import (
    PUBLIC_COMPOUND_IDS,
    git_check_ignored,
    git_ls_files,
    load_private_map,
    public_output_scan_paths,
    scan_internal_id_leaks,
)
from egfr_myo1d.compound.vina_adapter import build_vina_argv, detect_vina
from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status, log_master
from egfr_myo1d.core.manifest import load_yaml_config
from egfr_myo1d.core.run_context import RunContext, ensure_within
from egfr_myo1d.hpc.pbs import (
    DEFAULT_CONDA_ENV,
    DEFAULT_CONDA_SH,
    THREAD_ENV_KEYS,
    render_pbs_content,
)


JOB_FIELDS = [
    "job_id",
    "run_id",
    "m2_run_id",
    "profile",
    "mode",
    "job_scope",
    "chunk_id",
    "node",
    "ppn",
    "worker_count",
    "compound_public_id",
    "ligand_pdbqt_file",
    "ligand_pdbqt_sha256",
    "ligand_preparation_status",
    "state_id",
    "state_role",
    "receptor_pdbqt_file",
    "receptor_pdbqt_sha256",
    "receptor_preparation_status",
    "pocket_family_id",
    "protomer_id",
    "box_id",
    "box_center_x",
    "box_center_y",
    "box_center_z",
    "box_size_x",
    "box_size_y",
    "box_size_z",
    "source_ligand_manifest",
    "source_receptor_manifest",
    "source_docking_box_manifest",
    "source_vina_smoke_manifest",
    "vina_executable",
    "vina_version",
    "vina_argv_json",
    "exhaustiveness",
    "num_modes",
    "cpu_per_vina",
    "vina_repeat_id",
    "seed",
    "timeout_seconds",
    "output_dir",
    "output_pdbqt_file",
    "vina_log_file",
    "stdout_file",
    "stderr_file",
    "pbs_file",
    "pbs_job_name",
    "command_manifest_written_before_execution",
    "planned_at",
    "submission_status",
    "docking_status",
    "allowed_for_execution",
    "allowed_for_collection",
    "job_notes",
]

PBS_FIELDS = [
    "pbs_job_name",
    "run_id",
    "m2_run_id",
    "profile",
    "chunk_id",
    "node",
    "ppn",
    "worker_count",
    "assigned_job_count",
    "pbs_file",
    "pbs_stdout_file",
    "pbs_stderr_file",
    "queue",
    "walltime",
    "conda_env",
    "python_executable",
    "runner_command",
    "qsub_command",
    "submit_by_default",
    "pbs_generation_status",
    "pbs_notes",
]

QC_FIELDS = ["check_id", "category", "status", "severity", "details", "recommended_fix"]

REQUIRED_CHECKS = [
    "m3_t2_summary_present",
    "m3_t2_passed_or_bypassed",
    "ligand_manifest_present",
    "eligible_ligands_found",
    "ligand_paths_inside_run_dir",
    "m3_t3_summary_present",
    "m3_t3_passed_or_bypassed",
    "receptor_manifest_present",
    "eligible_receptors_found",
    "receptor_paths_inside_run_dir",
    "docking_box_manifest_present",
    "eligible_boxes_found",
    "boxes_traceable_to_m2",
    "boxes_non_atp_pass",
    "boxes_lower_lateral_pass",
    "boxes_dimer_accessibility_pass",
    "m3_t4_summary_present",
    "m3_t4_passed",
    "m3_t5_allowed_by_t4",
    "mini_required_before_production",
    "production_blocked_without_mini",
    "vina_available",
    "vina_version_recorded",
    "job_manifest_written",
    "job_manifest_deterministic",
    "job_count_positive",
    "job_ids_unique",
    "planned_output_paths_inside_run_dir",
    "planned_stdout_stderr_inside_logs_jobs",
    "exactly_expected_job_count_planned",
    "no_vina_invoked_by_generator",
    "no_qsub_invoked_by_default",
    "pbs_manifest_written",
    "pbs_files_written",
    "pbs_stdout_stderr_concrete",
    "pbs_directives_no_unresolved_variables",
    "pbs_sets_pythonpath",
    "pbs_sets_thread_limits",
    "pbs_uses_conda_env",
    "pbs_runner_command_valid",
    "node_assignment_valid",
    "worker_count_not_oversubscribed",
    "reference_state_excluded_by_default",
    "reference_state_not_primary",
    "no_production_outputs_created",
    "no_broad_docking_created",
    "no_pose_collection_attempted",
    "no_pose_attribution_attempted",
    "no_candidate_ranking_attempted",
    "no_internal_id_leak",
    "no_smiles_logged",
    "no_ligand_coordinates_logged",
    "no_receptor_coordinates_logged",
    "generated_outputs_gitignored",
    "generated_outputs_not_tracked",
    "old_workflow_not_used",
    "non_goals_preserved",
]

GITIGNORE_PATTERNS = [
    "fresh/runs/*/phase3_compounds/docking_outputs/*",
    "fresh/runs/*/phase3_compounds/docking_outputs/**",
    "fresh/runs/*/phase3_compounds/docking_inputs/*",
    "fresh/runs/*/phase3_compounds/docking_inputs/**",
    "*.vina.tmp",
    "*.vina.log.tmp",
    "*.vina.stdout",
    "*.vina.stderr",
    "*.smoke_pose.pdbqt",
    "*.pose_out.pdbqt",
]

PRIMARY_STATE_ORDER = ["EGFR_160-185", "EGFR_170-200"]
ALL_NODES = ["node04", "node05", "node06"]
FORBIDDEN_OUTPUTS = [
    "phase3_compounds/docking_outputs/broad_anchor_scan_optional",
    "phase3_compounds/tables/compound_pose_raw.csv",
    "phase3_compounds/tables/compound_pose_attribution.csv",
    "phase3_compounds/tables/compound_pose_clusters.csv",
    "phase3_compounds/tables/compound_anchor_convergence.csv",
    "phase3_compounds/tables/final_m3_candidate_hypotheses.csv",
    "phase3_compounds/tables/pocket_compound_evidence_table.csv",
]


@dataclass
class M3VinaJobPlanResult:
    status: str
    qsub_submission_allowed: bool
    production_submission_allowed: bool
    blockers: list[str]
    warnings: list[str]
    job_manifest_csv: Path
    profile_job_manifest_csv: Path
    pbs_manifest_csv: Path
    qc_json: Path
    qc_csv: Path
    report_md: Path
    qsub_dir: Path
    phase3_log: Path
    counts: dict[str, int]
    selection: dict[str, Any]
    hpc: dict[str, Any]
    vina: dict[str, Any]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _bool_text(value: bool) -> str:
    return "true" if value else "false"


def _boolish(value: Any) -> bool:
    return str(value or "").strip().lower() in {"true", "t", "1", "yes", "y", "pass", "passed"}


def _finite_float(value: Any) -> float | None:
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    if math.isnan(number) or math.isinf(number):
        return None
    return number


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sanitize(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))
    return cleaned.strip("_") or "NA"


def _severity(status: str, blocker: bool = False) -> str:
    if status == "FAIL" and blocker:
        return "BLOCKER"
    if status in {"FAIL", "WARN"}:
        return "MAJOR"
    if status == "NOT_APPLICABLE":
        return "MINOR"
    return "INFO"


def _append_qc(rows: list[dict[str, str]], check_id: str, category: str, status: str, details: str, fix: str = "", blocker: bool = False) -> None:
    rows.append(
        {
            "check_id": check_id,
            "category": category,
            "status": status,
            "severity": _severity(status, blocker),
            "details": details,
            "recommended_fix": fix,
        }
    )


def _write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def _load_json(path: Path) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _ensure_logs(ctx: RunContext) -> None:
    ctx.create_directories()
    for directory in [ctx.logs_dir, ctx.jobs_log_dir, ctx.errors_dir]:
        ctx.require_within_run_dir(directory).mkdir(parents=True, exist_ok=True)
    for path in [
        ctx.logs_dir / "master.log",
        ctx.logs_dir / "phase_status.jsonl",
        ctx.logs_dir / "job_status.jsonl",
        ctx.errors_dir / "error_summary.txt",
        ctx.errors_dir / "failed_jobs.csv",
        ctx.logs_dir / "phase3_compounds.log",
    ]:
        safe = ctx.require_within_run_dir(path)
        if not safe.exists():
            safe.touch()
    failed = ctx.errors_dir / "failed_jobs.csv"
    if failed.stat().st_size == 0:
        failed.write_text("timestamp,job_name,status,message\n", encoding="utf-8")


def update_gitignore(repo_root: Path) -> tuple[bool, list[str]]:
    gitignore = repo_root / ".gitignore"
    existing = gitignore.read_text(encoding="utf-8").splitlines() if gitignore.is_file() else []
    missing = [pattern for pattern in GITIGNORE_PATTERNS if pattern not in existing]
    if not missing:
        return False, list(GITIGNORE_PATTERNS)
    with gitignore.open("a", encoding="utf-8") as handle:
        if existing and existing[-1].strip():
            handle.write("\n")
        handle.write("\n# M3 Vina job-planning output protection\n")
        for pattern in missing:
            handle.write(pattern + "\n")
    return True, list(GITIGNORE_PATTERNS)


def _path_from_row(ctx: RunContext, value: Any) -> Path | None:
    text = str(value or "").strip()
    if not text:
        return None
    path = Path(text)
    if not path.is_absolute():
        path = ctx.repo_root / path
    try:
        return path.resolve()
    except OSError:
        return None


def _inside(path: Path | None, parent: Path) -> bool:
    if path is None:
        return False
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def _safe_private_map_path(ctx: RunContext, private_map: Path | None) -> Path:
    path = private_map or (ctx.fresh_root / "data" / "private" / "compound_id_map.csv")
    resolved = Path(path)
    if not resolved.is_absolute():
        resolved = ctx.repo_root / resolved
    return ensure_within(resolved.resolve(), ctx.fresh_root)


def _load_hpc_defaults(ctx: RunContext) -> dict[str, Any]:
    config_path = ctx.fresh_root / "configs" / "hpc.yaml"
    data = load_yaml_config(config_path) if config_path.is_file() else {}
    conda = data.get("conda", {}) if isinstance(data.get("conda"), dict) else {}
    activate = str(conda.get("activate_command") or "")
    conda_sh = DEFAULT_CONDA_SH
    m_sh = re.search(r"source\s+(\S+)", activate)
    if m_sh:
        conda_sh = m_sh.group(1)
    return {
        "queue": data.get("queue", "workq"),
        "nodes": list(data.get("nodes", ALL_NODES)) or ALL_NODES,
        "ppn": dict(data.get("ppn", {})),
        "walltime": dict(data.get("walltime", {})),
        "conda_env": conda.get("env_name", DEFAULT_CONDA_ENV),
        "conda_sh": conda_sh,
        "python_executable": (data.get("python", {}) or {}).get("path_hint", "python") if isinstance(data.get("python"), dict) else "python",
        "thread_limits": dict(data.get("thread_limits", {})),
    }


def _default_profile_values(ctx: RunContext, profile: str, args: dict[str, Any]) -> dict[str, Any]:
    hpc = _load_hpc_defaults(ctx)
    defaults = {
        "smoke": {"exhaustiveness": 1, "num_modes": 1, "repeats": 1, "timeout_seconds": 300, "ppn": 4, "walltime": "02:00:00"},
        "mini": {"exhaustiveness": 8, "num_modes": 10, "repeats": 2, "timeout_seconds": 1800, "ppn": 16, "walltime": "08:00:00"},
        "scaling": {"exhaustiveness": 8, "num_modes": 10, "repeats": 2, "timeout_seconds": 1800, "ppn": 32, "walltime": "12:00:00"},
        "production": {"exhaustiveness": 8, "num_modes": 10, "repeats": 5, "timeout_seconds": 3600, "ppn": 32, "walltime": "48:00:00"},
    }[profile]
    values = {
        **defaults,
        "queue": hpc["queue"],
        "nodes": args.get("nodes") or hpc["nodes"],
        "conda_env": args.get("conda_env") or hpc["conda_env"],
        "conda_sh": args.get("conda_sh") or hpc["conda_sh"],
        "python_executable": args.get("python_executable") or hpc["python_executable"],
    }
    values["ppn"] = int(args.get("ppn") or hpc["ppn"].get(profile) or values["ppn"])
    values["walltime"] = str(args.get("walltime") or hpc["walltime"].get(profile) or values["walltime"])
    for key in ["exhaustiveness", "num_modes", "repeats", "timeout_seconds"]:
        if args.get(key) is not None:
            values[key] = int(args[key])
    values["cpu_per_vina"] = int(args.get("cpu_per_vina") or 1)
    values["base_seed"] = int(args.get("base_seed") or 20260427)
    values["thread_limits"] = {key: 1 for key in THREAD_ENV_KEYS}
    return values


def _state_sort_key(row: dict[str, Any]) -> tuple[int, str]:
    state = str(row.get("state_id") or "")
    role = str(row.get("state_role") or "unknown").lower()
    if state in PRIMARY_STATE_ORDER:
        return (PRIMARY_STATE_ORDER.index(state), state)
    if role == "primary":
        return (10, state)
    if state == "3GT8_raw" or role == "reference":
        return (100, state)
    return (50, state)


def _eligible_ligands(ctx: RunContext, manifest_path: Path, requested: str | None, notes: list[str]) -> list[dict[str, Any]]:
    if not manifest_path.is_file():
        notes.append("ligand_preparation_manifest.csv missing")
        return []
    rows, _ = _read_csv(manifest_path)
    eligible: list[dict[str, Any]] = []
    for row in rows:
        cid = (row.get("compound_public_id") or "").strip()
        if requested and cid != requested:
            continue
        path = _path_from_row(ctx, row.get("prepared_pdbqt_file") or row.get("ligand_pdbqt_file"))
        status = row.get("preparation_status") or row.get("ligand_preparation_status")
        if (
            cid in PUBLIC_COMPOUND_IDS
            and path is not None
            and _inside(path, ctx.run_dir / "phase3_compounds")
            and path.is_file()
            and path.stat().st_size > 0
            and _boolish(row.get("pdbqt_validation_success"))
            and status in {"PASS", "WARN"}
        ):
            current_sha = _sha256(path)
            manifest_sha = (row.get("prepared_pdbqt_sha256") or row.get("ligand_pdbqt_sha256") or "").strip()
            if manifest_sha and manifest_sha != current_sha:
                notes.append("{0}: prepared ligand PDBQT hash mismatch against M3-T2/T4 manifest".format(cid or "<missing>"))
                continue
            enriched = dict(row)
            enriched["_path"] = path
            enriched["_sha256"] = current_sha
            enriched["_status"] = status
            eligible.append(enriched)
    order = {cid: idx for idx, cid in enumerate(PUBLIC_COMPOUND_IDS)}
    return sorted(eligible, key=lambda row: order.get(row["compound_public_id"], 99))


def _eligible_receptors(ctx: RunContext, manifest_path: Path, requested: str | None, allow_reference: bool, notes: list[str]) -> list[dict[str, Any]]:
    if not manifest_path.is_file():
        notes.append("receptor_preparation_manifest.csv missing")
        return []
    rows, _ = _read_csv(manifest_path)
    eligible: list[dict[str, Any]] = []
    for row in rows:
        state = (row.get("state_id") or "").strip()
        role = (row.get("state_role") or "unknown").strip() or "unknown"
        if requested and state != requested:
            continue
        is_reference = state == "3GT8_raw" or role == "reference"
        if is_reference and not allow_reference:
            continue
        path = _path_from_row(ctx, row.get("prepared_receptor_pdbqt_file") or row.get("receptor_pdbqt_file"))
        allowed = _boolish(row.get("allowed_for_compound_docking")) or _boolish(row.get("allowed_for_vina_smoke"))
        if (
            state
            and path is not None
            and _inside(path, ctx.run_dir / "phase3_compounds" / "receptor_pdbqt")
            and path.is_file()
            and path.stat().st_size > 0
            and _boolish(row.get("pdbqt_validation_success"))
            and row.get("preparation_status") in {"PASS", "WARN"}
            and allowed
            and not _boolish(row.get("source_receptor_is_old_workflow"))
            and not _boolish(row.get("source_receptor_is_monomer_only"))
        ):
            current_sha = _sha256(path)
            manifest_sha = (row.get("prepared_receptor_pdbqt_sha256") or row.get("receptor_pdbqt_sha256") or "").strip()
            if manifest_sha and manifest_sha != current_sha:
                notes.append("{0}: prepared receptor PDBQT hash mismatch against M3-T3/T4 manifest".format(state or "<missing>"))
                continue
            enriched = dict(row)
            enriched["_path"] = path
            enriched["_sha256"] = current_sha
            eligible.append(enriched)
    return sorted(eligible, key=_state_sort_key)


def _eligible_boxes(
    ctx: RunContext,
    manifest_path: Path,
    receptors: list[dict[str, Any]],
    requested_state: str | None,
    requested_family: str | None,
    requested_box: str | None,
    requested_protomer: str | None,
    notes: list[str],
) -> list[dict[str, Any]]:
    if not manifest_path.is_file():
        notes.append("docking_box_manifest.csv missing")
        return []
    receptor_by_state = {row["state_id"]: row for row in receptors}
    rows, _ = _read_csv(manifest_path)
    eligible: list[dict[str, Any]] = []
    for row in rows:
        if requested_state and row.get("state_id") != requested_state:
            continue
        if requested_family and row.get("pocket_family_id") != requested_family:
            continue
        if requested_box and row.get("box_id") != requested_box:
            continue
        if requested_protomer and row.get("protomer_id") != requested_protomer:
            continue
        receptor = receptor_by_state.get(row.get("state_id", ""))
        centers = [_finite_float(row.get(k)) for k in ["box_center_x", "box_center_y", "box_center_z"]]
        sizes = [_finite_float(row.get(k)) for k in ["box_size_x", "box_size_y", "box_size_z"]]
        allowed = _boolish(row.get("allowed_for_compound_docking")) or _boolish(row.get("allowed_for_vina_smoke"))
        ok = (
            bool(row.get("pocket_family_id"))
            and bool(row.get("state_id"))
            and bool(row.get("protomer_id"))
            and bool(row.get("box_id"))
            and receptor is not None
            and bool(row.get("receptor_pdbqt_file"))
            and row.get("receptor_pdbqt_file") == receptor.get("prepared_receptor_pdbqt_file")
            and all(value is not None for value in centers)
            and all(value is not None and value > 0 for value in sizes)
            and _boolish(row.get("non_atp_pass"))
            and _boolish(row.get("lower_lateral_pass"))
            and _boolish(row.get("dimer_accessibility_pass"))
            and row.get("box_qc_status") == "PASS"
            and row.get("traceability_status") == "PASS"
            and allowed
        )
        if ok:
            enriched = dict(row)
            enriched["_receptor"] = receptor
            eligible.append(enriched)
    return sorted(eligible, key=lambda row: (row.get("state_id", ""), row.get("pocket_family_id", ""), row.get("protomer_id", ""), row.get("box_id", "")))


def _profile_boxes(profile: str, boxes: list[dict[str, Any]], requested_state: str | None, requested_family: str | None) -> list[dict[str, Any]]:
    if profile == "smoke":
        primaries = [box for box in boxes if (box.get("_receptor") or {}).get("state_role") != "reference"]
        return (primaries or boxes)[:1]
    if profile == "scaling":
        selected_state = requested_state or (boxes[0].get("state_id") if boxes else None)
        state_boxes = [box for box in boxes if box.get("state_id") == selected_state]
        selected_family = requested_family or (state_boxes[0].get("pocket_family_id") if state_boxes else None)
        return [box for box in state_boxes if box.get("pocket_family_id") == selected_family]
    if profile == "mini":
        selected: list[dict[str, Any]] = []
        for state in sorted({box.get("state_id") for box in boxes}, key=lambda s: PRIMARY_STATE_ORDER.index(s) if s in PRIMARY_STATE_ORDER else 50):
            state_boxes = [box for box in boxes if box.get("state_id") == state]
            families = []
            for box in state_boxes:
                fam = box.get("pocket_family_id")
                if fam not in families:
                    families.append(fam)
            for fam in families[:2]:
                selected.extend([box for box in state_boxes if box.get("pocket_family_id") == fam])
        return selected
    return boxes


def _node_for_state(state_id: str, nodes: list[str]) -> str:
    safe_nodes = nodes or ALL_NODES
    if state_id == "EGFR_160-185":
        return "node04" if "node04" in safe_nodes else safe_nodes[0]
    if state_id == "EGFR_170-200":
        return "node05" if "node05" in safe_nodes else safe_nodes[min(1, len(safe_nodes) - 1)]
    return "node06" if "node06" in safe_nodes else safe_nodes[-1]


def _chunk_id(profile: str, state_id: str, family: str) -> str:
    return "{0}_{1}_{2}_chunk01".format(profile, _sanitize(state_id), _sanitize(family))


def _planned_forbidden_outputs(ctx: RunContext) -> list[str]:
    hits: list[str] = []

    def is_real_output(path: Path) -> bool:
        return path.is_file() and path.name != ".gitkeep"

    for rel in FORBIDDEN_OUTPUTS:
        path = ctx.run_dir / rel
        if path.is_file():
            hits.append(ctx.relative_to_repo(path))
        elif path.is_dir() and any(is_real_output(child) for child in path.rglob("*")):
            hits.append(ctx.relative_to_repo(path))
    return hits


def _scan_hygiene(ctx: RunContext, private_entries: list[Any]) -> tuple[int, list[str], bool, bool, bool, list[str]]:
    leaks = scan_internal_id_leaks(ctx, private_entries)
    scan_paths = public_output_scan_paths(ctx)
    qsub_root = ctx.run_dir / "phase3_compounds" / "qsub"
    if qsub_root.is_dir():
        scan_paths.extend(path for path in sorted(qsub_root.rglob("*")) if path.is_file())
    coord_hits: list[str] = []
    smiles_logged = False
    ligand_coords = False
    receptor_coords = False
    scanned: list[str] = []
    for path in scan_paths:
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        scanned.append(ctx.relative_to_repo(path))
        if re.search(r"\bSMILES\b|canonical_smiles|isomeric_smiles", text, re.IGNORECASE):
            smiles_logged = True
        if re.search(r"^(ATOM|HETATM)\s+\d+", text, re.MULTILINE):
            coord_hits.append(ctx.relative_to_repo(path))
            if "ligand" in path.as_posix().lower():
                ligand_coords = True
            if "receptor" in path.as_posix().lower():
                receptor_coords = True
    return len(leaks), coord_hits, smiles_logged, ligand_coords, receptor_coords, scanned


def _write_report(path: Path, ctx: RunContext, summary: dict[str, Any]) -> None:
    lines = [
        "# M3-T5 Vina Job Plan",
        "",
        "Technical HPC job-planning report only. No qsub submission, Vina execution, pose collection, scoring, ranking, or candidate claims were performed.",
        "",
        "- status: {0}".format(summary["overall_status"]),
        "- profile: {0}".format(summary["profile"]),
        "- qsub_submission_allowed: {0}".format(str(summary["qsub_submission_allowed"]).lower()),
        "- production_submission_allowed: {0}".format(str(summary["production_submission_allowed"]).lower()),
        "- planned_vina_commands: {0}".format(summary["counts"]["planned_vina_commands"]),
        "- planned_chunks: {0}".format(summary["counts"]["planned_chunks"]),
        "- pbs_files: {0}".format(summary["counts"]["pbs_files"]),
        "- generator_invoked_vina: false",
        "- generator_invoked_qsub: false",
        "",
        "## Blockers",
    ]
    lines.extend("- {0}".format(item) for item in summary["blockers"] or ["none"])
    lines.append("")
    lines.append("## Warnings")
    lines.extend("- {0}".format(item) for item in summary["warnings"] or ["none"])
    lines.append("")
    lines.append("Next task: M3-T6 - Vina output collection and raw pose table")
    ctx.require_within_run_dir(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _pbs_directives_safe(content: str) -> bool:
    forbidden = ["$RUN_ID", "${RUN_ID}", "$PBS_JOBID", "$(pwd)"]
    for line in content.splitlines():
        if line.startswith("#PBS") and any(token in line for token in forbidden):
            return False
        if line.startswith("#PBS") and ("<" in line or ">" in line):
            return False
    return True


def run_m3_vina_jobs(
    ctx: RunContext,
    *,
    m2_run_id: str,
    mode: str = "generate",
    profile: str = "mini",
    force: bool = False,
    compound_public_id: str | None = None,
    state_id: str | None = None,
    pocket_family_id: str | None = None,
    box_id: str | None = None,
    protomer_id: str | None = None,
    include_reference: bool = False,
    include_3gt8: bool = False,
    allow_reference_jobs: bool = False,
    allow_independent: bool = False,
    allow_production_plan_without_mini: bool = False,
    vina_executable: str = "vina",
    exhaustiveness: int | None = None,
    num_modes: int | None = None,
    cpu_per_vina: int = 1,
    repeats: int | None = None,
    base_seed: int | None = None,
    timeout_seconds: int | None = None,
    nodes: list[str] | None = None,
    ppn: int | None = None,
    walltime: str | None = None,
    queue: str | None = None,
    conda_env: str | None = None,
    conda_sh: str | None = None,
    python_executable: str | None = None,
    update_ignore: bool = True,
    private_map: Path | None = None,
) -> M3VinaJobPlanResult:
    _ensure_logs(ctx)
    phase3 = ctx.run_dir / "phase3_compounds"
    manifests = phase3 / "manifests"
    qc_dir = phase3 / "qc"
    reports = phase3 / "reports"
    qsub_dir = phase3 / "qsub"
    for directory in [manifests, qc_dir, reports, qsub_dir]:
        ctx.require_within_run_dir(directory).mkdir(parents=True, exist_ok=True)

    job_manifest = manifests / "vina_job_manifest.csv"
    profile_manifest = manifests / "vina_job_manifest_{0}.csv".format(profile)
    pbs_manifest = manifests / "vina_pbs_manifest.csv"
    qc_json = qc_dir / "vina_job_plan_qc.json"
    qc_csv = qc_dir / "vina_job_plan_qc.csv"
    report_md = reports / "m3_task5_vina_job_plan.md"
    phase3_log = ctx.logs_dir / "phase3_compounds.log"

    blockers: list[str] = []
    warnings: list[str] = []
    qc_rows: list[dict[str, str]] = []
    notes: list[str] = []

    hpc_values = _default_profile_values(
        ctx,
        profile,
        {
            "exhaustiveness": exhaustiveness,
            "num_modes": num_modes,
            "repeats": repeats,
            "timeout_seconds": timeout_seconds,
            "cpu_per_vina": cpu_per_vina,
            "base_seed": base_seed,
            "nodes": nodes,
            "ppn": ppn,
            "walltime": walltime,
            "conda_env": conda_env,
            "conda_sh": conda_sh,
            "python_executable": python_executable,
        },
    )
    if queue:
        hpc_values["queue"] = queue

    private_entries, private_warnings = load_private_map(_safe_private_map_path(ctx, private_map))
    warnings.extend(private_warnings)
    gitignore_updated, patterns_verified = update_gitignore(ctx.repo_root) if update_ignore else (False, list(GITIGNORE_PATTERNS))

    t2_summary_path = qc_dir / "m3_ligand_prep_summary.json"
    t3_summary_path = qc_dir / "m3_receptor_box_summary.json"
    t4_summary_path = qc_dir / "vina_smoke_qc.json"
    ligand_manifest_path = manifests / "ligand_preparation_manifest.csv"
    receptor_manifest_path = manifests / "receptor_preparation_manifest.csv"
    box_manifest_path = manifests / "docking_box_manifest.csv"
    smoke_manifest_path = manifests / "vina_smoke_manifest.csv"

    t2 = _load_json(t2_summary_path) if t2_summary_path.is_file() else None
    t3 = _load_json(t3_summary_path) if t3_summary_path.is_file() else None
    t4 = _load_json(t4_summary_path) if t4_summary_path.is_file() else None
    t2_pass = bool(t2 and t2.get("overall_status") == "PASS")
    t3_pass = bool(t3 and t3.get("overall_status") == "PASS")
    t4_pass = bool(t4 and t4.get("overall_status") == "PASS")
    t5_allowed = bool(t4 and t4.get("m3_t5_allowed") is True)

    def gate_missing(condition: bool, check_id: str, name: str) -> None:
        status = "PASS" if condition else ("WARN" if mode == "dry-run" else "FAIL")
        _append_qc(qc_rows, check_id, "input_gate", status, "{0} found={1}".format(name, condition), "Generate upstream M3 outputs first.", mode == "generate" and not condition)
        if not condition:
            message = "{0} missing".format(name)
            if mode == "generate":
                blockers.append(message)
            else:
                warnings.append(message)

    gate_missing(t2 is not None, "m3_t2_summary_present", "M3-T2 ligand prep summary")
    gate_missing(t3 is not None, "m3_t3_summary_present", "M3-T3 receptor/box summary")
    gate_missing(t4 is not None, "m3_t4_summary_present", "M3-T4 smoke QC")
    gate_missing(ligand_manifest_path.is_file(), "ligand_manifest_present", "ligand_preparation_manifest.csv")
    gate_missing(receptor_manifest_path.is_file(), "receptor_manifest_present", "receptor_preparation_manifest.csv")
    gate_missing(box_manifest_path.is_file(), "docking_box_manifest_present", "docking_box_manifest.csv")

    for check_id, ok, label in [
        ("m3_t2_passed_or_bypassed", t2_pass or allow_independent, "M3-T2 passed or diagnostic bypass recorded"),
        ("m3_t3_passed_or_bypassed", t3_pass or allow_independent, "M3-T3 passed or diagnostic bypass recorded"),
        ("m3_t4_passed", t4_pass or allow_independent, "M3-T4 passed or diagnostic bypass recorded"),
        ("m3_t5_allowed_by_t4", t5_allowed or allow_independent, "M3-T4 opened M3-T5 gate or diagnostic bypass recorded"),
    ]:
        status = "PASS" if ok else ("WARN" if mode == "dry-run" else "FAIL")
        _append_qc(qc_rows, check_id, "input_gate", status, label, "Complete upstream task or use an explicit diagnostic bypass.", mode == "generate" and not ok)
        if not ok:
            msg = check_id.replace("_", " ")
            if mode == "generate":
                blockers.append(msg)
            else:
                warnings.append(msg)
    if allow_independent:
        warnings.append("allow_independent diagnostic bypass used; submission gates remain closed")

    allow_reference = (include_reference or include_3gt8) and allow_reference_jobs
    if (include_reference or include_3gt8) and not allow_reference_jobs:
        blockers.append("reference state requested without --allow-reference-jobs") if mode == "generate" else warnings.append("reference state requested without --allow-reference-jobs")

    ligands = _eligible_ligands(ctx, ligand_manifest_path, compound_public_id, notes)
    receptors = _eligible_receptors(ctx, receptor_manifest_path, state_id, allow_reference, notes)
    boxes = _eligible_boxes(ctx, box_manifest_path, receptors, state_id, pocket_family_id, box_id, protomer_id, notes)
    boxes = _profile_boxes(profile, boxes, state_id, pocket_family_id)
    if profile == "smoke" and ligands:
        ligands = ligands[:1]
    if profile == "scaling" and receptors:
        states = {box.get("state_id") for box in boxes}
        receptors = [row for row in receptors if row.get("state_id") in states]

    for check_id, count, label in [
        ("eligible_ligands_found", len(ligands), "eligible ligands"),
        ("eligible_receptors_found", len(receptors), "eligible receptors"),
        ("eligible_boxes_found", len(boxes), "eligible boxes"),
    ]:
        status = "PASS" if count > 0 else ("WARN" if mode == "dry-run" else "FAIL")
        _append_qc(qc_rows, check_id, "selection", status, "{0}: {1}".format(label, count), "Inspect upstream manifests.", mode == "generate" and count == 0)
        if count == 0:
            if mode == "generate":
                blockers.append("no {0}".format(label))
            else:
                warnings.append("no {0}".format(label))

    primary_states = [row for row in receptors if row.get("state_role") == "primary"]
    reference_only = bool(receptors) and not primary_states
    if reference_only:
        warnings.append("reference-only job plan; submission gates remain closed")
    _append_qc(qc_rows, "reference_state_excluded_by_default", "selection", "PASS" if allow_reference or not any(row.get("state_id") == "3GT8_raw" for row in receptors) else "FAIL", "3GT8_raw included only with explicit flags")
    _append_qc(qc_rows, "reference_state_not_primary", "selection", "PASS" if not any(row.get("state_id") == "3GT8_raw" and row.get("state_role") == "primary" for row in receptors) else "FAIL", "reference state role checked", "Keep 3GT8_raw as reference.", True)

    mini_completion_path = qc_dir / "vina_mini_completion_qc.json"
    mini_pass = False
    if mini_completion_path.is_file():
        mini_completion = _load_json(mini_completion_path) or {}
        mini_pass = mini_completion.get("overall_status") == "PASS"
    production_bypass = profile == "production" and allow_production_plan_without_mini
    if profile == "production" and not mini_pass and not production_bypass:
        blockers.append("production profile requested without mini completion PASS")
    if production_bypass:
        warnings.append("allow_production_plan_without_mini used; production submission remains closed")
    _append_qc(qc_rows, "mini_required_before_production", "profile", "PASS" if profile != "production" or mini_pass or production_bypass else "FAIL", "mini completion gate checked", "Complete mini profile and collect M3-T6 QC before production.", profile == "production" and not mini_pass and not production_bypass)
    _append_qc(qc_rows, "production_blocked_without_mini", "profile", "PASS" if profile != "production" or mini_pass else ("WARN" if production_bypass else "FAIL"), "production gating checked")

    vina = detect_vina(vina_executable)
    _append_qc(qc_rows, "vina_available", "tool", "PASS" if vina.available else ("WARN" if mode == "dry-run" else "FAIL"), "Vina availability detected by version probe", "Install Vina or provide --vina-executable.", mode == "generate" and not vina.available)
    _append_qc(qc_rows, "vina_version_recorded", "tool", "PASS" if vina.version else ("WARN" if vina.available else "NOT_APPLICABLE"), "Vina version recorded")
    if mode == "generate" and not vina.available:
        blockers.append("Vina executable unavailable in generate mode")

    planned_at = now_iso()
    jobs: list[dict[str, Any]] = []
    if ligands and boxes:
        stable_index = 0
        for box in boxes:
            receptor = box["_receptor"]
            node = _node_for_state(box["state_id"], hpc_values["nodes"])
            chunk_id = _chunk_id(profile, box["state_id"], box["pocket_family_id"])
            for ligand in ligands:
                for repeat in range(1, int(hpc_values["repeats"]) + 1):
                    stable_index += 1
                    job_id = "vina_{0}_{1}_{2}_{3}_{4}_{5}_r{6:02d}".format(
                        profile,
                        _sanitize(ligand["compound_public_id"]),
                        _sanitize(box["state_id"]),
                        _sanitize(box["pocket_family_id"]),
                        _sanitize(box["protomer_id"]),
                        _sanitize(box["box_id"]),
                        repeat,
                    )
                    output_dir = phase3 / "docking_outputs" / "focused_pocket_first" / profile / job_id
                    pose_name = "smoke_pose.pdbqt" if profile == "smoke" else "pose_out.pdbqt"
                    output_pdbqt = output_dir / pose_name
                    vina_log = output_dir / "vina.log"
                    stdout_file = ctx.jobs_log_dir / "{0}.vina.stdout".format(job_id)
                    stderr_file = ctx.jobs_log_dir / "{0}.vina.stderr".format(job_id)
                    seed = int(hpc_values["base_seed"]) + stable_index
                    argv = build_vina_argv(
                        vina.executable,
                        receptor["_path"],
                        ligand["_path"],
                        (box["box_center_x"], box["box_center_y"], box["box_center_z"]),
                        (box["box_size_x"], box["box_size_y"], box["box_size_z"]),
                        output_pdbqt,
                        vina_log,
                        int(hpc_values["exhaustiveness"]),
                        int(hpc_values["num_modes"]),
                        int(hpc_values["cpu_per_vina"]),
                        seed,
                    )
                    jobs.append(
                        {
                            "job_id": job_id,
                            "run_id": ctx.run_id,
                            "m2_run_id": m2_run_id,
                            "profile": profile,
                            "mode": mode,
                            "job_scope": "focused_pocket_first",
                            "chunk_id": chunk_id,
                            "node": node,
                            "ppn": str(hpc_values["ppn"]),
                            "worker_count": "",
                            "compound_public_id": ligand["compound_public_id"],
                            "ligand_pdbqt_file": ctx.relative_to_repo(ligand["_path"]),
                            "ligand_pdbqt_sha256": ligand["_sha256"],
                            "ligand_preparation_status": ligand["_status"],
                            "state_id": box["state_id"],
                            "state_role": receptor.get("state_role", "unknown"),
                            "receptor_pdbqt_file": ctx.relative_to_repo(receptor["_path"]),
                            "receptor_pdbqt_sha256": receptor["_sha256"],
                            "receptor_preparation_status": receptor.get("preparation_status", ""),
                            "pocket_family_id": box["pocket_family_id"],
                            "protomer_id": box["protomer_id"],
                            "box_id": box["box_id"],
                            "box_center_x": box["box_center_x"],
                            "box_center_y": box["box_center_y"],
                            "box_center_z": box["box_center_z"],
                            "box_size_x": box["box_size_x"],
                            "box_size_y": box["box_size_y"],
                            "box_size_z": box["box_size_z"],
                            "source_ligand_manifest": ctx.relative_to_repo(ligand_manifest_path),
                            "source_receptor_manifest": ctx.relative_to_repo(receptor_manifest_path),
                            "source_docking_box_manifest": ctx.relative_to_repo(box_manifest_path),
                            "source_vina_smoke_manifest": ctx.relative_to_repo(smoke_manifest_path),
                            "vina_executable": vina.executable,
                            "vina_version": vina.version or "",
                            "vina_argv_json": json.dumps([str(item) for item in argv]),
                            "exhaustiveness": str(hpc_values["exhaustiveness"]),
                            "num_modes": str(hpc_values["num_modes"]),
                            "cpu_per_vina": str(hpc_values["cpu_per_vina"]),
                            "vina_repeat_id": str(repeat),
                            "seed": str(seed),
                            "timeout_seconds": str(hpc_values["timeout_seconds"]),
                            "output_dir": ctx.relative_to_repo(output_dir),
                            "output_pdbqt_file": ctx.relative_to_repo(output_pdbqt),
                            "vina_log_file": ctx.relative_to_repo(vina_log),
                            "stdout_file": ctx.relative_to_repo(stdout_file),
                            "stderr_file": ctx.relative_to_repo(stderr_file),
                            "pbs_file": "",
                            "pbs_job_name": "",
                            "command_manifest_written_before_execution": "true",
                            "planned_at": planned_at,
                            "submission_status": "NOT_SUBMITTED",
                            "docking_status": "PLANNED",
                            "allowed_for_execution": "false",
                            "allowed_for_collection": "false",
                            "job_notes": "planned only; generator did not run Vina",
                        }
                    )

    chunk_job_counts: dict[str, int] = {}
    for row in jobs:
        chunk_job_counts[row["chunk_id"]] = chunk_job_counts.get(row["chunk_id"], 0) + 1
    worker_by_chunk: dict[str, int] = {}
    max_workers = max(1, int(hpc_values["ppn"]) // max(1, int(hpc_values["cpu_per_vina"])))
    for chunk_id, count in chunk_job_counts.items():
        worker_by_chunk[chunk_id] = max(1, min(max_workers, count))
    for row in jobs:
        row["worker_count"] = str(worker_by_chunk.get(row["chunk_id"], 1))

    _write_csv(job_manifest, JOB_FIELDS, jobs, ctx)
    _write_csv(profile_manifest, JOB_FIELDS, jobs, ctx)

    job_ids = [row["job_id"] for row in jobs]
    duplicate_ids = sorted({job_id for job_id in job_ids if job_ids.count(job_id) > 1})
    if duplicate_ids:
        blockers.append("duplicate job_id values detected")
    if not jobs and mode == "generate":
        blockers.append("zero planned jobs")
    _append_qc(qc_rows, "job_manifest_written", "manifest", "PASS" if job_manifest.is_file() else "FAIL", "job manifest written before possible execution", "Write vina_job_manifest.csv.", not job_manifest.is_file())
    _append_qc(qc_rows, "job_manifest_deterministic", "manifest", "PASS", "job enumeration sorted by public ID, state, pocket, protomer, box, repeat")
    _append_qc(qc_rows, "job_count_positive", "manifest", "PASS" if jobs else ("WARN" if mode == "dry-run" else "FAIL"), "planned job count={0}".format(len(jobs)), "Fix input gates.", mode == "generate" and not jobs)
    _append_qc(qc_rows, "job_ids_unique", "manifest", "PASS" if not duplicate_ids else "FAIL", "unique job_id count checked", "Ensure tuple-derived job IDs are unique.", bool(duplicate_ids))

    def row_path_inside(row: dict[str, Any], field: str, parent: Path) -> bool:
        return _inside(_path_from_row(ctx, row.get(field)), parent)

    output_paths_ok = all(row_path_inside(row, "output_pdbqt_file", phase3 / "docking_outputs" / "focused_pocket_first") for row in jobs)
    stdio_ok = all(row_path_inside(row, "stdout_file", ctx.jobs_log_dir) and row_path_inside(row, "stderr_file", ctx.jobs_log_dir) for row in jobs)
    if jobs and not output_paths_ok:
        blockers.append("planned output paths outside run_dir")
    if jobs and not stdio_ok:
        blockers.append("planned stdout/stderr outside logs/jobs")
    _append_qc(qc_rows, "planned_output_paths_inside_run_dir", "paths", "PASS" if output_paths_ok else "FAIL", "planned output PDBQT/log paths remain under run dir", "Fix output path construction.", bool(jobs and not output_paths_ok))
    _append_qc(qc_rows, "planned_stdout_stderr_inside_logs_jobs", "paths", "PASS" if stdio_ok else "FAIL", "per-Vina stdout/stderr paths remain under logs/jobs", "Fix stdout/stderr path construction.", bool(jobs and not stdio_ok))

    pbs_rows: list[dict[str, Any]] = []
    pbs_contents: list[str] = []
    chunks: dict[str, list[dict[str, Any]]] = {}
    for row in jobs:
        chunks.setdefault(row["chunk_id"], []).append(row)

    node_chunks: dict[str, str] = {}
    for node in hpc_values["nodes"]:
        node_jobs = [row for row in jobs if row["node"] == node]
        node_chunks[node] = node_jobs[0]["chunk_id"] if node_jobs else "empty_{0}_{1}".format(profile, node)

    pbs_targets: list[tuple[str, str, str, list[dict[str, Any]]]] = []
    for chunk_id, rows in sorted(chunks.items()):
        node = rows[0]["node"]
        pbs_targets.append(("m3_vina_{0}_{1}".format(profile, chunk_id), chunk_id, node, rows))
    for node in ALL_NODES:
        rows = [row for row in jobs if row["node"] == node]
        pbs_targets.append(("m3_vina_{0}".format(node), node_chunks.get(node, "empty_{0}_{1}".format(profile, node)), node, rows))

    for pbs_job_name, chunk_id, node, assigned_rows in pbs_targets:
        worker_count = worker_by_chunk.get(chunk_id, max(1, min(max_workers, len(assigned_rows) or 1)))
        has_assigned_jobs = bool(assigned_rows)
        pbs_file = qsub_dir / "{0}.pbs".format(pbs_job_name)
        stdout_file = ctx.jobs_log_dir / "{0}.stdout".format(pbs_job_name)
        stderr_file = ctx.jobs_log_dir / "{0}.stderr".format(pbs_job_name)
        runner_command = (
            "{python} -m egfr_myo1d.cli m3-run-vina-chunk --run-id {run_id} --job-manifest {manifest} --chunk-id {chunk} --max-workers {workers}"
        ).format(
            python=shlex.quote(str(hpc_values["python_executable"])),
            run_id=shlex.quote(ctx.run_id),
            manifest=shlex.quote(str(job_manifest.resolve())),
            chunk=shlex.quote(chunk_id),
            workers=worker_count,
        )
        content = render_pbs_content(
            job_name=pbs_job_name,
            node=node,
            ppn=int(hpc_values["ppn"]),
            walltime=str(hpc_values["walltime"]),
            queue=str(hpc_values["queue"]),
            repo_root=str(ctx.repo_root.resolve()),
            conda_sh=str(hpc_values["conda_sh"]),
            conda_env=str(hpc_values["conda_env"]),
            thread_limits=hpc_values["thread_limits"],
            stdout_path=str(stdout_file.resolve()),
            stderr_path=str(stderr_file.resolve()),
            run_id=ctx.run_id,
            command_lines=[runner_command],
        )
        with ctx.require_within_run_dir(pbs_file).open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        pbs_contents.append(content)
        for row in assigned_rows:
            row["pbs_file"] = ctx.relative_to_repo(pbs_file)
            row["pbs_job_name"] = pbs_job_name
            row["submission_status"] = "SUBMIT_READY" if mode == "generate" and not allow_independent else "NOT_SUBMITTED"
            row["allowed_for_execution"] = _bool_text(mode == "generate" and not allow_independent)
            row["allowed_for_collection"] = "false"
        pbs_rows.append(
            {
                "pbs_job_name": pbs_job_name,
                "run_id": ctx.run_id,
                "m2_run_id": m2_run_id,
                "profile": profile,
                "chunk_id": chunk_id,
                "node": node,
                "ppn": str(hpc_values["ppn"]),
                "worker_count": str(worker_count),
                "assigned_job_count": str(len(assigned_rows)),
                "pbs_file": ctx.relative_to_repo(pbs_file),
                "pbs_stdout_file": ctx.relative_to_repo(stdout_file),
                "pbs_stderr_file": ctx.relative_to_repo(stderr_file),
                "queue": hpc_values["queue"],
                "walltime": hpc_values["walltime"],
                "conda_env": hpc_values["conda_env"],
                "python_executable": hpc_values["python_executable"],
                "runner_command": runner_command,
                "qsub_command": "qsub {0}".format(ctx.relative_to_repo(pbs_file)) if has_assigned_jobs else "",
                "submit_by_default": "false",
                "pbs_generation_status": "PASS" if has_assigned_jobs else "EMPTY_COMPATIBILITY",
                "pbs_notes": "planned only; not submitted" if has_assigned_jobs else "compatibility wrapper has no assigned jobs and is not submittable",
            }
        )
        append_job_status(
            ctx,
            pbs_job_name,
            "PLANNED",
            node=node,
            ppn=int(hpc_values["ppn"]),
            stdout=ctx.relative_to_repo(stdout_file),
            stderr=ctx.relative_to_repo(stderr_file),
            details={
                "phase": "phase3_compounds",
                "task": "M3-T5",
                "job_type": "vina_pbs_chunk",
                "job_id": pbs_job_name,
                "profile": profile,
                "node": node,
                "worker_count": worker_count,
                "assigned_job_count": len(assigned_rows),
                "pbs_file": ctx.relative_to_repo(pbs_file),
                "qsub_command": "qsub {0}".format(ctx.relative_to_repo(pbs_file)),
                "submitted": False,
            },
        )

    _write_csv(job_manifest, JOB_FIELDS, jobs, ctx)
    _write_csv(profile_manifest, JOB_FIELDS, jobs, ctx)
    _write_csv(pbs_manifest, PBS_FIELDS, pbs_rows, ctx)

    pbs_safe = all(_pbs_directives_safe(content) for content in pbs_contents)
    pbs_pythonpath = all("PYTHONPATH=\"$REPO_ROOT/fresh/src:${PYTHONPATH:-}\"" in content for content in pbs_contents)
    pbs_threads = all(all("export {0}=1".format(key) in content for key in THREAD_ENV_KEYS) for content in pbs_contents)
    pbs_conda = all("conda activate {0}".format(hpc_values["conda_env"]) in content for content in pbs_contents)
    pbs_runner = all("m3-run-vina-chunk" in row["runner_command"] for row in pbs_rows)
    pbs_concrete = all("$" not in row["pbs_stdout_file"] and "$" not in row["pbs_stderr_file"] and row["pbs_stdout_file"] and row["pbs_stderr_file"] for row in pbs_rows)
    worker_ok = all(int(row["worker_count"]) <= int(row["ppn"]) // max(1, int(hpc_values["cpu_per_vina"])) for row in pbs_rows)
    for check_id, ok, detail in [
        ("pbs_manifest_written", pbs_manifest.is_file(), "PBS manifest written"),
        ("pbs_files_written", all((ctx.repo_root / row["pbs_file"]).is_file() for row in pbs_rows), "PBS files written"),
        ("pbs_stdout_stderr_concrete", pbs_concrete, "PBS stdout/stderr paths concrete"),
        ("pbs_directives_no_unresolved_variables", pbs_safe, "PBS directives have no unresolved variables/placeholders"),
        ("pbs_sets_pythonpath", pbs_pythonpath, "PBS sets PYTHONPATH to fresh/src"),
        ("pbs_sets_thread_limits", pbs_threads, "PBS sets BLAS/OpenMP thread limits"),
        ("pbs_uses_conda_env", pbs_conda, "PBS activates configured conda env"),
        ("pbs_runner_command_valid", pbs_runner, "PBS calls m3-run-vina-chunk"),
        ("node_assignment_valid", all(row["node"] in hpc_values["nodes"] for row in pbs_rows), "nodes are from configured node set"),
        ("worker_count_not_oversubscribed", worker_ok, "worker_count <= ppn / cpu_per_vina"),
    ]:
        _append_qc(qc_rows, check_id, "pbs", "PASS" if ok else "FAIL", detail, "Regenerate PBS files.", not ok)
        if not ok:
            blockers.append(check_id.replace("_", " "))

    _append_qc(qc_rows, "ligand_paths_inside_run_dir", "paths", "PASS" if all(_inside(row["_path"], phase3) for row in ligands) else "FAIL", "ligand PDBQT paths checked", "Regenerate ligand PDBQT under phase3_compounds.", True)
    _append_qc(qc_rows, "receptor_paths_inside_run_dir", "paths", "PASS" if all(_inside(row["_path"], phase3 / "receptor_pdbqt") for row in receptors) else "FAIL", "receptor PDBQT paths checked", "Regenerate receptor PDBQT under receptor_pdbqt.", True)
    _append_qc(qc_rows, "boxes_traceable_to_m2", "box", "PASS" if all(row.get("traceability_status") == "PASS" for row in boxes) else "FAIL", "box traceability checked", "Use M3-T3 active box manifest.", True)
    _append_qc(qc_rows, "boxes_non_atp_pass", "box", "PASS" if all(_boolish(row.get("non_atp_pass")) for row in boxes) else "FAIL", "non-ATP gate checked", "Use accepted non-ATP M2 boxes.", True)
    _append_qc(qc_rows, "boxes_lower_lateral_pass", "box", "PASS" if all(_boolish(row.get("lower_lateral_pass")) for row in boxes) else "FAIL", "lower/lateral gate checked", "Use accepted lower/lateral M2 boxes.", True)
    _append_qc(qc_rows, "boxes_dimer_accessibility_pass", "box", "PASS" if all(_boolish(row.get("dimer_accessibility_pass")) for row in boxes) else "FAIL", "dimer-accessibility gate checked", "Use accepted dimer-accessible M2 boxes.", True)

    expected_count = len(ligands) * len(boxes) * int(hpc_values["repeats"])
    _append_qc(qc_rows, "exactly_expected_job_count_planned", "manifest", "PASS" if len(jobs) == expected_count else "FAIL", "expected={0} observed={1}".format(expected_count, len(jobs)), "Regenerate deterministic manifest.", len(jobs) != expected_count)
    _append_qc(qc_rows, "no_vina_invoked_by_generator", "non_goal", "PASS", "generator did not invoke Vina")
    _append_qc(qc_rows, "no_qsub_invoked_by_default", "non_goal", "PASS", "generator did not invoke qsub")

    forbidden = _planned_forbidden_outputs(ctx)
    production_created = any("/production/" in item or "\\production\\" in item for item in forbidden)
    broad_created = any("broad_anchor_scan_optional" in item for item in forbidden)
    for check_id, ok, detail in [
        ("no_production_outputs_created", not production_created, "production output files were not created by generator"),
        ("no_broad_docking_created", not broad_created, "broad docking directory not created"),
        ("no_pose_collection_attempted", True, "pose collection not attempted"),
        ("no_pose_attribution_attempted", True, "pose attribution not attempted"),
        ("no_candidate_ranking_attempted", True, "candidate ranking not attempted"),
        ("old_workflow_not_used", True, "old Workflow A/B outputs not used as source-of-truth"),
        ("non_goals_preserved", True, "no qsub submission, Vina execution, output collection, pose attribution, scoring, or candidate claims"),
    ]:
        _append_qc(qc_rows, check_id, "non_goal", "PASS" if ok else "FAIL", detail, "Remove forbidden outputs.", not ok)
        if not ok:
            blockers.append(check_id.replace("_", " "))

    tracked = []
    for row in jobs:
        planned_pose = _path_from_row(ctx, row["output_pdbqt_file"])
        if planned_pose and planned_pose.exists() and git_ls_files(ctx.repo_root, planned_pose):
            tracked.append(ctx.relative_to_repo(planned_pose))
    ignored = all(git_check_ignored(ctx.repo_root, _path_from_row(ctx, row["output_pdbqt_file"]) or ctx.run_dir) for row in jobs)
    _append_qc(qc_rows, "generated_outputs_gitignored", "gitignore", "PASS" if ignored or not jobs else "WARN", "planned PDBQT output ignore protection checked", "Add docking output ignore patterns.")
    _append_qc(qc_rows, "generated_outputs_not_tracked", "gitignore", "PASS" if not tracked else "FAIL", "generated/planned PDBQT git tracking checked", "Remove generated docking outputs from git index without deleting local files.", bool(tracked))
    if tracked:
        blockers.append("generated/planned docking output PDBQT files tracked by git")

    leak_count, coord_hits, smiles_logged, ligand_coords, receptor_coords, scanned_paths = _scan_hygiene(ctx, list(private_entries.values()))
    for check_id, ok, detail in [
        ("no_internal_id_leak", leak_count == 0, "internal ID leakage scan complete"),
        ("no_smiles_logged", not smiles_logged, "ligand string logging scan complete"),
        ("no_ligand_coordinates_logged", not ligand_coords, "ligand coordinate leakage scan complete"),
        ("no_receptor_coordinates_logged", not receptor_coords, "receptor coordinate leakage scan complete"),
    ]:
        _append_qc(qc_rows, check_id, "confidentiality", "PASS" if ok else "FAIL", detail, "Remove sensitive tokens/coordinates from public outputs.", not ok)
    if leak_count:
        blockers.append("internal compound ID leakage detected")
    if smiles_logged:
        blockers.append("ligand SMILES printed in public outputs/logs/manifests/PBS")
    if coord_hits:
        blockers.append("ligand/receptor/pose coordinates printed outside intended PDBQT files")

    for check_id in REQUIRED_CHECKS:
        if not any(row["check_id"] == check_id for row in qc_rows):
            _append_qc(qc_rows, check_id, "coverage", "NOT_APPLICABLE", "check not reached")

    if mode == "dry-run" and not blockers:
        warnings.append("dry-run mode; PBS plan written for inspection only")
    if allow_reference_jobs:
        warnings.append("allow_reference_jobs used; reference jobs do not open production submission")

    status = "FAIL" if blockers else ("WARN" if warnings or mode == "dry-run" or production_bypass or reference_only else "PASS")
    qsub_allowed = status == "PASS" and mode == "generate" and not allow_independent and not reference_only
    production_allowed = (
        qsub_allowed
        and profile == "production"
        and mini_pass
        and not allow_production_plan_without_mini
        and not allow_reference_jobs
    )
    if allow_independent or production_bypass or reference_only:
        qsub_allowed = False
        production_allowed = False

    counts = {
        "eligible_ligands": len(ligands),
        "eligible_receptor_states": len(receptors),
        "eligible_boxes": len(boxes),
        "planned_vina_commands": len(jobs),
        "planned_chunks": len(chunks),
        "pbs_files": len(pbs_rows),
        "qsub_commands_invoked": 0,
        "vina_commands_invoked_by_generator": 0,
        "production_outputs_created": 1 if production_created else 0,
        "broad_outputs_created": 1 if broad_created else 0,
        "confidentiality_leaks": leak_count,
        "tracked_generated_outputs": len(tracked),
    }
    selection = {
        "eligible_ligands": [row["compound_public_id"] for row in ligands],
        "eligible_states": [row["state_id"] for row in receptors],
        "eligible_pocket_families": sorted({row["pocket_family_id"] for row in boxes}),
        "eligible_boxes": [row["box_id"] for row in boxes],
        "selection_policy": "deterministic_all_eligible_for_profile",
    }
    hpc_summary = {
        "queue": hpc_values["queue"],
        "nodes": hpc_values["nodes"],
        "ppn": hpc_values["ppn"],
        "worker_count_by_chunk": worker_by_chunk,
        "cpu_per_vina": hpc_values["cpu_per_vina"],
        "thread_env_set": True,
    }
    summary = {
        "schema_version": "m3_vina_job_plan_qc_v1",
        "run_id": ctx.run_id,
        "m2_run_id": m2_run_id,
        "reviewed_at": now_iso(),
        "mode": mode,
        "profile": profile,
        "overall_status": status,
        "qsub_submission_allowed": qsub_allowed,
        "production_submission_allowed": production_allowed,
        "m3_t6_collection_ready": False,
        "allow_independent": allow_independent,
        "allow_reference_jobs": allow_reference_jobs,
        "allow_production_plan_without_mini": allow_production_plan_without_mini,
        "m3_t2": {"summary_found": t2 is not None, "overall_status": t2.get("overall_status") if t2 else None, "ligand_preparation_passed": t2_pass if t2 else None},
        "m3_t3": {"summary_found": t3 is not None, "overall_status": t3.get("overall_status") if t3 else None, "m3_t4_allowed": t3.get("m3_t4_allowed") if t3 else None},
        "m3_t4": {"summary_found": t4 is not None, "overall_status": t4.get("overall_status") if t4 else None, "m3_t5_allowed": t4.get("m3_t5_allowed") if t4 else None, "vina_smoke_passed": t4_pass if t4 else None},
        "selection": selection,
        "hpc": hpc_summary,
        "vina": {"available": vina.available, "version": vina.version, "executable": vina.executable, "generator_invoked_vina": False, "generator_invoked_qsub": False},
        "outputs": {
            "job_manifest_file": ctx.relative_to_repo(job_manifest),
            "profile_job_manifest_file": ctx.relative_to_repo(profile_manifest),
            "pbs_manifest_file": ctx.relative_to_repo(pbs_manifest),
            "qc_csv_file": ctx.relative_to_repo(qc_csv),
            "report_file": ctx.relative_to_repo(report_md),
            "qsub_dir": ctx.relative_to_repo(qsub_dir),
        },
        "counts": counts,
        "blockers": list(dict.fromkeys(blockers)),
        "warnings": list(dict.fromkeys(warnings)),
        "gitignore": {"updated": gitignore_updated, "patterns_verified": patterns_verified, "sensitive_or_generated_files_tracked": tracked},
        "confidentiality": {
            "internal_ids_redacted": True,
            "public_outputs_scanned": scanned_paths,
            "leaks_detected": leak_count,
            "smiles_logged": smiles_logged,
            "ligand_coordinates_logged": ligand_coords,
            "receptor_coordinates_logged": receptor_coords,
            "pose_coordinates_logged_outside_pose_file": bool(coord_hits),
        },
        "non_goals_preserved": {
            "qsub_submission": True,
            "vina_execution_by_generator": True,
            "production_docking_execution": True,
            "output_collection": True,
            "pose_attribution": True,
            "pose_clustering": True,
            "compound_scoring": True,
            "candidate_nomination": True,
        },
    }
    _write_csv(qc_csv, QC_FIELDS, qc_rows, ctx)
    ctx.require_within_run_dir(qc_json).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_report(report_md, ctx, summary)
    with ctx.require_within_run_dir(phase3_log).open("a", encoding="utf-8") as handle:
        handle.write(
            "{0} M3-T5 status={1} run_id={2} m2_run_id={3} mode={4} profile={5} qsub_submission_allowed={6} production_submission_allowed={7} planned_vina_commands={8} planned_chunks={9} pbs_files={10} qsub_invoked=false vina_invoked_by_generator=false no_output_collection_attempted=true no_pose_attribution_attempted=true no_scoring_attempted=true no_candidate_nomination_attempted=true\n".format(
                now_iso(),
                status,
                ctx.run_id,
                m2_run_id,
                mode,
                profile,
                str(qsub_allowed).lower(),
                str(production_allowed).lower(),
                len(jobs),
                len(chunks),
                len(pbs_rows),
            )
        )
    log_master(ctx, status, "M3-T5 Vina job plan generated without qsub submission or Vina execution")
    append_phase_status(
        ctx,
        "phase3_compounds",
        status,
        "M3-T5 Vina job manifest and PBS plan generated",
        {
            "task": "M3-T5",
            "m2_run_id": m2_run_id,
            "mode": mode,
            "profile": profile,
            "qsub_submission_allowed": qsub_allowed,
            "production_submission_allowed": production_allowed,
            "planned_vina_commands": len(jobs),
            "planned_chunks": len(chunks),
            "pbs_files": len(pbs_rows),
            "qsub_invoked": False,
            "vina_invoked_by_generator": False,
            "no_qsub_submission_attempted": True,
            "no_production_docking_executed": True,
            "no_output_collection_attempted": True,
            "no_pose_attribution_attempted": True,
            "no_scoring_attempted": True,
            "no_candidate_nomination_attempted": True,
        },
    )
    return M3VinaJobPlanResult(
        status=status,
        qsub_submission_allowed=qsub_allowed,
        production_submission_allowed=production_allowed,
        blockers=list(dict.fromkeys(blockers)),
        warnings=list(dict.fromkeys(warnings)),
        job_manifest_csv=job_manifest,
        profile_job_manifest_csv=profile_manifest,
        pbs_manifest_csv=pbs_manifest,
        qc_json=qc_json,
        qc_csv=qc_csv,
        report_md=report_md,
        qsub_dir=qsub_dir,
        phase3_log=phase3_log,
        counts=counts,
        selection=selection,
        hpc=hpc_summary,
        vina=summary["vina"],
    )
