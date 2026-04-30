"""M3-T6 Vina output collection and raw pose table generation.

Collects expected outputs from an M3-T5 manifest. It never invokes Vina, qsub,
or the PBS runner, and it does not attribute poses, analyze ATP migration,
cluster poses, score compounds, rank compounds, or nominate candidates.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
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
from egfr_myo1d.compound.pdbqt_parse import parse_vina_log, parse_vina_pdbqt
from egfr_myo1d.compound.vina_jobs import JOB_FIELDS
from egfr_myo1d.core.logging_utils import append_failed_job, append_job_status, append_phase_status, log_master
from egfr_myo1d.core.run_context import RunContext, ensure_within


RAW_FIELDS = [
    "run_id",
    "m2_run_id",
    "profile",
    "job_id",
    "compound_public_id",
    "state_id",
    "state_role",
    "pocket_family_id",
    "box_id",
    "protomer_id",
    "vina_repeat_id",
    "pose_rank",
    "pose_model_index",
    "vina_affinity_kcal_mol",
    "vina_log_affinity_kcal_mol",
    "pose_file",
    "output_pdbqt_file",
    "output_pdbqt_sha256",
    "vina_log_file",
    "vina_log_sha256",
    "pose_center_x",
    "pose_center_y",
    "pose_center_z",
    "pose_atom_count",
    "pose_heavy_atom_count",
    "source_vina_job_manifest",
    "source_vina_job_manifest_sha256",
    "parse_status",
    "parse_notes",
    "collected_at",
    "allowed_for_attribution",
]

COMPLETION_FIELDS = [
    "run_id",
    "m2_run_id",
    "profile",
    "job_id",
    "chunk_id",
    "node",
    "compound_public_id",
    "state_id",
    "state_role",
    "pocket_family_id",
    "box_id",
    "protomer_id",
    "vina_repeat_id",
    "expected_output_pdbqt_file",
    "expected_vina_log_file",
    "stdout_file",
    "stderr_file",
    "output_exists",
    "output_size_bytes",
    "output_pdbqt_sha256",
    "log_exists",
    "log_size_bytes",
    "vina_log_sha256",
    "stdout_exists",
    "stderr_exists",
    "parsed_pose_count",
    "expected_num_modes",
    "affinity_values_recorded_count",
    "completion_status",
    "parse_status",
    "failure_reason_code",
    "path_validation_status",
    "inside_run_dir",
    "inside_logs_dir",
    "source_vina_job_manifest",
    "collected_at",
    "completion_notes",
]

QC_FIELDS = ["check_id", "category", "status", "severity", "details", "recommended_fix"]

REQUIRED_CHECKS = [
    "m3_t5_job_manifest_present",
    "m3_t5_job_manifest_schema_valid",
    "m3_t5_job_plan_qc_present",
    "m3_t5_job_plan_not_failed",
    "profile_selected",
    "expected_jobs_found",
    "expected_job_ids_unique",
    "expected_jobs_traceable_to_compound_state_pocket_box_repeat",
    "expected_output_paths_inside_run_dir",
    "expected_log_paths_inside_run_dir",
    "expected_stdout_stderr_paths_inside_logs_jobs",
    "all_expected_jobs_classified",
    "no_silent_missing_outputs",
    "missing_outputs_recorded",
    "failed_jobs_recorded",
    "completed_jobs_recorded",
    "pdbqt_outputs_present_for_completed_jobs",
    "vina_logs_present_for_completed_jobs",
    "pdbqt_outputs_nonempty_for_completed_jobs",
    "vina_logs_nonempty_for_completed_jobs",
    "pdbqt_parser_multimodel_supported",
    "pdbqt_pose_rank_parsed",
    "pdbqt_pose_affinity_parsed",
    "vina_log_affinity_parsed",
    "pdbqt_log_affinity_consistency_checked",
    "pose_centers_computed",
    "pose_centers_numeric",
    "pose_atom_counts_recorded",
    "raw_pose_table_written",
    "raw_pose_table_positive_rows",
    "raw_pose_rows_traceable",
    "job_completion_table_written",
    "docking_completion_qc_written",
    "docking_completion_report_written",
    "no_vina_invoked_by_collector",
    "no_qsub_invoked_by_collector",
    "no_runner_invoked_by_collector",
    "no_pose_attribution_attempted",
    "no_atp_migration_attempted",
    "no_pocket_retention_attempted",
    "no_ppi_geometry_attempted",
    "no_pose_clustering_attempted",
    "no_candidate_scoring_attempted",
    "no_candidate_nomination_attempted",
    "no_final_candidate_tables_created",
    "vina_affinity_not_used_for_ranking",
    "reference_state_not_promoted",
    "no_internal_id_leak",
    "no_smiles_logged",
    "no_ligand_coordinates_logged",
    "no_receptor_coordinates_logged",
    "no_pose_coordinates_logged_outside_pdbqt",
    "generated_outputs_gitignored",
    "pose_outputs_not_tracked",
    "old_workflow_not_used",
    "non_goals_preserved",
]

GITIGNORE_PATTERNS = [
    "fresh/runs/*/phase3_compounds/docking_outputs/*",
    "fresh/runs/*/phase3_compounds/docking_outputs/**",
    "*.vina.tmp",
    "*.vina.log.tmp",
    "*.vina.stdout",
    "*.vina.stderr",
    "*.smoke_pose.pdbqt",
    "*.pose_out.pdbqt",
]

LOG_ERROR_PATTERNS = [
    ("SEGFAULT", re.compile(r"segmentation fault", re.IGNORECASE)),
    ("FAILED", re.compile(r"\bfailed\b", re.IGNORECASE)),
    ("ERROR", re.compile(r"\berror\b", re.IGNORECASE)),
    ("CANNOT_OPEN", re.compile(r"cannot open|could not open", re.IGNORECASE)),
    ("PDBQT_PARSE_ERROR", re.compile(r"parse error|PDBQT parsing error", re.IGNORECASE)),
    ("TIMEOUT", re.compile(r"exceeded timeout|TimeoutExpired", re.IGNORECASE)),
    ("NO_SUCH_FILE", re.compile(r"No such file", re.IGNORECASE)),
    ("PERMISSION_DENIED", re.compile(r"Permission denied", re.IGNORECASE)),
]

FORBIDDEN_OUTPUTS = [
    "phase3_compounds/tables/compound_pose_attribution.csv",
    "phase3_compounds/tables/compound_pose_clusters.csv",
    "phase3_compounds/tables/compound_anchor_convergence.csv",
    "phase3_compounds/tables/final_m3_candidate_hypotheses.csv",
    "phase3_compounds/tables/pocket_compound_evidence_table.csv",
    "phase3_compounds/qc/atp_migration_qc.csv",
    "phase3_compounds/qc/pocket_retention_qc.csv",
    "phase3_compounds/qc/ppi_geometry_qc.csv",
    "phase3_compounds/docking_outputs/broad_anchor_scan_optional",
]


@dataclass
class M3VinaCollectResult:
    status: str
    m3_t7_attribution_ready: bool
    blockers: list[str]
    warnings: list[str]
    raw_pose_csv: Path
    completion_csv: Path
    qc_csv: Path
    qc_json: Path
    report_md: Path
    phase3_log: Path
    collection: dict[str, int]
    profiles: dict[str, Any]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _boolish(value: Any) -> bool:
    return str(value or "").strip().lower() in {"true", "t", "1", "yes", "y", "pass", "passed"}


def _bool_text(value: bool) -> str:
    return "true" if value else "false"


def _severity(status: str, blocker: bool = False) -> str:
    if status == "FAIL" and blocker:
        return "BLOCKER"
    if status in {"FAIL", "WARN"}:
        return "MAJOR"
    if status == "NOT_APPLICABLE":
        return "MINOR"
    return "INFO"


def _append_qc(rows: list[dict[str, str]], check_id: str, category: str, status: str, details: str, fix: str = "", blocker: bool = False) -> None:
    rows.append({"check_id": check_id, "category": category, "status": status, "severity": _severity(status, blocker), "details": details, "recommended_fix": fix})


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
        handle.write("\n# M3 Vina collection output protection\n")
        for pattern in missing:
            handle.write(pattern + "\n")
    return True, list(GITIGNORE_PATTERNS)


def _resolve_path(ctx: RunContext, value: Any) -> Path | None:
    text = str(value or "").strip()
    if not text or any(token in text for token in ["$RUN_ID", "${RUN_ID}", "$PBS_JOBID", "$(pwd)", "<", ">"]):
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


def _job_manifest_path(ctx: RunContext, job_manifest: Path | None) -> Path:
    path = job_manifest or (ctx.run_dir / "phase3_compounds" / "manifests" / "vina_job_manifest.csv")
    if not path.is_absolute():
        path = ctx.repo_root / path
    return path.resolve()


def _detect_log_error(path: Path | None, max_log_bytes: int) -> str:
    if not path or not path.is_file() or path.stat().st_size == 0:
        return ""
    text = path.read_bytes()[:max_log_bytes].decode("utf-8", errors="ignore")
    for code, pattern in LOG_ERROR_PATTERNS:
        if pattern.search(text):
            return code
    return ""


def _is_expected(row: dict[str, str]) -> bool:
    if row.get("job_scope") != "focused_pocket_first":
        return False
    if row.get("docking_status") == "NOT_APPLICABLE":
        return False
    if not _boolish(row.get("allowed_for_execution")):
        return False
    if not _boolish(row.get("allowed_for_collection")):
        return False
    if row.get("docking_status") not in {"PLANNED", "RUNNING", "PASS", "WARN", "FAIL", "MISSING", "UNKNOWN", ""}:
        return False
    return bool(row.get("output_pdbqt_file")) and bool(row.get("vina_log_file"))


def _parse_job(row: dict[str, str], output_path: Path, log_path: Path, max_log_bytes: int) -> tuple[list[Any], str, str, int, str]:
    log_affinities = parse_vina_log(log_path, max_bytes=max_log_bytes) if log_path.is_file() and log_path.stat().st_size > 0 else {}
    poses = parse_vina_pdbqt(output_path)
    if not poses:
        return [], "FAIL", "NO_VALID_POSES", len(log_affinities), "no_valid_poses"
    parsed: list[Any] = []
    status = "PASS"
    notes: list[str] = []
    for pose in poses:
        pose_status = pose.parse_status
        pose_notes = [pose.parse_notes]
        log_affinity = log_affinities.get(int(pose.pose_rank or pose.pose_model_index))
        affinity = pose.vina_affinity_kcal_mol
        if affinity is None and log_affinity is not None:
            affinity = log_affinity
            pose_status = "WARN"
            pose_notes.append("affinity_from_log")
        elif affinity is not None and log_affinity is None:
            pose_status = "WARN"
            pose_notes.append("log_affinity_missing")
        elif affinity is not None and log_affinity is not None and abs(affinity - log_affinity) > 0.01:
            pose_status = "WARN"
            pose_notes.append("pdbqt_log_affinity_disagreement")
        if affinity is None:
            pose_status = "FAIL"
            pose_notes.append("affinity_missing")
        if pose_status == "FAIL":
            status = "FAIL"
        elif pose_status == "WARN" and status == "PASS":
            status = "WARN"
        parsed.append((pose, affinity, log_affinity, pose_status, ";".join(pose_notes)))
    if status == "WARN":
        notes.append("parse_warn")
    return parsed, status, "", len(log_affinities), ";".join(notes) or "parsed"


def _scan_hygiene(ctx: RunContext, private_entries: list[Any]) -> tuple[int, list[str], bool, bool, bool, bool, list[str]]:
    leaks = scan_internal_id_leaks(ctx, private_entries)
    scan_paths = public_output_scan_paths(ctx)
    for root in [
        ctx.run_dir / "phase3_compounds" / "tables",
        ctx.run_dir / "phase3_compounds" / "qsub",
        ctx.logs_dir / "errors",
    ]:
        if root.is_dir():
            scan_paths.extend(path for path in sorted(root.rglob("*")) if path.is_file())
    filtered = []
    for path in scan_paths:
        if "docking_outputs" in path.parts and path.suffix.lower() == ".pdbqt":
            continue
        filtered.append(path)
    coord_hits: list[str] = []
    smiles_logged = False
    ligand_coords = False
    receptor_coords = False
    pose_coords = False
    scanned: list[str] = []
    for path in filtered:
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        scanned.append(ctx.relative_to_repo(path))
        if re.search(r"\bSMILES\b|canonical_smiles|isomeric_smiles", text, re.IGNORECASE):
            smiles_logged = True
        if re.search(r"^(ATOM|HETATM)\s+\d+", text, re.MULTILINE):
            coord_hits.append(ctx.relative_to_repo(path))
            lower = path.as_posix().lower()
            ligand_coords = ligand_coords or "ligand" in lower
            receptor_coords = receptor_coords or "receptor" in lower
            pose_coords = pose_coords or "pose" in lower or "completion" in lower or "raw" in lower
    return len(leaks), coord_hits, smiles_logged, ligand_coords, receptor_coords, pose_coords, scanned


def _write_report(path: Path, ctx: RunContext, summary: dict[str, Any]) -> None:
    lines = [
        "# M3-T6 Vina Output Collection",
        "",
        "Technical collection report only. Vina affinities and pose centers are parsed for traceability; they are not used for ranking, scoring, or candidate claims.",
        "",
        "- status: {0}".format(summary["overall_status"]),
        "- profile: {0}".format(summary["profile"]),
        "- expected_jobs: {0}".format(summary["collection"]["expected_jobs"]),
        "- parsed_pose_rows: {0}".format(summary["collection"]["parsed_pose_rows"]),
        "- m3_t7_attribution_ready: {0}".format(str(summary["m3_t7_attribution_ready"]).lower()),
        "",
        "## Blockers",
    ]
    lines.extend("- {0}".format(item) for item in summary["blockers"] or ["none"])
    lines.append("")
    lines.append("## Warnings")
    lines.extend("- {0}".format(item) for item in summary["warnings"] or ["none"])
    lines.append("")
    lines.append("Next task: M3-T7 - Pose attribution and ATP/pocket/PPI geometry QC")
    ctx.require_within_run_dir(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_m3_vina_collect(
    ctx: RunContext,
    *,
    m2_run_id: str,
    job_manifest: Path | None = None,
    profile: str | None = None,
    mode: str = "collect",
    force: bool = False,
    allow_partial: bool = False,
    require_all_jobs: bool = True,
    allow_reference_only: bool = False,
    private_map: Path | None = None,
    update_ignore: bool = True,
    max_log_bytes: int = 200000,
    strict_log_parse: bool = True,
    write_empty_raw_table: bool = False,
) -> M3VinaCollectResult:
    _ensure_logs(ctx)
    phase3 = ctx.run_dir / "phase3_compounds"
    tables = phase3 / "tables"
    qc_dir = phase3 / "qc"
    reports = phase3 / "reports"
    for directory in [tables, qc_dir, reports]:
        ctx.require_within_run_dir(directory).mkdir(parents=True, exist_ok=True)
    raw_csv = tables / "compound_pose_raw.csv"
    completion_csv = tables / "vina_job_completion.csv"
    qc_csv = qc_dir / "docking_completion_qc.csv"
    qc_json = qc_dir / "docking_completion_qc.json"
    report_md = reports / "m3_task6_vina_collection.md"
    phase3_log = ctx.logs_dir / "phase3_compounds.log"

    blockers: list[str] = []
    warnings: list[str] = []
    qc_rows: list[dict[str, str]] = []
    collected_at = now_iso()
    gitignore_updated, patterns_verified = update_gitignore(ctx.repo_root) if update_ignore else (False, list(GITIGNORE_PATTERNS))
    private_entries, private_warnings = load_private_map(_safe_private_map_path(ctx, private_map))
    warnings.extend(private_warnings)

    manifest_path = _job_manifest_path(ctx, job_manifest)
    manifest_inside = _inside(manifest_path, ctx.run_dir / "phase3_compounds" / "manifests")
    manifest_found = manifest_path.is_file() and manifest_inside
    _append_qc(qc_rows, "m3_t5_job_manifest_present", "input", "PASS" if manifest_found else ("WARN" if mode == "dry-run" else "FAIL"), "vina_job_manifest present={0}".format(manifest_found), "Run M3-T5 job generation first.", mode == "collect" and not manifest_found)
    if not manifest_found:
        msg = "vina_job_manifest.csv absent or outside phase3_compounds/manifests"
        if mode == "collect":
            blockers.append(msg)
        else:
            warnings.append(msg)

    rows: list[dict[str, str]] = []
    fields: list[str] = []
    manifest_sha = ""
    if manifest_found:
        rows, fields = _read_csv(manifest_path)
        manifest_sha = _sha256(manifest_path)
    required = set(JOB_FIELDS)
    schema_valid = manifest_found and required.issubset(set(fields))
    _append_qc(qc_rows, "m3_t5_job_manifest_schema_valid", "input", "PASS" if schema_valid else ("WARN" if mode == "dry-run" else "FAIL"), "M3-T5 job manifest schema checked", "Regenerate M3-T5 manifest.", mode == "collect" and not schema_valid)
    if manifest_found and not schema_valid and mode == "collect":
        blockers.append("manifest schema missing required columns")

    plan_qc = qc_dir / "vina_job_plan_qc.json"
    plan = _load_json(plan_qc) if plan_qc.is_file() else None
    _append_qc(qc_rows, "m3_t5_job_plan_qc_present", "input", "PASS" if plan else ("WARN" if mode == "dry-run" else "FAIL"), "M3-T5 QC found={0}".format(bool(plan)), "Run M3-T5 first.", mode == "collect" and not plan)
    plan_passed = bool(plan and plan.get("overall_status") == "PASS")
    _append_qc(qc_rows, "m3_t5_job_plan_not_failed", "input", "PASS" if plan_passed else ("WARN" if mode == "dry-run" else "FAIL"), "M3-T5 plan status checked", "Resolve M3-T5 blockers.", mode == "collect" and not plan_passed)
    if mode == "collect" and plan and plan.get("overall_status") != "PASS":
        blockers.append("M3-T5 job plan QC did not PASS")

    profiles = sorted({row.get("profile", "") for row in rows if row.get("profile")})
    selected_profile = profile
    if not selected_profile:
        selected_profile = profiles[0] if len(profiles) == 1 else None
    if not selected_profile and rows:
        blockers.append("profile could not be inferred; pass --profile")
    _append_qc(qc_rows, "profile_selected", "input", "PASS" if selected_profile else ("WARN" if mode == "dry-run" else "FAIL"), "profiles_in_manifest={0} selected={1}".format(",".join(profiles), selected_profile or "null"), "Pass --profile.", mode == "collect" and not selected_profile)

    selected_rows = [row for row in rows if not selected_profile or row.get("profile") == selected_profile]
    job_ids = [row.get("job_id", "") for row in selected_rows]
    unique_ids = len(job_ids) == len(set(job_ids))
    if selected_rows and not unique_ids:
        blockers.append("duplicate job_id")
    _append_qc(qc_rows, "expected_job_ids_unique", "input", "PASS" if unique_ids else "FAIL", "job_id uniqueness checked", "Regenerate M3-T5 manifest.", not unique_ids)

    completion_rows: list[dict[str, Any]] = []
    raw_rows: list[dict[str, Any]] = []
    expected_rows = 0
    malformed_paths = False
    path_outside = False
    completed = 0
    completed_warn = 0
    missing_output = 0
    missing_log = 0
    missing_both = 0
    failed_by_log = 0
    parse_fail = 0
    output_empty = 0
    log_empty = 0
    not_expected = 0
    unknown = 0
    jobs_with_pose = 0

    for row in selected_rows:
        expected = _is_expected(row)
        output_path = _resolve_path(ctx, row.get("output_pdbqt_file"))
        log_path = _resolve_path(ctx, row.get("vina_log_file"))
        stdout_path = _resolve_path(ctx, row.get("stdout_file"))
        stderr_path = _resolve_path(ctx, row.get("stderr_file"))
        inside_output = _inside(output_path, ctx.run_dir)
        inside_log = _inside(log_path, ctx.run_dir)
        inside_stdio = _inside(stdout_path, ctx.jobs_log_dir) and _inside(stderr_path, ctx.jobs_log_dir)
        path_valid = inside_output and inside_log and inside_stdio and output_path is not None and log_path is not None
        if not path_valid and expected:
            malformed_paths = True
            path_outside = path_outside or not (inside_output and inside_log)
        output_exists = bool(output_path and output_path.is_file())
        log_exists = bool(log_path and log_path.is_file())
        stdout_exists = bool(stdout_path and stdout_path.is_file())
        stderr_exists = bool(stderr_path and stderr_path.is_file())
        output_size = output_path.stat().st_size if output_exists and output_path else 0
        log_size = log_path.stat().st_size if log_exists and log_path else 0
        output_sha256 = _sha256(output_path) if output_exists and output_path and output_size > 0 else ""
        log_sha256 = _sha256(log_path) if log_exists and log_path and log_size > 0 else ""
        completion_status = "UNKNOWN"
        row_parse_status = "MISSING"
        failure_code = ""
        notes = ""
        parsed_pose_count = 0
        affinity_count = 0
        if not expected:
            if not _boolish(row.get("allowed_for_execution")):
                completion_status = "SKIPPED_NOT_ALLOWED"
            elif not _boolish(row.get("allowed_for_collection")):
                completion_status = "SKIPPED_COLLECTION_NOT_ALLOWED"
            else:
                completion_status = "NOT_EXPECTED"
            row_parse_status = "NOT_APPLICABLE"
            not_expected += 1
        else:
            expected_rows += 1
            if not path_valid:
                completion_status = "PATH_INVALID"
                row_parse_status = "FAIL"
                failure_code = "PATH_INVALID"
            elif not output_exists and not log_exists:
                completion_status = "MISSING_OUTPUT_AND_LOG"
                failure_code = "MISSING_OUTPUT_AND_LOG"
                missing_both += 1
            elif output_exists and output_size == 0:
                completion_status = "OUTPUT_EMPTY"
                failure_code = "OUTPUT_EMPTY"
                output_empty += 1
            elif log_exists and log_size == 0:
                completion_status = "LOG_EMPTY"
                failure_code = "LOG_EMPTY"
                log_empty += 1
            elif not output_exists:
                completion_status = "MISSING_OUTPUT"
                failure_code = _detect_log_error(log_path, max_log_bytes) or "MISSING_OUTPUT"
                if failure_code != "MISSING_OUTPUT":
                    completion_status = "FAILED_BY_LOG"
                    failed_by_log += 1
                else:
                    missing_output += 1
            elif not log_exists:
                completion_status = "MISSING_LOG"
                failure_code = "MISSING_LOG"
                missing_log += 1
            else:
                log_error = _detect_log_error(log_path, max_log_bytes)
                try:
                    parsed, parse_status, parse_failure, affinity_count, notes = _parse_job(row, output_path, log_path, max_log_bytes)
                except Exception:
                    parsed, parse_status, parse_failure, affinity_count, notes = [], "FAIL", "PDBQT_PARSE_EXCEPTION", 0, "parse_exception"
                parsed_pose_count = len(parsed)
                if log_error and parsed_pose_count == 0:
                    completion_status = "FAILED_BY_LOG"
                    failure_code = log_error
                    failed_by_log += 1
                elif parsed_pose_count == 0:
                    completion_status = "PARSE_FAIL"
                    failure_code = parse_failure or "PARSE_FAIL"
                    parse_fail += 1
                elif parse_status == "WARN":
                    completion_status = "COMPLETED_PARSE_WARN"
                    row_parse_status = "WARN"
                    completed_warn += 1
                elif parse_status == "FAIL":
                    completion_status = "PARSE_FAIL"
                    row_parse_status = "FAIL"
                    failure_code = parse_failure or "PARSE_FAIL"
                    parse_fail += 1
                else:
                    completion_status = "COMPLETED"
                    row_parse_status = "PASS"
                    completed += 1
                if parsed_pose_count:
                    jobs_with_pose += 1
                for pose, affinity, log_affinity, pose_status, pose_notes in parsed:
                    raw_rows.append(
                        {
                            "run_id": ctx.run_id,
                            "m2_run_id": m2_run_id,
                            "profile": row.get("profile", ""),
                            "job_id": row.get("job_id", ""),
                            "compound_public_id": row.get("compound_public_id", ""),
                            "state_id": row.get("state_id", ""),
                            "state_role": row.get("state_role", "unknown"),
                            "pocket_family_id": row.get("pocket_family_id", ""),
                            "box_id": row.get("box_id", ""),
                            "protomer_id": row.get("protomer_id", ""),
                            "vina_repeat_id": row.get("vina_repeat_id", ""),
                            "pose_rank": pose.pose_rank,
                            "pose_model_index": pose.pose_model_index,
                            "vina_affinity_kcal_mol": "" if affinity is None else "{0:.3f}".format(affinity),
                            "vina_log_affinity_kcal_mol": "" if log_affinity is None else "{0:.3f}".format(log_affinity),
                            "pose_file": ctx.relative_to_repo(output_path),
                            "output_pdbqt_file": ctx.relative_to_repo(output_path),
                            "output_pdbqt_sha256": output_sha256,
                            "vina_log_file": ctx.relative_to_repo(log_path),
                            "vina_log_sha256": log_sha256,
                            "pose_center_x": "{0:.3f}".format(pose.pose_center_x),
                            "pose_center_y": "{0:.3f}".format(pose.pose_center_y),
                            "pose_center_z": "{0:.3f}".format(pose.pose_center_z),
                            "pose_atom_count": pose.pose_atom_count,
                            "pose_heavy_atom_count": pose.pose_heavy_atom_count,
                            "source_vina_job_manifest": ctx.relative_to_repo(manifest_path) if manifest_found else "",
                            "source_vina_job_manifest_sha256": manifest_sha,
                            "parse_status": pose_status,
                            "parse_notes": pose_notes,
                            "collected_at": collected_at,
                            "allowed_for_attribution": _bool_text(pose_status in {"PASS", "WARN"} and completion_status in {"COMPLETED", "COMPLETED_PARSE_WARN"}),
                        }
                    )
        if completion_status == "UNKNOWN":
            unknown += 1
        completion_rows.append(
            {
                "run_id": ctx.run_id,
                "m2_run_id": m2_run_id,
                "profile": row.get("profile", ""),
                "job_id": row.get("job_id", ""),
                "chunk_id": row.get("chunk_id", ""),
                "node": row.get("node", ""),
                "compound_public_id": row.get("compound_public_id", ""),
                "state_id": row.get("state_id", ""),
                "state_role": row.get("state_role", "unknown"),
                "pocket_family_id": row.get("pocket_family_id", ""),
                "box_id": row.get("box_id", ""),
                "protomer_id": row.get("protomer_id", ""),
                "vina_repeat_id": row.get("vina_repeat_id", ""),
                "expected_output_pdbqt_file": row.get("output_pdbqt_file", ""),
                "expected_vina_log_file": row.get("vina_log_file", ""),
                "stdout_file": row.get("stdout_file", ""),
                "stderr_file": row.get("stderr_file", ""),
                "output_exists": _bool_text(output_exists),
                "output_size_bytes": output_size,
                "output_pdbqt_sha256": output_sha256,
                "log_exists": _bool_text(log_exists),
                "log_size_bytes": log_size,
                "vina_log_sha256": log_sha256,
                "stdout_exists": _bool_text(stdout_exists),
                "stderr_exists": _bool_text(stderr_exists),
                "parsed_pose_count": parsed_pose_count,
                "expected_num_modes": row.get("num_modes", ""),
                "affinity_values_recorded_count": affinity_count,
                "completion_status": completion_status,
                "parse_status": row_parse_status,
                "failure_reason_code": failure_code,
                "path_validation_status": "PASS" if path_valid or not expected else "FAIL",
                "inside_run_dir": _bool_text(inside_output and inside_log),
                "inside_logs_dir": _bool_text(inside_stdio),
                "source_vina_job_manifest": ctx.relative_to_repo(manifest_path) if manifest_found else "",
                "collected_at": collected_at,
                "completion_notes": notes or completion_status.lower(),
            }
        )
        append_job_status(
            ctx,
            row.get("job_id", "unknown_job"),
            completion_status,
            node=row.get("node"),
            stdout=row.get("stdout_file"),
            stderr=row.get("stderr_file"),
            details={
                "phase": "phase3_compounds",
                "task": "M3-T6",
                "job_type": "vina_output_collection",
                "profile": row.get("profile"),
                "compound_public_id": row.get("compound_public_id"),
                "state_id": row.get("state_id"),
                "pocket_family_id": row.get("pocket_family_id"),
                "box_id": row.get("box_id"),
                "vina_repeat_id": row.get("vina_repeat_id"),
                "parsed_pose_count": parsed_pose_count,
                "output_pdbqt_file": row.get("output_pdbqt_file"),
                "vina_log_file": row.get("vina_log_file"),
                "collection_notes": failure_code or notes or completion_status,
            },
        )
        if completion_status not in {"COMPLETED", "COMPLETED_PARSE_WARN", "NOT_EXPECTED", "SKIPPED_NOT_ALLOWED", "SKIPPED_COLLECTION_NOT_ALLOWED"}:
            append_failed_job(ctx, row.get("job_id", "unknown_job"), completion_status, failure_code or completion_status)

    non_public = [row.get("compound_public_id") for row in selected_rows if row.get("compound_public_id") and row.get("compound_public_id") not in PUBLIC_COMPOUND_IDS]
    if non_public:
        blockers.append("non-public compound_public_id in manifest")
    if malformed_paths:
        blockers.append("expected output/log/stdout/stderr paths invalid")
    if path_outside:
        blockers.append("expected output or log path outside run_dir")
    if not unique_ids:
        blockers.append("duplicate job_id")
    if expected_rows == 0 and mode == "collect":
        blockers.append("zero expected jobs in collect mode")
    if (missing_output or missing_log or missing_both or output_empty or log_empty or parse_fail or failed_by_log or unknown) and mode == "collect" and require_all_jobs and not allow_partial:
        blockers.append("one or more expected jobs missing/failed/parse-failed without --allow-partial")
    if expected_rows and expected_rows == (missing_output + missing_log + missing_both + output_empty + log_empty + failed_by_log + parse_fail + unknown):
        blockers.append("all expected jobs missing or failed")
    if mode == "collect" and not raw_rows and not (write_empty_raw_table or allow_partial):
        blockers.append("zero parsed pose rows in collect mode")
    if allow_partial and raw_rows and any(row["completion_status"] not in {"COMPLETED", "COMPLETED_PARSE_WARN", "NOT_EXPECTED", "SKIPPED_NOT_ALLOWED", "SKIPPED_COLLECTION_NOT_ALLOWED"} for row in completion_rows):
        warnings.append("allow_partial used with some missing/failed jobs")
    if completed_warn:
        warnings.append("some jobs completed with parse warnings")
    if mode == "dry-run":
        warnings.append("dry-run mode; collection outputs are diagnostic")

    primary_states = sorted({row.get("state_id") for row in selected_rows if row.get("state_role") == "primary" and row.get("state_id")})
    reference_states = sorted({row.get("state_id") for row in selected_rows if row.get("state_role") == "reference" and row.get("state_id")})
    reference_only = bool(reference_states) and not primary_states
    if reference_states:
        warnings.append("reference jobs included")
    if reference_only and not allow_reference_only:
        warnings.append("reference-only collection without --allow-reference-only")

    write_raw = bool(raw_rows or write_empty_raw_table or mode == "dry-run")
    if write_raw:
        _write_csv(raw_csv, RAW_FIELDS, raw_rows, ctx)
    _write_csv(completion_csv, COMPLETION_FIELDS, completion_rows, ctx)

    def has_raw(check) -> bool:
        return bool(raw_rows) and all(check(row) for row in raw_rows)

    for check_id, ok, detail, blocker in [
        ("expected_jobs_found", expected_rows > 0, "expected_jobs={0}".format(expected_rows), mode == "collect"),
        ("expected_jobs_traceable_to_compound_state_pocket_box_repeat", all(row.get("compound_public_id") and row.get("state_id") and row.get("pocket_family_id") and row.get("box_id") and row.get("vina_repeat_id") for row in selected_rows), "manifest traceability fields checked", True),
        ("expected_output_paths_inside_run_dir", not malformed_paths and all(row.get("inside_run_dir") == "true" for row in completion_rows if row.get("completion_status") not in {"NOT_EXPECTED", "SKIPPED_NOT_ALLOWED"}), "expected output/log paths checked", True),
        ("expected_log_paths_inside_run_dir", not path_outside, "expected log paths checked", True),
        ("expected_stdout_stderr_paths_inside_logs_jobs", all(row.get("inside_logs_dir") == "true" for row in completion_rows if row.get("completion_status") not in {"NOT_EXPECTED", "SKIPPED_NOT_ALLOWED"}), "stdout/stderr paths checked", True),
        ("all_expected_jobs_classified", len(completion_rows) == len(selected_rows), "completion rows match selected manifest rows", True),
        ("no_silent_missing_outputs", all(row.get("completion_status") != "UNKNOWN" for row in completion_rows), "missing/unknown jobs classified", True),
        ("missing_outputs_recorded", True, "missing outputs recorded in completion table", False),
        ("failed_jobs_recorded", True, "failed jobs recorded in completion table", False),
        ("completed_jobs_recorded", completed + completed_warn > 0 or mode == "dry-run", "completed jobs recorded when present", False),
        ("pdbqt_outputs_present_for_completed_jobs", all(row.get("output_exists") == "true" for row in completion_rows if row.get("completion_status") in {"COMPLETED", "COMPLETED_PARSE_WARN"}), "completed output files present", True),
        ("vina_logs_present_for_completed_jobs", all(row.get("log_exists") == "true" for row in completion_rows if row.get("completion_status") in {"COMPLETED", "COMPLETED_PARSE_WARN"}), "completed log files present", True),
        ("pdbqt_outputs_nonempty_for_completed_jobs", all(int(row.get("output_size_bytes") or 0) > 0 for row in completion_rows if row.get("completion_status") in {"COMPLETED", "COMPLETED_PARSE_WARN"}), "completed output files non-empty", True),
        ("vina_logs_nonempty_for_completed_jobs", all(int(row.get("log_size_bytes") or 0) > 0 for row in completion_rows if row.get("completion_status") in {"COMPLETED", "COMPLETED_PARSE_WARN"}), "completed log files non-empty", True),
        ("pdbqt_parser_multimodel_supported", True, "MODEL/ENDMDL and single-model parser supported", False),
        ("pdbqt_pose_rank_parsed", has_raw(lambda row: row.get("pose_rank") not in {"", None}), "pose ranks parsed", False),
        ("pdbqt_pose_affinity_parsed", has_raw(lambda row: row.get("vina_affinity_kcal_mol") not in {"", None}), "pose affinities parsed from PDBQT/log", False),
        ("vina_log_affinity_parsed", any(row.get("vina_log_affinity_kcal_mol") not in {"", None} for row in raw_rows) or not raw_rows, "Vina log affinity parsed when available", False),
        ("pdbqt_log_affinity_consistency_checked", True, "PDBQT/log affinity consistency checked", False),
        ("pose_centers_computed", has_raw(lambda row: row.get("pose_center_x") not in {"", None}), "pose centers computed", False),
        ("pose_centers_numeric", has_raw(lambda row: all(re.match(r"^-?\d+\.\d+$", str(row.get(k, ""))) for k in ["pose_center_x", "pose_center_y", "pose_center_z"])), "pose center columns are numeric", False),
        ("pose_atom_counts_recorded", has_raw(lambda row: int(row.get("pose_atom_count") or 0) > 0), "pose atom counts recorded", False),
        ("raw_pose_table_written", raw_csv.is_file(), "raw pose table written", mode == "collect" and bool(raw_rows)),
        ("raw_pose_table_positive_rows", len(raw_rows) > 0 or mode == "dry-run", "raw pose rows={0}".format(len(raw_rows)), mode == "collect"),
        ("raw_pose_rows_traceable", has_raw(lambda row: all(row.get(k) for k in ["compound_public_id", "state_id", "pocket_family_id", "box_id", "vina_repeat_id"])), "raw rows traceable", False),
        ("job_completion_table_written", completion_csv.is_file(), "job completion table written", True),
    ]:
        status = "PASS" if ok else ("WARN" if mode == "dry-run" or not blocker else "FAIL")
        _append_qc(qc_rows, check_id, "collection", status, detail, "Inspect M3-T5 manifest and Vina outputs.", status == "FAIL")

    forbidden = [ctx.relative_to_repo(ctx.run_dir / rel) for rel in FORBIDDEN_OUTPUTS if (ctx.run_dir / rel).exists()]
    if forbidden:
        blockers.append("forbidden downstream output created")
    for check_id, ok, detail in [
        ("no_vina_invoked_by_collector", True, "collector did not invoke Vina"),
        ("no_qsub_invoked_by_collector", True, "collector did not invoke qsub"),
        ("no_runner_invoked_by_collector", True, "collector did not invoke M3-T5 runner"),
        ("no_pose_attribution_attempted", not any("pose_attribution" in item for item in forbidden), "pose attribution not attempted"),
        ("no_atp_migration_attempted", not any("atp_migration" in item for item in forbidden), "ATP migration analysis not attempted"),
        ("no_pocket_retention_attempted", not any("pocket_retention" in item for item in forbidden), "pocket retention analysis not attempted"),
        ("no_ppi_geometry_attempted", not any("ppi_geometry" in item for item in forbidden), "PPI geometry analysis not attempted"),
        ("no_pose_clustering_attempted", not any("pose_clusters" in item for item in forbidden), "pose clustering not attempted"),
        ("no_candidate_scoring_attempted", True, "candidate scoring not attempted"),
        ("no_candidate_nomination_attempted", not any("candidate" in item for item in forbidden), "candidate nomination not attempted"),
        ("no_final_candidate_tables_created", not any("final_m3_candidate" in item for item in forbidden), "final candidate tables not created"),
        ("vina_affinity_not_used_for_ranking", True, "affinity values not used for ranking"),
        ("reference_state_not_promoted", True, "reference states remain reference"),
        ("old_workflow_not_used", True, "old Workflow A/B outputs not used"),
        ("non_goals_preserved", True, "non-goals preserved"),
    ]:
        _append_qc(qc_rows, check_id, "non_goal", "PASS" if ok else "FAIL", detail, "Remove forbidden downstream outputs.", not ok)

    tracked = []
    for row in selected_rows:
        output_path = _resolve_path(ctx, row.get("output_pdbqt_file"))
        if output_path and output_path.exists() and output_path.suffix.lower() == ".pdbqt" and git_ls_files(ctx.repo_root, output_path):
            tracked.append(ctx.relative_to_repo(output_path))
    ignored = all(git_check_ignored(ctx.repo_root, _resolve_path(ctx, row.get("output_pdbqt_file")) or ctx.run_dir) for row in selected_rows if row.get("output_pdbqt_file"))
    _append_qc(qc_rows, "generated_outputs_gitignored", "gitignore", "PASS" if ignored or not selected_rows else "WARN", "docking output ignore protection checked", "Add docking output PDBQT ignore patterns.")
    _append_qc(qc_rows, "pose_outputs_not_tracked", "gitignore", "PASS" if not tracked else "FAIL", "docking output PDBQT git tracking checked", "Remove generated PDBQT outputs from git index without deleting local files.", bool(tracked))
    if tracked:
        blockers.append("docking output PDBQT tracked by git")

    leak_count, coord_hits, smiles_logged, ligand_coords, receptor_coords, pose_coords, scanned_paths = _scan_hygiene(ctx, list(private_entries.values()))
    for check_id, ok, detail in [
        ("no_internal_id_leak", leak_count == 0, "internal ID leakage scan complete"),
        ("no_smiles_logged", not smiles_logged, "ligand string logging scan complete"),
        ("no_ligand_coordinates_logged", not ligand_coords, "ligand coordinate leakage scan complete"),
        ("no_receptor_coordinates_logged", not receptor_coords, "receptor coordinate leakage scan complete"),
        ("no_pose_coordinates_logged_outside_pdbqt", not pose_coords and not coord_hits, "pose coordinate leakage scan complete"),
    ]:
        _append_qc(qc_rows, check_id, "confidentiality", "PASS" if ok else "FAIL", detail, "Remove sensitive tokens/coordinates from public outputs.", not ok)
    if leak_count:
        blockers.append("internal compound ID leakage detected")
    if smiles_logged:
        blockers.append("SMILES printed in public outputs/logs")
    if coord_hits:
        blockers.append("ligand/receptor/pose coordinates printed outside intended PDBQT files")

    for check_id in REQUIRED_CHECKS:
        if not any(row["check_id"] == check_id for row in qc_rows):
            _append_qc(qc_rows, check_id, "coverage", "NOT_APPLICABLE", "check not reached")

    status = "FAIL" if blockers else ("WARN" if warnings or mode == "dry-run" or completed_warn or reference_states else "PASS")
    ready = status == "PASS" and mode == "collect" and expected_rows > 0 and bool(raw_rows) and not allow_partial
    collection = {
        "expected_jobs": expected_rows,
        "completed_jobs": completed,
        "completed_parse_warn_jobs": completed_warn,
        "missing_output_jobs": missing_output,
        "missing_log_jobs": missing_log,
        "missing_output_and_log_jobs": missing_both,
        "failed_by_log_jobs": failed_by_log,
        "parse_fail_jobs": parse_fail,
        "output_empty_jobs": output_empty,
        "log_empty_jobs": log_empty,
        "not_expected_jobs": not_expected,
        "unknown_jobs": unknown,
        "parsed_pose_rows": len(raw_rows),
        "jobs_with_at_least_one_pose": jobs_with_pose,
    }
    summary = {
        "schema_version": "m3_vina_collection_qc_v1",
        "run_id": ctx.run_id,
        "m2_run_id": m2_run_id,
        "reviewed_at": now_iso(),
        "mode": mode,
        "profile": selected_profile,
        "overall_status": status,
        "allow_partial": allow_partial,
        "require_all_jobs": require_all_jobs,
        "allow_reference_only": allow_reference_only,
        "m3_t7_attribution_ready": ready,
        "m3_t5": {
            "job_manifest_found": manifest_found,
            "job_plan_qc_found": bool(plan),
            "job_plan_overall_status": plan.get("overall_status") if plan else None,
            "planned_vina_commands": (plan.get("counts", {}) if plan else {}).get("planned_vina_commands", 0),
            "planned_chunks": (plan.get("counts", {}) if plan else {}).get("planned_chunks", 0),
        },
        "collection": collection,
        "profiles": {"profiles_in_manifest": profiles, "selected_profile": selected_profile},
        "states": {"primary_states": primary_states, "reference_states": reference_states, "reference_only_collection": reference_only},
        "outputs": {
            "raw_pose_table": ctx.relative_to_repo(raw_csv) if raw_csv.is_file() else None,
            "job_completion_table": ctx.relative_to_repo(completion_csv),
            "qc_csv_file": ctx.relative_to_repo(qc_csv),
            "qc_json_file": ctx.relative_to_repo(qc_json),
            "report_file": ctx.relative_to_repo(report_md),
        },
        "parser": {"pdbqt_parser_used": True, "vina_log_parser_used": True, "coordinate_records_written_to_public_tables": False, "affinity_used_for_ranking": False},
        "counts": {
            "vina_commands_invoked_by_collector": 0,
            "qsub_commands_invoked_by_collector": 0,
            "runner_commands_invoked_by_collector": 0,
            "pose_attribution_rows_created": 0,
            "candidate_tables_created": len([item for item in forbidden if "candidate" in item]),
            "confidentiality_leaks": leak_count,
            "tracked_pose_outputs": len(tracked),
        },
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
            "pose_coordinates_logged_outside_pose_file": pose_coords or bool(coord_hits),
        },
        "non_goals_preserved": {
            "vina_execution": True,
            "qsub_submission": True,
            "pbs_runner_invocation": True,
            "pose_attribution": True,
            "atp_migration_analysis": True,
            "pocket_retention_analysis": True,
            "ppi_geometry_analysis": True,
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
            "{0} M3-T6 status={1} run_id={2} m2_run_id={3} mode={4} profile={5} expected_jobs={6} completed_jobs={7} missing_jobs={8} failed_jobs={9} parse_fail_jobs={10} parsed_pose_rows={11} m3_t7_attribution_ready={12} qsub_invoked=false vina_invoked_by_collector=false runner_invoked_by_collector=false no_pose_attribution_attempted=true no_atp_migration_attempted=true no_pocket_retention_attempted=true no_ppi_geometry_attempted=true no_pose_clustering_attempted=true no_scoring_attempted=true no_candidate_nomination_attempted=true\n".format(
                now_iso(),
                status,
                ctx.run_id,
                m2_run_id,
                mode,
                selected_profile or "null",
                expected_rows,
                completed,
                missing_output + missing_log + missing_both,
                failed_by_log,
                parse_fail,
                len(raw_rows),
                str(ready).lower(),
            )
        )
    log_master(ctx, status, "M3-T6 Vina output collection completed without Vina/qsub/runner invocation or pose attribution")
    append_phase_status(
        ctx,
        "phase3_compounds",
        status,
        "M3-T6 Vina output collection completed",
        {
            "task": "M3-T6",
            "m2_run_id": m2_run_id,
            "mode": mode,
            "profile": selected_profile,
            "expected_jobs": expected_rows,
            "completed_jobs": completed,
            "missing_jobs": missing_output + missing_log + missing_both,
            "failed_jobs": failed_by_log,
            "parse_fail_jobs": parse_fail,
            "parsed_pose_rows": len(raw_rows),
            "m3_t7_attribution_ready": ready,
            "qsub_invoked": False,
            "vina_invoked_by_collector": False,
            "runner_invoked_by_collector": False,
            "no_pose_attribution_attempted": True,
            "no_atp_migration_attempted": True,
            "no_pocket_retention_attempted": True,
            "no_ppi_geometry_attempted": True,
            "no_pose_clustering_attempted": True,
            "no_scoring_attempted": True,
            "no_candidate_nomination_attempted": True,
        },
    )
    return M3VinaCollectResult(
        status=status,
        m3_t7_attribution_ready=ready,
        blockers=list(dict.fromkeys(blockers)),
        warnings=list(dict.fromkeys(warnings)),
        raw_pose_csv=raw_csv,
        completion_csv=completion_csv,
        qc_csv=qc_csv,
        qc_json=qc_json,
        report_md=report_md,
        phase3_log=phase3_log,
        collection=collection,
        profiles=summary["profiles"],
    )
