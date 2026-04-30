"""M3-T4 Vina smoke adapter and focused docking smoke test.

This task plans at most one focused Vina smoke command and, in run mode,
invokes exactly one Vina command. It does not generate production jobs, qsub,
mini/scaling/broad docking, pose attribution, clustering, scoring, ranking, or
candidate claims.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
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
from egfr_myo1d.compound.vina_adapter import VinaAvailability, build_vina_argv, detect_vina, run_vina_once
from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status, log_master
from egfr_myo1d.core.run_context import RunContext, ensure_within


MANIFEST_FIELDS = [
    "smoke_id",
    "run_id",
    "m2_run_id",
    "mode",
    "command_manifest_written_before_execution",
    "command_manifest_written_at",
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
    "source_docking_box_manifest",
    "source_receptor_manifest",
    "source_ligand_manifest",
    "vina_executable",
    "vina_version",
    "vina_argv_json",
    "exhaustiveness",
    "num_modes",
    "cpu",
    "seed",
    "timeout_seconds",
    "stdout_file",
    "stderr_file",
    "vina_log_file",
    "output_pdbqt_file",
    "started_at",
    "finished_at",
    "duration_seconds",
    "vina_invoked",
    "vina_invocation_count",
    "exit_code",
    "timeout",
    "output_pdbqt_exists",
    "output_pdbqt_size_bytes",
    "output_pdbqt_sha256",
    "vina_log_exists",
    "vina_log_size_bytes",
    "minimal_pose_detected",
    "affinity_table_detected",
    "allowed_for_m3_t5",
    "smoke_status",
    "smoke_notes",
]

QC_FIELDS = ["smoke_id", "check_id", "category", "status", "severity", "details", "recommended_fix"]

REQUIRED_CHECKS = [
    "m3_t2_summary_present",
    "m3_t2_passed_or_bypassed",
    "ligand_manifest_present",
    "eligible_ligand_selected",
    "ligand_pdbqt_exists",
    "ligand_pdbqt_nonempty",
    "ligand_path_inside_run_dir",
    "m3_t3_summary_present",
    "m3_t3_passed_or_bypassed",
    "m3_t4_allowed_by_t3",
    "receptor_manifest_present",
    "eligible_receptor_selected",
    "receptor_pdbqt_exists",
    "receptor_pdbqt_nonempty",
    "receptor_path_inside_run_dir",
    "docking_box_manifest_present",
    "eligible_box_selected",
    "box_traceable_to_m2",
    "box_non_atp_pass",
    "box_lower_lateral_pass",
    "box_dimer_accessibility_pass",
    "box_numeric_dimensions",
    "box_positive_dimensions",
    "vina_available",
    "vina_version_recorded",
    "vina_command_manifest_written_before_execution",
    "vina_argv_safe",
    "exactly_one_vina_command_planned",
    "exactly_one_vina_command_invoked",
    "stdout_captured",
    "stderr_captured",
    "output_path_inside_run_dir",
    "vina_output_pdbqt_exists",
    "vina_output_pdbqt_nonempty",
    "vina_log_exists",
    "minimal_pose_detected",
    "no_production_docking_launched",
    "no_qsub_launched",
    "no_broad_docking_launched",
    "no_pose_attribution_attempted",
    "no_candidate_ranking_attempted",
    "no_internal_id_leak",
    "no_smiles_logged",
    "no_ligand_coordinates_logged",
    "no_receptor_coordinates_logged",
    "generated_outputs_gitignored",
    "generated_outputs_not_tracked",
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
]

PRIMARY_STATE_ORDER = ["EGFR_160-185", "EGFR_170-200"]
FORBIDDEN_OUTPUTS = [
    "phase3_compounds/qsub",
    "phase3_compounds/docking_inputs/production",
    "phase3_compounds/docking_outputs/focused_pocket_first/production",
    "phase3_compounds/docking_outputs/broad_anchor_scan_optional",
    "phase3_compounds/vina_raw/production",
    "phase3_compounds/tables/compound_pose_raw.csv",
    "phase3_compounds/tables/compound_pose_attribution.csv",
    "phase3_compounds/tables/compound_pose_clusters.csv",
    "phase3_compounds/tables/compound_anchor_convergence.csv",
    "phase3_compounds/tables/final_m3_candidate_hypotheses.csv",
    "phase3_compounds/tables/pocket_compound_evidence_table.csv",
]


@dataclass
class SmokeSelection:
    ligand: dict[str, Any] | None
    receptor: dict[str, Any] | None
    box: dict[str, Any] | None
    smoke_id: str


@dataclass
class M3VinaSmokeResult:
    status: str
    m3_t5_allowed: bool
    blockers: list[str]
    warnings: list[str]
    manifest_csv: Path
    qc_json: Path
    qc_csv: Path
    report_md: Path
    phase3_log: Path
    counts: dict[str, int]
    selection: dict[str, Any]
    vina: dict[str, Any]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def _severity(status: str, blocker: bool = False) -> str:
    if status == "FAIL" and blocker:
        return "BLOCKER"
    if status in {"FAIL", "WARN"}:
        return "MAJOR"
    if status == "NOT_APPLICABLE":
        return "MINOR"
    return "INFO"


def _append_qc(rows: list[dict[str, str]], smoke_id: str, check_id: str, category: str, status: str, details: str, fix: str = "", blocker: bool = False) -> None:
    rows.append(
        {
            "smoke_id": smoke_id,
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
    gitignore.parent.mkdir(parents=True, exist_ok=True)
    with gitignore.open("a", encoding="utf-8") as handle:
        if existing and existing[-1].strip():
            handle.write("\n")
        handle.write("\n# M3 Vina smoke/docking-output protection\n")
        for pattern in missing:
            handle.write(pattern + "\n")
    return True, list(GITIGNORE_PATTERNS)


def _safe_private_map_path(ctx: RunContext, private_map: Path | None) -> Path:
    path = private_map or (ctx.fresh_root / "data" / "private" / "compound_id_map.csv")
    resolved = Path(path)
    if not resolved.is_absolute():
        resolved = ctx.repo_root / resolved
    return ensure_within(resolved.resolve(), ctx.fresh_root)


def _run_path(ctx: RunContext, text: str) -> Path | None:
    if not text:
        return None
    path = Path(text)
    if not path.is_absolute():
        path = ctx.repo_root / path
    try:
        return ctx.require_within_run_dir(path.resolve())
    except Exception:
        return None


def _inside_phase3(ctx: RunContext, path: Path | None) -> bool:
    if not path:
        return False
    try:
        path.resolve().relative_to((ctx.run_dir / "phase3_compounds").resolve())
        return True
    except ValueError:
        return False


def _inside_receptor_pdbqt(ctx: RunContext, path: Path | None) -> bool:
    if not path:
        return False
    try:
        path.resolve().relative_to((ctx.run_dir / "phase3_compounds" / "receptor_pdbqt").resolve())
        return True
    except ValueError:
        return False


def _default_int(ctx: RunContext, keys: list[str], default: int) -> int:
    path = ctx.fresh_root / "configs" / "gates.yaml"
    if not path.is_file():
        return default
    try:
        import yaml

        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception:
        return default
    scopes = [data, data.get("compound", {}) if isinstance(data.get("compound"), dict) else {}, data.get("m3", {}) if isinstance(data.get("m3"), dict) else {}]
    for scope in scopes:
        for key in keys:
            value = scope.get(key) if isinstance(scope, dict) else None
            if isinstance(value, int):
                return value
    return default


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


def _eligible_ligands(ctx: RunContext, manifest_path: Path | None, requested: str | None) -> tuple[list[dict[str, Any]], list[str]]:
    notes: list[str] = []
    if not manifest_path or not manifest_path.is_file():
        return [], ["ligand_preparation_manifest.csv missing"]
    rows, _ = _read_csv(manifest_path)
    eligible: list[dict[str, Any]] = []
    for row in rows:
        cid = (row.get("compound_public_id") or "").strip()
        if requested and cid != requested:
            continue
        path = _run_path(ctx, row.get("prepared_pdbqt_file", ""))
        ok = (
            cid in PUBLIC_COMPOUND_IDS
            and path is not None
            and _inside_phase3(ctx, path)
            and path.is_file()
            and path.stat().st_size > 0
            and _boolish(row.get("pdbqt_validation_success"))
            and row.get("preparation_status") == "PASS"
        )
        if ok:
            current_sha = _sha256(path)
            manifest_sha = (row.get("prepared_pdbqt_sha256") or "").strip()
            if manifest_sha and manifest_sha != current_sha:
                notes.append(f"{cid or '<missing>'}: prepared ligand PDBQT hash mismatch against M3-T2 manifest")
                continue
            enriched = dict(row)
            enriched["_path"] = path
            enriched["_sha256"] = current_sha
            eligible.append(enriched)
    order = {cid: idx for idx, cid in enumerate(PUBLIC_COMPOUND_IDS)}
    eligible.sort(key=lambda row: order.get(row["compound_public_id"], 99))
    return eligible, notes


def _eligible_receptors(ctx: RunContext, manifest_path: Path | None, requested: str | None, allow_reference: bool) -> tuple[list[dict[str, Any]], list[str]]:
    if not manifest_path or not manifest_path.is_file():
        return [], ["receptor_preparation_manifest.csv missing"]
    rows, _ = _read_csv(manifest_path)
    notes: list[str] = []
    eligible: list[dict[str, Any]] = []
    for row in rows:
        state = (row.get("state_id") or "").strip()
        role = (row.get("state_role") or "unknown").strip()
        if requested and state != requested:
            continue
        if (state == "3GT8_raw" or role == "reference") and not allow_reference:
            continue
        path = _run_path(ctx, row.get("prepared_receptor_pdbqt_file", ""))
        status_ok = row.get("preparation_status") in {"PASS", "WARN"}
        ok = (
            state
            and path is not None
            and _inside_receptor_pdbqt(ctx, path)
            and path.is_file()
            and path.stat().st_size > 0
            and _boolish(row.get("pdbqt_validation_success"))
            and _boolish(row.get("allowed_for_vina_smoke"))
            and status_ok
            and not _boolish(row.get("source_receptor_is_old_workflow"))
            and not _boolish(row.get("source_receptor_is_monomer_only"))
        )
        if ok:
            current_sha = _sha256(path)
            manifest_sha = (row.get("prepared_receptor_pdbqt_sha256") or "").strip()
            if manifest_sha and manifest_sha != current_sha:
                notes.append(f"{state or '<missing>'}: prepared receptor PDBQT hash mismatch against M3-T3 manifest")
                continue
            enriched = dict(row)
            enriched["_path"] = path
            enriched["_sha256"] = current_sha
            eligible.append(enriched)
    eligible.sort(key=_state_sort_key)
    return eligible, notes


def _eligible_boxes(
    ctx: RunContext,
    manifest_path: Path | None,
    receptors: list[dict[str, Any]],
    state_id: str | None,
    pocket_family_id: str | None,
    box_id: str | None,
    protomer_id: str | None,
) -> tuple[list[dict[str, Any]], list[str]]:
    if not manifest_path or not manifest_path.is_file():
        return [], ["docking_box_manifest.csv missing"]
    receptor_by_state = {row["state_id"]: row for row in receptors}
    rows, _ = _read_csv(manifest_path)
    eligible: list[dict[str, Any]] = []
    for row in rows:
        if state_id and row.get("state_id") != state_id:
            continue
        if pocket_family_id and row.get("pocket_family_id") != pocket_family_id:
            continue
        if box_id and row.get("box_id") != box_id:
            continue
        if protomer_id and row.get("protomer_id") != protomer_id:
            continue
        receptor = receptor_by_state.get(row.get("state_id", ""))
        sizes = [_finite_float(row.get(k)) for k in ["box_size_x", "box_size_y", "box_size_z"]]
        centers = [_finite_float(row.get(k)) for k in ["box_center_x", "box_center_y", "box_center_z"]]
        ok = (
            bool(row.get("pocket_family_id"))
            and bool(row.get("state_id"))
            and bool(row.get("protomer_id"))
            and bool(row.get("box_id"))
            and receptor is not None
            and row.get("receptor_pdbqt_file") == receptor.get("prepared_receptor_pdbqt_file")
            and all(value is not None for value in centers)
            and all(value is not None and value > 0 for value in sizes)
            and _boolish(row.get("non_atp_pass"))
            and _boolish(row.get("lower_lateral_pass"))
            and _boolish(row.get("dimer_accessibility_pass"))
            and row.get("box_qc_status") == "PASS"
            and row.get("traceability_status") == "PASS"
            and _boolish(row.get("allowed_for_vina_smoke"))
        )
        if ok:
            enriched = dict(row)
            enriched["_receptor"] = receptor
            eligible.append(enriched)
    eligible.sort(key=lambda row: (row.get("state_id", ""), row.get("pocket_family_id", ""), row.get("protomer_id", ""), row.get("box_id", "")))
    return eligible, []


def _select_smoke(
    ctx: RunContext,
    ligand_manifest: Path | None,
    receptor_manifest: Path | None,
    box_manifest: Path | None,
    requested_compound: str | None,
    requested_state: str | None,
    requested_family: str | None,
    requested_box: str | None,
    requested_protomer: str | None,
    allow_reference: bool,
) -> tuple[SmokeSelection, dict[str, int], list[str]]:
    ligand_rows, ligand_notes = _eligible_ligands(ctx, ligand_manifest, requested_compound)
    receptor_rows, receptor_notes = _eligible_receptors(ctx, receptor_manifest, requested_state, allow_reference)
    box_rows, box_notes = _eligible_boxes(ctx, box_manifest, receptor_rows, requested_state, requested_family, requested_box, requested_protomer)
    box = box_rows[0] if box_rows else None
    receptor = box.get("_receptor") if box else (receptor_rows[0] if receptor_rows else None)
    ligand = ligand_rows[0] if ligand_rows else None
    if box and receptor and receptor.get("state_id") != box.get("state_id"):
        receptor = None
    parts = [
        ligand.get("compound_public_id") if ligand else "NO_LIGAND",
        receptor.get("state_id") if receptor else "NO_STATE",
        box.get("pocket_family_id") if box else "NO_FAMILY",
        box.get("box_id") if box else "NO_BOX",
    ]
    smoke_id = "smoke_" + "_".join(re.sub(r"[^A-Za-z0-9_.-]+", "_", part) for part in parts)
    counts = {"eligible_ligands": len(ligand_rows), "eligible_receptor_states": len(receptor_rows), "eligible_boxes": len(box_rows)}
    return SmokeSelection(ligand=ligand, receptor=receptor, box=box, smoke_id=smoke_id), counts, ligand_notes + receptor_notes + box_notes


def _pdbqt_has_pose(path: Path | None) -> str:
    if not path or not path.is_file():
        return "false"
    text = path.read_text(encoding="utf-8", errors="ignore")
    return "true" if any(line.startswith(("ATOM", "HETATM", "MODEL")) for line in text.splitlines()) else "false"


def _affinity_detected(vina_log: Path | None, stdout_file: Path | None) -> str:
    texts = []
    for path in [vina_log, stdout_file]:
        if path and path.is_file():
            texts.append(path.read_text(encoding="utf-8", errors="ignore"))
    text = "\n".join(texts).lower()
    return "true" if "affinity" in text or "mode" in text and "kcal" in text else "false"


def _forbidden_outputs(ctx: RunContext) -> list[str]:
    hits: list[str] = []
    for rel in FORBIDDEN_OUTPUTS:
        path = ctx.run_dir / rel
        if path.exists():
            hits.append(ctx.relative_to_repo(path))
    return hits


def _scan_hygiene(ctx: RunContext, private_entries: list[Any], pose_file: Path | None) -> tuple[int, list[str], bool, bool, bool, bool, list[str]]:
    leaks = scan_internal_id_leaks(ctx, private_entries)
    coord_hits: list[str] = []
    smiles_logged = False
    for path in public_output_scan_paths(ctx):
        if pose_file and path.resolve() == pose_file.resolve():
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        upper = text.upper()
        if ("SMILES=" in upper or "SMILES:" in upper) and re.search(r"(^|\s)(C|N|O|S|P|Cl|Br|F|I)[A-Za-z0-9@+\-\[\]\(\)=#$\\/]+(\s|$)", text):
            smiles_logged = True
        for line in text.splitlines():
            if line.startswith(("ATOM  ", "HETATM", "MODEL ")):
                coord_hits.append(ctx.relative_to_repo(path))
                break
    coord_leak = bool(coord_hits)
    scanned = [ctx.relative_to_repo(path) for path in public_output_scan_paths(ctx)]
    return len(leaks), coord_hits, smiles_logged, coord_leak, coord_leak, coord_leak, scanned


def _write_report(path: Path, ctx: RunContext, summary: dict[str, Any]) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# M3-T4 Vina Focused Docking Smoke Test",
        "",
        f"- run_id: {summary['run_id']}",
        f"- m2_run_id: {summary['m2_run_id']}",
        f"- mode: {summary['mode']}",
        f"- overall_status: {summary['overall_status']}",
        f"- m3_t5_allowed: {str(summary['m3_t5_allowed']).lower()}",
        "- technical_smoke_only: true",
        "- no_candidate_claims: true",
        "",
        "## Selection",
    ]
    for key, value in summary["selection"].items():
        lines.append(f"- {key}: {value}")
    lines.extend(["", "## Vina", ""])
    for key, value in summary["vina"].items():
        lines.append(f"- {key}: {value}")
    lines.extend(["", "## Blockers", ""])
    lines.extend([f"- {item}" for item in summary["blockers"]] or ["- none"])
    lines.extend(["", "## Warnings", ""])
    lines.extend([f"- {item}" for item in summary["warnings"]] or ["- none"])
    lines.extend(["", "## Next Task", "", "M3-T5 - PBS generator and compound docking job manifest"])
    safe.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_m3_vina_smoke(
    ctx: RunContext,
    m2_run_id: str,
    mode: str = "run",
    force: bool = False,
    compound_public_id: str | None = None,
    state_id: str | None = None,
    pocket_family_id: str | None = None,
    box_id: str | None = None,
    protomer_id: str | None = None,
    vina_executable: str = "vina",
    exhaustiveness: int | None = None,
    num_modes: int | None = None,
    cpu: int = 1,
    seed: int | None = None,
    timeout_seconds: int | None = None,
    require_m3_t2_pass: bool = True,
    require_m3_t3_pass: bool = True,
    allow_independent: bool = False,
    allow_reference_smoke: bool = False,
    update_ignore: bool = True,
    private_map: Path | None = None,
) -> M3VinaSmokeResult:
    if mode not in {"dry-run", "run"}:
        raise ValueError("mode must be dry-run or run")
    if compound_public_id and compound_public_id not in PUBLIC_COMPOUND_IDS:
        raise ValueError("compound_public_id must be one of Cpd-A, Cpd-B, Cpd-C")
    ctx.create_directories()
    _ensure_logs(ctx)
    phase3 = ctx.require_within_run_dir(ctx.run_dir / "phase3_compounds")
    manifest_dir = ctx.require_within_run_dir(phase3 / "manifests")
    qc_dir = ctx.require_within_run_dir(phase3 / "qc")
    report_dir = ctx.require_within_run_dir(phase3 / "reports")
    for directory in [manifest_dir, qc_dir, report_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    manifest_csv = manifest_dir / "vina_smoke_manifest.csv"
    qc_json = qc_dir / "vina_smoke_qc.json"
    qc_csv = qc_dir / "vina_smoke_qc.csv"
    report_md = report_dir / "m3_task4_vina_smoke.md"
    phase3_log = ctx.logs_dir / "phase3_compounds.log"
    ligand_manifest = manifest_dir / "ligand_preparation_manifest.csv"
    receptor_manifest = manifest_dir / "receptor_preparation_manifest.csv"
    box_manifest = manifest_dir / "docking_box_manifest.csv"
    t2_summary_path = qc_dir / "m3_ligand_prep_summary.json"
    t3_summary_path = qc_dir / "m3_receptor_box_summary.json"
    exhaustiveness = exhaustiveness if exhaustiveness is not None else _default_int(ctx, ["exhaustiveness", "vina_exhaustiveness"], 1)
    num_modes = num_modes if num_modes is not None else _default_int(ctx, ["num_modes", "vina_num_modes"], 1)
    seed = seed if seed is not None else _default_int(ctx, ["seed", "vina_seed"], 20260427)
    timeout_seconds = timeout_seconds if timeout_seconds is not None else _default_int(ctx, ["timeout_seconds", "vina_timeout_seconds"], 300)

    blockers: list[str] = []
    warnings: list[str] = []
    qc_rows: list[dict[str, str]] = []
    if update_ignore:
        gitignore_updated, patterns_verified = update_gitignore(ctx.repo_root)
    else:
        gitignore_updated, patterns_verified = False, list(GITIGNORE_PATTERNS)
    private_entries, private_warnings = load_private_map(_safe_private_map_path(ctx, private_map))
    warnings.extend(private_warnings)

    t2 = _load_json(t2_summary_path) if t2_summary_path.is_file() else None
    t3 = _load_json(t3_summary_path) if t3_summary_path.is_file() else None
    smoke_id = "smoke_unselected"
    _append_qc(qc_rows, smoke_id, "m3_t2_summary_present", "precondition", "PASS" if t2 else ("FAIL" if mode == "run" else "WARN"), "M3-T2 summary checked", "Run M3-T2 ligand preparation.", mode == "run" and not t2)
    if not t2 and mode == "run":
        blockers.append("run mode with missing M3-T2 ligand prep summary")
    elif not t2:
        warnings.append("M3-T2 ligand prep summary absent in dry-run")
    t2_passed = bool(t2 and t2.get("overall_status") == "PASS")
    t2_allows_t3 = bool(t2 and t2.get("m3_t3_allowed") is True)
    t2_ok = t2_passed and t2_allows_t3
    if require_m3_t2_pass and not t2_passed and not allow_independent:
        if mode == "run":
            blockers.append("M3-T2 did not PASS")
        else:
            warnings.append("M3-T2 did not PASS")
    if require_m3_t2_pass and not t2_allows_t3 and not allow_independent:
        if mode == "run":
            blockers.append("m3_ligand_prep_summary.json has m3_t3_allowed != true")
        else:
            warnings.append("M3-T2 m3_t3_allowed is not true")
    _append_qc(qc_rows, smoke_id, "m3_t2_passed_or_bypassed", "precondition", "PASS" if t2_ok else ("WARN" if allow_independent or mode == "dry-run" else "FAIL"), "M3-T2 pass and m3_t3_allowed handoff checked", "Run or fix M3-T2.", mode == "run" and require_m3_t2_pass and not t2_ok and not allow_independent)

    _append_qc(qc_rows, smoke_id, "m3_t3_summary_present", "precondition", "PASS" if t3 else ("FAIL" if mode == "run" else "WARN"), "M3-T3 summary checked", "Run M3-T3 receptor/box preparation.", mode == "run" and not t3)
    if not t3 and mode == "run":
        blockers.append("run mode with missing M3-T3 receptor/box summary")
    elif not t3:
        warnings.append("M3-T3 receptor/box summary absent in dry-run")
    t3_passed = bool(t3 and t3.get("overall_status") == "PASS")
    t4_allowed = bool(t3 and t3.get("m3_t4_allowed") is True)
    if require_m3_t3_pass and not t3_passed and not allow_independent:
        if mode == "run":
            blockers.append("M3-T3 did not PASS")
        else:
            warnings.append("M3-T3 did not PASS")
    if not t4_allowed and not allow_independent:
        if mode == "run":
            blockers.append("m3_receptor_box_summary.json has m3_t4_allowed != true")
        else:
            warnings.append("M3-T3 m3_t4_allowed is not true")
    if allow_independent:
        warnings.append("allow_independent diagnostic bypass used; production gates remain closed")
    _append_qc(qc_rows, smoke_id, "m3_t3_passed_or_bypassed", "precondition", "PASS" if t3_passed else ("WARN" if allow_independent or mode == "dry-run" else "FAIL"), "M3-T3 pass/bypass checked", "Run or fix M3-T3.", mode == "run" and require_m3_t3_pass and not t3_passed and not allow_independent)
    _append_qc(qc_rows, smoke_id, "m3_t4_allowed_by_t3", "precondition", "PASS" if t4_allowed else ("WARN" if allow_independent or mode == "dry-run" else "FAIL"), "M3-T3 m3_t4_allowed checked", "Fix M3-T3 gates.", mode == "run" and not t4_allowed and not allow_independent)

    _append_qc(qc_rows, smoke_id, "ligand_manifest_present", "input", "PASS" if ligand_manifest.is_file() else ("FAIL" if mode == "run" else "WARN"), "ligand preparation manifest checked", "Run M3-T2.", mode == "run" and not ligand_manifest.is_file())
    _append_qc(qc_rows, smoke_id, "receptor_manifest_present", "input", "PASS" if receptor_manifest.is_file() else ("FAIL" if mode == "run" else "WARN"), "receptor preparation manifest checked", "Run M3-T3.", mode == "run" and not receptor_manifest.is_file())
    _append_qc(qc_rows, smoke_id, "docking_box_manifest_present", "input", "PASS" if box_manifest.is_file() else ("FAIL" if mode == "run" else "WARN"), "docking box manifest checked", "Run M3-T3.", mode == "run" and not box_manifest.is_file())
    for path, label in [(ligand_manifest, "missing ligand_preparation_manifest.csv"), (receptor_manifest, "missing receptor_preparation_manifest.csv"), (box_manifest, "missing docking_box_manifest.csv")]:
        if not path.is_file():
            if mode == "run":
                blockers.append(label)
            else:
                warnings.append(label)

    selection, selection_counts, selection_notes = _select_smoke(
        ctx,
        ligand_manifest if ligand_manifest.is_file() else None,
        receptor_manifest if receptor_manifest.is_file() else None,
        box_manifest if box_manifest.is_file() else None,
        compound_public_id,
        state_id,
        pocket_family_id,
        box_id,
        protomer_id,
        allow_reference_smoke,
    )
    smoke_id = selection.smoke_id
    warnings.extend(selection_notes if mode == "dry-run" else [])
    for label, obj in [("eligible_ligand_selected", selection.ligand), ("eligible_receptor_selected", selection.receptor), ("eligible_box_selected", selection.box)]:
        ok = obj is not None
        _append_qc(qc_rows, smoke_id, label, "selection", "PASS" if ok else ("FAIL" if mode == "run" else "WARN"), label.replace("_", " "), "Fix input manifests or selection filters.", mode == "run" and not ok)
        if not ok and mode == "run":
            blockers.append("no " + label.replace("eligible_", "").replace("_", " "))

    ligand_path = selection.ligand.get("_path") if selection.ligand else None
    receptor_path = selection.receptor.get("_path") if selection.receptor else None
    box = selection.box or {}
    state_role = (selection.receptor or {}).get("state_role") or box.get("state_role") or None
    if state_role == "reference" or (selection.receptor and selection.receptor.get("state_id") == "3GT8_raw"):
        if not allow_reference_smoke and mode == "run":
            blockers.append("selected state is 3GT8_raw/reference without --allow-reference-smoke")
        elif allow_reference_smoke:
            warnings.append("allow_reference_smoke used; m3_t5_allowed remains false")
    for path, prefix in [(ligand_path, "ligand"), (receptor_path, "receptor")]:
        _append_qc(qc_rows, smoke_id, f"{prefix}_pdbqt_exists", "input", "PASS" if path and path.is_file() else ("FAIL" if mode == "run" else "WARN"), f"{prefix} PDBQT exists checked", f"Prepare valid {prefix} PDBQT.", mode == "run" and not (path and path.is_file()))
        _append_qc(qc_rows, smoke_id, f"{prefix}_pdbqt_nonempty", "input", "PASS" if path and path.is_file() and path.stat().st_size > 0 else ("FAIL" if mode == "run" else "WARN"), f"{prefix} PDBQT non-empty checked", f"Prepare non-empty {prefix} PDBQT.", mode == "run" and not (path and path.is_file() and path.stat().st_size > 0))
    _append_qc(qc_rows, smoke_id, "ligand_path_inside_run_dir", "path_safety", "PASS" if _inside_phase3(ctx, ligand_path) else ("FAIL" if mode == "run" else "WARN"), "ligand path inside phase3 checked", "Use run-local prepared ligand PDBQT.", mode == "run" and not _inside_phase3(ctx, ligand_path))
    _append_qc(qc_rows, smoke_id, "receptor_path_inside_run_dir", "path_safety", "PASS" if _inside_receptor_pdbqt(ctx, receptor_path) else ("FAIL" if mode == "run" else "WARN"), "receptor path inside receptor_pdbqt checked", "Use run-local receptor PDBQT.", mode == "run" and not _inside_receptor_pdbqt(ctx, receptor_path))

    centers = tuple(str(box.get(k, "")) for k in ["box_center_x", "box_center_y", "box_center_z"])
    sizes = tuple(str(box.get(k, "")) for k in ["box_size_x", "box_size_y", "box_size_z"])
    center_values = [_finite_float(v) for v in centers]
    size_values = [_finite_float(v) for v in sizes]
    box_numeric = all(v is not None for v in center_values + size_values)
    box_positive = all(v is not None and v > 0 for v in size_values)
    for check, ok, detail in [
        ("box_traceable_to_m2", box.get("traceability_status") == "PASS", "box traceability checked"),
        ("box_non_atp_pass", _boolish(box.get("non_atp_pass")), "non-ATP checked"),
        ("box_lower_lateral_pass", _boolish(box.get("lower_lateral_pass")), "lower/lateral checked"),
        ("box_dimer_accessibility_pass", _boolish(box.get("dimer_accessibility_pass")), "dimer accessibility checked"),
        ("box_numeric_dimensions", box_numeric, "numeric dimensions checked"),
        ("box_positive_dimensions", box_positive, "positive dimensions checked"),
    ]:
        _append_qc(qc_rows, smoke_id, check, "box", "PASS" if ok else ("FAIL" if mode == "run" else "WARN"), detail, "Use an eligible M3-T3 PASS smoke box.", mode == "run" and not ok)
        if not ok and mode == "run":
            blockers.append("selected " + check.replace("_", " "))

    vina = detect_vina(vina_executable)
    if mode == "run" and not vina.available:
        blockers.append("Vina executable unavailable in run mode")
    elif mode == "dry-run" and not vina.available:
        warnings.append("Vina executable unavailable in dry-run")
    _append_qc(qc_rows, smoke_id, "vina_available", "tooling", "PASS" if vina.available else ("FAIL" if mode == "run" else "WARN"), "Vina availability checked", "Install Vina or pass --vina-executable.", mode == "run" and not vina.available)
    _append_qc(qc_rows, smoke_id, "vina_version_recorded", "tooling", "PASS" if vina.version else ("WARN" if vina.available else "NOT_APPLICABLE"), "Vina version probe recorded", "Check Vina --version.", False)
    if vina.available and not vina.version:
        warnings.append("Vina available but version string could not be parsed")

    smoke_root = ctx.require_within_run_dir(phase3 / "docking_outputs" / "focused_pocket_first" / "smoke" / smoke_id)
    output_pdbqt = ctx.require_within_run_dir(smoke_root / "smoke_pose.pdbqt")
    vina_log = ctx.require_within_run_dir(smoke_root / "vina.log")
    stdout_file = ctx.require_within_run_dir(ctx.jobs_log_dir / f"{smoke_id}.vina.stdout")
    stderr_file = ctx.require_within_run_dir(ctx.jobs_log_dir / f"{smoke_id}.vina.stderr")
    output_inside = True
    try:
        ctx.require_within_run_dir(output_pdbqt)
    except Exception:
        output_inside = False
    _append_qc(qc_rows, smoke_id, "output_path_inside_run_dir", "path_safety", "PASS" if output_inside else "FAIL", "Vina output path checked", "Write smoke output under run_dir.", not output_inside)
    if not output_inside:
        blockers.append("output PDBQT path outside run_dir")
    if output_pdbqt.exists() and not force and mode == "run":
        blockers.append("smoke output already exists and force=false")

    planned = bool(selection.ligand and selection.receptor and selection.box)
    argv: list[str] = []
    if planned and ligand_path and receptor_path:
        argv = build_vina_argv(
            vina.executable,
            receptor_path,
            ligand_path,
            centers,
            sizes,
            output_pdbqt,
            vina_log,
            exhaustiveness,
            num_modes,
            cpu,
            seed,
        )
    _append_qc(qc_rows, smoke_id, "vina_argv_safe", "command", "PASS" if argv and all(isinstance(x, str) and x for x in argv) else ("WARN" if mode == "dry-run" else "FAIL"), "Vina argv list constructed without shell", "Build a complete argv list.", mode == "run" and not argv)
    _append_qc(qc_rows, smoke_id, "exactly_one_vina_command_planned", "command", "PASS" if len([argv] if argv else []) == 1 else ("FAIL" if mode == "run" else "WARN"), "exactly one Vina command planned", "Plan exactly one smoke command.", mode == "run" and not argv)

    manifest_written_at = now_iso()
    row = {
        "smoke_id": smoke_id,
        "run_id": ctx.run_id,
        "m2_run_id": m2_run_id,
        "mode": mode,
        "command_manifest_written_before_execution": "true",
        "command_manifest_written_at": manifest_written_at,
        "compound_public_id": (selection.ligand or {}).get("compound_public_id", ""),
        "ligand_pdbqt_file": ctx.relative_to_repo(ligand_path) if ligand_path else "",
        "ligand_pdbqt_sha256": (selection.ligand or {}).get("_sha256", ""),
        "ligand_preparation_status": (selection.ligand or {}).get("preparation_status", ""),
        "state_id": (selection.receptor or {}).get("state_id", ""),
        "state_role": state_role or "",
        "receptor_pdbqt_file": ctx.relative_to_repo(receptor_path) if receptor_path else "",
        "receptor_pdbqt_sha256": (selection.receptor or {}).get("_sha256", ""),
        "receptor_preparation_status": (selection.receptor or {}).get("preparation_status", ""),
        "pocket_family_id": box.get("pocket_family_id", ""),
        "protomer_id": box.get("protomer_id", ""),
        "box_id": box.get("box_id", ""),
        "box_center_x": centers[0],
        "box_center_y": centers[1],
        "box_center_z": centers[2],
        "box_size_x": sizes[0],
        "box_size_y": sizes[1],
        "box_size_z": sizes[2],
        "source_docking_box_manifest": ctx.relative_to_repo(box_manifest),
        "source_receptor_manifest": ctx.relative_to_repo(receptor_manifest),
        "source_ligand_manifest": ctx.relative_to_repo(ligand_manifest),
        "vina_executable": vina.executable,
        "vina_version": vina.version or "",
        "vina_argv_json": json.dumps(argv),
        "exhaustiveness": str(exhaustiveness),
        "num_modes": str(num_modes),
        "cpu": str(cpu),
        "seed": str(seed),
        "timeout_seconds": str(timeout_seconds),
        "stdout_file": ctx.relative_to_repo(stdout_file),
        "stderr_file": ctx.relative_to_repo(stderr_file),
        "vina_log_file": ctx.relative_to_repo(vina_log),
        "output_pdbqt_file": ctx.relative_to_repo(output_pdbqt),
        "started_at": "",
        "finished_at": "",
        "duration_seconds": "",
        "vina_invoked": "false",
        "vina_invocation_count": "0",
        "exit_code": "",
        "timeout": "false",
        "output_pdbqt_exists": "false",
        "output_pdbqt_size_bytes": "",
        "output_pdbqt_sha256": "",
        "vina_log_exists": "false",
        "vina_log_size_bytes": "",
        "minimal_pose_detected": "unknown",
        "affinity_table_detected": "unknown",
        "allowed_for_m3_t5": "false",
        "smoke_status": "NOT_APPLICABLE" if mode == "dry-run" else "FAIL",
        "smoke_notes": "manifest written before execution",
    }
    _write_csv(manifest_csv, MANIFEST_FIELDS, [row], ctx)
    _append_qc(qc_rows, smoke_id, "vina_command_manifest_written_before_execution", "command", "PASS", "command manifest written before execution")

    invoked = False
    invocation_count = 0
    exit_code = None
    timeout = False
    duration = None
    started_at = None
    finished_at = None
    if mode == "run" and not blockers and argv:
        smoke_root.mkdir(parents=True, exist_ok=True)
        run_result = run_vina_once(argv, stdout_file, stderr_file, timeout_seconds, now_iso)
        invoked = run_result.invoked
        invocation_count = 1 if invoked else 0
        exit_code = run_result.exit_code
        timeout = run_result.timeout
        duration = run_result.duration_seconds
        started_at = run_result.started_at
        finished_at = run_result.finished_at
        if timeout:
            blockers.append("Vina timed out in run mode")
        elif exit_code != 0:
            blockers.append("Vina exited non-zero in run mode")
    _append_qc(qc_rows, smoke_id, "exactly_one_vina_command_invoked", "command", "PASS" if (mode == "run" and invocation_count == 1) or mode == "dry-run" else "FAIL", "Vina invocation count checked", "Invoke exactly one command in run mode.", mode == "run" and invocation_count != 1)
    _append_qc(qc_rows, smoke_id, "stdout_captured", "logging", "PASS" if stdout_file.is_file() or mode == "dry-run" else "FAIL", "stdout captured under logs/jobs", "Capture stdout.", mode == "run" and not stdout_file.is_file())
    _append_qc(qc_rows, smoke_id, "stderr_captured", "logging", "PASS" if stderr_file.is_file() or mode == "dry-run" else "FAIL", "stderr captured under logs/jobs", "Capture stderr.", mode == "run" and not stderr_file.is_file())

    output_exists = output_pdbqt.is_file()
    output_nonempty = output_exists and output_pdbqt.stat().st_size > 0
    output_sha = _sha256(output_pdbqt) if output_exists else ""
    log_exists = vina_log.is_file()
    pose_detected = _pdbqt_has_pose(output_pdbqt)
    affinity_detected = _affinity_detected(vina_log, stdout_file)
    if mode == "run" and not output_nonempty:
        blockers.append("output PDBQT missing or empty in run mode")
    if mode == "run" and pose_detected != "true":
        blockers.append("minimal pose records missing from Vina output PDBQT")
    if mode == "run" and not log_exists:
        blockers.append("Vina log missing in run mode")
    if mode == "run" and not stdout_file.is_file():
        blockers.append("Vina stdout capture missing in run mode")
    if mode == "run" and not stderr_file.is_file():
        blockers.append("Vina stderr capture missing in run mode")
    _append_qc(qc_rows, smoke_id, "vina_output_pdbqt_exists", "output", "PASS" if output_exists else ("NOT_APPLICABLE" if mode == "dry-run" else "FAIL"), "output PDBQT exists checked", "Inspect Vina stderr/log.", mode == "run" and not output_exists)
    _append_qc(qc_rows, smoke_id, "vina_output_pdbqt_nonempty", "output", "PASS" if output_nonempty else ("NOT_APPLICABLE" if mode == "dry-run" else "FAIL"), "output PDBQT non-empty checked", "Inspect Vina stderr/log.", mode == "run" and not output_nonempty)
    _append_qc(qc_rows, smoke_id, "vina_log_exists", "output", "PASS" if log_exists else ("NOT_APPLICABLE" if mode == "dry-run" else "FAIL"), "Vina log existence checked", "Vina log is required for a successful M3-T4 run.", mode == "run" and not log_exists)
    _append_qc(qc_rows, smoke_id, "minimal_pose_detected", "output", "PASS" if pose_detected == "true" else ("NOT_APPLICABLE" if mode == "dry-run" else "FAIL"), "minimal pose marker checked", "Ensure Vina output PDBQT contains pose records.", mode == "run" and pose_detected != "true")

    forbidden = _forbidden_outputs(ctx)
    qsub_count = 0
    for check, ok, detail in [
        ("no_production_docking_launched", not forbidden, "no production docking outputs created"),
        ("no_qsub_launched", qsub_count == 0, "qsub/PBS submission not launched"),
        ("no_broad_docking_launched", not any("broad_anchor" in item for item in forbidden), "broad docking not launched"),
        ("no_pose_attribution_attempted", not any("pose_attribution" in item for item in forbidden), "pose attribution not attempted"),
        ("no_candidate_ranking_attempted", not any("candidate" in item for item in forbidden), "candidate ranking not attempted"),
    ]:
        _append_qc(qc_rows, smoke_id, check, "non_goal", "PASS" if ok else "FAIL", detail, "Remove forbidden outputs and rerun.", not ok)
    if forbidden:
        blockers.append("production/broad/pose/candidate outputs created")

    tracked = []
    for generated in [output_pdbqt, vina_log, stdout_file, stderr_file]:
        if generated.exists() and git_ls_files(ctx.repo_root, generated):
            tracked.append(ctx.relative_to_repo(generated))
    output_ignored = (not output_pdbqt.exists()) or git_check_ignored(ctx.repo_root, output_pdbqt)
    generated_ignored = output_ignored and ((not stdout_file.exists()) or git_check_ignored(ctx.repo_root, stdout_file) or stdout_file.is_relative_to(ctx.logs_dir))
    _append_qc(qc_rows, smoke_id, "generated_outputs_gitignored", "gitignore", "PASS" if generated_ignored else "WARN", "generated output ignore protection checked", "Add smoke output ignore patterns.")
    _append_qc(qc_rows, smoke_id, "generated_outputs_not_tracked", "gitignore", "PASS" if not tracked else "FAIL", "generated output git tracking checked", "Remove generated outputs from git index without deleting local files.", bool(tracked))
    if tracked:
        blockers.append("generated docking output PDBQT/log files tracked by git")

    leak_count, coord_hits, smiles_logged, ligand_coords, receptor_coords, pose_coords, scanned_paths = _scan_hygiene(ctx, list(private_entries.values()), output_pdbqt)
    for check, ok, detail in [
        ("no_internal_id_leak", leak_count == 0, "internal ID leakage scan complete"),
        ("no_smiles_logged", not smiles_logged, "SMILES logging scan complete"),
        ("no_ligand_coordinates_logged", not ligand_coords, "ligand coordinate leakage scan complete"),
        ("no_receptor_coordinates_logged", not receptor_coords, "receptor coordinate leakage scan complete"),
    ]:
        _append_qc(qc_rows, smoke_id, check, "confidentiality", "PASS" if ok else "FAIL", detail, "Remove sensitive tokens/coordinates from public outputs.", not ok)
    if leak_count:
        blockers.append("internal compound ID leakage detected")
    if smiles_logged:
        blockers.append("ligand SMILES printed in public outputs/logs")
    if coord_hits:
        blockers.append("ligand/receptor/pose coordinates printed outside intended PDBQT output")
    _append_qc(qc_rows, smoke_id, "non_goals_preserved", "scope", "PASS", "no qsub, PBS generation, production/mini/scaling/broad docking, pose attribution, clustering, scoring, or candidate claims")

    for check_id in REQUIRED_CHECKS:
        if not any(row["check_id"] == check_id for row in qc_rows):
            _append_qc(qc_rows, smoke_id, check_id, "coverage", "NOT_APPLICABLE", "check not reached")

    if mode == "dry-run" and not blockers:
        warnings.append("dry-run mode; Vina was not invoked")
    if allow_reference_smoke:
        warnings.append("reference smoke allowed diagnostically; M3-T5 remains closed")
    status = "FAIL" if blockers else ("WARN" if warnings or mode == "dry-run" else "PASS")
    m3_t5_allowed = status == "PASS" and mode == "run" and not allow_independent and not allow_reference_smoke and t2_ok and t3_passed and t4_allowed and invocation_count == 1
    row.update(
        {
            "started_at": started_at or "",
            "finished_at": finished_at or "",
            "duration_seconds": "" if duration is None else str(duration),
            "vina_invoked": _bool_text(invoked),
            "vina_invocation_count": str(invocation_count),
            "exit_code": "" if exit_code is None else str(exit_code),
            "timeout": _bool_text(timeout),
            "output_pdbqt_exists": _bool_text(output_exists),
            "output_pdbqt_size_bytes": str(output_pdbqt.stat().st_size) if output_exists else "",
            "output_pdbqt_sha256": output_sha,
            "vina_log_exists": _bool_text(log_exists),
            "vina_log_size_bytes": str(vina_log.stat().st_size) if log_exists else "",
            "minimal_pose_detected": pose_detected if output_exists else ("unknown" if mode == "dry-run" else "false"),
            "affinity_table_detected": affinity_detected if (log_exists or stdout_file.exists()) else "unknown",
            "allowed_for_m3_t5": _bool_text(m3_t5_allowed),
            "smoke_status": status if mode == "run" else "NOT_APPLICABLE",
            "smoke_notes": "; ".join(dict.fromkeys(blockers + warnings)) or "technical smoke command completed",
        }
    )
    _write_csv(manifest_csv, MANIFEST_FIELDS, [row], ctx)
    _write_csv(qc_csv, QC_FIELDS, qc_rows, ctx)
    counts = {
        **selection_counts,
        "vina_commands_planned": 1 if argv else 0,
        "vina_commands_invoked": invocation_count,
        "output_pose_files": 1 if output_exists else 0,
        "production_outputs_created": len(forbidden),
        "qsub_commands_invoked": qsub_count,
        "confidentiality_leaks": leak_count,
        "tracked_generated_outputs": len(tracked),
    }
    summary = {
        "schema_version": "m3_vina_smoke_qc_v1",
        "run_id": ctx.run_id,
        "m2_run_id": m2_run_id,
        "reviewed_at": now_iso(),
        "mode": mode,
        "overall_status": status,
        "m3_t5_allowed": m3_t5_allowed,
        "allow_independent": allow_independent,
        "allow_reference_smoke": allow_reference_smoke,
        "require_m3_t2_pass": require_m3_t2_pass,
        "require_m3_t3_pass": require_m3_t3_pass,
        "m3_t2": {"summary_found": t2 is not None, "overall_status": t2.get("overall_status") if t2 else None, "ligand_preparation_passed": t2_ok if t2 else None},
        "m3_t3": {"summary_found": t3 is not None, "overall_status": t3.get("overall_status") if t3 else None, "m3_t4_allowed": t3.get("m3_t4_allowed") if t3 else None},
        "selection": {
            "compound_public_id": row["compound_public_id"] or None,
            "state_id": row["state_id"] or None,
            "state_role": row["state_role"] or None,
            "pocket_family_id": row["pocket_family_id"] or None,
            "protomer_id": row["protomer_id"] or None,
            "box_id": row["box_id"] or None,
            "selection_policy": "deterministic_first_eligible",
        },
        "vina": {
            "available": vina.available,
            "version": vina.version,
            "executable": vina.executable,
            "invoked": invoked,
            "invocation_count": invocation_count,
            "exit_code": exit_code,
            "timeout": timeout,
            "duration_seconds": duration,
        },
        "outputs": {
            "manifest_file": ctx.relative_to_repo(manifest_csv),
            "qc_csv_file": ctx.relative_to_repo(qc_csv),
            "report_file": ctx.relative_to_repo(report_md),
            "stdout_file": ctx.relative_to_repo(stdout_file),
            "stderr_file": ctx.relative_to_repo(stderr_file),
            "vina_log_file": ctx.relative_to_repo(vina_log),
            "output_pdbqt_file": ctx.relative_to_repo(output_pdbqt),
            "output_pdbqt_sha256": output_sha or None,
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
            "pose_coordinates_logged_outside_pose_file": pose_coords,
        },
        "non_goals_preserved": {
            "production_docking": True,
            "pbs_generation": True,
            "qsub_submission": True,
            "mini_docking": True,
            "scaling_docking": True,
            "broad_docking": True,
            "pose_attribution": True,
            "pose_clustering": True,
            "compound_scoring": True,
            "candidate_nomination": True,
        },
    }
    ctx.require_within_run_dir(qc_json).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_report(report_md, ctx, summary)
    with ctx.require_within_run_dir(phase3_log).open("a", encoding="utf-8") as handle:
        handle.write(
            "{0} M3-T4 status={1} run_id={2} m2_run_id={3} mode={4} m3_t5_allowed={5} compound_public_id={6} state_id={7} pocket_family_id={8} box_id={9} vina_invoked={10} vina_invocation_count={11} exit_code={12} no_qsub_submission_attempted=true no_production_docking_attempted=true no_broad_docking_attempted=true no_pose_attribution_attempted=true no_scoring_attempted=true no_candidate_nomination_attempted=true\n".format(
                now_iso(),
                status,
                ctx.run_id,
                m2_run_id,
                mode,
                str(m3_t5_allowed).lower(),
                row["compound_public_id"] or "null",
                row["state_id"] or "null",
                row["pocket_family_id"] or "null",
                row["box_id"] or "null",
                str(invoked).lower(),
                invocation_count,
                "null" if exit_code is None else exit_code,
            )
        )
    log_master(ctx, status, "M3-T4 Vina smoke completed without production docking/qsub/scoring/candidate claims")
    append_phase_status(
        ctx,
        "phase3_compounds",
        status,
        "M3-T4 Vina focused docking smoke completed",
        {
            "task": "M3-T4",
            "m2_run_id": m2_run_id,
            "mode": mode,
            "m3_t5_allowed": m3_t5_allowed,
            "compound_public_id": row["compound_public_id"] or None,
            "state_id": row["state_id"] or None,
            "pocket_family_id": row["pocket_family_id"] or None,
            "box_id": row["box_id"] or None,
            "vina_invoked": invoked,
            "vina_invocation_count": invocation_count,
            "exit_code": exit_code,
            "no_qsub_submission_attempted": True,
            "no_production_docking_attempted": True,
            "no_broad_docking_attempted": True,
            "no_pose_attribution_attempted": True,
            "no_scoring_attempted": True,
            "no_candidate_nomination_attempted": True,
        },
    )
    if invoked:
        append_job_status(
            ctx,
            smoke_id,
            status,
            stdout=ctx.relative_to_repo(stdout_file),
            stderr=ctx.relative_to_repo(stderr_file),
            details={
                "phase": "phase3_compounds",
                "task": "M3-T4",
                "job_type": "vina_smoke",
                "job_id": smoke_id,
                "compound_public_id": row["compound_public_id"],
                "state_id": row["state_id"],
                "pocket_family_id": row["pocket_family_id"],
                "box_id": row["box_id"],
                "output_pdbqt_file": ctx.relative_to_repo(output_pdbqt),
                "exit_code": exit_code,
            },
        )
    return M3VinaSmokeResult(
        status=status,
        m3_t5_allowed=m3_t5_allowed,
        blockers=list(dict.fromkeys(blockers)),
        warnings=list(dict.fromkeys(warnings)),
        manifest_csv=manifest_csv,
        qc_json=qc_json,
        qc_csv=qc_csv,
        report_md=report_md,
        phase3_log=phase3_log,
        counts=counts,
        selection=summary["selection"],
        vina=summary["vina"],
    )
