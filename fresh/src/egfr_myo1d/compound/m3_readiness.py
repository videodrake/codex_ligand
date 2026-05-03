"""Milestone 3 Task 0 readiness audit and directory skeleton.

This module validates the M2 accepted-pocket handoff package and M1 receptor
inputs before M3 begins. It deliberately avoids ligand/receptor preparation,
docking job creation, Vina execution, scoring, and candidate nomination.
"""

from __future__ import annotations

import csv
import glob
import hashlib
import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import (
    append_job_status,
    append_phase_status,
    log_master,
)
from egfr_myo1d.core.run_context import RunContext, RunContextError, ensure_within, validate_run_id


COMPOUND_IDS = ["Cpd-A", "Cpd-B", "Cpd-C"]
PRIMARY_STATES = ["EGFR_160-185", "EGFR_170-200"]
REFERENCE_STATE = "3GT8_raw"
RECEPTOR_STATES = [*PRIMARY_STATES, REFERENCE_STATE]
ALLOWED_LIGAND_EXTENSIONS = {".sdf", ".mol2", ".pdb", ".pdbqt", ".smi", ".smiles"}

READINESS_FIELDS = [
    "check_id",
    "category",
    "requirement",
    "expected_path",
    "observed_path",
    "status",
    "severity",
    "m3_blocking",
    "details",
    "recommended_fix",
]

PATH_RESOLUTION_FIELDS = [
    "input_id",
    "preferred_paths",
    "observed_path",
    "status",
    "details",
]

SKELETON_QC_FIELDS = [
    "path",
    "status",
    "created_or_present",
    "details",
]


@dataclass
class M3ReadinessResult:
    status: str
    m3_docking_allowed: bool
    manual_review_required: bool
    blockers: list[str]
    warnings: list[str]
    report_csv: Path
    summary_json: Path
    report_md: Path
    path_resolution_csv: Path
    skeleton_qc_csv: Path
    phase3_log: Path
    counts: dict[str, int]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _as_repo_path(ctx: RunContext, path: Path | None) -> str | None:
    if path is None:
        return None
    return ctx.relative_to_repo(path)


def _boolish(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "pass"}


def _finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _safe_write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]], ctx: RunContext) -> None:
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    with safe_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _read_csv_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _severity_for(status: str, blocking: bool) -> str:
    if status == "FAIL" and blocking:
        return "BLOCKER"
    if status == "WARN":
        return "MAJOR"
    if status == "FAIL":
        return "MAJOR"
    return "INFO"


def _add_check(
    checks: list[dict[str, Any]],
    check_id: str,
    category: str,
    requirement: str,
    expected_path: str,
    observed_path: str | None,
    status: str,
    blocking: bool,
    details: str,
    recommended_fix: str = "",
) -> None:
    checks.append(
        {
            "check_id": check_id,
            "category": category,
            "requirement": requirement,
            "expected_path": expected_path,
            "observed_path": observed_path or "",
            "status": status,
            "severity": _severity_for(status, blocking),
            "m3_blocking": str(bool(blocking)),
            "details": details,
            "recommended_fix": recommended_fix,
        }
    )


def _safe_m2_dir(ctx: RunContext, m2_run_id: str) -> Path:
    final_m2_run_id = validate_run_id(m2_run_id)
    runs_root = ensure_within(ctx.fresh_root / "runs", ctx.fresh_root)
    return ensure_within(runs_root / final_m2_run_id, runs_root)


def _first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.is_file():
            return path.resolve()
    return None


def _resolve_inside_fresh(ctx: RunContext, base: Path, relatives: list[str]) -> tuple[Path | None, list[str]]:
    preferred = []
    for relative in relatives:
        candidate = ensure_within(base / relative, ctx.fresh_root)
        preferred.append(ctx.relative_to_repo(candidate))
        if candidate.is_file():
            return candidate, preferred
    return None, preferred


def _phase3_dirs(ctx: RunContext, create_skeleton: bool) -> list[Path]:
    root = ctx.run_dir / "phase3_compounds"
    required_output_dirs = [
        root,
        root / "manifests",
        root / "qc",
        root / "reports",
        ctx.logs_dir,
        ctx.jobs_log_dir,
        ctx.errors_dir,
    ]
    if not create_skeleton:
        return required_output_dirs
    return [
        *required_output_dirs,
        root / "prepared_ligands",
        root / "prepared_ligands" / "Cpd-A",
        root / "prepared_ligands" / "Cpd-B",
        root / "prepared_ligands" / "Cpd-C",
        root / "receptor_pdbqt",
        root / "receptor_pdbqt" / "EGFR_160-185",
        root / "receptor_pdbqt" / "EGFR_170-200",
        root / "receptor_pdbqt" / "3GT8_raw",
        root / "docking_inputs",
        root / "docking_inputs" / "focused_pocket_first",
        root / "docking_inputs" / "broad_anchor_scan_optional",
        root / "docking_outputs",
        root / "docking_outputs" / "focused_pocket_first",
        root / "docking_outputs" / "broad_anchor_scan_optional",
        root / "tables",
    ]


def _create_skeleton(ctx: RunContext, create_skeleton: bool) -> list[dict[str, Any]]:
    rows = []
    for directory in _phase3_dirs(ctx, create_skeleton):
        safe_dir = ctx.require_within_run_dir(directory)
        existed = safe_dir.is_dir()
        safe_dir.mkdir(parents=True, exist_ok=True)
        rows.append(
            {
                "path": ctx.relative_to_repo(safe_dir),
                "status": "PASS",
                "created_or_present": "present" if existed else "created",
                "details": "M3-T0 directory skeleton path",
            }
        )
        gitkeep = safe_dir / ".gitkeep"
        ctx.require_within_run_dir(gitkeep)
        if not gitkeep.exists():
            gitkeep.write_text("", encoding="utf-8")
    return rows


def _touch_logs(ctx: RunContext) -> None:
    for path in [
        ctx.logs_dir / "master.log",
        ctx.logs_dir / "phase_status.jsonl",
        ctx.logs_dir / "job_status.jsonl",
        ctx.errors_dir / "error_summary.txt",
        ctx.errors_dir / "failed_jobs.csv",
        ctx.logs_dir / "phase3_compounds.log",
    ]:
        safe_path = ctx.require_within_run_dir(path)
        safe_path.parent.mkdir(parents=True, exist_ok=True)
        if not safe_path.exists():
            safe_path.touch()
    failed_jobs = ctx.errors_dir / "failed_jobs.csv"
    if failed_jobs.stat().st_size == 0:
        failed_jobs.write_text("timestamp,job_name,status,message\n", encoding="utf-8")


def _find_ligand(ctx: RunContext, compound_id: str) -> Path | None:
    patterns = [
        ctx.fresh_root / "data" / "raw" / "ligands" / f"{compound_id}.*",
        ctx.run_dir / "phase3_compounds" / "ligands" / f"{compound_id}.*",
    ]
    for pattern in patterns:
        for match in sorted(glob.glob(str(pattern))):
            path = Path(match).resolve()
            if path.suffix.lower() not in ALLOWED_LIGAND_EXTENSIONS:
                continue
            if not path.is_file():
                continue
            ensure_within(path, ctx.fresh_root)
            if path.stat().st_size > 0:
                return path
    return None


def _check_m2_review(
    ctx: RunContext,
    m2_dir: Path,
    mode: str,
    checks: list[dict[str, Any]],
    warnings: list[str],
    blockers: list[str],
) -> None:
    review_path = m2_dir / "reviews" / "m2_review_status.json"
    if not review_path.is_file():
        final_status_path = _first_existing(
            [
                m2_dir / "phase2_pockets" / "final" / "m2_final_status.json",
                m2_dir / "phase2_pockets" / "export_for_m3" / "m2_final_status.json",
                m2_dir / "manifest" / "m2_final_status.json",
            ]
        )
        if final_status_path is not None:
            try:
                payload = _read_json(final_status_path)
            except (OSError, json.JSONDecodeError) as exc:
                blocker = f"M2 final status is unreadable: {exc}"
                blockers.append(blocker)
                _add_check(
                    checks,
                    "M2_REVIEW_STATUS",
                    "m2_review",
                    "M2.8 final status must be parseable when dedicated review status is absent",
                    ctx.relative_to_repo(final_status_path),
                    ctx.relative_to_repo(final_status_path),
                    "FAIL",
                    True,
                    blocker,
                    "Regenerate M2.8 export/final status JSON.",
                )
                return
            review_blocks = []
            if payload.get("m3_docking_allowed") is not True:
                review_blocks.append("m3_docking_allowed_not_true")
            if payload.get("blockers"):
                review_blocks.append("blockers_present")
            if payload.get("synthetic_fixture") is True:
                review_blocks.append("synthetic_fixture=true")
            if review_blocks:
                blocker = "M2 final status blocks M3 readiness: " + "; ".join(review_blocks)
                blockers.append(blocker)
                _add_check(
                    checks,
                    "M2_REVIEW_STATUS",
                    "m2_review",
                    "M2.8 final status must allow M3 before readiness can pass",
                    ctx.relative_to_repo(final_status_path),
                    ctx.relative_to_repo(final_status_path),
                    "FAIL",
                    True,
                    blocker,
                    "Fix M2 final-status blockers or obtain an explicit documented override.",
                )
            else:
                _add_check(
                    checks,
                    "M2_REVIEW_STATUS",
                    "m2_review",
                    "M2.8 final status permits M3 start when dedicated review status is absent",
                    ctx.relative_to_repo(final_status_path),
                    ctx.relative_to_repo(final_status_path),
                    "PASS",
                    False,
                    "M2.8 final status has m3_docking_allowed=true and no blockers.",
                )
            return
        message = "M2 review status file not found"
        status = "FAIL" if mode == "docking" else "WARN"
        blocking = mode == "docking"
        if blocking:
            blockers.append(message)
        else:
            warnings.append(message)
        _add_check(
            checks,
            "M2_REVIEW_STATUS",
            "m2_review",
            "M2 review GO/NO-GO must be available before docking mode",
            ctx.relative_to_repo(review_path),
            None,
            status,
            blocking,
            message,
            "Run the M2 completion audit before production M3 docking.",
        )
        return
    try:
        payload = _read_json(review_path)
    except (OSError, json.JSONDecodeError) as exc:
        blocker = f"M2 review status is unreadable: {exc}"
        blockers.append(blocker)
        _add_check(
            checks,
            "M2_REVIEW_STATUS",
            "m2_review",
            "M2 review GO/NO-GO must be parseable",
            ctx.relative_to_repo(review_path),
            ctx.relative_to_repo(review_path),
            "FAIL",
            True,
            blocker,
            "Regenerate M2 review status JSON.",
        )
        return

    review_blocks = []
    if payload.get("m3_go_no_go") == "NO_GO":
        review_blocks.append("m3_go_no_go=NO_GO")
    if payload.get("m3_docking_allowed") is False:
        review_blocks.append("m3_docking_allowed=false")
    if payload.get("blockers"):
        review_blocks.append("blockers_present")
    if review_blocks:
        blocker = "M2 review blocks M3 readiness: " + "; ".join(review_blocks)
        blockers.append(blocker)
        _add_check(
            checks,
            "M2_REVIEW_STATUS",
            "m2_review",
            "M2 review must allow M3 before readiness can pass",
            ctx.relative_to_repo(review_path),
            ctx.relative_to_repo(review_path),
            "FAIL",
            True,
            blocker,
            "Fix M2 review blockers or obtain an explicit documented override.",
        )
    else:
        _add_check(
            checks,
            "M2_REVIEW_STATUS",
            "m2_review",
            "M2 review permits M3 start",
            ctx.relative_to_repo(review_path),
            ctx.relative_to_repo(review_path),
            "PASS",
            False,
            "M2 review does not block M3.",
        )


def _path_inputs(ctx: RunContext, m2_dir: Path) -> tuple[dict[str, Path | None], list[dict[str, Any]]]:
    specs = {
        "accepted_pockets_for_m3": [
            "phase2_pockets/export_for_m3/accepted_pockets_for_m3.csv",
            "phase2_pockets/tables/accepted_pockets_for_m3.csv",
        ],
        "accepted_pocket_boxes": [
            "phase2_pockets/export_for_m3/zone_composite_boxes.csv",
            "phase2_pockets/export_for_m3/accepted_pocket_boxes.csv",
            "phase2_pockets/tables/accepted_pocket_boxes.csv",
        ],
        "accepted_pocket_residue_sets": [
            "phase2_pockets/export_for_m3/accepted_pocket_residue_sets.json",
        ],
        "ppi_consensus_patch": [
            "phase2_pockets/export_for_m3/ppi_consensus_patch.csv",
            "phase1_ppi/tables/ppi_consensus_patch.csv",
            "phase1_ppi/consensus/ppi_consensus_patch_merged.csv",
        ],
        "pocket_gate_qc": [
            "phase2_pockets/export_for_m3/pocket_gate_qc.csv",
            "phase2_pockets/tables/pocket_gate_qc.csv",
            "phase2_pockets/gated/pocket_gate_qc.csv",
        ],
        "m2_final_report": [
            "phase2_pockets/export_for_m3/m2_final_report.md",
            "reports/milestone2_summary.md",
            "reports/m2_final_report.md",
        ],
        "m2_export_manifest": [
            "phase2_pockets/export_for_m3/m2_export_manifest.json",
        ],
        "atp_reference": [
            "phase2_pockets/export_for_m3/atp_site_reference.csv",
            "phase2_pockets/atp_reference/atp_site_reference.csv",
            "phase2_pockets/tables/atp_site_reference.csv",
        ],
        "pocket_candidates_merged": [
            "phase2_pockets/merged/pocket_candidates_merged.csv",
            "phase2_pockets/tables/pocket_candidates_merged.csv",
        ],
    }
    resolved: dict[str, Path | None] = {}
    rows = []
    for input_id, relatives in specs.items():
        path, preferred = _resolve_inside_fresh(ctx, m2_dir, relatives)
        resolved[input_id] = path
        rows.append(
            {
                "input_id": input_id,
                "preferred_paths": ";".join(preferred),
                "observed_path": ctx.relative_to_repo(path) if path else "",
                "status": "PASS" if path else "FAIL",
                "details": "resolved from documented M2 export/fallback paths" if path else "not found",
            }
        )
    return resolved, rows


def _check_parseable_csv(
    ctx: RunContext,
    path: Path | None,
    checks: list[dict[str, Any]],
    blockers: list[str],
    check_id: str,
    category: str,
    requirement: str,
    expected_path: str,
    blocking: bool = True,
) -> tuple[list[dict[str, str]], list[str]]:
    if path is None:
        message = f"missing {requirement}"
        if blocking:
            blockers.append(message)
        _add_check(checks, check_id, category, requirement, expected_path, None, "FAIL", blocking, message)
        return [], []
    try:
        rows, header = _read_csv_rows(path)
    except (OSError, csv.Error, UnicodeDecodeError) as exc:
        message = f"{requirement} is not parseable CSV: {exc}"
        if blocking:
            blockers.append(message)
        _add_check(
            checks,
            check_id,
            category,
            requirement,
            expected_path,
            ctx.relative_to_repo(path),
            "FAIL",
            blocking,
            message,
        )
        return [], []
    _add_check(
        checks,
        check_id,
        category,
        requirement,
        expected_path,
        ctx.relative_to_repo(path),
        "PASS",
        False,
        f"parseable CSV with {len(rows)} rows",
    )
    return rows, header


def _row_value(row: dict[str, str], *names: str) -> str:
    for name in names:
        if name in row and str(row[name]).strip() != "":
            return str(row[name]).strip()
    return ""


def _is_reference_only(row: dict[str, str]) -> bool:
    if _boolish(row.get("reference_only")):
        return True
    if _row_value(row, "export_tier", "acceptance_tier") == "reference_only":
        return True
    state = _row_value(row, "representative_state_id", "state_id")
    return state == REFERENCE_STATE and int(_row_value(row, "primary_state_count") or 0) == 0


def _is_primary_or_secondary(row: dict[str, str]) -> bool:
    gate_class = _row_value(row, "gate_class")
    export_tier = _row_value(row, "export_tier", "acceptance_tier")
    if gate_class in {"accepted_primary_pocket", "accepted_secondary_cryptic"}:
        return True
    return export_tier in {"primary", "secondary_review"}


def _docking_ready(row: dict[str, str]) -> bool:
    if _is_reference_only(row) or not _is_primary_or_secondary(row):
        return False
    gate_class = _row_value(row, "gate_class")
    if gate_class and gate_class not in {"accepted_primary_pocket", "accepted_secondary_cryptic"}:
        return False
    export_tier = _row_value(row, "export_tier", "acceptance_tier")
    if gate_class == "accepted_primary_pocket" or export_tier == "primary":
        if "m3_primary_allowed" in row and not _boolish(row.get("m3_primary_allowed")):
            return False
    if gate_class == "accepted_secondary_cryptic" or export_tier == "secondary_review":
        if "m3_review_allowed" in row and not _boolish(row.get("m3_review_allowed")):
            return False
    hard_gate = row.get("hard_gate_pass")
    if hard_gate is not None and hard_gate != "" and not _boolish(hard_gate):
        return False
    for column in ["non_atp_pass", "lower_lateral_pass", "dimer_accessibility_pass"]:
        if column in row and not _boolish(row[column]):
            return False
    state = _row_value(row, "representative_state_id", "state_id")
    if state == REFERENCE_STATE:
        return False
    return True


def _accepted_counts(rows: list[dict[str, str]]) -> dict[str, int]:
    counts = {
        "accepted_pocket_rows": len(rows),
        "docking_ready_pocket_rows": 0,
        "primary_state_pocket_rows": 0,
        "secondary_pocket_rows": 0,
        "reference_only_rows": 0,
    }
    for row in rows:
        if _docking_ready(row):
            counts["docking_ready_pocket_rows"] += 1
        gate_class = _row_value(row, "gate_class")
        tier = _row_value(row, "export_tier", "acceptance_tier")
        if gate_class == "accepted_primary_pocket" or tier == "primary":
            counts["primary_state_pocket_rows"] += 1
        if gate_class == "accepted_secondary_cryptic" or tier == "secondary_review":
            counts["secondary_pocket_rows"] += 1
        if _is_reference_only(row):
            counts["reference_only_rows"] += 1
    return counts


def _check_accepted_pockets(
    ctx: RunContext,
    path: Path | None,
    checks: list[dict[str, Any]],
    warnings: list[str],
    blockers: list[str],
) -> tuple[list[dict[str, str]], dict[str, int], set[str]]:
    rows, header = _check_parseable_csv(
        ctx,
        path,
        checks,
        blockers,
        "M2_ACCEPTED_POCKETS_PARSE",
        "m2_handoff",
        "accepted_pockets_for_m3.csv",
        "phase2_pockets/export_for_m3/accepted_pockets_for_m3.csv",
    )
    counts = _accepted_counts(rows)
    states = {_row_value(row, "representative_state_id", "state_id") for row in rows if _row_value(row, "representative_state_id", "state_id")}
    if rows and not {"pocket_family_id", "gate_class"}.issubset(set(header)):
        warning = "accepted_pockets_for_m3.csv uses a noncanonical but parseable schema"
        warnings.append(warning)
        _add_check(
            checks,
            "M2_ACCEPTED_POCKETS_SCHEMA",
            "m2_handoff",
            "accepted pocket table should expose pocket family and gate class",
            "pocket_family_id,gate_class or documented equivalents",
            ctx.relative_to_repo(path) if path else None,
            "WARN",
            False,
            warning,
        )
    if path and not rows:
        blocker = "accepted_pockets_for_m3.csv has zero accepted rows"
        blockers.append(blocker)
        _add_check(checks, "M2_ACCEPTED_POCKETS_ROWS", "m2_handoff", "at least one accepted pocket row", ctx.relative_to_repo(path), ctx.relative_to_repo(path), "FAIL", True, blocker)
    elif counts["docking_ready_pocket_rows"] == 0:
        blocker = "no accepted primary or secondary pocket is M3 docking-ready"
        blockers.append(blocker)
        _add_check(checks, "M2_ACCEPTED_POCKETS_READY", "m2_handoff", "at least one primary or secondary M3-ready pocket", ctx.relative_to_repo(path) if path else "", ctx.relative_to_repo(path) if path else None, "FAIL", True, blocker)
    else:
        _add_check(checks, "M2_ACCEPTED_POCKETS_READY", "m2_handoff", "at least one primary or secondary M3-ready pocket", ctx.relative_to_repo(path) if path else "", ctx.relative_to_repo(path) if path else None, "PASS", False, f"{counts['docking_ready_pocket_rows']} docking-ready pocket rows")
    invalid_states = sorted(state for state in states if state not in RECEPTOR_STATES)
    if invalid_states:
        blocker = "accepted pockets contain unsupported state_id values: " + ",".join(invalid_states)
        blockers.append(blocker)
        _add_check(checks, "M2_ACCEPTED_POCKETS_STATES", "m2_handoff", "accepted pocket states must be documented EGFR states", ",".join(RECEPTOR_STATES), ctx.relative_to_repo(path) if path else None, "FAIL", True, blocker)
    return rows, counts, states


def _check_boxes(
    ctx: RunContext,
    path: Path | None,
    accepted_rows: list[dict[str, str]],
    checks: list[dict[str, Any]],
    blockers: list[str],
) -> int:
    rows, header = _check_parseable_csv(
        ctx,
        path,
        checks,
        blockers,
        "M2_BOXES_PARSE",
        "m2_handoff",
        "accepted_pocket_boxes.csv",
        "phase2_pockets/export_for_m3/accepted_pocket_boxes.csv",
    )
    if path is None:
        return 0
    required = {"pocket_family_id", "box_center_x", "box_center_y", "box_center_z", "box_size_x", "box_size_y", "box_size_z"}
    if not required.issubset(set(header)):
        blocker = "accepted_pocket_boxes.csv is missing required box columns"
        blockers.append(blocker)
        _add_check(checks, "M2_BOXES_SCHEMA", "m2_handoff", "box table must include finite center and positive size columns", ",".join(sorted(required)), ctx.relative_to_repo(path), "FAIL", True, blocker)
        return 0
    valid = 0
    bad_messages = []
    for row in rows:
        center_ok = all(_finite_float(row.get(col)) is not None for col in ["box_center_x", "box_center_y", "box_center_z"])
        sizes = [_finite_float(row.get(col)) for col in ["box_size_x", "box_size_y", "box_size_z"]]
        size_ok = all(value is not None and value > 0 for value in sizes)
        state = _row_value(row, "state_id", "representative_state_id")
        protomer = _row_value(row, "protomer_id", "chain_id")
        if center_ok and size_ok and state and protomer:
            valid += 1
        else:
            bad_messages.append(row.get("pocket_family_id", "<unknown>"))
    accepted_ids = {row.get("pocket_family_id", "") for row in accepted_rows if _docking_ready(row)}
    box_ids = {row.get("pocket_family_id", "") for row in rows}
    missing_boxes = sorted(accepted_ids - box_ids)
    if bad_messages or missing_boxes:
        blocker = "malformed or missing accepted pocket boxes"
        if bad_messages:
            blocker += "; invalid=" + ",".join(bad_messages)
        if missing_boxes:
            blocker += "; missing=" + ",".join(missing_boxes)
        blockers.append(blocker)
        _add_check(checks, "M2_BOXES_VALID", "m2_handoff", "all active accepted pockets must have valid boxes", ctx.relative_to_repo(path), ctx.relative_to_repo(path), "FAIL", True, blocker)
    else:
        _add_check(checks, "M2_BOXES_VALID", "m2_handoff", "all active accepted pockets must have valid boxes", ctx.relative_to_repo(path), ctx.relative_to_repo(path), "PASS", False, f"{valid} valid box rows")
    return valid


def _check_reference_csv(
    ctx: RunContext,
    path: Path | None,
    input_name: str,
    check_id: str,
    checks: list[dict[str, Any]],
    blockers: list[str],
    required_columns_any: list[str],
) -> None:
    rows, header = _check_parseable_csv(
        ctx,
        path,
        checks,
        blockers,
        f"{check_id}_PARSE",
        "m2_handoff",
        input_name,
        f"documented M2 path for {input_name}",
    )
    if path is None:
        return
    if not rows:
        blocker = f"{input_name} contains no evidence rows"
        blockers.append(blocker)
        _add_check(checks, f"{check_id}_ROWS", "m2_handoff", f"{input_name} must contain evidence", ctx.relative_to_repo(path), ctx.relative_to_repo(path), "FAIL", True, blocker)
        return
    if not set(required_columns_any).intersection(set(header)):
        blocker = f"{input_name} lacks required/equivalent mapping columns"
        blockers.append(blocker)
        _add_check(checks, f"{check_id}_SCHEMA", "m2_handoff", f"{input_name} must expose residue/protomer or centroid evidence", ",".join(required_columns_any), ctx.relative_to_repo(path), "FAIL", True, blocker)


def _check_receptor_inputs(
    ctx: RunContext,
    m2_dir: Path,
    accepted_states: set[str],
    checks: list[dict[str, Any]],
    blockers: list[str],
    warnings: list[str],
) -> dict[str, dict[str, Any]]:
    receptor_inputs: dict[str, dict[str, Any]] = {}
    for state in RECEPTOR_STATES:
        data_state_dir = ctx.fresh_root / "data" / "normalized" / "receptors" / state
        prepared_receptor_dir = m2_dir / "prepared" / "m2_1_ppi_inputs" / state / "receptor"
        mapping = _first_existing(
            [
                data_state_dir / "receptor_mapping.csv",
                m2_dir / "qc" / f"{state}_receptor_mapping.csv",
                m2_dir / "normalized" / "receptors" / f"{state}_receptor_mapping.csv",
                m2_dir / "normalized" / "receptors" / state / "receptor_mapping.csv",
            ]
        )
        m2_runtime = _first_existing(
            [
                prepared_receptor_dir / f"{state}_runtime_offset_receptor_only.pdb",
                m2_dir / "normalized" / "receptors" / f"{state}_runtime_offset_receptor_only.pdb",
                m2_dir / "normalized" / "receptors" / state / "runtime_offset_receptor_only.pdb",
            ]
        )
        runtime_candidates = [data_state_dir / "runtime_offset_receptor_only.pdb"]
        if m2_runtime is not None:
            runtime_candidates.append(m2_runtime)
        runtime = _first_existing(runtime_candidates)
        if state == REFERENCE_STATE:
            dockable_candidates = [
                data_state_dir / "dockable_reference_explicit_AB.pdb",
                data_state_dir / "dockable_669_1014_explicit_AB.pdb",
                prepared_receptor_dir / f"{state}_dockable_reference_explicit_AB.pdb",
                prepared_receptor_dir / f"{state}_dockable_669_1014_explicit_AB.pdb",
            ]
        else:
            dockable_candidates = [
                data_state_dir / "dockable_669_1014_explicit_AB.pdb",
                prepared_receptor_dir / f"{state}_dockable_669_1014_explicit_AB.pdb",
            ]
        if m2_runtime is not None:
            dockable_candidates.append(m2_runtime)
        dockable = _first_existing(dockable_candidates)
        record = {
            "mapping": {"status": "PASS" if mapping else "FAIL", "path": ctx.relative_to_repo(mapping) if mapping else ""},
            "dockable": {"status": "PASS" if dockable else "FAIL", "path": ctx.relative_to_repo(dockable) if dockable else ""},
            "runtime": {"status": "PASS" if runtime else "FAIL", "path": ctx.relative_to_repo(runtime) if runtime else ""},
        }
        receptor_inputs[state] = record
        if state in accepted_states:
            missing_required = []
            if not mapping:
                missing_required.append("receptor mapping")
            if not dockable:
                missing_required.append("dockable receptor PDB")
            if not runtime:
                missing_required.append("runtime-offset receptor PDB")
            if missing_required:
                blocker = "missing {0} for accepted pocket state {1}".format(", ".join(missing_required), state)
                blockers.append(blocker)
                _add_check(
                    checks,
                    f"M1_RECEPTOR_INPUTS_{state}",
                    "m1_inputs",
                    "mapping, dockable receptor PDB, and runtime-offset receptor PDB must exist for exported states",
                    ";".join(
                        [
                            ctx.relative_to_repo(data_state_dir / "receptor_mapping.csv"),
                            ctx.relative_to_repo(data_state_dir / ("dockable_reference_explicit_AB.pdb" if state == REFERENCE_STATE else "dockable_669_1014_explicit_AB.pdb")),
                            ctx.relative_to_repo(data_state_dir / "runtime_offset_receptor_only.pdb"),
                        ]
                    ),
                    None,
                    "FAIL",
                    True,
                    blocker,
                )
            else:
                _add_check(
                    checks,
                    f"M1_RECEPTOR_INPUTS_{state}",
                    "m1_inputs",
                    "mapping, dockable receptor PDB, and runtime-offset receptor PDB must exist for exported states",
                    ";".join([record["mapping"]["path"], record["dockable"]["path"], record["runtime"]["path"]]),
                    ";".join([record["mapping"]["path"], record["dockable"]["path"], record["runtime"]["path"]]),
                    "PASS",
                    False,
                    "accepted-state receptor inputs present",
                )
        elif not mapping:
            _add_check(
                checks,
                f"M1_MAPPING_{state}",
                "m1_inputs",
                "receptor mapping exists for non-exported provenance states when available",
                ctx.relative_to_repo(data_state_dir / "receptor_mapping.csv"),
                None,
                "WARN",
                False,
                f"missing receptor mapping for non-exported state {state}",
            )
        else:
            _add_check(checks, f"M1_MAPPING_{state}", "m1_inputs", "receptor mapping exists", ctx.relative_to_repo(mapping), ctx.relative_to_repo(mapping), "PASS", False, "mapping present")
    geometry_root = ctx.fresh_root / "data" / "normalized" / "geometry"
    membrane = _first_existing(
        [
            geometry_root / "membrane_frame.json",
            m2_dir / "manifest" / "membrane_frame.json",
            m2_dir / "qc" / "membrane_frame.json",
            m2_dir / "phase0_inputs" / "membrane_frame.json",
        ]
    )
    if not membrane:
        blocker = "missing membrane_frame.json"
        blockers.append(blocker)
        _add_check(checks, "M1_MEMBRANE_FRAME", "m1_inputs", "membrane_frame.json required for M3 readiness", ctx.relative_to_repo(geometry_root / "membrane_frame.json"), None, "FAIL", True, blocker)
    else:
        _add_check(checks, "M1_MEMBRANE_FRAME", "m1_inputs", "membrane_frame.json required for M3 readiness", ctx.relative_to_repo(membrane), ctx.relative_to_repo(membrane), "PASS", False, "membrane frame present")
    for optional_name in ["frame_qc.csv", "dimer_geometry_qc.csv"]:
        optional_path = geometry_root / optional_name
        if not optional_path.is_file():
            _add_check(checks, f"M1_{optional_name}", "m1_inputs", f"{optional_name} should be available for provenance", ctx.relative_to_repo(optional_path), None, "WARN", False, "optional geometry QC missing")
    return receptor_inputs


def _check_ligands(
    ctx: RunContext,
    mode: str,
    checks: list[dict[str, Any]],
    blockers: list[str],
    warnings: list[str],
) -> tuple[dict[str, dict[str, Any]], int, int]:
    ligand_status: dict[str, dict[str, Any]] = {}
    found = 0
    missing = 0
    for compound_id in COMPOUND_IDS:
        path = _find_ligand(ctx, compound_id)
        if path is None:
            missing += 1
            message = f"ligand input missing for {compound_id}"
            status = "FAIL" if mode == "docking" else "WARN"
            blocking = mode == "docking"
            if blocking:
                blockers.append(message)
            else:
                warnings.append(message)
            ligand_status[compound_id] = {"status": status, "path": None, "sha256": None}
            _add_check(
                checks,
                f"LIGAND_{compound_id}",
                "ligands",
                f"{compound_id} public ligand input present",
                f"fresh/data/raw/ligands/{compound_id}.* or run-local ligand input",
                None,
                status,
                blocking,
                message,
                "Provide public anonymized ligand input before docking mode.",
            )
        else:
            found += 1
            ligand_status[compound_id] = {
                "status": "PASS",
                "path": ctx.relative_to_repo(path),
                "sha256": _sha256(path),
            }
            _add_check(
                checks,
                f"LIGAND_{compound_id}",
                "ligands",
                f"{compound_id} public ligand input present",
                f"fresh/data/raw/ligands/{compound_id}.* or run-local ligand input",
                ctx.relative_to_repo(path),
                "PASS",
                False,
                "ligand file present and hashed",
            )
    return ligand_status, found, missing


def _write_markdown_report(
    path: Path,
    ctx: RunContext,
    mode: str,
    m2_run_id: str,
    status: str,
    m3_docking_allowed: bool,
    blockers: list[str],
    warnings: list[str],
    counts: dict[str, int],
    outputs: dict[str, Path],
) -> None:
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# M3-T0 Readiness Audit",
        "",
        f"- run_id: {ctx.run_id}",
        f"- m2_run_id: {m2_run_id}",
        f"- mode: {mode}",
        f"- overall_status: {status}",
        f"- m3_docking_allowed: {m3_docking_allowed}",
        f"- no_docking_attempted: true",
        "",
        "## Counts",
        "",
    ]
    for key in sorted(counts):
        lines.append(f"- {key}: {counts[key]}")
    lines.extend(["", "## Blockers", ""])
    lines.extend([f"- {item}" for item in blockers] or ["- none"])
    lines.extend(["", "## Warnings", ""])
    lines.extend([f"- {item}" for item in warnings] or ["- none"])
    lines.extend(
        [
            "",
            "## Non-goals Preserved",
            "",
            "- no ligand preparation",
            "- no ligand PDBQT generation",
            "- no receptor PDBQT generation",
            "- no Vina execution",
            "- no docking job generation",
            "- no compound scoring",
            "- no candidate nomination",
            "",
            "## Outputs",
            "",
        ]
    )
    for label, output_path in outputs.items():
        lines.append(f"- {label}: {ctx.relative_to_repo(output_path)}")
    lines.extend(["", "## Next Task", "", "M3-T1 - ligand manifest, confidentiality audit, ligand QC"])
    safe_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_m3_readiness(
    ctx: RunContext,
    m2_run_id: str,
    mode: str = "dry-run",
    create_skeleton: bool = True,
) -> M3ReadinessResult:
    if mode not in {"dry-run", "docking"}:
        raise ValueError("mode must be dry-run or docking")

    checks: list[dict[str, Any]] = []
    warnings: list[str] = []
    blockers: list[str] = []

    ctx.create_directories()
    skeleton_rows = _create_skeleton(ctx, create_skeleton=create_skeleton)
    _touch_logs(ctx)

    phase3_root = ctx.run_dir / "phase3_compounds"
    report_csv = phase3_root / "manifests" / "m3_readiness_report.csv"
    path_resolution_csv = phase3_root / "manifests" / "m3_input_path_resolution.csv"
    skeleton_qc_csv = phase3_root / "qc" / "m3_directory_skeleton_qc.csv"
    summary_json = phase3_root / "qc" / "m3_readiness_summary.json"
    report_md = phase3_root / "reports" / "m3_task0_readiness.md"
    phase3_log = ctx.logs_dir / "phase3_compounds.log"

    try:
        m2_dir = _safe_m2_dir(ctx, m2_run_id)
    except RunContextError as exc:
        m2_dir = ctx.fresh_root / "runs" / "invalid"
        blockers.append(str(exc))
        _add_check(checks, "M2_RUN_ID", "path_safety", "m2_run_id must resolve inside fresh/runs", f"fresh/runs/{m2_run_id}", None, "FAIL", True, str(exc))

    if not m2_dir.is_dir():
        blocker = f"M2 run directory not found: {m2_dir}"
        blockers.append(blocker)
        _add_check(checks, "M2_RUN_DIR", "m2_handoff", "M2 source run directory exists", ctx.relative_to_repo(m2_dir), None, "FAIL", True, blocker)
        resolved_inputs: dict[str, Path | None] = {}
        resolution_rows: list[dict[str, Any]] = []
    else:
        _add_check(checks, "M2_RUN_DIR", "m2_handoff", "M2 source run directory exists", ctx.relative_to_repo(m2_dir), ctx.relative_to_repo(m2_dir), "PASS", False, "M2 source run directory exists")
        _check_m2_review(ctx, m2_dir, mode, checks, warnings, blockers)
        resolved_inputs, resolution_rows = _path_inputs(ctx, m2_dir)

    accepted_rows, accepted_counts, accepted_states = _check_accepted_pockets(
        ctx,
        resolved_inputs.get("accepted_pockets_for_m3") if resolved_inputs else None,
        checks,
        warnings,
        blockers,
    )
    valid_box_rows = _check_boxes(
        ctx,
        resolved_inputs.get("accepted_pocket_boxes") if resolved_inputs else None,
        accepted_rows,
        checks,
        blockers,
    )
    _check_reference_csv(
        ctx,
        resolved_inputs.get("atp_reference") if resolved_inputs else None,
        "atp_site_reference.csv",
        "M2_ATP_REFERENCE",
        checks,
        blockers,
        ["uniprot_residue_number", "atp_centroid_x", "ligand_resname", "residue_name"],
    )
    _check_reference_csv(
        ctx,
        resolved_inputs.get("ppi_consensus_patch") if resolved_inputs else None,
        "ppi_consensus_patch.csv",
        "M2_PPI_CONSENSUS",
        checks,
        blockers,
        ["egfr_protomer_id", "egfr_chain_id", "egfr_residue_number", "state_id", "uniprot_residue_number"],
    )
    _check_parseable_csv(
        ctx,
        resolved_inputs.get("pocket_gate_qc") if resolved_inputs else None,
        checks,
        blockers,
        "M2_POCKET_GATE_QC",
        "m2_handoff",
        "pocket_gate_qc.csv",
        "phase2_pockets/export_for_m3/pocket_gate_qc.csv",
    )

    receptor_inputs = _check_receptor_inputs(ctx, m2_dir, accepted_states, checks, blockers, warnings)
    membrane_path = _first_existing(
        [
            ctx.fresh_root / "data" / "normalized" / "geometry" / "membrane_frame.json",
            m2_dir / "manifest" / "membrane_frame.json",
            m2_dir / "qc" / "membrane_frame.json",
            m2_dir / "phase0_inputs" / "membrane_frame.json",
        ]
    )
    if membrane_path is None:
        membrane_path = ctx.fresh_root / "data" / "normalized" / "geometry" / "membrane_frame.json"
    ligand_inputs, ligands_found, ligands_missing = _check_ligands(ctx, mode, checks, blockers, warnings)

    counts = {
        **accepted_counts,
        "valid_box_rows": valid_box_rows,
        "ligands_expected": len(COMPOUND_IDS),
        "ligands_found": ligands_found,
        "ligands_missing": ligands_missing,
    }
    status = "FAIL" if blockers else ("WARN" if warnings else "PASS")
    m3_docking_allowed = status == "PASS"
    if mode == "dry-run":
        m3_docking_allowed = False
    manual_review_required = status != "PASS" or accepted_counts.get("secondary_pocket_rows", 0) > 0

    inputs_summary = {
        "accepted_pockets_for_m3": {
            "status": "PASS" if resolved_inputs.get("accepted_pockets_for_m3") else "FAIL",
            "path": _as_repo_path(ctx, resolved_inputs.get("accepted_pockets_for_m3")),
        },
        "accepted_pocket_boxes": {
            "status": "PASS" if resolved_inputs.get("accepted_pocket_boxes") and valid_box_rows > 0 else "FAIL",
            "path": _as_repo_path(ctx, resolved_inputs.get("accepted_pocket_boxes")),
        },
        "atp_reference": {
            "status": "PASS" if resolved_inputs.get("atp_reference") else "FAIL",
            "path": _as_repo_path(ctx, resolved_inputs.get("atp_reference")),
        },
        "ppi_consensus_patch": {
            "status": "PASS" if resolved_inputs.get("ppi_consensus_patch") else "FAIL",
            "path": _as_repo_path(ctx, resolved_inputs.get("ppi_consensus_patch")),
        },
        "membrane_frame": {
            "status": "PASS" if membrane_path.is_file() else "FAIL",
            "path": ctx.relative_to_repo(membrane_path) if membrane_path.is_file() else None,
        },
        "receptor_mappings": receptor_inputs,
        "ligands": ligand_inputs,
    }

    summary = {
        "schema_version": "m3_readiness_summary_v1",
        "run_id": ctx.run_id,
        "m2_run_id": m2_run_id,
        "reviewed_at": now_iso(),
        "mode": mode,
        "create_skeleton": create_skeleton,
        "overall_status": status,
        "m3_docking_allowed": m3_docking_allowed,
        "manual_review_required": manual_review_required,
        "no_docking_attempted": True,
        "blockers": blockers,
        "warnings": warnings,
        "counts": counts,
        "inputs": inputs_summary,
        "non_goals_preserved": {
            "ligand_preparation": True,
            "ligand_pdbqt_generation": True,
            "receptor_pdbqt_generation": True,
            "vina_execution": True,
            "docking_job_generation": True,
            "compound_scoring": True,
            "candidate_nomination": True,
        },
    }

    _safe_write_csv(report_csv, READINESS_FIELDS, checks, ctx)
    _safe_write_csv(path_resolution_csv, PATH_RESOLUTION_FIELDS, resolution_rows, ctx)
    _safe_write_csv(skeleton_qc_csv, SKELETON_QC_FIELDS, skeleton_rows, ctx)
    summary_path = ctx.require_within_run_dir(summary_json)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_markdown_report(
        report_md,
        ctx,
        mode,
        m2_run_id,
        status,
        m3_docking_allowed,
        blockers,
        warnings,
        counts,
        {
            "m3_readiness_report": report_csv,
            "m3_readiness_summary": summary_json,
            "m3_task0_readiness": report_md,
            "input_path_resolution": path_resolution_csv,
            "directory_skeleton_qc": skeleton_qc_csv,
            "phase3_log": phase3_log,
        },
    )

    with ctx.require_within_run_dir(phase3_log).open("a", encoding="utf-8") as handle:
        handle.write(
            "{0} M3-T0 status={1} run_id={2} m2_run_id={3} mode={4} no_docking_attempted=true blockers={5} warnings={6}\n".format(
                now_iso(),
                status,
                ctx.run_id,
                m2_run_id,
                mode,
                len(blockers),
                len(warnings),
            )
        )
    log_master(ctx, status, "M3-T0 readiness audit completed")
    append_phase_status(
        ctx,
        "phase3_compounds",
        status if status != "WARN" else "WARN",
        "M3-T0 readiness audit and directory skeleton completed",
        {
            "task": "M3-T0",
            "run_id": ctx.run_id,
            "m2_run_id": m2_run_id,
            "mode": mode,
            "m3_docking_allowed": m3_docking_allowed,
            "no_docking_attempted": True,
        },
    )
    append_job_status(
        ctx,
        "M3-T0-readiness-audit",
        status if status != "WARN" else "WARN",
        details={
            "message": "readiness audit only; no external engines invoked",
            "m3_docking_allowed": m3_docking_allowed,
            "no_docking_attempted": True,
        },
    )

    return M3ReadinessResult(
        status=status,
        m3_docking_allowed=m3_docking_allowed,
        manual_review_required=manual_review_required,
        blockers=blockers,
        warnings=warnings,
        report_csv=report_csv,
        summary_json=summary_json,
        report_md=report_md,
        path_resolution_csv=path_resolution_csv,
        skeleton_qc_csv=skeleton_qc_csv,
        phase3_log=phase3_log,
        counts=counts,
    )
