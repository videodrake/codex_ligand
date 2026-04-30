"""M2.3 PPI pose collection and residue mapping restoration scaffold.

This stage collects M2.2 dry-run job status records and optionally restores a
run-local raw contact table through the M1 receptor mapping CSV. It does not
execute docking, parse new scientific outputs from engines, score poses, or
classify accepted poses.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError, ensure_within
from egfr_myo1d.ppi.restore_residue_mapping import (
    MAPPING_QC_FIELDS,
    POSE_CONTACT_FIELDS,
    UNMAPPED_CONTACT_FIELDS,
    MappingRestorationResult,
    restore_contact_table_mapping,
)


M2_3_PHASE = "m2.3-ppi-pose-collection-mapping-restoration"
HARNESS_ROOT = "phase1_ppi/pyrosetta_adapter"
TABLE_ROOT = "phase1_ppi/tables"
RAW_CONTACT_REQUIRED_FIELDS = [
    "state_id",
    "seed",
    "pose_id",
    "protomer_id",
    "egfr_runtime_residue",
    "myo1d_residue",
    "contact_type",
]
RAW_POSE_FIELDS = [
    "run_id",
    "state_id",
    "state_role",
    "method",
    "seed",
    "job_name",
    "pose_id",
    "score",
    "pose_path",
    "source_status_path",
    "status",
    "notes",
]
NON_GOALS_CONFIRMED = [
    "no_pyrosetta_import",
    "no_pyrosetta_docking_or_relaxation",
    "no_lightdock_execution",
    "no_vina_execution",
    "no_fpocket_or_p2rank_runtime",
    "no_compound_docking_or_scoring",
    "no_pose_acceptance_classification",
    "no_candidate_nomination",
    "no_scheduler_submission",
    "no_cleanup_deletion",
]


@dataclass
class M2PpiCollectionReport:
    run_id: str
    mode: str
    profile: str
    status: str
    raw_pose_count: int = 0
    restored_contact_count: int = 0
    unmapped_contact_count: int = 0
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    manifest_path: Path | None = None
    raw_pose_table: Path | None = None
    pose_contacts_csv: Path | None = None
    mapping_qc_csv: Path | None = None
    unmapped_contacts_csv: Path | None = None
    summary_report_path: Path | None = None


def _repo_path(ctx: RunContext, path: Path) -> str:
    return ctx.relative_to_repo(path)


def _resolve_repo_or_abs(ctx: RunContext, text: str | Path) -> Path:
    path = Path(text)
    if path.is_absolute():
        return path.resolve()
    return (ctx.repo_root / path).resolve()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            encoded = {}
            for key, value in row.items():
                if isinstance(value, (dict, list)):
                    encoded[key] = json.dumps(value, sort_keys=True)
                else:
                    encoded[key] = value
            writer.writerow(encoded)


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def _read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _default_job_manifest(ctx: RunContext) -> Path:
    return ctx.run_dir / HARNESS_ROOT / "pyrosetta_job_manifest.jsonl"


def _load_job_manifest(ctx: RunContext, job_manifest: Path | None) -> tuple[list[dict[str, Any]], Path, list[str]]:
    path = job_manifest or _default_job_manifest(ctx)
    path = _resolve_repo_or_abs(ctx, path)
    ensure_within(path, ctx.run_dir)
    if not path.is_file():
        return [], path, ["missing_m2_2_job_manifest:{0}".format(path)]
    return _read_jsonl(path), path, []


def _state_roles_from_m2_1_manifest(ctx: RunContext) -> dict[str, str]:
    path = ctx.manifest_dir / "m2_1_ppi_input_manifest.json"
    if not path.is_file():
        return {}
    payload = _load_json(path)
    return {
        str(pack.get("state_id", "")): str(pack.get("role", ""))
        for pack in payload.get("packs", [])
    }


def _mapping_paths_by_state(ctx: RunContext, jobs: list[dict[str, Any]]) -> tuple[dict[str, Path], list[str]]:
    mapping_paths: dict[str, Path] = {}
    warnings: list[str] = []
    for job in jobs:
        state_id = str(job.get("state_id", ""))
        mapping_text = str(job.get("mapping_csv", ""))
        if not state_id or not mapping_text:
            warnings.append("job_missing_state_or_mapping:{0}".format(job.get("job_name", "")))
            continue
        mapping_path = _resolve_repo_or_abs(ctx, mapping_text)
        try:
            ensure_within(mapping_path, ctx.run_dir)
        except RunContextError:
            warnings.append("mapping_csv_escapes_run_dir:{0}:{1}".format(state_id, mapping_path))
            continue
        if not mapping_path.is_file():
            warnings.append("mapping_csv_missing:{0}:{1}".format(state_id, mapping_path))
            continue
        mapping_paths[state_id] = mapping_path
    return mapping_paths, warnings


def _raw_pose_rows(ctx: RunContext, jobs: list[dict[str, Any]], state_roles: dict[str, str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for job in jobs:
        state_id = str(job.get("state_id", ""))
        output_dir = _resolve_repo_or_abs(ctx, str(job.get("output_dir", "")))
        try:
            output_dir = ctx.require_within_run_dir(output_dir)
        except RunContextError:
            rows.append(
                {
                    "run_id": ctx.run_id,
                    "state_id": state_id,
                    "state_role": state_roles.get(state_id, ""),
                    "method": str(job.get("method", "")),
                    "seed": str(job.get("seed", "")),
                    "job_name": str(job.get("job_name", "")),
                    "pose_id": "",
                    "score": "",
                    "pose_path": "",
                    "source_status_path": "",
                    "status": "WARN",
                    "notes": "job_output_dir_escapes_run_dir",
                }
            )
            continue
        dry_run_status = output_dir / "dry_run_status.json"
        if dry_run_status.is_file():
            status = "DRY_RUN_ONLY"
            status_path = _repo_path(ctx, dry_run_status)
            notes = "dry_run_status_only_no_pose"
        else:
            status = "NO_POSE_OUTPUT"
            status_path = ""
            notes = "no_pose_or_dry_run_status_found"
        rows.append(
            {
                "run_id": ctx.run_id,
                "state_id": state_id,
                "state_role": state_roles.get(state_id, ""),
                "method": str(job.get("method", "")),
                "seed": str(job.get("seed", "")),
                "job_name": str(job.get("job_name", "")),
                "pose_id": "",
                "score": "",
                "pose_path": "",
                "source_status_path": status_path,
                "status": status,
                "notes": notes,
            }
        )
    return rows


def _pose_rows_from_contacts(
    ctx: RunContext,
    raw_contact_rows: list[dict[str, Any]],
    state_roles: dict[str, str],
) -> list[dict[str, Any]]:
    seen: set[tuple[str, str, str]] = set()
    rows: list[dict[str, Any]] = []
    for row in raw_contact_rows:
        key = (
            str(row.get("state_id", "")),
            str(row.get("seed", "")),
            str(row.get("pose_id", "")),
        )
        if key in seen:
            continue
        seen.add(key)
        state_id, seed, pose_id = key
        rows.append(
            {
                "run_id": str(row.get("run_id") or ctx.run_id),
                "state_id": state_id,
                "state_role": str(row.get("state_role") or state_roles.get(state_id, "")),
                "method": str(row.get("method", "raw_contact_table")),
                "seed": seed,
                "job_name": str(row.get("job_name", "")),
                "pose_id": pose_id,
                "score": str(row.get("score", "")),
                "pose_path": str(row.get("pose_path", "")),
                "source_status_path": "",
                "status": "CONTACT_TABLE_ONLY",
                "notes": "contact_table_supplied_for_mapping_restoration_only",
            }
        )
    return rows


def _validate_raw_contact_schema(rows: list[dict[str, Any]], path: Path) -> list[str]:
    if not rows:
        return []
    present = set(rows[0].keys())
    missing = [field for field in RAW_CONTACT_REQUIRED_FIELDS if field not in present]
    if missing:
        return ["raw_contact_table_missing_columns:{0}:{1}".format(path, ";".join(missing))]
    return []


def _load_raw_contact_rows(ctx: RunContext, raw_contact_table: Path | None) -> tuple[list[dict[str, Any]], Path | None, list[str]]:
    if raw_contact_table is None:
        return [], None, []
    path = _resolve_repo_or_abs(ctx, raw_contact_table)
    ensure_within(path, ctx.run_dir)
    if not path.is_file():
        return [], path, ["raw_contact_table_missing:{0}".format(path)]
    rows = _read_csv(path)
    return rows, path, _validate_raw_contact_schema(rows, path)


def _write_summary(ctx: RunContext, report: M2PpiCollectionReport) -> Path:
    target = ctx.require_within_run_dir(ctx.reports_dir / "m2_3_ppi_collection_mapping.md")
    lines = [
        "# M2.3 PPI Pose Collection and Mapping Restoration",
        "",
        "Status: {0}".format(report.status),
        "",
        "This phase records M2.2 dry-run job status and restores run-local raw contact tables through the M1 receptor mapping CSV.",
        "No docking, relaxation, pose scoring, pose acceptance classification, pocket discovery, or candidate nomination was run.",
        "",
        "Raw pose/status rows: {0}".format(report.raw_pose_count),
        "Restored contact rows: {0}".format(report.restored_contact_count),
        "Unmapped contact rows: {0}".format(report.unmapped_contact_count),
        "",
        "## Non-goals Confirmed",
    ]
    for item in NON_GOALS_CONFIRMED:
        lines.append("- {0}".format(item))
    if report.warnings:
        lines.extend(["", "## Warnings"])
        for warning in report.warnings:
            lines.append("- {0}".format(warning))
    if report.blockers:
        lines.extend(["", "## Blockers"])
        for blocker in report.blockers:
            lines.append("- {0}".format(blocker))
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def _status_from_parts(blockers: list[str], warnings: list[str], restoration: MappingRestorationResult) -> str:
    if blockers or restoration.status == "FAIL":
        return "FAIL"
    if warnings or restoration.status == "WARN":
        return "PASS_WITH_WARNINGS"
    return "PASS"


def collect_ppi_outputs(
    ctx: RunContext,
    mode: str = "smoke_input",
    profile: str = "codex_dev",
    job_manifest: Path | None = None,
    raw_contact_table: Path | None = None,
    unmapped_fraction_fail_threshold: float = 0.0,
) -> M2PpiCollectionReport:
    """Collect M2.2 dry-run outputs and restore optional raw contacts."""

    ctx.create_directories()
    tables_dir = ctx.require_within_run_dir(ctx.run_dir / TABLE_ROOT)
    tables_dir.mkdir(parents=True, exist_ok=True)

    jobs, job_manifest_path, job_blockers = _load_job_manifest(ctx, job_manifest)
    raw_contacts, raw_contact_path, raw_contact_blockers = _load_raw_contact_rows(
        ctx, raw_contact_table
    )
    state_roles = _state_roles_from_m2_1_manifest(ctx)
    mapping_paths, mapping_warnings = _mapping_paths_by_state(ctx, jobs)
    for row in raw_contacts:
        if not row.get("state_role"):
            row["state_role"] = state_roles.get(str(row.get("state_id", "")), "")

    blockers = list(job_blockers) + list(raw_contact_blockers)
    warnings = list(mapping_warnings)
    if raw_contact_table is None:
        warnings.append("no_raw_contact_table_supplied_mapping_restoration_header_only")

    raw_pose_rows = _raw_pose_rows(ctx, jobs, state_roles)
    raw_pose_rows.extend(_pose_rows_from_contacts(ctx, raw_contacts, state_roles))

    if blockers:
        restoration = MappingRestorationResult(
            status="FAIL",
            restored_rows=[],
            unmapped_rows=[],
            qc_rows=[
                {
                    "state_id": "",
                    "mapping_csv": "",
                    "input_contact_count": len(raw_contacts),
                    "restored_contact_count": 0,
                    "unmapped_contact_count": 0,
                    "unmapped_fraction": "0.000000",
                    "status": "FAIL",
                    "notes": "collection_blocked_before_mapping_restoration",
                }
            ],
        )
    else:
        restoration = restore_contact_table_mapping(
            run_id=ctx.run_id,
            raw_contact_rows=raw_contacts,
            mapping_csv_by_state=mapping_paths,
            unmapped_fraction_fail_threshold=unmapped_fraction_fail_threshold,
        )
        warnings.extend(restoration.warnings)
        blockers.extend(restoration.blockers)

    raw_pose_table = tables_dir / "ppi_raw_pose_table.csv"
    pose_contacts_csv = tables_dir / "ppi_pose_contacts.csv"
    mapping_qc_csv = tables_dir / "ppi_mapping_restoration_qc.csv"
    unmapped_contacts_csv = tables_dir / "unmapped_contacts.csv"
    _write_csv(raw_pose_table, raw_pose_rows, RAW_POSE_FIELDS, ctx)
    _write_csv(pose_contacts_csv, restoration.restored_rows, POSE_CONTACT_FIELDS, ctx)
    _write_csv(mapping_qc_csv, restoration.qc_rows, MAPPING_QC_FIELDS, ctx)
    _write_csv(
        unmapped_contacts_csv,
        restoration.unmapped_rows,
        UNMAPPED_CONTACT_FIELDS,
        ctx,
    )

    status = _status_from_parts(blockers, warnings, restoration)
    report = M2PpiCollectionReport(
        run_id=ctx.run_id,
        mode=mode,
        profile=profile,
        status=status,
        raw_pose_count=len(raw_pose_rows),
        restored_contact_count=len(restoration.restored_rows),
        unmapped_contact_count=len(restoration.unmapped_rows),
        warnings=warnings,
        blockers=blockers,
        raw_pose_table=raw_pose_table,
        pose_contacts_csv=pose_contacts_csv,
        mapping_qc_csv=mapping_qc_csv,
        unmapped_contacts_csv=unmapped_contacts_csv,
    )
    summary = _write_summary(ctx, report)
    report.summary_report_path = summary

    manifest = ctx.manifest_dir / "m2_3_ppi_collection_manifest.json"
    payload = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": M2_3_PHASE,
        "mode": mode,
        "profile": profile,
        "status": status,
        "source_job_manifest": _repo_path(ctx, job_manifest_path),
        "source_raw_contact_table": _repo_path(ctx, raw_contact_path) if raw_contact_path else None,
        "execution_allowed": False,
        "pyrosetta_imported": False,
        "docking_executed": False,
        "relaxation_executed": False,
        "contact_extraction_from_pose_pdb_executed": False,
        "pose_acceptance_classification_executed": False,
        "score_bonus_allowed": False,
        "raw_pose_count": len(raw_pose_rows),
        "restored_contact_count": len(restoration.restored_rows),
        "unmapped_contact_count": len(restoration.unmapped_rows),
        "unmapped_fraction_fail_threshold": unmapped_fraction_fail_threshold,
        "outputs": {
            "raw_pose_table": _repo_path(ctx, raw_pose_table),
            "pose_contacts_csv": _repo_path(ctx, pose_contacts_csv),
            "mapping_qc_csv": _repo_path(ctx, mapping_qc_csv),
            "unmapped_contacts_csv": _repo_path(ctx, unmapped_contacts_csv),
            "summary_report": _repo_path(ctx, summary),
        },
        "warnings": warnings,
        "blockers": blockers,
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
    }
    write_json(manifest, payload, ctx)
    report.manifest_path = manifest

    phase_status = "FAIL" if status == "FAIL" else (
        "WARN" if status == "PASS_WITH_WARNINGS" else "PASS"
    )
    append_phase_status(
        ctx,
        M2_3_PHASE,
        phase_status,
        "M2.3 PPI collection/mapping restoration {0}".format(status),
        {
            "raw_pose_count": len(raw_pose_rows),
            "restored_contact_count": len(restoration.restored_rows),
            "unmapped_contact_count": len(restoration.unmapped_rows),
            "docking_executed": False,
        },
    )
    append_job_status(
        ctx,
        "m2_3_ppi_collection_mapping_local",
        phase_status,
        details={
            "message": "local dry-run collection and mapping restoration only",
            "execution_allowed": False,
        },
    )
    return report
