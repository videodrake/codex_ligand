"""Task 8 PPI-guided pocket-discovery planning/intake writer."""

import csv
import json
from pathlib import Path

from egfr_myo1d.analysis.pocket_selection import (
    CONSENSUS_REQUIRED_COLUMNS,
    POCKET_FIELDS,
    build_pocket_records,
    load_consensus_patch_table,
    validate_consensus_schema,
)
from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import initialize_manifests, now_iso, write_json
from egfr_myo1d.core.run_context import RunContextError
from egfr_myo1d.io.hashing import describe_file


SELECTION_FIELDS = [
    "pocket_plan_id",
    "ppi_patch_id",
    "source_evidence_class",
    "accepted_for_pocket_planning",
    "future_discovery_allowed",
    "future_compound_docking_planned",
    "status",
    "notes",
]
ATP_FIELDS = ["pocket_plan_id", "ppi_patch_id", "atp_overlap_fraction", "atp_relation_class", "status", "notes"]
MEMBRANE_FIELDS = ["pocket_plan_id", "ppi_patch_id", "membrane_proximal_fraction", "membrane_accessibility_class", "status", "notes"]
DIMER_FIELDS = ["pocket_plan_id", "ppi_patch_id", "protomer_id", "dimer_accessibility_class", "status", "notes"]

NON_GOALS_CONFIRMED = [
    "no_docking_execution",
    "no_vina_pyrosetta_lightdock_execution",
    "no_fpocket_or_p2rank_runtime_execution",
    "no_qsub_pbs_sbatch_generation",
    "no_compound_scoring",
    "no_candidate_nomination",
    "no_egfr_mutation_or_residue_renumbering",
]


def _safe_consensus_path(ctx, input_root, consensus_path):
    root = input_root if input_root.is_absolute() else (ctx.repo_root / input_root)
    root = root.resolve()
    if consensus_path is None:
        path = root / "ppi_consensus_patch.csv"
    else:
        candidate = Path(consensus_path)
        path = candidate if candidate.is_absolute() else root / candidate
    path = path.resolve()
    try:
        path.relative_to(root)
    except ValueError:
        raise RunContextError("PPI consensus patch path escapes input root: {0}".format(path))
    return root, path


def _write_csv(path, rows, fieldnames, ctx):
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    with safe_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_text(path, text, ctx):
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    safe_path.write_text(text, encoding="utf-8")


def _input_record(ctx, path):
    record = {"path": ctx.relative_to_repo(path)}
    record.update(describe_file(path))
    return record


def _summary_text(report):
    lines = [
        "# Task 8 Pocket Discovery Planning Summary",
        "",
        "- status: {0}".format(report["status"]),
        "- input consensus: {0}".format(report["ppi_consensus_patch"]["path"]),
        "- accepted PPI patches: {0}".format(report["accepted_ppi_patch_count"]),
        "- planned pocket records: {0}".format(report["planned_pocket_count"]),
        "- record semantics: {0}".format(report["record_semantics"]),
        "- detected pocket records: {0}".format(report["detected_pocket_record"]),
        "- schema audit: {0}".format(report["outputs"].get("task7_consensus_schema_audit", "")),
        "- planned docking jobs: 0",
        "",
        "Task 8 is pocket-selection planning/intake only.",
        "No pocket detector, docking engine, scoring, qsub/PBS job, or candidate nomination was run.",
        "",
    ]
    if report.get("warnings"):
        lines.append("Warnings:")
        for item in report["warnings"]:
            lines.append("- {0}: {1}".format(item.get("ppi_patch_id", "plan"), item.get("message", item.get("warning", ""))))
    if report.get("blockers"):
        lines.append("Blockers:")
        for item in report["blockers"]:
            lines.append("- {0}: {1}".format(item.get("code", ""), item.get("message", "")))
    return "\n".join(lines) + "\n"


def plan_pocket_discovery(ctx, mode, profile, input_root, consensus_path=None, input_kind="task7_consensus"):
    initialize_manifests(ctx, mode)
    root, consensus_file = _safe_consensus_path(ctx, input_root, consensus_path)
    pockets_dir = ctx.require_within_run_dir(ctx.run_dir / "pockets")
    pockets_dir.mkdir(parents=True, exist_ok=True)

    pocket_csv = pockets_dir / "egfr_myo1d_ppi_adjacent_pockets.csv"
    pocket_alias_csv = pockets_dir / "egfr_myo1d_ppi_guided_pocket_plan_records.csv"
    plan_json = pockets_dir / "pocket_discovery_plan.json"
    schema_audit_path = ctx.qc_dir / "task7_consensus_schema_audit.csv"
    selection_audit_path = ctx.qc_dir / "pocket_selection_audit.csv"
    atp_audit_path = ctx.qc_dir / "atp_overlap_audit.csv"
    membrane_audit_path = ctx.qc_dir / "membrane_accessibility_audit.csv"
    dimer_audit_path = ctx.qc_dir / "dimer_accessibility_audit.csv"
    manifest_path = ctx.manifest_dir / "pocket_discovery_manifest.json"
    report_path = ctx.reports_dir / "pocket_discovery_summary.md"

    blockers = []
    warnings = []
    consensus_rows = []
    schema_audit = []
    if not consensus_file.exists():
        blockers.append({"code": "missing_task7_consensus_patch", "message": "Task 7 ppi_consensus_patch.csv is missing"})
        schema_audit = validate_consensus_schema([])[1]
        schema_audit.insert(0, {
            "check": "input_file",
            "field": "ppi_consensus_patch.csv",
            "status": "FAIL",
            "details": "missing Task 7 consensus patch table",
        })
        status = "NO_GO"
    else:
        missing, consensus_rows, fieldnames = load_consensus_patch_table(consensus_file)
        missing_schema, schema_audit = validate_consensus_schema(fieldnames)
        if missing:
            blockers.append({"code": "invalid_task7_consensus_schema", "message": ";".join(missing)})
            status = "NO_GO"
        else:
            status = "PASS"

    pocket_rows = []
    selection_audit = []
    atp_audit = []
    membrane_audit = []
    dimer_audit = []
    accepted_count = 0
    if status != "NO_GO":
        pocket_rows, selection_audit, atp_audit, membrane_audit, dimer_audit, accepted_count = build_pocket_records(consensus_rows)
        if accepted_count == 0:
            status = "NO_GO"
            warnings.append({"ppi_patch_id": "all", "message": "no accepted Task 7 PPI patch available for pocket planning"})
        for row in pocket_rows:
            if row.get("warnings"):
                warnings.append({"ppi_patch_id": row.get("ppi_patch_id", ""), "message": row.get("warnings", "")})
    else:
        selection_audit = [{
            "pocket_plan_id": "",
            "ppi_patch_id": "",
            "source_evidence_class": "SCHEMA_INVALID",
            "accepted_for_pocket_planning": False,
            "future_discovery_allowed": False,
            "future_compound_docking_planned": False,
            "status": "FAIL",
            "notes": ";".join([item["message"] for item in blockers]),
        }]

    _write_csv(pocket_csv, pocket_rows if status != "NO_GO" else [], POCKET_FIELDS, ctx)
    _write_csv(pocket_alias_csv, pocket_rows if status != "NO_GO" else [], POCKET_FIELDS, ctx)
    _write_csv(schema_audit_path, schema_audit, ["check", "field", "status", "details"], ctx)
    _write_csv(selection_audit_path, selection_audit, SELECTION_FIELDS, ctx)
    _write_csv(atp_audit_path, atp_audit, ATP_FIELDS, ctx)
    _write_csv(membrane_audit_path, membrane_audit, MEMBRANE_FIELDS, ctx)
    _write_csv(dimer_audit_path, dimer_audit, DIMER_FIELDS, ctx)

    outputs = {
        "pockets_csv": ctx.relative_to_repo(pocket_csv),
        "legacy_pockets_csv": ctx.relative_to_repo(pocket_csv),
        "pocket_plan_records_csv": ctx.relative_to_repo(pocket_alias_csv),
        "pocket_discovery_plan": ctx.relative_to_repo(plan_json),
        "task7_consensus_schema_audit": ctx.relative_to_repo(schema_audit_path),
        "pocket_selection_audit": ctx.relative_to_repo(selection_audit_path),
        "atp_overlap_audit": ctx.relative_to_repo(atp_audit_path),
        "membrane_accessibility_audit": ctx.relative_to_repo(membrane_audit_path),
        "dimer_accessibility_audit": ctx.relative_to_repo(dimer_audit_path),
        "manifest": ctx.relative_to_repo(manifest_path),
        "summary_report": ctx.relative_to_repo(report_path),
    }
    plan = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": "plan-pocket-discovery",
        "status": status,
        "mode": mode,
        "profile": profile,
        "input_kind": input_kind,
        "ppi_consensus_patch": _input_record(ctx, consensus_file),
        "task7_consensus_schema_audit": schema_audit,
        "record_semantics": "ppi_guided_pocket_plan",
        "detected_pocket_record": False,
        "accepted_ppi_patch_count": accepted_count,
        "planned_pocket_count": len(pocket_rows) if status != "NO_GO" else 0,
        "planned_docking_job_count": 0,
        "pocket_detector_runtime_executed": False,
        "compound_docking_or_scoring_executed": False,
        "candidate_nomination_executed": False,
        "outputs": outputs,
        "warnings": warnings,
        "blockers": blockers,
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
    }
    manifest = dict(plan)
    manifest["input_root"] = str(root)
    write_json(plan_json, plan, ctx)
    write_json(manifest_path, manifest, ctx)
    _write_text(report_path, _summary_text(plan), ctx)
    phase_status = "WARN" if status == "NO_GO" or warnings else "PASS"
    append_phase_status(
        ctx,
        "plan-pocket-discovery",
        phase_status,
        "PPI-guided pocket selection plan completed with {0}".format(status),
        {"mode": mode, "profile": profile, "accepted_ppi_patch_count": accepted_count},
    )
    append_job_status(
        ctx,
        "plan_pocket_discovery_local",
        phase_status,
        details={"message": "local pure-Python pocket-selection planning only"},
    )
    return plan
