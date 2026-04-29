"""Task 7 PPI consensus patch intake and report writer."""

import csv
import json
from pathlib import Path

from egfr_myo1d.analysis.ppi_consensus import (
    CONSENSUS_FIELDS,
    REQUIRED_COLUMNS,
    build_consensus_patches,
    load_contact_table,
    normalize_contact_rows,
)
from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import initialize_manifests, now_iso, write_json
from egfr_myo1d.core.run_context import RunContextError
from egfr_myo1d.io.hashing import describe_file


SCHEMA_FIELDS = ["column", "required", "present", "status", "notes"]
PARSING_FIELDS = ["row_index", "pose_id", "column", "status", "parsed_count", "errors"]
TAIL_FIELDS = ["pose_id", "tail_noise_contact_count", "terminal_artifact_flag", "quarantine_recommended", "status", "notes"]
ACTIVE_FIELDS = ["pose_id", "active_face_contact_count", "sheet12_support_count", "score_bonus_allowed", "key_residue_bonus_weight", "status", "notes"]
MEMBRANE_FIELDS = ["pose_id", "atp_overlap_flag", "membrane_proximal_flag", "membrane_side_class", "status", "notes"]
CONVERGENCE_FIELDS = [
    "ppi_patch_id",
    "receptor_id",
    "receptor_state",
    "protomer_id",
    "n_accepted_poses",
    "unique_egfr_residue_count",
    "max_residue_occupancy",
    "centroid_spread_angstrom",
    "evidence_class",
    "status",
    "notes",
]

NON_GOALS_CONFIRMED = [
    "no_docking_execution",
    "no_pyrosetta_execution",
    "no_lightdock_vina_fpocket_p2rank_execution",
    "no_scoring_execution",
    "no_qsub_pbs_sbatch_generation",
    "no_candidate_nomination",
    "no_egfr_mutation_or_receptor_redesign",
]


def _safe_input_path(ctx, input_root, relative_or_absolute):
    root = input_root if input_root.is_absolute() else (ctx.repo_root / input_root)
    root = root.resolve()
    if relative_or_absolute is None:
        path = root / "accepted_ppi_contacts.csv"
    else:
        candidate = Path(relative_or_absolute)
        path = candidate if candidate.is_absolute() else root / candidate
    path = path.resolve()
    try:
        path.relative_to(root)
    except ValueError:
        raise RunContextError("PPI consensus contact table escapes input root: {0}".format(path))
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


def _output_record(ctx, path):
    record = {"path": ctx.relative_to_repo(path)}
    record.update(describe_file(path))
    return record


def _status_from_classes(patch_rows):
    if not patch_rows:
        return "FAIL"
    warn_classes = set(["DISPERSED", "TAIL_DOMINATED", "ATP_OVERLAPPING", "MEMBRANE_PROXIMAL", "INSUFFICIENT_EVIDENCE"])
    for row in patch_rows:
        if row.get("evidence_class") in warn_classes:
            return "WARN"
    return "PASS"


def _report_text(report, patch_rows):
    lines = [
        "Task 7 PPI consensus patch report",
        "run_id: {0}".format(report["run_id"]),
        "status: {0}".format(report["status"]),
        "input_kind: {0}".format(report["input_kind"]),
        "contact_table: {0}".format(report["contact_table"]["path"]),
        "patch_count: {0}".format(len(patch_rows)),
        "",
        "This is consensus evidence from supplied PPI contact records only.",
        "It is not a docking result, not a pocket-discovery result, and not a candidate nomination.",
        "All consensus evidence requires visual and scientific review.",
        "",
    ]
    for row in patch_rows:
        lines.append(
            "{0}: {1} {2} protomer={3} class={4} residues={5}".format(
                row["ppi_patch_id"],
                row["receptor_id"],
                row["receptor_state"],
                row["protomer_id"],
                row["evidence_class"],
                row["egfr_consensus_residues"],
            )
        )
    return "\n".join(lines) + "\n"


def summarize_ppi_consensus(ctx, mode, profile, input_root, contact_table=None, input_kind="synthetic_fixture"):
    initialize_manifests(ctx, mode)
    root, table_path = _safe_input_path(ctx, input_root, contact_table)
    manifest_path = ctx.manifest_dir / "ppi_consensus_input_manifest.json"
    report_path = ctx.manifest_dir / "ppi_consensus_qc_report.json"
    schema_path = ctx.qc_dir / "ppi_consensus_schema_audit.csv"
    parsing_path = ctx.qc_dir / "ppi_residue_parsing_audit.csv"
    tail_path = ctx.qc_dir / "ppi_tail_noise_audit.csv"
    active_path = ctx.qc_dir / "ppi_active_face_audit.csv"
    membrane_path = ctx.qc_dir / "ppi_membrane_atp_audit.csv"
    convergence_path = ctx.qc_dir / "ppi_convergence_audit.csv"
    output_path = ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv"
    text_report_path = ctx.reports_dir / "ppi_consensus_report.txt"

    missing_schema = []
    schema_rows = []
    raw_rows = []
    parse_errors = []
    normalized = []
    parsing_rows = []
    active_rows = []
    tail_rows = []
    membrane_rows = []
    patch_rows = []
    convergence_rows = []

    if not table_path.exists():
        missing_schema = list(REQUIRED_COLUMNS)
        schema_rows = [{
            "column": column,
            "required": True,
            "present": False,
            "status": "FAIL",
            "notes": "contact table missing",
        } for column in REQUIRED_COLUMNS]
        status = "FAIL"
    else:
        missing_schema, schema_rows, raw_rows = load_contact_table(table_path)
        if missing_schema:
            status = "FAIL"
        else:
            normalized, parsing_rows, active_rows, tail_rows, membrane_rows, parse_errors = normalize_contact_rows(raw_rows)
            if parse_errors:
                status = "FAIL"
            else:
                patch_rows, convergence_rows = build_consensus_patches(normalized)
                status = _status_from_classes(patch_rows)

    _write_csv(schema_path, schema_rows, SCHEMA_FIELDS, ctx)
    _write_csv(parsing_path, parsing_rows, PARSING_FIELDS, ctx)
    _write_csv(tail_path, tail_rows, TAIL_FIELDS, ctx)
    _write_csv(active_path, active_rows, ACTIVE_FIELDS, ctx)
    _write_csv(membrane_path, membrane_rows, MEMBRANE_FIELDS, ctx)
    _write_csv(convergence_path, convergence_rows, CONVERGENCE_FIELDS, ctx)
    if status != "FAIL":
        _write_csv(output_path, patch_rows, CONSENSUS_FIELDS, ctx)

    outputs = {
        "schema_audit": ctx.relative_to_repo(schema_path),
        "residue_parsing_audit": ctx.relative_to_repo(parsing_path),
        "tail_noise_audit": ctx.relative_to_repo(tail_path),
        "active_face_audit": ctx.relative_to_repo(active_path),
        "membrane_atp_audit": ctx.relative_to_repo(membrane_path),
        "convergence_audit": ctx.relative_to_repo(convergence_path),
        "consensus_patch_csv": ctx.relative_to_repo(output_path) if status != "FAIL" else None,
        "text_report": ctx.relative_to_repo(text_report_path),
    }
    warnings = []
    for row in patch_rows:
        if row.get("warnings"):
            warnings.append({"ppi_patch_id": row["ppi_patch_id"], "warnings": row["warnings"], "evidence_class": row["evidence_class"]})
    blockers = []
    if missing_schema:
        blockers.append({"code": "missing_required_schema", "message": ";".join(missing_schema)})
    if parse_errors:
        blockers.append({"code": "residue_or_coordinate_parse_error", "message": ";".join(parse_errors)})

    manifest = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": "summarize-ppi-consensus",
        "mode": mode,
        "profile": profile,
        "input_root": str(root),
        "input_kind": input_kind,
        "contact_table": _input_record(ctx, table_path),
        "required_columns": REQUIRED_COLUMNS,
        "outputs": outputs,
        "no_docking_executed": True,
        "no_scoring_executed": True,
        "score_bonus_allowed": False,
    }
    report = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "status": status,
        "mode": mode,
        "profile": profile,
        "input_kind": input_kind,
        "contact_table": _input_record(ctx, table_path),
        "row_count": len(raw_rows),
        "patch_count": len(patch_rows),
        "warnings": warnings,
        "blockers": blockers,
        "outputs": outputs,
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
        "score_bonus_allowed": False,
        "active_face_and_sheet12_are_annotation_only": True,
        "tail_noise_residues_are_not_primary_binding_evidence": True,
    }

    write_json(manifest_path, manifest, ctx)
    write_json(report_path, report, ctx)
    _write_text(text_report_path, _report_text(report, patch_rows), ctx)
    phase_status = "FAIL" if status == "FAIL" else ("WARN" if status == "WARN" else "PASS")
    append_phase_status(
        ctx,
        "summarize-ppi-consensus",
        phase_status,
        "PPI consensus patch summary completed with {0}".format(status),
        {"mode": mode, "profile": profile, "input_kind": input_kind, "patch_count": len(patch_rows)},
    )
    append_job_status(
        ctx,
        "summarize_ppi_consensus_local",
        phase_status,
        details={"message": "local pure-Python PPI contact-table post-processing"},
    )
    return report

