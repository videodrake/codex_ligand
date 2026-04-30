"""M2.8 accepted pocket export and Milestone 2 report generation."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError


M2_8_PHASE = "m2.8-accepted-pocket-export-and-report"
FINAL_ROOT = "phase2_pockets/final"
EXPORT_ROOT = "phase2_pockets/export_for_m3"

ACCEPTED_M3_FIELDS = [
    "pocket_family_id",
    "export_tier",
    "gate_class",
    "acceptance_tier",
    "representative_candidate_pocket_id",
    "representative_state_id",
    "representative_state_role",
    "representative_protomer_id",
    "member_candidate_pocket_ids",
    "member_state_ids",
    "member_state_roles",
    "member_protomer_ids",
    "primary_state_count",
    "reference_state_count",
    "hard_gate_pass",
    "non_atp_pass",
    "ppi_relationship_pass",
    "lower_lateral_pass",
    "dimer_accessibility_pass",
    "primary_state_support_or_review_flag",
    "ppi_relationship_class",
    "atp_overlap_fraction",
    "atp_centroid_distance",
    "z_percentile_within_receptor_surface",
    "lateral_score",
    "interface_overlap_fraction",
    "sasa_ratio_dimer_vs_isolated_protomer",
    "family_center_x",
    "family_center_y",
    "family_center_z",
    "residue_union_uniprot",
    "residue_union_state_residue_ids",
    "box_export_status",
    "m3_primary_allowed",
    "m3_review_allowed",
    "reference_only",
    "acceptance_reason",
    "review_notes",
    "warnings",
]

BOX_FIELDS = [
    "pocket_family_id",
    "box_id",
    "export_tier",
    "state_id",
    "state_role",
    "protomer_id",
    "chain_id",
    "receptor_pdb",
    "receptor_mapping_csv",
    "box_center_x",
    "box_center_y",
    "box_center_z",
    "box_size_x",
    "box_size_y",
    "box_size_z",
    "box_margin_angstrom",
    "box_source",
    "box_qc_status",
    "box_qc_notes",
    "non_atp_pass",
    "ppi_relationship_class",
    "lower_lateral_pass",
    "dimer_accessibility_pass",
    "reference_only",
    "m3_primary_allowed",
    "warnings",
]

STATE_MAP_FIELDS = [
    "pocket_family_id",
    "state_id",
    "state_role",
    "protomer_id",
    "chain_id",
    "receptor_pdb",
    "receptor_mapping_csv",
    "member_candidate_pocket_ids",
    "primary_state_support",
    "reference_state_support",
    "allowed_for_primary_m3",
    "allowed_for_secondary_review_m3",
    "mapping_status",
    "warnings",
]

EVIDENCE_FIELDS = [
    "pocket_family_id",
    "gate_class",
    "export_tier",
    "primary_state_count",
    "reference_state_count",
    "ppi_patch_ids",
    "ppi_relationship_class",
    "residue_overlap_count",
    "minimum_heavy_atom_distance_to_ppi",
    "pocket_centroid_to_ppi_centroid_distance",
    "atp_overlap_fraction",
    "atp_centroid_distance",
    "lower_gate_class",
    "z_percentile_within_receptor_surface",
    "lateral_score",
    "interface_overlap_fraction",
    "sasa_ratio_dimer_vs_isolated_protomer",
    "hard_gate_pass",
    "m3_primary_allowed",
    "m3_review_allowed",
    "evidence_status",
    "reject_or_review_reason",
    "warnings",
]

TRANSITION_GATE_FIELDS = [
    "gate_id",
    "gate_name",
    "required_for_m3",
    "status",
    "evidence_file",
    "evidence_count",
    "failure_reason",
    "warnings",
]

NON_GOALS_CONFIRMED = [
    "no_milestone3_started",
    "no_vina_execution",
    "no_ligand_preparation",
    "no_receptor_pdbqt_preparation",
    "no_pyrosetta_docking_or_relaxation",
    "no_lightdock_execution",
    "no_p2rank_execution",
    "no_compound_scoring",
    "no_candidate_nomination",
    "no_scheduler_submission",
    "no_cleanup_deletion",
]


@dataclass
class M2ExportReport:
    run_id: str
    mode: str
    profile: str
    status: str
    transition_gate_status: str
    total_family_count: int = 0
    exported_primary_count: int = 0
    exported_secondary_count: int = 0
    reference_only_count: int = 0
    rejected_count: int = 0
    boxes_generated_count: int = 0
    boxes_failed_count: int = 0
    m3_docking_allowed: bool = False
    manual_review_required: bool = False
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    final_summary_md: Path | None = None
    final_evidence_csv: Path | None = None
    final_accepted_csv: Path | None = None
    transition_gate_csv: Path | None = None
    transition_gate_status_json: Path | None = None
    final_status_json: Path | None = None
    export_accepted_csv: Path | None = None
    export_boxes_csv: Path | None = None
    export_residue_sets_json: Path | None = None
    export_state_map_csv: Path | None = None
    export_gate_qc_csv: Path | None = None
    export_ppi_consensus_csv: Path | None = None
    export_report_md: Path | None = None
    export_manifest_json: Path | None = None
    cleanup_report_json: Path | None = None
    reports_summary_md: Path | None = None


def export_m2_results(
    ctx: RunContext,
    mode: str = "smoke_input",
    profile: str = "codex_dev",
    ppi_consensus: Path | None = None,
    atp_reference: Path | None = None,
    pocket_families: Path | None = None,
    pocket_gate_qc: Path | None = None,
    accepted_families: Path | None = None,
    rejected_families: Path | None = None,
    synthetic_fixture: bool = False,
) -> M2ExportReport:
    """Build M2.8 final aggregation and M3 handoff metadata."""

    ctx.create_directories()
    production = mode == "production" and not synthetic_fixture
    final_dir = ctx.require_within_run_dir(ctx.run_dir / FINAL_ROOT)
    export_dir = ctx.require_within_run_dir(ctx.run_dir / EXPORT_ROOT)
    final_dir.mkdir(parents=True, exist_ok=True)
    export_dir.mkdir(parents=True, exist_ok=True)

    ppi_path = _resolve_optional_path(ctx, ppi_consensus) or _discover_ppi_consensus(ctx)
    atp_path = _resolve_optional_path(ctx, atp_reference) or _discover_atp_reference(ctx)
    families_path = _resolve_optional_path(ctx, pocket_families) or _discover_pocket_families(ctx)
    gate_qc_path = _resolve_optional_path(ctx, pocket_gate_qc) or _discover_gate_qc(ctx)
    accepted_path = _resolve_optional_path(ctx, accepted_families) or _discover_accepted_families(ctx)
    rejected_path = _resolve_optional_path(ctx, rejected_families) or _discover_rejected_families(ctx)
    membrane_frame_path = _discover_membrane_frame(ctx)

    warnings: list[str] = []
    blockers: list[str] = []
    if synthetic_fixture:
        warnings.append("synthetic_fixture_not_production_evidence_m3_docking_not_allowed")
    _required_input("missing_ppi_consensus_patch_blocks_m2_export", ppi_path, warnings, blockers, production)
    _required_input("missing_atp_reference_blocks_m2_export", atp_path, warnings, blockers, production)
    _required_input("missing_merged_pocket_candidates_blocks_m2_export", families_path, warnings, blockers, production)
    _required_input("missing_pocket_gate_qc_blocks_m2_export", gate_qc_path, warnings, blockers, production)
    _required_input("missing_accepted_pocket_family_table_blocks_m2_export", accepted_path, warnings, blockers, production)
    _required_input("missing_m1_membrane_frame_blocks_m2_export", membrane_frame_path, warnings, blockers, production)
    if production and (rejected_path is None or not rejected_path.is_file()):
        blockers.append("missing_rejected_pocket_family_table_blocks_m2_export")
    elif rejected_path is None or not rejected_path.is_file():
        warnings.append("missing_rejected_pocket_family_table_synthetic_or_smoke_not_production_evidence")

    ppi_rows = _read_csv(ppi_path) if ppi_path and ppi_path.is_file() else []
    atp_rows = _read_csv(atp_path) if atp_path and atp_path.is_file() else []
    family_rows = _read_csv(families_path) if families_path and families_path.is_file() else []
    gate_rows = _read_csv(gate_qc_path) if gate_qc_path and gate_qc_path.is_file() else []
    accepted_rows = _read_csv(accepted_path) if accepted_path and accepted_path.is_file() else []
    rejected_rows = _read_csv(rejected_path) if rejected_path and rejected_path.is_file() else []

    gate_by_family = {str(row.get("pocket_family_id", "")): row for row in gate_rows}
    class_rows = list(accepted_rows)
    export_rows = _build_export_rows(class_rows, gate_by_family)
    if not any(row.get("export_tier") in {"primary", "secondary_review"} for row in export_rows):
        warnings.append("no_accepted_primary_or_secondary_pockets_blocks_m3_docking_pending_manual_review")

    state_map_rows, state_map_warnings, state_map_blockers = _build_state_map_rows(ctx, export_rows, production)
    warnings.extend(state_map_warnings)
    blockers.extend(state_map_blockers)
    box_rows = [_build_box_row(ctx, row) for row in export_rows]
    box_by_family = {row.get("pocket_family_id", ""): row for row in box_rows}
    for row in export_rows:
        box = box_by_family.get(row.get("pocket_family_id", ""))
        row["box_export_status"] = box.get("box_qc_status", "FAIL") if box else "FAIL"

    residue_sets = _build_residue_sets(ctx, export_rows, state_map_rows)
    evidence_rows = _build_evidence_rows(gate_rows, export_rows, ppi_rows)
    transition_rows = _build_transition_gates(
        ppi_path,
        ppi_rows,
        atp_path,
        atp_rows,
        families_path,
        family_rows,
        gate_qc_path,
        gate_rows,
        accepted_path,
        accepted_rows,
        export_rows,
        production,
    )
    transition_status = _transition_status(transition_rows)
    m3_docking_allowed = (
        transition_status != "FAIL"
        and not blockers
        and not synthetic_fixture
        and any(row.get("export_tier") in {"primary", "secondary_review"} for row in export_rows)
    )
    manual_review_required = not m3_docking_allowed or any(row.get("export_tier") == "secondary_review" for row in export_rows)

    _write_csv(final_dir / "ppi_to_pocket_evidence_table.csv", evidence_rows, EVIDENCE_FIELDS, ctx)
    _write_csv(final_dir / "accepted_pocket_families_for_compound_docking.csv", export_rows, ACCEPTED_M3_FIELDS, ctx)
    _write_csv(final_dir / "m2_transition_gate_qc.csv", transition_rows, TRANSITION_GATE_FIELDS, ctx)
    _write_csv(export_dir / "accepted_pockets_for_m3.csv", export_rows, ACCEPTED_M3_FIELDS, ctx)
    _write_csv(export_dir / "accepted_pocket_boxes.csv", box_rows, BOX_FIELDS, ctx)
    _write_csv(export_dir / "accepted_pocket_receptor_state_map.csv", state_map_rows, STATE_MAP_FIELDS, ctx)
    _write_csv(export_dir / "pocket_gate_qc.csv", gate_rows, _fieldnames_or_default(gate_rows, []), ctx)
    _write_csv(export_dir / "ppi_consensus_patch.csv", ppi_rows, _fieldnames_or_default(ppi_rows, []), ctx)
    write_json(export_dir / "accepted_pocket_residue_sets.json", residue_sets, ctx)

    cleanup_payload = _cleanup_payload(ctx)
    write_json(export_dir / "m2_cleanup_report.json", cleanup_payload, ctx)

    counts = {
        "total_pocket_families": len(gate_rows) if gate_rows else len(family_rows),
        "accepted_primary_pocket": _count_export(export_rows, "primary"),
        "accepted_secondary_cryptic": _count_export(export_rows, "secondary_review"),
        "reference_only": _count_export(export_rows, "reference_only"),
        "rejected": max(0, (len(gate_rows) if gate_rows else len(family_rows)) - len(export_rows)),
        "exported_primary_for_m3": _count_export(export_rows, "primary"),
        "exported_secondary_review_for_m3": _count_export(export_rows, "secondary_review"),
    }
    status = _overall_status(blockers, warnings, transition_status, m3_docking_allowed)
    source_files = _source_files(
        ctx,
        ppi_path=ppi_path,
        atp_path=atp_path,
        families_path=families_path,
        gate_qc_path=gate_qc_path,
        accepted_path=accepted_path,
        rejected_path=rejected_path,
        membrane_frame_path=membrane_frame_path,
    )
    output_files = {
        "final_summary_md": _repo_path(ctx, final_dir / "milestone2_summary.md"),
        "final_evidence_csv": _repo_path(ctx, final_dir / "ppi_to_pocket_evidence_table.csv"),
        "final_accepted_csv": _repo_path(ctx, final_dir / "accepted_pocket_families_for_compound_docking.csv"),
        "transition_gate_qc": _repo_path(ctx, final_dir / "m2_transition_gate_qc.csv"),
        "transition_gate_status": _repo_path(ctx, final_dir / "m2_transition_gate_status.json"),
        "m2_final_status": _repo_path(ctx, final_dir / "m2_final_status.json"),
        "export_accepted_pockets": _repo_path(ctx, export_dir / "accepted_pockets_for_m3.csv"),
        "export_boxes": _repo_path(ctx, export_dir / "accepted_pocket_boxes.csv"),
        "export_residue_sets": _repo_path(ctx, export_dir / "accepted_pocket_residue_sets.json"),
        "export_receptor_state_map": _repo_path(ctx, export_dir / "accepted_pocket_receptor_state_map.csv"),
        "export_gate_qc": _repo_path(ctx, export_dir / "pocket_gate_qc.csv"),
        "export_ppi_consensus": _repo_path(ctx, export_dir / "ppi_consensus_patch.csv"),
        "export_report": _repo_path(ctx, export_dir / "m2_final_report.md"),
        "export_manifest": _repo_path(ctx, export_dir / "m2_export_manifest.json"),
        "cleanup_report": _repo_path(ctx, export_dir / "m2_cleanup_report.json"),
    }
    manifest_payload = {
        "schema_version": "m2_8_export_manifest_v1",
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "source_files": source_files,
        "output_files": output_files,
        "counts": counts,
        "transition_gate_status": transition_status,
        "m3_docking_allowed": m3_docking_allowed,
        "warnings": sorted(set(warnings)),
    }
    transition_payload = {
        "schema_version": "m2_8_transition_gate_status_v1",
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "status": transition_status,
        "m3_docking_allowed": m3_docking_allowed,
        "manual_review_required": manual_review_required,
        "gates": transition_rows,
        "warnings": sorted(set(warnings)),
        "blockers": blockers,
    }
    final_status_payload = {
        "schema_version": "m2_8_final_status_v1",
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": M2_8_PHASE,
        "mode": mode,
        "profile": profile,
        "synthetic_fixture": synthetic_fixture,
        "status": status,
        "transition_gate_status": transition_status,
        "counts": counts,
        "m3_docking_allowed": m3_docking_allowed,
        "manual_review_required": manual_review_required,
        "vina_executed": False,
        "ligand_preparation_executed": False,
        "receptor_pdbqt_preparation_executed": False,
        "pyrosetta_docking_or_relaxation_executed": False,
        "lightdock_executed": False,
        "p2rank_executed": False,
        "compound_scoring_executed": False,
        "candidate_nomination_executed": False,
        "milestone3_started": False,
        "cleanup_mode": "report_only",
        "warnings": sorted(set(warnings)),
        "blockers": blockers,
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
        "source_files": source_files,
        "output_files": output_files,
    }
    write_json(final_dir / "m2_transition_gate_status.json", transition_payload, ctx)
    write_json(final_dir / "m2_final_status.json", final_status_payload, ctx)
    write_json(export_dir / "m2_export_manifest.json", manifest_payload, ctx)

    report = M2ExportReport(
        run_id=ctx.run_id,
        mode=mode,
        profile=profile,
        status=status,
        transition_gate_status=transition_status,
        total_family_count=counts["total_pocket_families"],
        exported_primary_count=counts["exported_primary_for_m3"],
        exported_secondary_count=counts["exported_secondary_review_for_m3"],
        reference_only_count=counts["reference_only"],
        rejected_count=counts["rejected"],
        boxes_generated_count=sum(1 for row in box_rows if row.get("box_qc_status") == "PASS"),
        boxes_failed_count=sum(1 for row in box_rows if row.get("box_qc_status") != "PASS"),
        m3_docking_allowed=m3_docking_allowed,
        manual_review_required=manual_review_required,
        warnings=sorted(set(warnings)),
        blockers=blockers,
        final_summary_md=final_dir / "milestone2_summary.md",
        final_evidence_csv=final_dir / "ppi_to_pocket_evidence_table.csv",
        final_accepted_csv=final_dir / "accepted_pocket_families_for_compound_docking.csv",
        transition_gate_csv=final_dir / "m2_transition_gate_qc.csv",
        transition_gate_status_json=final_dir / "m2_transition_gate_status.json",
        final_status_json=final_dir / "m2_final_status.json",
        export_accepted_csv=export_dir / "accepted_pockets_for_m3.csv",
        export_boxes_csv=export_dir / "accepted_pocket_boxes.csv",
        export_residue_sets_json=export_dir / "accepted_pocket_residue_sets.json",
        export_state_map_csv=export_dir / "accepted_pocket_receptor_state_map.csv",
        export_gate_qc_csv=export_dir / "pocket_gate_qc.csv",
        export_ppi_consensus_csv=export_dir / "ppi_consensus_patch.csv",
        export_report_md=export_dir / "m2_final_report.md",
        export_manifest_json=export_dir / "m2_export_manifest.json",
        cleanup_report_json=export_dir / "m2_cleanup_report.json",
        reports_summary_md=ctx.reports_dir / "milestone2_summary.md",
    )
    _write_report(ctx, report, source_files, output_files, counts, transition_rows, export_rows, gate_rows)
    write_json(ctx.manifest_dir / "m2_8_export_manifest.json", manifest_payload, ctx)

    phase_status = "FAIL" if status == "FAIL" else ("WARN" if status == "PASS_WITH_WARNINGS" else "PASS")
    append_phase_status(
        ctx,
        M2_8_PHASE,
        phase_status,
        "M2.8 export/report {0}: exported_primary={1}, exported_secondary={2}, m3_allowed={3}".format(
            status,
            report.exported_primary_count,
            report.exported_secondary_count,
            m3_docking_allowed,
        ),
        {
            "transition_gate_status": transition_status,
            "warnings": report.warnings,
            "blockers": blockers,
        },
    )
    append_job_status(
        ctx,
        "m2_8_export_report",
        phase_status,
        details={
            "status": status,
            "output_manifest": _repo_path(ctx, export_dir / "m2_export_manifest.json"),
            "m3_docking_allowed": m3_docking_allowed,
        },
    )
    return report


def validate_box_row(row: dict[str, Any]) -> tuple[str, str]:
    """Validate numeric M2.8 box metadata."""

    center_fields = ["box_center_x", "box_center_y", "box_center_z"]
    size_fields = ["box_size_x", "box_size_y", "box_size_z"]
    for field in center_fields:
        value = _safe_float(row.get(field))
        if value is None:
            return "FAIL", "box_center_not_numeric:{0}".format(field)
    for field in size_fields:
        value = _safe_float(row.get(field))
        if value is None:
            return "FAIL", "box_size_not_numeric:{0}".format(field)
        if value <= 0.0:
            return "FAIL", "box_size_not_positive:{0}".format(field)
    if not row.get("pocket_family_id"):
        return "FAIL", "box_missing_pocket_family_id"
    if not row.get("state_id"):
        return "FAIL", "box_missing_state_id"
    if not row.get("protomer_id"):
        return "FAIL", "box_missing_protomer_id"
    if _as_bool(row.get("m3_primary_allowed")) and row.get("state_id") == "3GT8_raw":
        return "FAIL", "3GT8_raw_box_not_primary_allowed"
    if not _as_bool(row.get("non_atp_pass")) and not _as_bool(row.get("reference_only")):
        return "FAIL", "atp_failed_box_not_allowed_for_primary_export"
    return "PASS", "box_metadata_valid"


def _repo_path(ctx: RunContext, path: Path) -> str:
    return ctx.relative_to_repo(path)


def _is_within(child: Path, parent: Path) -> bool:
    try:
        Path(child).resolve().relative_to(Path(parent).resolve())
        return True
    except ValueError:
        return False


def _resolve_optional_path(ctx: RunContext, value: Path | None) -> Path | None:
    if value is None:
        return None
    raw = Path(value)
    path = raw.resolve() if raw.is_absolute() else (ctx.repo_root / raw).resolve()
    if not (_is_within(path, ctx.fresh_root) or _is_within(path, ctx.run_dir)):
        path = (ctx.run_dir / raw).resolve()
    if not (_is_within(path, ctx.fresh_root) or _is_within(path, ctx.run_dir)):
        raise RunContextError("input path escapes fresh/ and run dir: {0}".format(path))
    return path


def _discover_ppi_consensus(ctx: RunContext) -> Path | None:
    candidates = [
        ctx.run_dir / "phase1_ppi" / "consensus" / "ppi_consensus_patch_merged.csv",
        ctx.run_dir / "phase1_ppi" / "tables" / "ppi_consensus_patch.csv",
        ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _discover_atp_reference(ctx: RunContext) -> Path | None:
    candidates = [
        ctx.run_dir / "phase2_pockets" / "atp_reference" / "atp_site_reference.csv",
        ctx.run_dir / "phase2_pockets" / "tables" / "atp_site_reference.csv",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _discover_pocket_families(ctx: RunContext) -> Path | None:
    candidates = [
        ctx.run_dir / "phase2_pockets" / "merged" / "pocket_candidates_merged.csv",
        ctx.run_dir / "phase2_pockets" / "tables" / "pocket_candidates_merged.csv",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _discover_gate_qc(ctx: RunContext) -> Path | None:
    candidates = [
        ctx.run_dir / "phase2_pockets" / "gated" / "pocket_gate_qc.csv",
        ctx.run_dir / "phase2_pockets" / "tables" / "pocket_gate_qc.csv",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _discover_accepted_families(ctx: RunContext) -> Path | None:
    candidates = [
        ctx.run_dir / "phase2_pockets" / "gated" / "accepted_pocket_families.csv",
        ctx.run_dir / "phase2_pockets" / "tables" / "accepted_pockets_for_m3.csv",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _discover_rejected_families(ctx: RunContext) -> Path | None:
    candidates = [
        ctx.run_dir / "phase2_pockets" / "gated" / "rejected_pocket_families.csv",
        ctx.run_dir / "phase2_pockets" / "tables" / "pocket_rejects.csv",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _discover_membrane_frame(ctx: RunContext) -> Path | None:
    candidates = [
        ctx.run_dir / "manifest" / "membrane_frame.json",
        ctx.run_dir / "normalized" / "geometry" / "membrane_frame.json",
        ctx.fresh_root / "data" / "normalized" / "geometry" / "membrane_frame.json",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _required_input(message: str, path: Path | None, warnings: list[str], blockers: list[str], production: bool) -> None:
    if path is not None and path.is_file():
        return
    if production:
        blockers.append(message)
    else:
        warnings.append(message.replace("_blocks_m2_export", "_synthetic_or_smoke_not_production_evidence"))


def _read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    if not fieldnames:
        fieldnames = ["status"]
        rows = [{"status": ""}] if rows else []
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _fieldnames_or_default(rows: list[dict[str, Any]], default: list[str]) -> list[str]:
    if rows:
        return list(rows[0].keys())
    return default or ["status"]


def _build_export_rows(class_rows: list[dict[str, Any]], gate_by_family: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for source in class_rows:
        family_id = str(source.get("pocket_family_id", ""))
        gate = gate_by_family.get(family_id, {})
        gate_class = str(gate.get("gate_class") or source.get("gate_class", ""))
        hard_gate_pass = _as_bool(gate.get("hard_gate_pass"))
        primary_count = _int(gate.get("primary_state_count") or source.get("primary_state_count"))
        reference_count = _int(gate.get("reference_state_count") or source.get("reference_state_count"))
        export_tier = _export_tier(gate_class, hard_gate_pass, primary_count, reference_count)
        if export_tier is None:
            continue
        reference_only = export_tier == "reference_only"
        m3_primary_allowed = export_tier == "primary" and hard_gate_pass and not reference_only
        m3_review_allowed = export_tier in {"primary", "secondary_review"} and hard_gate_pass and not reference_only
        row = {
            "pocket_family_id": family_id,
            "export_tier": export_tier,
            "gate_class": gate_class,
            "acceptance_tier": export_tier,
            "representative_candidate_pocket_id": source.get("representative_candidate_pocket_id", ""),
            "representative_state_id": source.get("representative_state_id", ""),
            "representative_state_role": source.get("representative_state_role", ""),
            "representative_protomer_id": source.get("representative_protomer_id", ""),
            "member_candidate_pocket_ids": source.get("member_candidate_pocket_ids", ""),
            "member_state_ids": source.get("member_state_ids", ""),
            "member_state_roles": source.get("member_state_roles", ""),
            "member_protomer_ids": source.get("member_protomer_ids", ""),
            "primary_state_count": primary_count,
            "reference_state_count": reference_count,
            "hard_gate_pass": hard_gate_pass,
            "non_atp_pass": gate.get("non_atp_pass", source.get("non_atp_pass", "")),
            "ppi_relationship_pass": gate.get("ppi_relationship_pass", source.get("ppi_relationship_pass", "")),
            "lower_lateral_pass": gate.get("lower_lateral_pass", source.get("lower_lateral_pass", "")),
            "dimer_accessibility_pass": gate.get("dimer_accessibility_pass", source.get("dimer_accessibility_pass", "")),
            "primary_state_support_or_review_flag": gate.get("primary_state_support_or_review_flag", source.get("primary_state_support_or_review_flag", "")),
            "ppi_relationship_class": gate.get("ppi_relationship_class", source.get("ppi_relationship_class", "")),
            "atp_overlap_fraction": gate.get("atp_overlap_fraction", ""),
            "atp_centroid_distance": gate.get("atp_centroid_distance", ""),
            "z_percentile_within_receptor_surface": gate.get("z_percentile_within_receptor_surface", ""),
            "lateral_score": gate.get("lateral_score", ""),
            "interface_overlap_fraction": gate.get("interface_overlap_fraction", ""),
            "sasa_ratio_dimer_vs_isolated_protomer": gate.get("sasa_ratio_dimer_vs_isolated_protomer", ""),
            "family_center_x": source.get("family_center_x", gate.get("family_center_x", "")),
            "family_center_y": source.get("family_center_y", gate.get("family_center_y", "")),
            "family_center_z": source.get("family_center_z", gate.get("family_center_z", "")),
            "residue_union_uniprot": source.get("residue_union_uniprot", ""),
            "residue_union_state_residue_ids": source.get("residue_union_uniprot", ""),
            "box_export_status": "",
            "m3_primary_allowed": m3_primary_allowed,
            "m3_review_allowed": m3_review_allowed,
            "reference_only": reference_only,
            "acceptance_reason": source.get("acceptance_reason", gate.get("review_reason", "")),
            "review_notes": source.get("review_notes", gate.get("reject_reason", "")),
            "warnings": ";".join(_dedupe(_split_values(source.get("warnings", "")) + _split_values(gate.get("warnings", "")))),
        }
        rows.append(row)
    return sorted(rows, key=lambda item: (str(item.get("export_tier")), str(item.get("pocket_family_id"))))


def _export_tier(gate_class: str, hard_gate_pass: bool, primary_count: int, reference_count: int) -> str | None:
    if gate_class == "accepted_primary_pocket" and hard_gate_pass and primary_count >= 2:
        return "primary"
    if gate_class == "accepted_secondary_cryptic" and hard_gate_pass and primary_count >= 1:
        return "secondary_review"
    return None


def _build_state_map_rows(
    ctx: RunContext,
    export_rows: list[dict[str, Any]],
    production: bool,
) -> tuple[list[dict[str, Any]], list[str], list[str]]:
    rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    blockers: list[str] = []
    for export in export_rows:
        states = _split_values(export.get("member_state_ids", "")) or [str(export.get("representative_state_id", ""))]
        roles = _split_values(export.get("member_state_roles", "")) or [str(export.get("representative_state_role", ""))]
        protomers = _split_values(export.get("member_protomer_ids", "")) or [str(export.get("representative_protomer_id", ""))]
        for index, state_id in enumerate(states):
            state_role = roles[index] if index < len(roles) else (roles[0] if roles else "")
            protomer_id = protomers[index] if index < len(protomers) else (protomers[0] if protomers else "")
            receptor_pdb = _state_receptor_path(ctx, state_id)
            mapping_csv = _mapping_path(ctx, state_id)
            row_warnings: list[str] = []
            mapping_status = "PASS"
            if receptor_pdb is None:
                row_warnings.append("missing_receptor_pdb:{0}".format(state_id))
            if mapping_csv is None:
                row_warnings.append("missing_receptor_mapping:{0}".format(state_id))
                mapping_status = "MISSING"
            if production and row_warnings:
                blockers.extend("{0}_blocks_m2_export".format(item) for item in row_warnings)
            warnings.extend(row_warnings)
            rows.append(
                {
                    "pocket_family_id": export.get("pocket_family_id", ""),
                    "state_id": state_id,
                    "state_role": state_role,
                    "protomer_id": protomer_id,
                    "chain_id": protomer_id,
                    "receptor_pdb": _repo_path(ctx, receptor_pdb) if receptor_pdb else "",
                    "receptor_mapping_csv": _repo_path(ctx, mapping_csv) if mapping_csv else "",
                    "member_candidate_pocket_ids": export.get("member_candidate_pocket_ids", ""),
                    "primary_state_support": str(state_role).lower() == "primary",
                    "reference_state_support": state_id == "3GT8_raw" or str(state_role).lower() == "reference",
                    "allowed_for_primary_m3": _as_bool(export.get("m3_primary_allowed")) and state_id != "3GT8_raw",
                    "allowed_for_secondary_review_m3": _as_bool(export.get("m3_review_allowed")) and state_id != "3GT8_raw",
                    "mapping_status": mapping_status,
                    "warnings": ";".join(row_warnings),
                }
            )
    return rows, sorted(set(warnings)), sorted(set(blockers))


def _build_box_row(ctx: RunContext, export: dict[str, Any]) -> dict[str, Any]:
    center = [
        _safe_float(export.get("family_center_x")),
        _safe_float(export.get("family_center_y")),
        _safe_float(export.get("family_center_z")),
    ]
    state_id = str(export.get("representative_state_id", ""))
    protomer_id = str(export.get("representative_protomer_id", ""))
    receptor_pdb = _state_receptor_path(ctx, state_id)
    mapping_csv = _mapping_path(ctx, state_id)
    margin = 4.0
    size = 12.0
    row = {
        "pocket_family_id": export.get("pocket_family_id", ""),
        "box_id": "box_{0}".format(export.get("pocket_family_id", "")),
        "export_tier": export.get("export_tier", ""),
        "state_id": state_id,
        "state_role": export.get("representative_state_role", ""),
        "protomer_id": protomer_id,
        "chain_id": protomer_id,
        "receptor_pdb": _repo_path(ctx, receptor_pdb) if receptor_pdb else "",
        "receptor_mapping_csv": _repo_path(ctx, mapping_csv) if mapping_csv else "",
        "box_center_x": _format_float(center[0]),
        "box_center_y": _format_float(center[1]),
        "box_center_z": _format_float(center[2]),
        "box_size_x": _format_float(size),
        "box_size_y": _format_float(size),
        "box_size_z": _format_float(size),
        "box_margin_angstrom": _format_float(margin),
        "box_source": "family_center_minimum_box_m2_8_metadata_only",
        "box_qc_status": "",
        "box_qc_notes": "",
        "non_atp_pass": export.get("non_atp_pass", ""),
        "ppi_relationship_class": export.get("ppi_relationship_class", ""),
        "lower_lateral_pass": export.get("lower_lateral_pass", ""),
        "dimer_accessibility_pass": export.get("dimer_accessibility_pass", ""),
        "reference_only": export.get("reference_only", ""),
        "m3_primary_allowed": export.get("m3_primary_allowed", ""),
        "warnings": export.get("warnings", ""),
    }
    status, notes = validate_box_row(row)
    row["box_qc_status"] = status
    row["box_qc_notes"] = notes
    return row


def _build_residue_sets(
    ctx: RunContext,
    export_rows: list[dict[str, Any]],
    state_map_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    state_map_by_family: dict[str, list[dict[str, Any]]] = {}
    for row in state_map_rows:
        state_map_by_family.setdefault(str(row.get("pocket_family_id", "")), []).append(row)
    pockets: list[dict[str, Any]] = []
    for export in export_rows:
        residues = _split_values(export.get("residue_union_uniprot", ""))
        residue_sets_by_state: dict[str, Any] = {}
        for state_row in state_map_by_family.get(str(export.get("pocket_family_id", "")), []):
            state_id = str(state_row.get("state_id", ""))
            protomer_id = str(state_row.get("protomer_id", ""))
            chain_id = str(state_row.get("chain_id", ""))
            residue_sets_by_state.setdefault(
                state_id,
                {
                    "state_role": state_row.get("state_role", ""),
                    "protomer_residue_sets": {},
                },
            )
            residue_sets_by_state[state_id]["protomer_residue_sets"][protomer_id] = [
                {
                    "chain_id": chain_id,
                    "residue_id": residue,
                    "uniprot_residue": residue,
                    "resname": "",
                    "mapping_status": state_row.get("mapping_status", ""),
                }
                for residue in residues
            ]
        pockets.append(
            {
                "pocket_family_id": export.get("pocket_family_id", ""),
                "export_tier": export.get("export_tier", ""),
                "gate_class": export.get("gate_class", ""),
                "reference_only": _as_bool(export.get("reference_only")),
                "residue_union_uniprot": residues,
                "residue_sets_by_state": residue_sets_by_state,
                "ppi_relationship_class": export.get("ppi_relationship_class", ""),
                "non_atp_pass": _as_bool(export.get("non_atp_pass")),
                "lower_lateral_pass": _as_bool(export.get("lower_lateral_pass")),
                "dimer_accessibility_pass": _as_bool(export.get("dimer_accessibility_pass")),
                "warnings": _split_values(export.get("warnings", "")),
            }
        )
    return {
        "schema_version": "m2_8_accepted_pocket_residue_sets_v1",
        "run_id": ctx.run_id,
        "pockets": pockets,
    }


def _build_evidence_rows(
    gate_rows: list[dict[str, Any]],
    export_rows: list[dict[str, Any]],
    ppi_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    export_by_family = {str(row.get("pocket_family_id", "")): row for row in export_rows}
    ppi_patch_ids = ";".join(_dedupe(str(row.get("ppi_patch_id", "")) for row in ppi_rows if row.get("ppi_patch_id")))
    rows: list[dict[str, Any]] = []
    for gate in gate_rows:
        export = export_by_family.get(str(gate.get("pocket_family_id", "")), {})
        export_tier = export.get("export_tier", "not_exported")
        reason = gate.get("review_reason") or gate.get("reject_reason") or export.get("review_notes", "")
        rows.append(
            {
                "pocket_family_id": gate.get("pocket_family_id", ""),
                "gate_class": gate.get("gate_class", ""),
                "export_tier": export_tier,
                "primary_state_count": gate.get("primary_state_count", ""),
                "reference_state_count": gate.get("reference_state_count", ""),
                "ppi_patch_ids": ppi_patch_ids,
                "ppi_relationship_class": gate.get("ppi_relationship_class", ""),
                "residue_overlap_count": gate.get("residue_overlap_count", ""),
                "minimum_heavy_atom_distance_to_ppi": gate.get("minimum_heavy_atom_distance_to_ppi", ""),
                "pocket_centroid_to_ppi_centroid_distance": gate.get("pocket_centroid_to_ppi_centroid_distance", ""),
                "atp_overlap_fraction": gate.get("atp_overlap_fraction", ""),
                "atp_centroid_distance": gate.get("atp_centroid_distance", ""),
                "lower_gate_class": gate.get("lower_gate_class", ""),
                "z_percentile_within_receptor_surface": gate.get("z_percentile_within_receptor_surface", ""),
                "lateral_score": gate.get("lateral_score", ""),
                "interface_overlap_fraction": gate.get("interface_overlap_fraction", ""),
                "sasa_ratio_dimer_vs_isolated_protomer": gate.get("sasa_ratio_dimer_vs_isolated_protomer", ""),
                "hard_gate_pass": gate.get("hard_gate_pass", ""),
                "m3_primary_allowed": export.get("m3_primary_allowed", False),
                "m3_review_allowed": export.get("m3_review_allowed", False),
                "evidence_status": gate.get("evidence_status", "GATED_RAW_EVIDENCE"),
                "reject_or_review_reason": reason,
                "warnings": ";".join(_dedupe(_split_values(gate.get("warnings", "")) + _split_values(export.get("warnings", "")))),
            }
        )
    return rows


def _build_transition_gates(
    ppi_path: Path | None,
    ppi_rows: list[dict[str, Any]],
    atp_path: Path | None,
    atp_rows: list[dict[str, Any]],
    families_path: Path | None,
    family_rows: list[dict[str, Any]],
    gate_qc_path: Path | None,
    gate_rows: list[dict[str, Any]],
    accepted_path: Path | None,
    accepted_rows: list[dict[str, Any]],
    export_rows: list[dict[str, Any]],
    production: bool,
) -> list[dict[str, Any]]:
    return [
        _transition_row("G1", "PPI pose QC complete", True, ppi_path, len(ppi_rows), production, allow_warn_empty=True),
        _transition_row("G2", "accepted/rejected PPI pose reasons recorded", True, ppi_path, len(ppi_rows), production, allow_warn_empty=True),
        _transition_row("G3", "PPI consensus patch generated", True, ppi_path, len(ppi_rows), production, allow_warn_empty=True),
        _transition_row("G4", "ATP reference generated", True, atp_path, len(atp_rows), production, allow_warn_empty=True),
        _transition_row("G5", "raw pocket discovery complete", True, families_path, len(family_rows), production, allow_warn_empty=True),
        _transition_row("G6", "pocket hard gates applied", True, gate_qc_path, len(gate_rows), production, allow_warn_empty=True),
        _transition_row("G7", "accepted pocket family table generated", True, accepted_path, len(accepted_rows), production, allow_warn_empty=True),
        _m3_docking_transition_row(export_rows),
    ]


def _transition_row(
    gate_id: str,
    gate_name: str,
    required: bool,
    path: Path | None,
    count: int,
    production: bool,
    allow_warn_empty: bool = False,
) -> dict[str, Any]:
    if path is None or not path.is_file():
        status = "FAIL" if production else "PASS_WITH_WARNINGS"
        failure = "missing_required_evidence_file"
        warnings = "" if production else "synthetic_or_smoke_missing_evidence"
        evidence = ""
    elif count <= 0 and allow_warn_empty:
        status = "PASS_WITH_WARNINGS"
        failure = ""
        warnings = "evidence_file_present_but_empty"
        evidence = str(path)
    else:
        status = "PASS"
        failure = ""
        warnings = ""
        evidence = str(path)
    return {
        "gate_id": gate_id,
        "gate_name": gate_name,
        "required_for_m3": required,
        "status": status,
        "evidence_file": evidence,
        "evidence_count": count,
        "failure_reason": failure,
        "warnings": warnings,
    }


def _m3_docking_transition_row(export_rows: list[dict[str, Any]]) -> dict[str, Any]:
    accepted_count = sum(1 for row in export_rows if row.get("export_tier") in {"primary", "secondary_review"})
    if accepted_count:
        status = "PASS"
        warnings = ""
        failure = ""
    else:
        status = "PASS_WITH_WARNINGS"
        warnings = "no_accepted_primary_or_secondary_pockets_m3_docking_blocked_pending_manual_review"
        failure = ""
    return {
        "gate_id": "G8",
        "gate_name": "if no accepted pockets, M3 docking is blocked pending manual review",
        "required_for_m3": True,
        "status": status,
        "evidence_file": "phase2_pockets/export_for_m3/accepted_pockets_for_m3.csv",
        "evidence_count": accepted_count,
        "failure_reason": failure,
        "warnings": warnings,
    }


def _transition_status(rows: list[dict[str, Any]]) -> str:
    if any(row.get("status") == "FAIL" for row in rows):
        return "FAIL"
    if any(row.get("status") == "PASS_WITH_WARNINGS" for row in rows):
        return "PASS_WITH_WARNINGS"
    return "PASS"


def _cleanup_payload(ctx: RunContext) -> dict[str, Any]:
    preserved = [
        "manifest/",
        "logs/",
        "phase1_ppi/",
        "phase2_pockets/atp_reference/",
        "phase2_pockets/merged/",
        "phase2_pockets/gated/",
        "phase2_pockets/tables/",
        "phase2_pockets/final/",
        "phase2_pockets/export_for_m3/",
        "reports/",
        "qc/",
    ]
    return {
        "schema_version": "m2_8_cleanup_report_v1",
        "run_id": ctx.run_id,
        "cleanup_mode": "report_only",
        "deleted_files": [],
        "preserved_files": preserved,
        "warnings": ["production_m2_8_deletes_nothing"],
    }


def _source_files(ctx: RunContext, **paths: Path | None) -> dict[str, str | None]:
    return {key: (_repo_path(ctx, path) if path is not None else None) for key, path in paths.items()}


def _overall_status(blockers: list[str], warnings: list[str], transition_status: str, m3_docking_allowed: bool) -> str:
    if blockers or transition_status == "FAIL":
        return "FAIL"
    if warnings or transition_status == "PASS_WITH_WARNINGS" or not m3_docking_allowed:
        return "PASS_WITH_WARNINGS"
    return "PASS"


def _write_report(
    ctx: RunContext,
    report: M2ExportReport,
    source_files: dict[str, str | None],
    output_files: dict[str, str],
    counts: dict[str, int],
    transition_rows: list[dict[str, Any]],
    export_rows: list[dict[str, Any]],
    gate_rows: list[dict[str, Any]],
) -> None:
    lines = [
        "# Milestone 2 Summary",
        "",
        "## 1. Run identity and input provenance",
        "- run_id: {0}".format(ctx.run_id),
        "- status: {0}".format(report.status),
        "- transition_gate_status: {0}".format(report.transition_gate_status),
        "- m3_docking_allowed: {0}".format(report.m3_docking_allowed),
        "",
        "## 2. Receptor state summary",
        "- Primary states: EGFR_160-185, EGFR_170-200",
        "- Reference/control state: 3GT8_raw",
        "- 3GT8_raw-only pockets were not promoted to primary.",
        "",
        "## 3. PPI consensus patch summary",
        "- Source: {0}".format(source_files.get("ppi_path") or "missing"),
        "",
        "## 4. ATP reference summary",
        "- Source: {0}".format(source_files.get("atp_path") or "missing"),
        "- ATP-rejected pockets were not promoted.",
        "",
        "## 5. Pocket discovery summary",
        "- M2.6 merged pocket families: {0}".format(counts["total_pocket_families"]),
        "",
        "## 6. Pocket merge/family summary",
        "- Merged family source: {0}".format(source_files.get("families_path") or "missing"),
        "",
        "## 7. M2.7 hard gate summary",
        "- Gate QC source: {0}".format(source_files.get("gate_qc_path") or "missing"),
        "- Soft scores did not rescue hard-gate failures.",
        "",
        "## 8. Accepted primary pocket families",
    ]
    primary = [row for row in export_rows if row.get("export_tier") == "primary"]
    secondary = [row for row in export_rows if row.get("export_tier") == "secondary_review"]
    reference = [row for row in export_rows if row.get("export_tier") == "reference_only"]
    rejected = [row for row in gate_rows if row.get("gate_class") not in {"accepted_primary_pocket", "accepted_secondary_cryptic", "reference_only"}]
    lines.extend(_family_lines(primary))
    lines.extend(["", "## 9. Accepted secondary/review pocket families"])
    lines.extend(_family_lines(secondary))
    lines.extend(["", "## 10. Reference-only 3GT8_raw pockets"])
    lines.extend(_family_lines(reference))
    lines.extend(["", "## 11. Rejected pocket families and reasons"])
    if rejected:
        for row in rejected:
            lines.append("- {0}: {1}".format(row.get("pocket_family_id", ""), row.get("reject_reason") or row.get("gate_class", "")))
    else:
        lines.append("- none")
    lines.extend(["", "## 12. M2 -> M3 transition gate status"])
    for row in transition_rows:
        lines.append("- {0} {1}: {2}".format(row.get("gate_id"), row.get("gate_name"), row.get("status")))
    lines.extend(["", "## 13. Export package files"])
    for key in sorted(output_files):
        lines.append("- {0}: {1}".format(key, output_files[key]))
    lines.extend(["", "## 14. Limitations and manual review items"])
    if report.warnings:
        for warning in report.warnings:
            lines.append("- {0}".format(warning))
    else:
        lines.append("- none")
    if not report.m3_docking_allowed:
        lines.append("- No accepted primary pocket family is available for unrestricted Milestone 3 docking.")
        lines.append("- Milestone 3 docking should not proceed without manual review or gate-threshold review.")
    lines.extend(["", "## 15. Explicit non-goals preserved"])
    for item in NON_GOALS_CONFIRMED:
        lines.append("- {0}".format(item))
    lines.extend(
        [
            "- No compound docking was run.",
            "- No final compound candidate claims are made.",
            "",
            "## 16. Next recommended step",
            "- Milestone 3 / M3 Codex Task 0 - readiness audit and M3 directory skeleton",
        ]
    )
    text = "\n".join(lines) + "\n"
    for path in [report.final_summary_md, report.export_report_md, report.reports_summary_md]:
        if path is None:
            continue
        safe = ctx.require_within_run_dir(path)
        safe.parent.mkdir(parents=True, exist_ok=True)
        safe.write_text(text, encoding="utf-8")


def _family_lines(rows: list[dict[str, Any]]) -> list[str]:
    if not rows:
        return ["- none"]
    return [
        "- {0}: tier={1}, state={2}, protomer={3}".format(
            row.get("pocket_family_id", ""),
            row.get("export_tier", ""),
            row.get("representative_state_id", ""),
            row.get("representative_protomer_id", ""),
        )
        for row in rows
    ]


def _state_receptor_path(ctx: RunContext, state_id: str) -> Path | None:
    candidates = [
        ctx.run_dir / "normalized" / "receptors" / "{0}_dockable_669_1014_explicit_AB.pdb".format(state_id),
        ctx.run_dir / "normalized" / "receptors" / "{0}_dockable_reference_explicit_AB.pdb".format(state_id),
        ctx.run_dir / "normalized" / "receptors" / "{0}_full_frame_explicit_AB.pdb".format(state_id),
        ctx.run_dir / "normalized" / "receptors" / "{0}_runtime_offset_receptor_only.pdb".format(state_id),
    ]
    return next((ctx.require_within_run_dir(path) for path in candidates if path.is_file()), None)


def _mapping_path(ctx: RunContext, state_id: str) -> Path | None:
    candidates = [
        ctx.qc_dir / "{0}_receptor_mapping.csv".format(state_id),
        ctx.run_dir / "normalized" / "receptors" / "{0}_receptor_mapping.csv".format(state_id),
    ]
    return next((ctx.require_within_run_dir(path) for path in candidates if path.is_file()), None)


def _count_export(rows: list[dict[str, Any]], tier: str) -> int:
    return sum(1 for row in rows if row.get("export_tier") == tier)


def _safe_float(value: Any) -> float | None:
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _format_float(value: float | None) -> str:
    if value is None:
        return ""
    return "{0:.3f}".format(value)


def _int(value: Any) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return 0


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "pass"}


def _split_values(value: Any) -> list[str]:
    text = str(value or "").strip()
    if not text:
        return []
    values: list[str] = []
    for item in text.replace(",", ";").split(";"):
        item = item.strip()
        if item:
            values.append(item)
    return values


def _dedupe(values: Any) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        item = str(value).strip()
        if item and item not in seen:
            seen.add(item)
            result.append(item)
    return result
