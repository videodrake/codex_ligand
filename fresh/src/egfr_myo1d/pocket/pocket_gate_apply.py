"""M2.7 hard-gate application for raw pocket families."""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import load_yaml_config, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError, ensure_within
from egfr_myo1d.io.residue_mapping import read_residue_mapping
from egfr_myo1d.pocket.dimer_accessibility import (
    DIMER_ACCESSIBILITY_FIELDS,
    evaluate_dimer_accessibility,
)
from egfr_myo1d.pocket.membrane_gate import (
    MEMBRANE_GEOMETRY_FIELDS,
    evaluate_membrane_geometry,
)
from egfr_myo1d.pocket.ppi_relationship import (
    PPI_RELATIONSHIP_FIELDS,
    classify_ppi_relationship,
    distance,
    format_float,
    point,
    residue_set,
    safe_float,
    split_values,
)
from egfr_myo1d.structure.pdb_parser import AtomRecord, PDBParseError, parse_pdb


M2_7_PHASE = "m2.7-ppi-membrane-dimer-pocket-gates"
GATED_ROOT = "phase2_pockets/gated"
TABLE_ROOT = "phase2_pockets/tables"
GATE_QC_FIELDS = [
    "pocket_family_id",
    "representative_candidate_pocket_id",
    "member_candidate_pocket_ids",
    "member_state_ids",
    "member_state_roles",
    "member_protomer_ids",
    "primary_state_count",
    "reference_state_count",
    "raw_candidate_count",
    "dimer_origin_pass",
    "non_atp_pass",
    "ppi_relationship_pass",
    "lower_lateral_pass",
    "dimer_accessibility_pass",
    "primary_state_support_or_review_flag",
    "atp_overlap_fraction",
    "atp_centroid_distance",
    "ppi_relationship_class",
    "residue_overlap_count",
    "minimum_heavy_atom_distance_to_ppi",
    "pocket_centroid_to_ppi_centroid_distance",
    "lower_gate_class",
    "z_percentile_within_receptor_surface",
    "lateral_score",
    "interface_overlap_fraction",
    "sasa_ratio_dimer_vs_isolated_protomer",
    "gate_class",
    "hard_gate_pass",
    "reject_reason",
    "review_reason",
    "soft_priority_score",
    "soft_score_status",
    "evidence_status",
    "warnings",
]
ACCEPTED_FIELDS = [
    "pocket_family_id",
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
    "family_center_x",
    "family_center_y",
    "family_center_z",
    "residue_union_uniprot",
    "ppi_relationship_class",
    "non_atp_pass",
    "ppi_relationship_pass",
    "lower_lateral_pass",
    "dimer_accessibility_pass",
    "primary_state_support_or_review_flag",
    "acceptance_reason",
    "review_notes",
    "warnings",
]
NON_GOALS_CONFIRMED = [
    "no_pyrosetta_import",
    "no_pyrosetta_docking_or_relaxation",
    "no_lightdock_execution",
    "no_p2rank_execution",
    "no_vina_execution",
    "no_ligand_preparation",
    "no_compound_docking_or_scoring",
    "no_candidate_nomination",
    "no_m2_8_final_export_package",
    "no_milestone3_outputs",
    "no_scheduler_submission",
    "no_cleanup_deletion",
]
FORBIDDEN_EXPORT_OUTPUTS = [
    "phase2_pockets/export_for_m3/accepted_pocket_boxes.csv",
    "phase2_pockets/export_for_m3/accepted_pocket_residue_sets.json",
    "phase2_pockets/export_for_m3/accepted_pocket_receptor_state_map.csv",
    "phase2_pockets/export_for_m3/m2_final_report.md",
    "phase2_pockets/export_for_m3/m2_cleanup_report.json",
    "phase3_compounds",
]


@dataclass
class PocketGateReport:
    run_id: str
    mode: str
    profile: str
    status: str
    family_count: int = 0
    accepted_primary_count: int = 0
    accepted_secondary_count: int = 0
    manual_review_count: int = 0
    reference_only_count: int = 0
    atp_reject_count: int = 0
    ppi_reject_count: int = 0
    membrane_reject_count: int = 0
    dimer_reject_count: int = 0
    mapping_or_origin_reject_count: int = 0
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    gate_qc_csv: Path | None = None
    accepted_csv: Path | None = None
    rejected_csv: Path | None = None
    ppi_relationship_csv: Path | None = None
    membrane_geometry_csv: Path | None = None
    dimer_accessibility_csv: Path | None = None
    status_json: Path | None = None
    manifest_path: Path | None = None
    summary_report_path: Path | None = None


def apply_pocket_gates(
    ctx: RunContext,
    mode: str = "smoke_input",
    profile: str = "codex_dev",
    config_path: Path | None = None,
    ppi_consensus: Path | None = None,
    atp_reference: Path | None = None,
    atp_centroids: Path | None = None,
    pocket_families: Path | None = None,
    raw_candidates: Path | None = None,
    membrane_frame: Path | None = None,
    synthetic_fixture: bool = False,
) -> PocketGateReport:
    """Apply M2.7 hard gates to M2.6 raw pocket families."""

    ctx.create_directories()
    production = mode == "production" and not synthetic_fixture
    config = _load_config(ctx, config_path)
    thresholds = _thresholds(config)
    ppi_path = _resolve_optional_path(ctx, ppi_consensus) or _discover_ppi_consensus(ctx)
    atp_path = _resolve_optional_path(ctx, atp_reference) or _discover_atp_reference(ctx)
    centroid_path = _resolve_optional_path(ctx, atp_centroids) or _discover_atp_centroids(ctx)
    families_path = _resolve_optional_path(ctx, pocket_families) or _discover_pocket_families(ctx)
    raw_path = _resolve_optional_path(ctx, raw_candidates) or _discover_raw_candidates(ctx)
    membrane_path = _resolve_optional_path(ctx, membrane_frame) or _discover_membrane_frame(ctx)

    warnings: list[str] = []
    blockers: list[str] = []
    if ppi_path is None or not ppi_path.is_file():
        _missing("missing_m2_4_ppi_consensus_patch", warnings, blockers, production)
    if atp_path is None or not atp_path.is_file():
        _missing("missing_m2_5_atp_reference", warnings, blockers, production)
    if centroid_path is None or not centroid_path.is_file():
        _missing("missing_m2_5_atp_centroid_by_state", warnings, blockers, production)
    if families_path is None or not families_path.is_file():
        _missing("missing_m2_6_pocket_candidates_merged", warnings, blockers, production)
    if membrane_path is None or not membrane_path.is_file():
        _missing("missing_m1_membrane_frame", warnings, blockers, production)

    ppi_rows = _read_csv(ppi_path) if ppi_path and ppi_path.is_file() else []
    atp_rows = _read_csv(atp_path) if atp_path and atp_path.is_file() else []
    atp_centroid_rows = _read_csv(centroid_path) if centroid_path and centroid_path.is_file() else []
    families = _read_csv(families_path) if families_path and families_path.is_file() else []
    raw_rows = _read_csv(raw_path) if raw_path and raw_path.is_file() else []
    membrane_payload = _read_json(membrane_path) if membrane_path and membrane_path.is_file() else None

    states = _states_from_families(families)
    residue_geometry, receptor_z_ranges, interface_residues, geometry_warnings = _build_receptor_geometry(ctx, states, thresholds)
    warnings.extend(geometry_warnings)
    if production and not residue_geometry:
        blockers.append("missing_receptor_mapping_or_receptor_pdb_blocks_production_gating")
    raw_by_id = {str(row.get("candidate_pocket_id", "")): row for row in raw_rows}
    atp_residues = {
        str(row.get("uniprot_residue_number", "")).strip()
        for row in atp_rows
        if str(row.get("uniprot_residue_number", "")).strip()
    }

    ppi_relationship_rows: list[dict[str, Any]] = []
    membrane_rows: list[dict[str, Any]] = []
    dimer_rows: list[dict[str, Any]] = []
    gate_rows: list[dict[str, Any]] = []
    accepted_rows: list[dict[str, Any]] = []
    rejected_rows: list[dict[str, Any]] = []
    for family in families:
        mapping_status, mapping_warnings = _family_mapping_status(family, raw_by_id)
        warnings.extend(mapping_warnings)
        atp_metrics = _evaluate_atp_gate(family, atp_residues, atp_centroid_rows, thresholds)
        membrane_eval = evaluate_membrane_geometry(
            family,
            membrane_payload,
            receptor_z_ranges,
            thresholds,
            production=production,
        )
        ppi_eval = classify_ppi_relationship(
            family,
            ppi_rows,
            residue_geometry,
            thresholds,
            same_lower_lateral_neighborhood_flag=_as_bool(membrane_eval.get("lower_lateral_pass")),
            mapping_status=mapping_status,
        )
        dimer_eval = evaluate_dimer_accessibility(
            family,
            interface_residues,
            threshold=thresholds["interface_overlap_fraction_reject_threshold"],
        )
        dimer_origin_pass, dimer_origin_notes = _dimer_origin(family)
        gate_class, reject_reason, review_reason, support_flag = _gate_class(
            family,
            dimer_origin_pass,
            mapping_status,
            atp_metrics["non_atp_pass"],
            _as_bool(ppi_eval["ppi_relationship_pass"]),
            _as_bool(membrane_eval["lower_lateral_pass"]),
            _as_bool(dimer_eval["dimer_accessibility_pass"]),
            production_missing=bool(blockers),
            thresholds=thresholds,
        )
        hard_gate_pass = gate_class in {"accepted_primary_pocket", "accepted_secondary_cryptic"}
        soft_score, soft_status = _soft_score(
            hard_gate_pass,
            family,
            ppi_eval,
            atp_metrics,
            membrane_eval,
            dimer_eval,
        )
        row_warnings = sorted(
            set(
                split_values(family.get("warnings", ""))
                + mapping_warnings
                + split_values(ppi_eval.get("warnings", ""))
                + split_values(membrane_eval.get("warnings", ""))
                + split_values(dimer_eval.get("warnings", ""))
                + split_values(atp_metrics.get("warnings", ""))
                + ([dimer_origin_notes] if dimer_origin_notes else [])
            )
        )
        gate_row = {
            "pocket_family_id": family.get("pocket_family_id", ""),
            "representative_candidate_pocket_id": family.get("representative_candidate_pocket_id", ""),
            "member_candidate_pocket_ids": family.get("member_candidate_pocket_ids", ""),
            "member_state_ids": family.get("member_state_ids", ""),
            "member_state_roles": family.get("member_state_roles", ""),
            "member_protomer_ids": family.get("member_protomer_ids", ""),
            "primary_state_count": family.get("primary_state_count", ""),
            "reference_state_count": family.get("reference_state_count", ""),
            "raw_candidate_count": family.get("raw_candidate_count", ""),
            "dimer_origin_pass": dimer_origin_pass,
            "non_atp_pass": atp_metrics["non_atp_pass"],
            "ppi_relationship_pass": ppi_eval["ppi_relationship_pass"],
            "lower_lateral_pass": membrane_eval["lower_lateral_pass"],
            "dimer_accessibility_pass": dimer_eval["dimer_accessibility_pass"],
            "primary_state_support_or_review_flag": support_flag,
            "atp_overlap_fraction": atp_metrics["atp_overlap_fraction"],
            "atp_centroid_distance": atp_metrics["atp_centroid_distance"],
            "ppi_relationship_class": ppi_eval["ppi_relationship_class"],
            "residue_overlap_count": ppi_eval["residue_overlap_count"],
            "minimum_heavy_atom_distance_to_ppi": ppi_eval["minimum_heavy_atom_distance_to_ppi"],
            "pocket_centroid_to_ppi_centroid_distance": ppi_eval["pocket_centroid_to_ppi_centroid_distance"],
            "lower_gate_class": membrane_eval["lower_gate_class"],
            "z_percentile_within_receptor_surface": membrane_eval["z_percentile_within_receptor_surface"],
            "lateral_score": membrane_eval["lateral_score"],
            "interface_overlap_fraction": dimer_eval["interface_overlap_fraction"],
            "sasa_ratio_dimer_vs_isolated_protomer": dimer_eval["sasa_ratio_dimer_vs_isolated_protomer"],
            "gate_class": gate_class,
            "hard_gate_pass": hard_gate_pass,
            "reject_reason": reject_reason,
            "review_reason": review_reason,
            "soft_priority_score": soft_score,
            "soft_score_status": soft_status,
            "evidence_status": "GATED_RAW_EVIDENCE",
            "warnings": ";".join(row_warnings),
        }
        ppi_relationship_rows.append(ppi_eval)
        membrane_rows.append(membrane_eval)
        dimer_rows.append(dimer_eval)
        gate_rows.append(gate_row)
        if gate_class in {"accepted_primary_pocket", "accepted_secondary_cryptic"}:
            accepted_rows.append(_accepted_row(family, gate_row))
        else:
            rejected_rows.append(_accepted_row(family, gate_row))

    export_warnings = _forbidden_export_warnings(ctx)
    warnings.extend(export_warnings)
    status = _status(blockers, warnings, gate_rows)
    gated_dir = ctx.require_within_run_dir(ctx.run_dir / GATED_ROOT)
    tables_dir = ctx.require_within_run_dir(ctx.run_dir / TABLE_ROOT)
    gated_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    gate_qc = gated_dir / "pocket_gate_qc.csv"
    accepted = gated_dir / "accepted_pocket_families.csv"
    rejected = gated_dir / "rejected_pocket_families.csv"
    ppi_csv = gated_dir / "pocket_ppi_relationship.csv"
    membrane_csv = gated_dir / "pocket_membrane_geometry.csv"
    dimer_csv = gated_dir / "pocket_dimer_accessibility.csv"
    status_json = gated_dir / "pocket_gate_status.json"
    _write_csv(gate_qc, gate_rows, GATE_QC_FIELDS, ctx)
    _write_csv(accepted, accepted_rows, ACCEPTED_FIELDS, ctx)
    _write_csv(rejected, rejected_rows, ACCEPTED_FIELDS, ctx)
    _write_csv(ppi_csv, ppi_relationship_rows, PPI_RELATIONSHIP_FIELDS, ctx)
    _write_csv(membrane_csv, membrane_rows, MEMBRANE_GEOMETRY_FIELDS, ctx)
    _write_csv(dimer_csv, dimer_rows, DIMER_ACCESSIBILITY_FIELDS, ctx)

    _write_csv(tables_dir / "pocket_gate_qc.csv", gate_rows, GATE_QC_FIELDS, ctx)
    _write_csv(tables_dir / "accepted_pockets_for_m3.csv", accepted_rows, ACCEPTED_FIELDS, ctx)
    _write_csv(tables_dir / "pocket_rejects.csv", rejected_rows, ACCEPTED_FIELDS, ctx)
    _write_csv(tables_dir / "pocket_ppi_relationship.csv", ppi_relationship_rows, PPI_RELATIONSHIP_FIELDS, ctx)
    _write_csv(tables_dir / "pocket_membrane_geometry.csv", membrane_rows, MEMBRANE_GEOMETRY_FIELDS, ctx)
    _write_csv(tables_dir / "pocket_dimer_accessibility.csv", dimer_rows, DIMER_ACCESSIBILITY_FIELDS, ctx)

    report = PocketGateReport(
        run_id=ctx.run_id,
        mode=mode,
        profile=profile,
        status=status,
        family_count=len(gate_rows),
        accepted_primary_count=_count_gate(gate_rows, "accepted_primary_pocket"),
        accepted_secondary_count=_count_gate(gate_rows, "accepted_secondary_cryptic"),
        manual_review_count=_count_gate(gate_rows, "manual_review"),
        reference_only_count=_count_gate(gate_rows, "reference_only"),
        atp_reject_count=_count_gate(gate_rows, "atp_reject"),
        ppi_reject_count=_count_gate(gate_rows, "ppi_low_relevance_reject"),
        membrane_reject_count=_count_gate(gate_rows, "membrane_geometry_reject"),
        dimer_reject_count=_count_gate(gate_rows, "dimer_buried_reject"),
        mapping_or_origin_reject_count=sum(
            1 for row in gate_rows if row.get("gate_class") in {"mapping_reject", "dimer_origin_reject", "not_evaluable"}
        ),
        warnings=sorted(set(warnings)),
        blockers=blockers,
        gate_qc_csv=gate_qc,
        accepted_csv=accepted,
        rejected_csv=rejected,
        ppi_relationship_csv=ppi_csv,
        membrane_geometry_csv=membrane_csv,
        dimer_accessibility_csv=dimer_csv,
    )
    summary = _write_summary(ctx, report)
    report.summary_report_path = summary
    payload = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": M2_7_PHASE,
        "mode": mode,
        "profile": profile,
        "synthetic_fixture": synthetic_fixture,
        "status": status,
        "family_count": len(gate_rows),
        "gate_summary": {
            "accepted_primary_pocket": report.accepted_primary_count,
            "accepted_secondary_cryptic": report.accepted_secondary_count,
            "manual_review": report.manual_review_count,
            "reference_only": report.reference_only_count,
            "atp_reject": report.atp_reject_count,
            "ppi_low_relevance_reject": report.ppi_reject_count,
            "membrane_geometry_reject": report.membrane_reject_count,
            "dimer_buried_reject": report.dimer_reject_count,
            "mapping_or_origin_reject": report.mapping_or_origin_reject_count,
        },
        "inputs": {
            "ppi_consensus": _repo_path(ctx, ppi_path) if ppi_path else None,
            "atp_reference": _repo_path(ctx, atp_path) if atp_path else None,
            "atp_centroids": _repo_path(ctx, centroid_path) if centroid_path else None,
            "pocket_families": _repo_path(ctx, families_path) if families_path else None,
            "raw_candidates": _repo_path(ctx, raw_path) if raw_path else None,
            "membrane_frame": _repo_path(ctx, membrane_path) if membrane_path else None,
        },
        "outputs": {
            "gated_pocket_gate_qc": _repo_path(ctx, gate_qc),
            "gated_accepted_pocket_families": _repo_path(ctx, accepted),
            "gated_rejected_pocket_families": _repo_path(ctx, rejected),
            "gated_ppi_relationship": _repo_path(ctx, ppi_csv),
            "gated_membrane_geometry": _repo_path(ctx, membrane_csv),
            "gated_dimer_accessibility": _repo_path(ctx, dimer_csv),
            "tables_pocket_gate_qc": _repo_path(ctx, tables_dir / "pocket_gate_qc.csv"),
            "tables_accepted_pockets_for_m3": _repo_path(ctx, tables_dir / "accepted_pockets_for_m3.csv"),
            "tables_pocket_rejects": _repo_path(ctx, tables_dir / "pocket_rejects.csv"),
            "status_json": _repo_path(ctx, status_json),
            "summary_report": _repo_path(ctx, summary),
        },
        "thresholds": thresholds,
        "vina_executed": False,
        "ligand_preparation_executed": False,
        "pyrosetta_docking_or_relaxation_executed": False,
        "lightdock_executed": False,
        "p2rank_executed": False,
        "compound_scoring_executed": False,
        "candidate_nomination_executed": False,
        "m2_8_export_package_created": False,
        "milestone3_started": False,
        "warnings": sorted(set(warnings)),
        "blockers": blockers,
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
    }
    write_json(status_json, payload, ctx)
    write_json(tables_dir / "pocket_gate_status.json", payload, ctx)
    manifest = ctx.manifest_dir / "m2_7_pocket_gate_manifest.json"
    write_json(manifest, payload, ctx)
    report.status_json = status_json
    report.manifest_path = manifest

    phase_status = "FAIL" if status == "FAIL" else ("WARN" if status == "PASS_WITH_WARNINGS" else "PASS")
    append_phase_status(
        ctx,
        M2_7_PHASE,
        phase_status,
        "M2.7 pocket gates {0} for {1} family/families".format(status, len(gate_rows)),
        {"family_count": len(gate_rows), "m2_8_export_package_created": False},
    )
    append_job_status(
        ctx,
        "m2_7_pocket_gates_local",
        phase_status,
        details={"message": "local gate classification only"},
    )
    return report


def _load_config(ctx: RunContext, config_path: Path | None) -> dict[str, Any]:
    path = _resolve_optional_path(ctx, config_path) or (ctx.repo_root / "fresh" / "configs" / "pocket.yaml")
    if path.is_file():
        return load_yaml_config(path)
    return {}


def _thresholds(config: dict[str, Any]) -> dict[str, float]:
    gates = config.get("gates", {}) or config.get("pocket_gates", {}) or {}
    return {
        "atp_overlap_fraction_reject_threshold": float(gates.get("atp_overlap_fraction_reject_threshold", gates.get("atp_overlap_fraction_cutoff", 0.15))),
        "atp_centroid_distance_reject_angstrom": float(gates.get("atp_centroid_distance_reject_angstrom", gates.get("atp_centroid_distance_cutoff_angstrom", 10.0))),
        "ppi_rim_heavy_atom_cutoff_angstrom": float(gates.get("ppi_rim_heavy_atom_cutoff_angstrom", gates.get("ppi_rim_min_distance_cutoff_angstrom", 6.0))),
        "ppi_rim_ca_cutoff_angstrom": float(gates.get("ppi_rim_ca_cutoff_angstrom", 10.0)),
        "ppi_allosteric_heavy_atom_cutoff_angstrom": float(gates.get("ppi_allosteric_heavy_atom_cutoff_angstrom", 12.0)),
        "ppi_allosteric_ca_cutoff_angstrom": float(gates.get("ppi_allosteric_ca_cutoff_angstrom", 18.0)),
        "ppi_allosteric_centroid_cutoff_angstrom": float(gates.get("ppi_allosteric_near_centroid_cutoff_angstrom", 18.0)),
        "lateral_score_pass_threshold": float(gates.get("lateral_score_pass_threshold", gates.get("lateral_score_cutoff", 0.15))),
        "lower_z_minimum": float(gates.get("lower_z_minimum", -5.0)),
        "interface_heavy_atom_cutoff_angstrom": float(gates.get("interface_heavy_atom_cutoff_angstrom", 5.0)),
        "interface_overlap_fraction_reject_threshold": float(gates.get("interface_overlap_fraction_reject_threshold", gates.get("interface_overlap_fraction_cutoff", 0.25))),
        "required_primary_state_count_for_primary": float(gates.get("required_primary_state_count_for_primary", 2)),
        "allow_secondary_if_primary_state_count": float(gates.get("allow_secondary_if_primary_state_count", 1)),
    }


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
        ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv",
        ctx.run_dir / "phase1_ppi" / "tables" / "ppi_consensus_patch.csv",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _discover_atp_reference(ctx: RunContext) -> Path | None:
    candidates = [
        ctx.run_dir / "phase2_pockets" / "atp_reference" / "atp_site_reference.csv",
        ctx.run_dir / "phase2_pockets" / "tables" / "atp_site_reference.csv",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _discover_atp_centroids(ctx: RunContext) -> Path | None:
    path = ctx.run_dir / "phase2_pockets" / "atp_reference" / "atp_site_centroid_by_state.csv"
    return path if path.is_file() else None


def _discover_pocket_families(ctx: RunContext) -> Path | None:
    candidates = [
        ctx.run_dir / "phase2_pockets" / "merged" / "pocket_candidates_merged.csv",
        ctx.run_dir / "phase2_pockets" / "tables" / "pocket_candidates_merged.csv",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _discover_raw_candidates(ctx: RunContext) -> Path | None:
    path = ctx.run_dir / "phase2_pockets" / "merged" / "pocket_candidates_raw.csv"
    return path if path.is_file() else None


def _discover_membrane_frame(ctx: RunContext) -> Path | None:
    candidates = [
        ctx.run_dir / "manifest" / "membrane_frame.json",
        ctx.run_dir / "normalized" / "geometry" / "membrane_frame.json",
        ctx.repo_root / "fresh" / "data" / "normalized" / "geometry" / "membrane_frame.json",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _missing(message: str, warnings: list[str], blockers: list[str], production: bool) -> None:
    if production:
        blockers.append("{0}_blocks_production_gating".format(message))
    else:
        warnings.append("{0}_synthetic_or_smoke_not_production_evidence".format(message))


def _read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _states_from_families(families: list[dict[str, Any]]) -> list[str]:
    states: list[str] = []
    for family in families:
        for state in split_values(family.get("member_state_ids", "")):
            if state not in states:
                states.append(state)
    return states


def _state_receptor_path(ctx: RunContext, state_id: str) -> Path | None:
    candidates = [
        ctx.run_dir / "normalized" / "receptors" / "{0}_dockable_669_1014_explicit_AB.pdb".format(state_id),
        ctx.run_dir / "normalized" / "receptors" / "{0}_dockable_reference_explicit_AB.pdb".format(state_id),
        ctx.run_dir / "normalized" / "receptors" / "{0}_full_frame_explicit_AB.pdb".format(state_id),
        ctx.run_dir / "normalized" / "receptors" / "{0}_runtime_offset_receptor_only.pdb".format(state_id),
    ]
    return next((ctx.require_within_run_dir(path) for path in candidates if path.is_file()), None)


def _mapping_path(ctx: RunContext, state_id: str) -> Path | None:
    path = ctx.qc_dir / "{0}_receptor_mapping.csv".format(state_id)
    return ctx.require_within_run_dir(path) if path.is_file() else None


def _build_receptor_geometry(
    ctx: RunContext,
    states: list[str],
    thresholds: dict[str, float],
) -> tuple[dict[tuple[str, str, str], dict[str, Any]], dict[str, tuple[float, float]], set[str], list[str]]:
    geometry: dict[tuple[str, str, str], dict[str, Any]] = {}
    z_ranges: dict[str, tuple[float, float]] = {}
    interface_residues: set[str] = set()
    warnings: list[str] = []
    for state_id in states:
        receptor = _state_receptor_path(ctx, state_id)
        mapping = _mapping_path(ctx, state_id)
        if receptor is None:
            warnings.append("missing_receptor_pdb_for_gate_geometry:{0}".format(state_id))
            continue
        if mapping is None:
            warnings.append("missing_receptor_mapping_for_gate_geometry:{0}".format(state_id))
            continue
        try:
            structure = parse_pdb(receptor)
        except PDBParseError as exc:
            warnings.append("receptor_parse_failed_for_gate_geometry:{0}:{1}".format(state_id, exc))
            continue
        mapping_rows = read_residue_mapping(mapping)
        by_source = {
            (str(row.get("source_chain", "")), _safe_int(row.get("source_resseq"))): row
            for row in mapping_rows
            if _safe_int(row.get("source_resseq")) is not None
        }
        residue_atoms: dict[tuple[str, str, str], list[AtomRecord]] = {}
        projected_z = [atom.z for atom in structure.atoms if _is_heavy(atom)]
        if projected_z:
            z_ranges[state_id] = (min(projected_z), max(projected_z))
        for atom in structure.atoms:
            if not _is_heavy(atom):
                continue
            row = by_source.get((atom.chain_id, atom.residue_number))
            if row is None:
                continue
            protomer_id = str(row.get("protomer_id", ""))
            residue_number = str(row.get("source_resseq", ""))
            key = (state_id, protomer_id, residue_number)
            residue_atoms.setdefault(key, []).append(atom)
        for key, atoms in residue_atoms.items():
            heavy_points = [(atom.x, atom.y, atom.z) for atom in atoms]
            ca_points = [(atom.x, atom.y, atom.z) for atom in atoms if atom.atom_name == "CA"]
            geometry[key] = {"heavy": heavy_points, "ca": ca_points or heavy_points[:1]}
        interface_residues.update(_interface_residues(residue_atoms, thresholds["interface_heavy_atom_cutoff_angstrom"]))
    return geometry, z_ranges, interface_residues, warnings


def _interface_residues(residue_atoms: dict[tuple[str, str, str], list[AtomRecord]], cutoff: float) -> set[str]:
    residues: set[str] = set()
    items = sorted(residue_atoms.items())
    for i, (key_a, atoms_a) in enumerate(items):
        _state_a, protomer_a, residue_a = key_a
        for key_b, atoms_b in items[i + 1 :]:
            _state_b, protomer_b, residue_b = key_b
            if protomer_a == protomer_b:
                continue
            if _min_atom_distance(atoms_a, atoms_b) <= cutoff:
                residues.add(residue_a)
                residues.add(residue_b)
    return residues


def _min_atom_distance(left: list[AtomRecord], right: list[AtomRecord]) -> float:
    best = math.inf
    for atom_a in left:
        for atom_b in right:
            best = min(best, distance((atom_a.x, atom_a.y, atom_a.z), (atom_b.x, atom_b.y, atom_b.z)))
    return best


def _is_heavy(atom: AtomRecord) -> bool:
    element = (atom.element or atom.atom_name[:1]).upper()
    return element != "H"


def _safe_int(value: Any) -> int | None:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def _family_mapping_status(family: dict[str, Any], raw_by_id: dict[str, dict[str, Any]]) -> tuple[str, list[str]]:
    statuses = []
    warnings = []
    for candidate_id in split_values(family.get("member_candidate_pocket_ids", "")):
        raw = raw_by_id.get(candidate_id)
        if raw is None:
            warnings.append("missing_raw_candidate_for_mapping_status:{0}".format(candidate_id))
            continue
        status = str(raw.get("mapping_status", "")).upper()
        if status:
            statuses.append(status)
        warnings.extend(split_values(raw.get("warnings", "")))
    if not residue_set(family.get("residue_union_uniprot", "")):
        return "missing", warnings + ["family_residue_union_missing"]
    if any(status in {"MISSING", "NO_RESIDUES", "FAIL"} for status in statuses):
        return "missing", warnings
    if any(status in {"PARTIAL", "WARN"} for status in statuses):
        return "partial", warnings
    return "PASS", warnings


def _evaluate_atp_gate(
    family: dict[str, Any],
    atp_residues: set[str],
    atp_centroid_rows: list[dict[str, Any]],
    thresholds: dict[str, float],
) -> dict[str, Any]:
    warnings: list[str] = []
    pocket_residues = residue_set(family.get("residue_union_uniprot", ""))
    overlap = pocket_residues & atp_residues
    overlap_fraction = len(overlap) / float(len(pocket_residues)) if pocket_residues else 0.0
    centroid_distance = _nearest_atp_centroid_distance(family, atp_centroid_rows)
    if not atp_residues:
        warnings.append("atp_reference_residues_missing")
        non_atp_pass = False
        reason = "missing_atp_reference_blocks_non_atp_gate"
    elif overlap_fraction >= thresholds["atp_overlap_fraction_reject_threshold"]:
        non_atp_pass = False
        reason = "atp_overlap_fraction_exceeds_threshold"
    elif centroid_distance is not None and centroid_distance <= thresholds["atp_centroid_distance_reject_angstrom"]:
        non_atp_pass = False
        reason = "atp_centroid_distance_below_threshold"
    else:
        non_atp_pass = True
        reason = "non_atp_gate_pass"
    return {
        "atp_overlap_count": len(overlap),
        "atp_overlap_fraction": format_float(overlap_fraction),
        "atp_centroid_distance": format_float(centroid_distance),
        "atp_gate_reason": reason,
        "non_atp_pass": non_atp_pass,
        "warnings": ";".join(sorted(set(warnings))),
    }


def _nearest_atp_centroid_distance(family: dict[str, Any], centroid_rows: list[dict[str, Any]]) -> float | None:
    center = point(family, "family_center_x", "family_center_y", "family_center_z")
    if center is None:
        return None
    states = set(split_values(family.get("member_state_ids", "")))
    protomers = set(split_values(family.get("member_protomer_ids", "")))
    distances: list[float] = []
    for row in centroid_rows:
        state_id = str(row.get("state_id", ""))
        protomer_id = str(row.get("egfr_protomer_id", ""))
        if states and state_id and state_id not in states:
            continue
        if protomers and protomer_id and protomer_id not in protomers:
            continue
        atp_point = point(row, "atp_centroid_x", "atp_centroid_y", "atp_centroid_z")
        if atp_point:
            distances.append(distance(center, atp_point))
    return min(distances) if distances else None


def _dimer_origin(family: dict[str, Any]) -> tuple[bool, str]:
    primary_count = _int(family.get("primary_state_count"))
    reference_count = _int(family.get("reference_state_count"))
    protomers = split_values(family.get("member_protomer_ids", ""))
    if primary_count + reference_count <= 0:
        return False, "state_support_missing"
    if not protomers:
        return False, "protomer_identity_missing"
    return True, ""


def _gate_class(
    family: dict[str, Any],
    dimer_origin_pass: bool,
    mapping_status: str,
    non_atp_pass: bool,
    ppi_pass: bool,
    membrane_pass: bool,
    dimer_pass: bool,
    production_missing: bool,
    thresholds: dict[str, float],
) -> tuple[str, str, str, str]:
    primary_count = _int(family.get("primary_state_count"))
    reference_count = _int(family.get("reference_state_count"))
    if production_missing:
        return "not_evaluable", "missing_required_production_input", "", "not_evaluable"
    if mapping_status.lower() == "missing":
        return "mapping_reject", "mapping_failure_prevents_gate_acceptance", "", "mapping_reject"
    if not dimer_origin_pass:
        return "dimer_origin_reject", "dimer_origin_or_protomer_identity_failed", "", "dimer_origin_reject"
    if reference_count > 0 and primary_count == 0:
        return "reference_only", "", "3GT8_raw_or_reference_only_not_primary", "reference_only"
    if not non_atp_pass:
        return "atp_reject", "non_atp_hard_gate_failed", "", "primary_support" if primary_count else "not_evaluable"
    if not ppi_pass:
        return "ppi_low_relevance_reject", "ppi_relationship_hard_gate_failed", "", "primary_support" if primary_count else "not_evaluable"
    if not membrane_pass:
        return "membrane_geometry_reject", "lower_lateral_hard_gate_failed", "", "primary_support" if primary_count else "not_evaluable"
    if not dimer_pass:
        return "dimer_buried_reject", "dimer_accessibility_hard_gate_failed", "", "primary_support" if primary_count else "not_evaluable"
    if primary_count >= int(thresholds["required_primary_state_count_for_primary"]):
        return "accepted_primary_pocket", "", "all_hard_gates_pass_with_primary_state_recurrence", "primary_support"
    if primary_count >= int(thresholds["allow_secondary_if_primary_state_count"]):
        return "accepted_secondary_cryptic", "", "all_hard_gates_pass_single_primary_state_review", "secondary_review"
    return "manual_review", "", "hard_gates_pass_but_primary_support_ambiguous", "manual_review"


def _soft_score(
    hard_gate_pass: bool,
    family: dict[str, Any],
    ppi_eval: dict[str, Any],
    atp_metrics: dict[str, Any],
    membrane_eval: dict[str, Any],
    dimer_eval: dict[str, Any],
) -> tuple[str, str]:
    if not hard_gate_pass:
        return "", "not_applied_hard_gate_failed"
    score = 0.0
    score += {"orthosteric": 3.0, "rim": 2.0, "allosteric_near": 1.0}.get(str(ppi_eval.get("ppi_relationship_class")), 0.0)
    score += min(_int(family.get("primary_state_count")), 2)
    lateral = safe_float(membrane_eval.get("lateral_score"))
    if lateral is not None:
        score += max(0.0, min(1.0, lateral))
    return format_float(score), "applied_after_hard_gate_pass"


def _accepted_row(family: dict[str, Any], gate_row: dict[str, Any]) -> dict[str, Any]:
    gate_class = str(gate_row.get("gate_class", ""))
    if gate_class == "accepted_primary_pocket":
        tier = "primary"
    elif gate_class == "accepted_secondary_cryptic":
        tier = "secondary_review"
    elif gate_class == "reference_only":
        tier = "reference_only"
    else:
        tier = "rejected"
    return {
        "pocket_family_id": family.get("pocket_family_id", ""),
        "gate_class": gate_class,
        "acceptance_tier": tier,
        "representative_candidate_pocket_id": family.get("representative_candidate_pocket_id", ""),
        "representative_state_id": family.get("representative_state_id", ""),
        "representative_state_role": family.get("representative_state_role", ""),
        "representative_protomer_id": family.get("representative_protomer_id", ""),
        "member_candidate_pocket_ids": family.get("member_candidate_pocket_ids", ""),
        "member_state_ids": family.get("member_state_ids", ""),
        "member_state_roles": family.get("member_state_roles", ""),
        "member_protomer_ids": family.get("member_protomer_ids", ""),
        "primary_state_count": family.get("primary_state_count", ""),
        "reference_state_count": family.get("reference_state_count", ""),
        "family_center_x": family.get("family_center_x", ""),
        "family_center_y": family.get("family_center_y", ""),
        "family_center_z": family.get("family_center_z", ""),
        "residue_union_uniprot": family.get("residue_union_uniprot", ""),
        "ppi_relationship_class": gate_row.get("ppi_relationship_class", ""),
        "non_atp_pass": gate_row.get("non_atp_pass", ""),
        "ppi_relationship_pass": gate_row.get("ppi_relationship_pass", ""),
        "lower_lateral_pass": gate_row.get("lower_lateral_pass", ""),
        "dimer_accessibility_pass": gate_row.get("dimer_accessibility_pass", ""),
        "primary_state_support_or_review_flag": gate_row.get("primary_state_support_or_review_flag", ""),
        "acceptance_reason": gate_row.get("review_reason", ""),
        "review_notes": gate_row.get("reject_reason", ""),
        "warnings": gate_row.get("warnings", ""),
    }


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).lower() in {"true", "1", "yes", "pass"}


def _int(value: Any) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return 0


def _count_gate(rows: list[dict[str, Any]], gate_class: str) -> int:
    return sum(1 for row in rows if row.get("gate_class") == gate_class)


def _status(blockers: list[str], warnings: list[str], gate_rows: list[dict[str, Any]]) -> str:
    if blockers:
        return "FAIL"
    if warnings or not gate_rows:
        return "PASS_WITH_WARNINGS"
    return "PASS"


def _forbidden_export_warnings(ctx: RunContext) -> list[str]:
    warnings: list[str] = []
    for relative in FORBIDDEN_EXPORT_OUTPUTS:
        if (ctx.run_dir / relative).exists():
            warnings.append("preexisting_m2_8_or_m3_output_not_modified:{0}".format(relative))
    return warnings


def _write_summary(ctx: RunContext, report: PocketGateReport) -> Path:
    target = ctx.require_within_run_dir(ctx.reports_dir / "m2_7_ppi_membrane_dimer_pocket_gates.md")
    lines = [
        "# M2.7 PPI/membrane/dimer Pocket Gates",
        "",
        "Status: {0}".format(report.status),
        "Pocket families evaluated: {0}".format(report.family_count),
        "accepted_primary_pocket: {0}".format(report.accepted_primary_count),
        "accepted_secondary_cryptic: {0}".format(report.accepted_secondary_count),
        "manual_review: {0}".format(report.manual_review_count),
        "reference_only: {0}".format(report.reference_only_count),
        "atp_reject: {0}".format(report.atp_reject_count),
        "ppi_low_relevance_reject: {0}".format(report.ppi_reject_count),
        "membrane_geometry_reject: {0}".format(report.membrane_reject_count),
        "dimer_buried_reject: {0}".format(report.dimer_reject_count),
        "mapping_or_origin_reject: {0}".format(report.mapping_or_origin_reject_count),
        "",
        "Hard-gate failures are not rescued by soft scores. The M2.8 export package was not created.",
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
