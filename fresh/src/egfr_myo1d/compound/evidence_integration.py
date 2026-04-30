"""M3-T10 integrated evidence table and candidate hypothesis tiering.

This module integrates already-produced M2/M3 evidence. It does not run Vina,
submit qsub, invoke PBS runners, rerun pose attribution/clustering/anchor
convergence, create broad scans, or expose private compound data.
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
from statistics import mean, median
from typing import Any

from egfr_myo1d.compound.candidate_tiering import (
    TIER_RANK,
    as_float,
    assign_candidate_tier,
    bool_text,
    boolish,
    compute_soft_scores,
    hard_gate_failures,
)
from egfr_myo1d.compound.confidentiality import PUBLIC_COMPOUND_IDS, load_private_map, public_output_scan_paths, scan_internal_id_leaks
from egfr_myo1d.core.logging_utils import append_failed_job, append_job_status, append_phase_status
from egfr_myo1d.core.run_context import RunContext, ensure_within


DEFAULT_CONFIG: dict[str, Any] = {
    "schema_version": "m3_evidence_tiering_config_v1",
    "input": {
        "require_m2_accepted_pockets": True,
        "require_m2_pocket_gate_qc": True,
        "require_m2_ppi_trace": True,
        "require_m3_anchor_convergence": True,
        "require_m3_compound_pocket_support": True,
        "require_anchor_qc_ready": True,
        "allow_partial": False,
    },
    "compounds": {"expected_public_compound_ids": ["Cpd-A", "Cpd-B", "Cpd-C"], "require_all_expected_compounds_tested_for_tier1": False},
    "states": {
        "primary_state_ids": ["EGFR_160-185", "EGFR_170-200"],
        "reference_state_ids": ["3GT8_raw"],
        "reference_state_role_values": ["reference", "crystallographic_reference_control", "reference_control"],
        "do_not_count_ab_symmetry_as_independent_state_support": True,
    },
    "hard_gates": {
        "reject_if": {
            "atp_migration_present": True,
            "outside_accepted_pocket": True,
            "membrane_or_dimer_conflict": True,
            "missing_mapping": True,
            "confidentiality_violation": True,
            "internal_id_leak": True,
            "smiles_leak": True,
        }
    },
    "anchor_classes": {
        "tier1_allowed": ["multi_compound_primary_state_supported"],
        "tier2_allowed": ["single_compound_state_robust", "single_state_multi_compound"],
        "tier3_allowed": ["weak_or_scattered", "reference_only", "no_valid_cluster_support"],
        "reject_anchor_classes": ["fail_missing_input"],
    },
    "soft_scores": {
        "enabled": True,
        "max_score": 100.0,
        "weights": {
            "S_ppi_patch_support": 15.0,
            "S_pocket_geometry": 15.0,
            "S_state_support": 15.0,
            "S_symmetry_support": 5.0,
            "S_compound_anchor_convergence": 20.0,
            "S_pose_convergence": 15.0,
            "S_ppi_disruption_geometry": 10.0,
            "S_chemical_qc": 5.0,
            "S_affinity_normalized": 0.0,
        },
        "affinity": {
            "compute_descriptive_affinity_score": True,
            "use_affinity_for_tier_assignment": False,
            "use_affinity_for_candidate_promotion": False,
            "use_affinity_for_best_compound_selection": False,
            "use_affinity_for_primary_ranking": False,
        },
    },
    "broad_scan": {"read_if_present": False, "allow_only_with_cli_flag": True, "exploratory_only": True, "cannot_promote_to_tier1": True, "cannot_override_m2_hard_gates": True},
    "confidentiality": {
        "allowed_public_compound_ids": ["Cpd-A", "Cpd-B", "Cpd-C"],
        "scan_public_outputs": True,
        "fail_on_internal_id_leak": True,
        "fail_on_smiles_leak": True,
        "fail_on_candidate_overclaim": True,
    },
}

ANCHOR_REQUIRED = [
    "run_id", "m2_run_id", "profile", "anchor_id", "anchor_scope", "pocket_family_id",
    "anchor_convergence_class", "anchor_convergence_status", "allowed_for_evidence_integration",
    "affinity_used_for_ranking", "affinity_used_for_best_compound_selection",
]

SUPPORT_REQUIRED = [
    "run_id", "m2_run_id", "profile", "support_id", "compound_public_id", "pocket_family_id",
    "state_id", "state_role", "allowed_for_anchor_convergence", "allowed_for_evidence_integration",
    "support_status", "support_reject_reason", "affinity_used_for_ranking",
]

POCKET_REQUIRED = ["pocket_family_id"]

EVIDENCE_FIELDS = [
    "run_id", "m2_run_id", "profile", "evidence_id", "candidate_hypothesis_id", "candidate_scope",
    "compound_public_id", "pocket_family_id", "state_id", "state_role", "primary_state_flag",
    "reference_state_flag", "protomer_ids", "box_ids", "m2_pocket_accepted", "m2_pocket_acceptance_status",
    "non_atp_pass", "ppi_relationship_pass", "lower_lateral_pass", "dimer_accessibility_pass",
    "membrane_geometry_pass", "pose_retention_pass", "atp_migration_absent", "compound_convergence_pass",
    "primary_state_support", "reference_only_support", "confidentiality_pass", "no_membrane_or_dimer_conflict",
    "accepted_pocket_trace_present", "compound_trace_present", "anchor_trace_present", "anchor_convergence_class",
    "anchor_convergence_status", "anchor_scope", "support_status", "support_reject_reason",
    "candidate_hard_gate_status", "candidate_hard_gate_failures", "candidate_reject_reason",
    "S_ppi_patch_support", "S_pocket_geometry", "S_state_support", "S_symmetry_support",
    "S_compound_anchor_convergence", "S_pose_convergence", "S_affinity_normalized", "S_ppi_disruption_geometry",
    "S_chemical_qc", "candidate_soft_score_total", "candidate_soft_score_non_affinity", "candidate_priority_score",
    "affinity_used_for_tier_assignment", "affinity_used_for_candidate_promotion",
    "affinity_used_for_best_compound_selection", "vina_affinity_median_kcal_mol", "vina_affinity_mean_kcal_mol",
    "vina_affinity_best_kcal_mol", "vina_affinity_worst_kcal_mol", "mechanism_class_counts",
    "dominant_mechanism_class", "cluster_ids_supporting", "support_ids_supporting", "primary_state_ids_supported",
    "reference_state_ids_supported", "symmetry_support", "protomer_support_summary",
    "compound_public_ids_supported_for_pocket", "num_compounds_supported_for_pocket", "num_primary_states_supported",
    "num_supporting_clusters", "max_cluster_fraction", "median_cluster_fraction", "min_within_pocket_fraction",
    "max_atp_migration_fraction", "max_membrane_penetration_fraction", "max_dimer_interface_clash_fraction",
    "ppi_hotspot_contact_median_max", "ppi_hotspot_contact_mean_max", "nearest_ppi_hotspot_distance_median_min",
    "ppi_rim_distance_median_min", "evidence_notes", "source_accepted_pockets",
    "source_accepted_pockets_sha256", "source_pocket_gate_qc", "source_pocket_gate_qc_sha256",
    "source_ppi_consensus", "source_ppi_consensus_sha256", "source_compound_pocket_support",
    "source_compound_pocket_support_sha256", "source_compound_anchor_convergence",
    "source_compound_anchor_convergence_sha256", "source_compound_pose_clusters",
    "source_compound_pose_clusters_sha256", "evaluated_at",
]

CANDIDATE_FIELDS = [
    "run_id", "m2_run_id", "profile", "candidate_hypothesis_id", "candidate_scope", "candidate_tier",
    "candidate_tier_code", "tier_rank", "compound_public_id", "pocket_family_id", "candidate_name_public",
    "candidate_hypothesis_statement", "m2_pocket_accepted", "non_atp_pass", "ppi_relationship_pass",
    "lower_lateral_pass", "dimer_accessibility_pass", "pose_retention_pass", "atp_migration_absent",
    "compound_convergence_pass", "primary_state_support", "reference_only_support", "confidentiality_pass",
    "no_membrane_or_dimer_conflict", "all_tier1_hard_gates_pass", "all_tier2_hard_gates_pass",
    "candidate_reject_reason", "candidate_accept_reason", "candidate_risk_notes", "anchor_convergence_class",
    "anchor_convergence_status", "anchor_scope", "primary_state_ids_supported", "reference_state_ids_supported",
    "symmetry_support", "protomer_support_summary", "compound_public_ids_supported_for_pocket",
    "num_compounds_supported_for_pocket", "num_primary_states_supported", "best_compound_public_id",
    "best_compound_selection_rule", "affinity_used_for_best_compound_selection", "S_ppi_patch_support",
    "S_pocket_geometry", "S_state_support", "S_symmetry_support", "S_compound_anchor_convergence",
    "S_pose_convergence", "S_affinity_normalized", "S_ppi_disruption_geometry", "S_chemical_qc",
    "candidate_soft_score_total", "candidate_soft_score_non_affinity", "candidate_priority_score",
    "affinity_used_for_tier_assignment", "affinity_used_for_candidate_promotion",
    "median_affinity_across_supported", "mean_affinity_across_supported", "best_affinity_across_supported",
    "worst_affinity_across_supported", "dominant_mechanism_class", "mechanism_class_counts",
    "cluster_ids_supporting", "support_ids_supporting", "evidence_ids_supporting",
    "source_pocket_compound_evidence_table", "source_pocket_compound_evidence_table_sha256",
    "source_compound_anchor_convergence", "source_compound_anchor_convergence_sha256",
    "source_accepted_pockets", "source_accepted_pockets_sha256", "evaluated_at",
]

REJECT_FIELDS = [
    "run_id", "m2_run_id", "profile", "candidate_hypothesis_id", "candidate_scope", "compound_public_id",
    "pocket_family_id", "state_ids", "candidate_tier", "candidate_reject_reason", "hard_gate_failures",
    "evidence_notes", "recommended_fix", "evaluated_at",
]

QC_FIELDS = ["check_id", "category", "status", "severity", "details", "recommended_fix"]

REQUIRED_CHECKS = [
    "m2_accepted_pockets_present", "m2_accepted_pockets_schema_valid", "m2_pocket_gate_qc_present",
    "m2_ppi_trace_present", "m3_compound_anchor_convergence_present", "m3_compound_pocket_support_present",
    "m3_anchor_qc_present", "m3_anchor_ready_or_partial_allowed", "profile_selected",
    "compound_public_ids_valid", "accepted_pocket_trace_valid", "compound_support_trace_valid",
    "anchor_support_trace_valid", "primary_states_identified", "reference_states_identified",
    "ab_symmetry_not_counted_as_independent_state", "reference_only_not_promoted",
    "broad_scan_not_used_unless_explicit", "broad_scan_not_promoted_to_tier1", "hard_gates_computed",
    "soft_scores_computed", "hard_gates_separated_from_soft_scores", "tier1_requires_m2_acceptance",
    "tier1_requires_compound_convergence", "tier1_requires_no_atp_migration", "tier1_requires_primary_support",
    "tier2_requires_no_atp_migration", "atp_confounded_support_quarantined",
    "outside_pocket_support_rejected", "membrane_or_dimer_conflict_rejected", "missing_mapping_rejected",
    "candidate_reject_reasons_explicit", "all_evidence_rows_represented", "no_silent_evidence_drop",
    "pocket_compound_evidence_table_written", "final_m3_candidate_hypotheses_written",
    "rejected_candidate_reasons_written", "final_candidate_gate_qc_written", "evidence_integration_qc_written",
    "candidate_tiering_qc_written", "task_report_written", "raw_pose_affinity_not_used_for_promotion",
    "best_compound_not_selected_by_affinity", "multi_compound_convergence_reported_separately_from_affinity",
    "no_vina_invoked_by_evidence_tiering", "no_qsub_invoked_by_evidence_tiering",
    "no_runner_invoked_by_evidence_tiering", "no_pose_attribution_rerun_attempted",
    "no_pose_clustering_rerun_attempted", "no_anchor_convergence_rerun_attempted",
    "no_broad_anchor_scan_attempted", "no_ligand_preparation_attempted", "no_receptor_preparation_attempted",
    "no_internal_id_leak", "no_smiles_logged", "no_ligand_coordinates_logged",
    "no_receptor_coordinates_logged", "no_pose_coordinates_logged_outside_pdbqt",
    "candidate_claims_are_hypothesis_only", "old_workflow_not_used", "non_goals_preserved",
]


@dataclass
class M3EvidenceTieringResult:
    status: str
    m3_t11_report_ready: bool
    blockers: list[str]
    warnings: list[str]
    evidence_csv: Path
    candidates_csv: Path
    reject_csv: Path
    qc_csv: Path
    qc_json: Path
    evidence_qc_csv: Path
    tiering_qc_csv: Path
    report_md: Path
    phase3_log: Path
    inputs: dict[str, bool]
    evidence_context: dict[str, Any]
    evidence_integration: dict[str, int]
    candidate_tiering: dict[str, int]
    hard_gates: dict[str, Any]
    soft_scores: dict[str, Any]
    mechanism_classes_by_candidate: dict[str, int]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _sha256(path: Path | None) -> str:
    if not path or not path.is_file():
        return ""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fmt(value: float | None, digits: int = 6) -> str:
    return "" if value is None else f"{value:.{digits}f}"


def _join(values: set[str] | list[str], sep: str = ";") -> str:
    return sep.join(sorted({str(value) for value in values if str(value)}))


def _stat(values: list[float], which: str) -> float | None:
    if not values:
        return None
    if which == "median":
        return median(values)
    if which == "mean":
        return mean(values)
    if which == "best":
        return min(values)
    if which == "worst":
        return max(values)
    return None


def _write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def _safe_private_map_path(ctx: RunContext, private_map: Path | None) -> Path:
    path = private_map or (ctx.fresh_root / "data" / "private" / "compound_id_map.csv")
    if not path.is_absolute():
        path = ctx.repo_root / path
    return ensure_within(path.resolve(), ctx.fresh_root)


def _table_path(ctx: RunContext, given: Path | None, default_rel: str) -> Path:
    path = given or (ctx.run_dir / default_rel)
    if not path.is_absolute():
        path = ctx.repo_root / path
    return path.resolve()


def _first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.is_file():
            return path.resolve()
    return None


def _load_json(path: Path | None) -> dict[str, Any] | None:
    if not path or not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _load_config(ctx: RunContext, config_path: Path | None) -> dict[str, Any]:
    config = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (ctx.fresh_root / "configs" / "compound_evidence_tiering.yaml")
    if path and path.is_file():
        try:
            import yaml  # type: ignore

            loaded = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        except Exception:
            loaded = {}
        for section, values in loaded.items():
            if isinstance(values, dict) and isinstance(config.get(section), dict):
                config[section].update(values)
            else:
                config[section] = values
    snapshot = ctx.run_dir / "phase3_compounds" / "config_snapshots" / "compound_evidence_tiering.resolved.yaml"
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    try:
        import yaml  # type: ignore

        snapshot.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    except Exception:
        snapshot.write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")
    return config


def _filter_profile(rows: list[dict[str, str]], profile: str | None) -> tuple[list[dict[str, str]], str | None, list[str]]:
    profiles = sorted({row.get("profile", "") for row in rows if row.get("profile")})
    selected = profile or (profiles[0] if len(profiles) == 1 else None)
    return [row for row in rows if not selected or row.get("profile") == selected], selected, profiles


def _ensure_logs(ctx: RunContext) -> None:
    ctx.create_directories()
    for path in [
        ctx.logs_dir / "master.log", ctx.logs_dir / "phase_status.jsonl", ctx.logs_dir / "job_status.jsonl",
        ctx.logs_dir / "phase3_compounds.log", ctx.errors_dir / "error_summary.txt", ctx.errors_dir / "failed_jobs.csv",
    ]:
        safe = ctx.require_within_run_dir(path)
        safe.parent.mkdir(parents=True, exist_ok=True)
        if not safe.exists():
            safe.touch()
    if (ctx.errors_dir / "failed_jobs.csv").stat().st_size == 0:
        (ctx.errors_dir / "failed_jobs.csv").write_text("timestamp,job_name,status,message\n", encoding="utf-8")


def _m2_bool(row: dict[str, str] | None, aliases: list[str], default: bool = True) -> bool:
    if not row:
        return False
    lowered = {key.lower(): value for key, value in row.items()}
    for alias in aliases:
        if alias.lower() in lowered and str(lowered[alias.lower()]).strip() != "":
            value = str(lowered[alias.lower()]).strip().lower()
            if value in {"false", "f", "0", "no", "n", "fail", "failed", "reject", "rejected"}:
                return False
            return True
    return default


def _is_primary(state_id: str, config: dict[str, Any]) -> bool:
    return state_id in set(config["states"]["primary_state_ids"])


def _is_reference(state_id: str, state_role: str, config: dict[str, Any]) -> bool:
    return state_id in set(config["states"]["reference_state_ids"]) or state_role in set(config["states"]["reference_state_role_values"])


def _make_pocket_maps(rows: list[dict[str, str]]) -> tuple[dict[str, dict[str, str]], dict[tuple[str, str], dict[str, str]], set[str]]:
    by_family: dict[str, dict[str, str]] = {}
    by_pair: dict[tuple[str, str], dict[str, str]] = {}
    states: set[str] = set()
    for row in rows:
        family = row.get("pocket_family_id") or row.get("family_id") or row.get("pocket_id") or ""
        state = row.get("state_id", "")
        if family and family not in by_family:
            by_family[family] = row
        if family and state:
            by_pair[(family, state)] = row
            states.add(state)
    return by_family, by_pair, states


def _support_sort(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (row.get("pocket_family_id", ""), row.get("state_id", ""), row.get("compound_public_id", ""), row.get("support_id", ""))


def _anchor_sort(row: dict[str, str]) -> tuple[int, str, str, str]:
    scope_rank = {"primary_cross_state_pocket": 0, "state_pocket": 1, "reference_pocket": 2, "all_states_pocket": 3}
    return (scope_rank.get(row.get("anchor_scope", ""), 9), row.get("pocket_family_id", ""), row.get("state_id", ""), row.get("anchor_id", ""))


def _support_has_atp(row: dict[str, str]) -> bool:
    return (as_float(row.get("atp_migration_fraction_max")) or 0.0) > 0.0 or row.get("support_reject_reason") == "ATP_migration" or "ATP_like_reject" in (row.get("mechanism_class_counts") or "")


def _support_has_membrane_dimer(row: dict[str, str]) -> bool:
    return (as_float(row.get("membrane_penetration_fraction_max")) or 0.0) > 0.0 or (as_float(row.get("dimer_interface_clash_fraction_max")) or 0.0) > 0.0


def _anchor_for_support(anchor_rows: list[dict[str, str]], support: dict[str, str]) -> dict[str, str] | None:
    support_id = support.get("support_id", "")
    pocket = support.get("pocket_family_id", "")
    state = support.get("state_id", "")
    compound = support.get("compound_public_id", "")
    candidates = [
        row for row in anchor_rows
        if row.get("pocket_family_id") == pocket
        and (support_id in (row.get("support_ids_supporting") or "").split("|") or support_id in (row.get("support_ids_supporting") or "").split(";"))
    ]
    if not candidates:
        candidates = [
            row for row in anchor_rows
            if row.get("pocket_family_id") == pocket and row.get("anchor_scope") in {"primary_cross_state_pocket", "state_pocket", "all_states_pocket"}
            and (not row.get("state_id") or row.get("state_id") == state)
            and (not compound or compound in (row.get("compound_public_ids_supported") or row.get("compound_public_ids_tested") or compound))
        ]
    if not candidates:
        return None
    return sorted(candidates, key=_anchor_sort)[0]


def _build_evidence_row(
    *,
    ctx: RunContext,
    m2_run_id: str,
    profile: str | None,
    index: int,
    support: dict[str, str] | None,
    pocket: dict[str, str] | None,
    anchor: dict[str, str] | None,
    cfg: dict[str, Any],
    timestamp: str,
    sources: dict[str, str],
    source_hashes: dict[str, str],
    confidentiality_pass: bool,
    candidate_scope: str = "compound_pocket",
) -> dict[str, Any]:
    raw_compound = (support or {}).get("compound_public_id", "")
    compound = raw_compound if raw_compound in PUBLIC_COMPOUND_IDS else ("NON_PUBLIC_COMPOUND_ID_REDACTED" if raw_compound else "")
    pocket_family = (support or {}).get("pocket_family_id") or (pocket or {}).get("pocket_family_id", "")
    state_id = (support or {}).get("state_id") or (pocket or {}).get("state_id", "")
    state_role = (support or {}).get("state_role") or ("reference" if _is_reference(state_id, "", cfg) else ("primary" if _is_primary(state_id, cfg) else ""))
    primary = _is_primary(state_id, cfg) or (support or {}).get("primary_state_flag") == "true"
    reference = _is_reference(state_id, state_role, cfg) or (support or {}).get("reference_state_flag") == "true"
    m2_trace = pocket is not None
    compound_trace = support is not None and compound in PUBLIC_COMPOUND_IDS
    anchor_trace = anchor is not None
    accepted = m2_trace
    non_atp = _m2_bool(pocket, ["non_atp_pass", "non_atp_gate_pass", "atp_excluded", "non_atp"], default=True) and not _support_has_atp(support or {})
    ppi_pass = _m2_bool(pocket, ["ppi_relationship_pass", "ppi_adjacent_pass", "ppi_gate_pass"], default=True)
    if pocket and not any(key.lower() in {k.lower() for k in pocket} for key in ["ppi_relationship_pass", "ppi_adjacent_pass", "ppi_gate_pass"]):
        relationship = " ".join(str(pocket.get(key, "")) for key in pocket)
        ppi_pass = bool(re.search(r"ppi|adjacent|rim|hotspot|contact", relationship, re.IGNORECASE)) or ppi_pass
    lower_lateral = _m2_bool(pocket, ["lower_lateral_pass", "membrane_lower_lateral_pass", "lateral_pass"], default=True)
    dimer_access = _m2_bool(pocket, ["dimer_accessibility_pass", "dimer_accessible", "dimer_access_pass"], default=True)
    membrane_geometry = _m2_bool(pocket, ["membrane_geometry_pass", "membrane_compatible", "lower_lateral_pass"], default=True)
    pose_retention = boolish((support or {}).get("allowed_for_anchor_convergence")) or (as_float((support or {}).get("within_pocket_fraction_min")) or 0.0) >= 0.95
    atp_absent = not _support_has_atp(support or {})
    membrane_absent = not _support_has_membrane_dimer(support or {})
    anchor_class = (anchor or {}).get("anchor_convergence_class") or ("no_valid_cluster_support" if support else "not_applicable")
    compound_conv = anchor_class in {"multi_compound_primary_state_supported", "single_compound_state_robust", "single_state_multi_compound"} and boolish((anchor or {}).get("allowed_for_evidence_integration"))
    primary_support = primary or (as_float((anchor or {}).get("primary_state_support")) or 0.0) > 0.0
    reference_only = reference and not primary_support or anchor_class == "reference_only"
    flags = {
        "m2_pocket_accepted": accepted,
        "non_atp_pass": non_atp,
        "ppi_relationship_pass": ppi_pass,
        "lower_lateral_pass": lower_lateral,
        "dimer_accessibility_pass": dimer_access,
        "pose_retention_pass": pose_retention,
        "atp_migration_absent": atp_absent,
        "compound_convergence_pass": compound_conv,
        "primary_state_support": primary_support and not reference_only,
        "reference_only_support": reference_only,
        "confidentiality_pass": confidentiality_pass,
        "no_membrane_or_dimer_conflict": membrane_absent,
        "accepted_pocket_trace_present": m2_trace,
        "compound_trace_present": compound_trace,
        "anchor_trace_present": anchor_trace,
        "symmetry_support": (support or {}).get("protomer_symmetry_support") or (anchor or {}).get("symmetry_support"),
        "dominant_mechanism_class": (support or {}).get("mechanism_class_majority") or (anchor or {}).get("dominant_mechanism_class"),
        "mechanism_class_counts": (support or {}).get("mechanism_class_counts") or (anchor or {}).get("mechanism_class_counts"),
    }
    tier, _code, _rank, reject_reason, _accept = assign_candidate_tier(flags, anchor_class, (support or {}).get("support_reject_reason", "none"))
    failures = hard_gate_failures(flags)
    evidence_id = f"evidence_{index:04d}"
    anchor_compounds_public = _join({item for item in str((anchor or {}).get("compound_public_ids_supported", compound)).replace("|", ";").split(";") if item in PUBLIC_COMPOUND_IDS})
    candidate_hypothesis_id = "candidate_{0}_{1}".format(compound or "pocket", pocket_family or "unknown")
    row: dict[str, Any] = {
        "run_id": ctx.run_id, "m2_run_id": m2_run_id, "profile": profile, "evidence_id": evidence_id,
        "candidate_hypothesis_id": candidate_hypothesis_id, "candidate_scope": candidate_scope,
        "compound_public_id": compound, "pocket_family_id": pocket_family, "state_id": state_id,
        "state_role": state_role, "primary_state_flag": bool_text(primary), "reference_state_flag": bool_text(reference),
        "protomer_ids": (support or {}).get("protomer_ids_with_support", ""), "box_ids": (support or {}).get("box_ids_with_support", ""),
        "m2_pocket_accepted": bool_text(accepted), "m2_pocket_acceptance_status": "accepted" if accepted else "missing_or_rejected",
        "non_atp_pass": bool_text(non_atp), "ppi_relationship_pass": bool_text(ppi_pass),
        "lower_lateral_pass": bool_text(lower_lateral), "dimer_accessibility_pass": bool_text(dimer_access),
        "membrane_geometry_pass": bool_text(membrane_geometry), "pose_retention_pass": bool_text(pose_retention),
        "atp_migration_absent": bool_text(atp_absent), "compound_convergence_pass": bool_text(compound_conv),
        "primary_state_support": bool_text(primary_support and not reference_only), "reference_only_support": bool_text(reference_only),
        "confidentiality_pass": bool_text(confidentiality_pass), "no_membrane_or_dimer_conflict": bool_text(membrane_absent),
        "accepted_pocket_trace_present": bool_text(m2_trace), "compound_trace_present": bool_text(compound_trace),
        "anchor_trace_present": bool_text(anchor_trace), "anchor_convergence_class": anchor_class,
        "anchor_convergence_status": (anchor or {}).get("anchor_convergence_status", "NOT_APPLICABLE"),
        "anchor_scope": (anchor or {}).get("anchor_scope", ""), "support_status": (support or {}).get("support_status", "NOT_APPLICABLE"),
        "support_reject_reason": (support or {}).get("support_reject_reason", "missing_compound_trace" if not support else "none"),
        "candidate_hard_gate_status": "PASS" if not failures else ("WARN" if tier == "Tier 3" else "REJECT"),
        "candidate_hard_gate_failures": ";".join(failures), "candidate_reject_reason": reject_reason,
        "affinity_used_for_tier_assignment": "false", "affinity_used_for_candidate_promotion": "false",
        "affinity_used_for_best_compound_selection": "false",
        "vina_affinity_median_kcal_mol": (support or {}).get("vina_affinity_median_kcal_mol", "") or (anchor or {}).get("median_affinity_across_supported", ""),
        "vina_affinity_mean_kcal_mol": (support or {}).get("vina_affinity_mean_kcal_mol", "") or (anchor or {}).get("mean_affinity_across_supported", ""),
        "vina_affinity_best_kcal_mol": (support or {}).get("vina_affinity_best_kcal_mol", "") or (anchor or {}).get("best_affinity_across_supported", ""),
        "vina_affinity_worst_kcal_mol": (support or {}).get("vina_affinity_worst_kcal_mol", "") or (anchor or {}).get("worst_affinity_across_supported", ""),
        "mechanism_class_counts": flags["mechanism_class_counts"] or "", "dominant_mechanism_class": flags["dominant_mechanism_class"] or "",
        "cluster_ids_supporting": (support or {}).get("cluster_ids_supporting", "") or (anchor or {}).get("cluster_ids_supporting", ""),
        "support_ids_supporting": (support or {}).get("support_id", "") or (anchor or {}).get("support_ids_supporting", ""),
        "primary_state_ids_supported": (anchor or {}).get("primary_state_ids_supported", state_id if primary else ""),
        "reference_state_ids_supported": (anchor or {}).get("reference_state_ids_supported", state_id if reference else ""),
        "symmetry_support": bool_text(boolish(flags.get("symmetry_support"))),
        "protomer_support_summary": (anchor or {}).get("protomer_support_summary", ""),
        "compound_public_ids_supported_for_pocket": anchor_compounds_public or compound,
        "num_compounds_supported_for_pocket": len([item for item in (anchor_compounds_public or compound).split(";") if item]),
        "num_primary_states_supported": (anchor or {}).get("primary_state_support", "1" if primary_support and not reference_only else "0"),
        "num_supporting_clusters": (support or {}).get("num_anchor_allowed_clusters", "0"),
        "max_cluster_fraction": (support or {}).get("max_cluster_fraction", "") or (anchor or {}).get("max_cluster_fraction_across_supported", ""),
        "median_cluster_fraction": (support or {}).get("median_cluster_fraction", "") or (anchor or {}).get("median_cluster_fraction_across_supported", ""),
        "min_within_pocket_fraction": (support or {}).get("within_pocket_fraction_min", "") or (anchor or {}).get("min_within_pocket_fraction_across_supported", ""),
        "max_atp_migration_fraction": (support or {}).get("atp_migration_fraction_max", "") or (anchor or {}).get("max_atp_migration_fraction_across_supported", ""),
        "max_membrane_penetration_fraction": (support or {}).get("membrane_penetration_fraction_max", "") or (anchor or {}).get("max_membrane_penetration_fraction_across_supported", ""),
        "max_dimer_interface_clash_fraction": (support or {}).get("dimer_interface_clash_fraction_max", "") or (anchor or {}).get("max_dimer_interface_clash_fraction_across_supported", ""),
        "ppi_hotspot_contact_median_max": (support or {}).get("ppi_hotspot_contact_median_max", ""),
        "ppi_hotspot_contact_mean_max": (support or {}).get("ppi_hotspot_contact_mean_max", ""),
        "nearest_ppi_hotspot_distance_median_min": (support or {}).get("nearest_ppi_hotspot_distance_median_min", ""),
        "ppi_rim_distance_median_min": (support or {}).get("ppi_rim_distance_median_min", ""),
        "evidence_notes": "hard gates evaluated separately from soft scores; affinity descriptive only",
        "source_accepted_pockets": sources.get("accepted_pockets", ""), "source_accepted_pockets_sha256": source_hashes.get("accepted_pockets", ""),
        "source_pocket_gate_qc": sources.get("pocket_gate_qc", ""), "source_pocket_gate_qc_sha256": source_hashes.get("pocket_gate_qc", ""),
        "source_ppi_consensus": sources.get("ppi_consensus", ""), "source_ppi_consensus_sha256": source_hashes.get("ppi_consensus", ""),
        "source_compound_pocket_support": sources.get("compound_pocket_support", ""), "source_compound_pocket_support_sha256": source_hashes.get("compound_pocket_support", ""),
        "source_compound_anchor_convergence": sources.get("compound_anchor_convergence", ""), "source_compound_anchor_convergence_sha256": source_hashes.get("compound_anchor_convergence", ""),
        "source_compound_pose_clusters": sources.get("compound_pose_clusters", ""), "source_compound_pose_clusters_sha256": source_hashes.get("compound_pose_clusters", ""),
        "evaluated_at": timestamp,
    }
    row.update({key: _fmt(value) for key, value in compute_soft_scores(row, cfg).items()})
    return row


def _aggregate_candidate(
    *,
    ctx: RunContext,
    m2_run_id: str,
    profile: str | None,
    candidate_id: str,
    rows: list[dict[str, Any]],
    cfg: dict[str, Any],
    timestamp: str,
    sources: dict[str, str],
    source_hashes: dict[str, str],
    evidence_sha: str,
    strict_tier1: bool,
) -> dict[str, Any]:
    best = sorted(rows, key=lambda row: (TIER_RANK.get(_tier_for_evidence(row, cfg, strict_tier1)[0], 9), -float(row.get("candidate_priority_score") or 0), row.get("evidence_id", "")))[0]
    flags = {key: all(boolish(row.get(key)) for row in rows) for key in [
        "m2_pocket_accepted", "non_atp_pass", "ppi_relationship_pass", "lower_lateral_pass",
        "dimer_accessibility_pass", "pose_retention_pass", "atp_migration_absent",
        "compound_convergence_pass", "confidentiality_pass", "no_membrane_or_dimer_conflict",
    ]}
    flags["primary_state_support"] = any(boolish(row.get("primary_state_support")) for row in rows)
    flags["reference_only_support"] = all(boolish(row.get("reference_only_support")) for row in rows) if rows else False
    flags["accepted_pocket_trace_present"] = any(boolish(row.get("accepted_pocket_trace_present")) for row in rows)
    flags["compound_trace_present"] = any(boolish(row.get("compound_trace_present")) for row in rows)
    flags["anchor_trace_present"] = any(boolish(row.get("anchor_trace_present")) for row in rows)
    flags["symmetry_support"] = any(boolish(row.get("symmetry_support")) for row in rows)
    flags["dominant_mechanism_class"] = best.get("dominant_mechanism_class", "")
    flags["mechanism_class_counts"] = best.get("mechanism_class_counts", "")
    anchor_class = best.get("anchor_convergence_class", "")
    support_reject = best.get("support_reject_reason", "none")
    tier, code, rank, reject_reason, accept_reason = assign_candidate_tier(flags, anchor_class, support_reject, strict_tier1)
    all_t1 = tier == "Tier 1"
    all_t2 = tier in {"Tier 1", "Tier 2"}
    compounds = set()
    for row in rows:
        compounds.update(item for item in str(row.get("compound_public_ids_supported_for_pocket") or row.get("compound_public_id") or "").split(";") if item)
    aff = [value for row in rows for value in [as_float(row.get("vina_affinity_median_kcal_mol"))] if value is not None]
    scores: dict[str, float] = {}
    for key in [
        "S_ppi_patch_support", "S_pocket_geometry", "S_state_support", "S_symmetry_support",
        "S_compound_anchor_convergence", "S_pose_convergence", "S_affinity_normalized",
        "S_ppi_disruption_geometry", "S_chemical_qc", "candidate_soft_score_total",
        "candidate_soft_score_non_affinity", "candidate_priority_score",
    ]:
        vals = [as_float(row.get(key)) for row in rows]
        vals = [value for value in vals if value is not None]
        scores[key] = max(vals) if vals else 0.0
    compound = best.get("compound_public_id", "")
    pocket = best.get("pocket_family_id", "")
    statement = "Computational hypothesis for review: {0} support near accepted pocket {1}; requires experimental follow-up.".format(compound or "pocket-level anchor", pocket)
    row = {
        "run_id": ctx.run_id, "m2_run_id": m2_run_id, "profile": profile, "candidate_hypothesis_id": candidate_id,
        "candidate_scope": best.get("candidate_scope", ""), "candidate_tier": tier, "candidate_tier_code": code, "tier_rank": rank,
        "compound_public_id": compound, "pocket_family_id": pocket,
        "candidate_name_public": "{0}::{1}".format(compound or "Pocket", pocket),
        "candidate_hypothesis_statement": statement, "m2_pocket_accepted": bool_text(flags["m2_pocket_accepted"]),
        "non_atp_pass": bool_text(flags["non_atp_pass"]), "ppi_relationship_pass": bool_text(flags["ppi_relationship_pass"]),
        "lower_lateral_pass": bool_text(flags["lower_lateral_pass"]), "dimer_accessibility_pass": bool_text(flags["dimer_accessibility_pass"]),
        "pose_retention_pass": bool_text(flags["pose_retention_pass"]), "atp_migration_absent": bool_text(flags["atp_migration_absent"]),
        "compound_convergence_pass": bool_text(flags["compound_convergence_pass"]), "primary_state_support": bool_text(flags["primary_state_support"]),
        "reference_only_support": bool_text(flags["reference_only_support"]), "confidentiality_pass": bool_text(flags["confidentiality_pass"]),
        "no_membrane_or_dimer_conflict": bool_text(flags["no_membrane_or_dimer_conflict"]),
        "all_tier1_hard_gates_pass": bool_text(all_t1), "all_tier2_hard_gates_pass": bool_text(all_t2),
        "candidate_reject_reason": reject_reason, "candidate_accept_reason": accept_reason,
        "candidate_risk_notes": "hypothesis-only; no validation claim; affinity not used for promotion",
        "anchor_convergence_class": anchor_class, "anchor_convergence_status": best.get("anchor_convergence_status", ""),
        "anchor_scope": best.get("anchor_scope", ""), "primary_state_ids_supported": _join({str(row.get("state_id", "")) for row in rows if boolish(row.get("primary_state_flag"))}),
        "reference_state_ids_supported": _join({str(row.get("state_id", "")) for row in rows if boolish(row.get("reference_state_flag"))}),
        "symmetry_support": bool_text(flags["symmetry_support"]), "protomer_support_summary": best.get("protomer_support_summary", ""),
        "compound_public_ids_supported_for_pocket": _join(compounds), "num_compounds_supported_for_pocket": len(compounds),
        "num_primary_states_supported": len({str(row.get("state_id", "")) for row in rows if boolish(row.get("primary_state_flag")) and boolish(row.get("primary_state_support"))}),
        "best_compound_public_id": sorted(compounds, key=lambda cid: PUBLIC_COMPOUND_IDS.index(cid) if cid in PUBLIC_COMPOUND_IDS else 99)[0] if compounds else compound,
        "best_compound_selection_rule": "deterministic_public_id_order", "affinity_used_for_best_compound_selection": "false",
        "median_affinity_across_supported": _fmt(_stat(aff, "median")), "mean_affinity_across_supported": _fmt(_stat(aff, "mean")),
        "best_affinity_across_supported": _fmt(_stat(aff, "best")), "worst_affinity_across_supported": _fmt(_stat(aff, "worst")),
        "dominant_mechanism_class": best.get("dominant_mechanism_class", ""), "mechanism_class_counts": best.get("mechanism_class_counts", ""),
        "cluster_ids_supporting": _join({item for row in rows for item in str(row.get("cluster_ids_supporting") or "").replace("|", ";").split(";") if item}),
        "support_ids_supporting": _join({item for row in rows for item in str(row.get("support_ids_supporting") or "").replace("|", ";").split(";") if item}),
        "evidence_ids_supporting": _join({str(row.get("evidence_id", "")) for row in rows}),
        "source_pocket_compound_evidence_table": sources.get("pocket_compound_evidence_table", ""),
        "source_pocket_compound_evidence_table_sha256": evidence_sha,
        "source_compound_anchor_convergence": sources.get("compound_anchor_convergence", ""),
        "source_compound_anchor_convergence_sha256": source_hashes.get("compound_anchor_convergence", ""),
        "source_accepted_pockets": sources.get("accepted_pockets", ""),
        "source_accepted_pockets_sha256": source_hashes.get("accepted_pockets", ""),
        "evaluated_at": timestamp,
    }
    row.update({key: _fmt(value) for key, value in scores.items() if key in CANDIDATE_FIELDS})
    row["affinity_used_for_tier_assignment"] = "false"
    row["affinity_used_for_candidate_promotion"] = "false"
    return row


def _tier_for_evidence(row: dict[str, Any], cfg: dict[str, Any], strict_tier1: bool) -> tuple[str, str, int, str, str]:
    flags = {key: row.get(key) for key in [
        "m2_pocket_accepted", "non_atp_pass", "ppi_relationship_pass", "lower_lateral_pass",
        "dimer_accessibility_pass", "pose_retention_pass", "atp_migration_absent",
        "compound_convergence_pass", "primary_state_support", "reference_only_support",
        "confidentiality_pass", "no_membrane_or_dimer_conflict", "accepted_pocket_trace_present",
        "compound_trace_present", "anchor_trace_present", "symmetry_support", "dominant_mechanism_class",
        "mechanism_class_counts",
    ]}
    return assign_candidate_tier(flags, str(row.get("anchor_convergence_class") or ""), str(row.get("support_reject_reason") or "none"), strict_tier1)


def _scan_hygiene(ctx: RunContext, private_entries: list[Any]) -> tuple[int, list[str], bool, bool, bool, bool, bool, list[str]]:
    leaks = scan_internal_id_leaks(ctx, private_entries)
    scan_paths = list(public_output_scan_paths(ctx))
    for root in [
        ctx.run_dir / "phase3_compounds" / "tables",
        ctx.run_dir / "phase3_compounds" / "qsub",
        ctx.logs_dir / "errors",
    ]:
        if root.is_dir():
            scan_paths.extend(path for path in sorted(root.rglob("*")) if path.is_file())
    coord_hits: list[str] = []
    smiles_logged = False
    ligand_coords = False
    receptor_coords = False
    pose_coords = False
    overclaim = False
    scanned: list[str] = []
    for path in sorted(set(scan_paths)):
        if "docking_outputs" in path.parts and path.suffix.lower() == ".pdbqt":
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        scanned.append(ctx.relative_to_repo(path))
        if re.search(r"\bSMILES\b|canonical_smiles|isomeric_smiles|[A-Za-z0-9@+\-\[\]\(\)=#$\\/]{8,}", text) and "compound_id_map" not in path.name:
            if re.search(r"\bSMILES\b|canonical_smiles|isomeric_smiles", text, re.IGNORECASE):
                smiles_logged = True
        if re.search(r"^(ATOM|HETATM)\s+\d+", text, re.MULTILINE):
            coord_hits.append(ctx.relative_to_repo(path))
            lower = path.as_posix().lower()
            ligand_coords = ligand_coords or "ligand" in lower
            receptor_coords = receptor_coords or "receptor" in lower
            pose_coords = pose_coords or "pose" in lower or "candidate" in lower or "evidence" in lower or "qc" in lower
        if re.search(r"validated inhibitor|proven binder|proven PPI inhibitor|clinically relevant drug candidate|\bdrug candidate\b|clinically validated", text, re.IGNORECASE):
            overclaim = True
    return len(leaks), coord_hits, smiles_logged, ligand_coords, receptor_coords, pose_coords, overclaim, scanned


def _severity(status: str) -> str:
    if status == "FAIL":
        return "BLOCKER"
    if status == "WARN":
        return "MAJOR"
    if status == "NOT_APPLICABLE":
        return "MINOR"
    return "INFO"


def _qc_row(check_id: str, status: str, details: str = "") -> dict[str, str]:
    return {"check_id": check_id, "category": "m3_t10_evidence_tiering", "status": status, "severity": _severity(status), "details": details or check_id, "recommended_fix": ""}


def _write_report(path: Path, ctx: RunContext, summary: dict[str, Any]) -> None:
    lines = [
        "# M3-T10 Evidence Tiering",
        "",
        "This report summarizes computational hypothesis tiers only. It does not assert binding, inhibition, PPI disruption, cellular activity, or clinical relevance.",
        "",
        f"- status: {summary['overall_status']}",
        f"- profile: {summary['profile']}",
        f"- evidence_rows: {summary['evidence_integration']['pocket_compound_evidence_rows']}",
        f"- candidate_rows: {summary['candidate_tiering']['candidate_rows']}",
        f"- tier1: {summary['candidate_tiering']['tier1']}",
        f"- tier2: {summary['candidate_tiering']['tier2']}",
        f"- tier3: {summary['candidate_tiering']['tier3']}",
        f"- reject: {summary['candidate_tiering']['reject']}",
        f"- m3_t11_report_ready: {str(summary['m3_t11_report_ready']).lower()}",
        "",
        "Affinity is recorded as descriptive metadata and is not used for promotion or best-compound selection.",
        "",
        "## Blockers",
    ]
    lines.extend(f"- {item}" for item in summary["blockers"] or ["none"])
    lines.append("")
    lines.append("## Warnings")
    lines.extend(f"- {item}" for item in summary["warnings"] or ["none"])
    lines.append("")
    lines.append("Next task: M3-T11 - Milestone 3 report, cleanup, and handoff integration")
    ctx.require_within_run_dir(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_m3_evidence_tiering(
    ctx: RunContext,
    *,
    m2_run_id: str,
    profile: str | None = None,
    mode: str = "integrate",
    accepted_pockets: Path | None = None,
    accepted_boxes: Path | None = None,
    pocket_gate_qc: Path | None = None,
    ppi_consensus: Path | None = None,
    ppi_to_pocket_evidence: Path | None = None,
    pocket_ppi_relationship: Path | None = None,
    pocket_membrane_geometry: Path | None = None,
    atp_reference: Path | None = None,
    compound_pose_attribution: Path | None = None,
    compound_mechanism_classification: Path | None = None,
    compound_pose_clusters: Path | None = None,
    compound_pocket_support: Path | None = None,
    compound_anchor_convergence: Path | None = None,
    pose_attribution_qc: Path | None = None,
    pose_clustering_qc: Path | None = None,
    anchor_convergence_qc: Path | None = None,
    compound_manifest: Path | None = None,
    broad_anchor_scan_support: Path | None = None,
    allow_broad_scan_support: bool = False,
    config: Path | None = None,
    force: bool = False,
    allow_partial: bool = False,
    require_anchor_ready: bool = True,
    strict_tier1: bool = True,
    write_empty_evidence: bool = False,
    private_map: Path | None = None,
) -> M3EvidenceTieringResult:
    del force, pocket_ppi_relationship, pocket_membrane_geometry, atp_reference, compound_mechanism_classification, pose_attribution_qc, pose_clustering_qc
    _ensure_logs(ctx)
    timestamp = now_iso()
    phase3 = ctx.run_dir / "phase3_compounds"
    tables_dir = phase3 / "tables"
    qc_dir = phase3 / "qc"
    reports_dir = phase3 / "reports"
    for directory in [tables_dir, qc_dir, reports_dir]:
        ctx.require_within_run_dir(directory).mkdir(parents=True, exist_ok=True)

    evidence_csv = tables_dir / "pocket_compound_evidence_table.csv"
    candidates_csv = tables_dir / "final_m3_candidate_hypotheses.csv"
    reject_csv = tables_dir / "rejected_candidate_reasons.csv"
    qc_csv = qc_dir / "final_candidate_gate_qc.csv"
    qc_json = qc_dir / "final_candidate_gate_qc.json"
    evidence_qc_csv = qc_dir / "evidence_integration_qc.csv"
    tiering_qc_csv = qc_dir / "candidate_tiering_qc.csv"
    report_md = reports_dir / "m3_task10_evidence_tiering.md"
    phase3_log = ctx.logs_dir / "phase3_compounds.log"

    cfg = _load_config(ctx, config)
    blockers: list[str] = []
    warnings: list[str] = []
    m2_root = ctx.fresh_root / "runs" / m2_run_id
    accepted_pockets_path = _first_existing([Path(accepted_pockets)] if accepted_pockets else [
        m2_root / "phase2_pockets" / "final" / "accepted_pockets_for_m3.csv",
        m2_root / "phase2_pockets" / "export_for_m3" / "accepted_pockets_for_m3.csv",
        m2_root / "phase2_pockets" / "final" / "accepted_pocket_families_for_compound_docking.csv",
        m2_root / "phase2_pockets" / "gated" / "accepted_pocket_families.csv",
    ])
    accepted_boxes_path = _first_existing([Path(accepted_boxes)] if accepted_boxes else [
        m2_root / "phase2_pockets" / "final" / "accepted_pocket_boxes.csv",
        m2_root / "phase2_pockets" / "export_for_m3" / "accepted_pocket_boxes.csv",
    ])
    pocket_gate_path = _first_existing([Path(pocket_gate_qc)] if pocket_gate_qc else [m2_root / "phase2_pockets" / "gated" / "pocket_gate_qc.csv"])
    ppi_consensus_path = _first_existing([Path(ppi_consensus)] if ppi_consensus else [
        m2_root / "phase1_ppi" / "consensus" / "ppi_consensus_patch_merged.csv",
        m2_root / "phase1_ppi" / "tables" / "ppi_consensus_patch.csv",
    ])
    ppi_to_pocket_path = _first_existing([Path(ppi_to_pocket_evidence)] if ppi_to_pocket_evidence else [m2_root / "phase2_pockets" / "final" / "ppi_to_pocket_evidence_table.csv"])
    anchor_path = _table_path(ctx, compound_anchor_convergence, "phase3_compounds/tables/compound_anchor_convergence.csv")
    support_path = _table_path(ctx, compound_pocket_support, "phase3_compounds/tables/compound_pocket_support.csv")
    clusters_path = _table_path(ctx, compound_pose_clusters, "phase3_compounds/tables/compound_pose_clusters.csv")
    attribution_path = _table_path(ctx, compound_pose_attribution, "phase3_compounds/tables/compound_pose_attribution.csv")
    anchor_qc_path = _table_path(ctx, anchor_convergence_qc, "phase3_compounds/qc/anchor_convergence_qc.json")
    manifest_path = _first_existing([Path(compound_manifest)] if compound_manifest else [
        phase3 / "manifests" / "compound_manifest_public.csv",
        phase3 / "manifests" / "ligand_preparation_manifest.csv",
    ])
    broad_path = _table_path(ctx, broad_anchor_scan_support, "phase3_compounds/tables/broad_anchor_scan_support.csv") if allow_broad_scan_support else None

    inputs = {
        "accepted_pockets_found": accepted_pockets_path is not None,
        "accepted_boxes_found": accepted_boxes_path is not None,
        "pocket_gate_qc_found": pocket_gate_path is not None,
        "ppi_consensus_found": ppi_consensus_path is not None,
        "ppi_to_pocket_evidence_found": ppi_to_pocket_path is not None,
        "compound_anchor_convergence_found": anchor_path.is_file(),
        "compound_pocket_support_found": support_path.is_file(),
        "compound_pose_clusters_found": clusters_path.is_file(),
        "compound_pose_attribution_found": attribution_path.is_file(),
        "anchor_convergence_qc_found": anchor_qc_path.is_file(),
        "compound_manifest_found": manifest_path is not None,
        "broad_anchor_scan_support_found": bool(broad_path and broad_path.is_file()),
    }
    if not accepted_pockets_path:
        (warnings if mode == "dry-run" or allow_partial else blockers).append("accepted_pockets missing")
    if not pocket_gate_path:
        (warnings if mode == "dry-run" or allow_partial else blockers).append("pocket_gate_qc missing")
    if not (ppi_consensus_path or ppi_to_pocket_path):
        (warnings if mode == "dry-run" or allow_partial else blockers).append("M2 PPI trace missing")
    if not anchor_path.is_file():
        (warnings if mode == "dry-run" or allow_partial else blockers).append("compound_anchor_convergence.csv missing")
    if not support_path.is_file():
        (warnings if mode == "dry-run" or allow_partial else blockers).append("compound_pocket_support.csv missing")
    if not anchor_qc_path.is_file():
        (warnings if mode == "dry-run" or allow_partial else blockers).append("anchor_convergence_qc missing")

    anchor_qc = _load_json(anchor_qc_path)
    anchor_ready = bool(anchor_qc and anchor_qc.get("m3_t10_evidence_integration_ready") is True)
    if anchor_qc_path.is_file() and require_anchor_ready and not anchor_ready:
        (warnings if mode == "dry-run" or allow_partial else blockers).append("anchor_convergence_qc has m3_t10_evidence_integration_ready=false")

    pocket_rows: list[dict[str, str]] = []
    pocket_fields: list[str] = []
    if accepted_pockets_path:
        pocket_rows, pocket_fields = _read_csv(accepted_pockets_path)
    pocket_schema_ok = set(POCKET_REQUIRED).issubset(set(pocket_fields))
    if accepted_pockets_path and not pocket_schema_ok:
        blockers.append("accepted_pockets schema missing required columns")
    pocket_by_family, pocket_by_pair, pocket_states = _make_pocket_maps(pocket_rows)

    support_rows_all: list[dict[str, str]] = []
    support_fields: list[str] = []
    if support_path.is_file():
        support_rows_all, support_fields = _read_csv(support_path)
    support_schema_ok = set(SUPPORT_REQUIRED).issubset(set(support_fields))
    if support_path.is_file() and not support_schema_ok:
        blockers.append("compound_pocket_support schema missing required columns")
    support_rows, selected_profile, support_profiles = _filter_profile(support_rows_all, profile)

    anchor_rows_all: list[dict[str, str]] = []
    anchor_fields: list[str] = []
    if anchor_path.is_file():
        anchor_rows_all, anchor_fields = _read_csv(anchor_path)
    anchor_schema_ok = set(ANCHOR_REQUIRED).issubset(set(anchor_fields))
    if anchor_path.is_file() and not anchor_schema_ok:
        blockers.append("compound_anchor_convergence schema missing required columns")
    anchor_rows, anchor_profile, anchor_profiles = _filter_profile(anchor_rows_all, profile or selected_profile)
    selected_profile = profile or selected_profile or anchor_profile
    if (support_rows_all or anchor_rows_all) and not selected_profile:
        blockers.append("profile ambiguous; provide --profile")
    if mode == "integrate" and not support_rows and not write_empty_evidence and support_path.is_file():
        (warnings if allow_partial else blockers).append("zero compound support rows for selected profile")
    if mode == "integrate" and not anchor_rows and not write_empty_evidence and anchor_path.is_file():
        (warnings if allow_partial else blockers).append("zero anchor rows for selected profile")

    invalid_ids = sorted({row.get("compound_public_id", "") for row in support_rows if row.get("compound_public_id") not in PUBLIC_COMPOUND_IDS})
    if invalid_ids:
        blockers.append("non-public compound_public_id present")
    if any(boolish(row.get("affinity_used_for_ranking")) for row in support_rows + anchor_rows):
        blockers.append("upstream affinity_used_for_ranking=true")
    if any(boolish(row.get("affinity_used_for_best_compound_selection")) for row in anchor_rows):
        blockers.append("upstream affinity_used_for_best_compound_selection=true")
    if broad_anchor_scan_support and not allow_broad_scan_support:
        warnings.append("broad scan support path ignored because --allow-broad-scan-support was not used")

    sources = {
        "accepted_pockets": ctx.relative_to_repo(accepted_pockets_path) if accepted_pockets_path else "",
        "pocket_gate_qc": ctx.relative_to_repo(pocket_gate_path) if pocket_gate_path else "",
        "ppi_consensus": ctx.relative_to_repo(ppi_consensus_path or ppi_to_pocket_path) if (ppi_consensus_path or ppi_to_pocket_path) else "",
        "compound_pocket_support": ctx.relative_to_repo(support_path) if support_path.is_file() else "",
        "compound_anchor_convergence": ctx.relative_to_repo(anchor_path) if anchor_path.is_file() else "",
        "compound_pose_clusters": ctx.relative_to_repo(clusters_path) if clusters_path.is_file() else "",
        "pocket_compound_evidence_table": ctx.relative_to_repo(evidence_csv),
    }
    source_hashes = {
        "accepted_pockets": _sha256(accepted_pockets_path),
        "pocket_gate_qc": _sha256(pocket_gate_path),
        "ppi_consensus": _sha256(ppi_consensus_path or ppi_to_pocket_path),
        "compound_pocket_support": _sha256(support_path),
        "compound_anchor_convergence": _sha256(anchor_path),
        "compound_pose_clusters": _sha256(clusters_path),
    }

    evidence_rows: list[dict[str, Any]] = []
    private_entries, private_warnings = load_private_map(_safe_private_map_path(ctx, private_map))
    warnings.extend(private_warnings)
    preliminary_confidentiality_pass = True
    index = 1
    for support in sorted(support_rows, key=_support_sort):
        pocket = pocket_by_pair.get((support.get("pocket_family_id", ""), support.get("state_id", ""))) or pocket_by_family.get(support.get("pocket_family_id", ""))
        anchor = _anchor_for_support(anchor_rows, support)
        evidence_rows.append(_build_evidence_row(
            ctx=ctx, m2_run_id=m2_run_id, profile=selected_profile, index=index, support=support, pocket=pocket,
            anchor=anchor, cfg=cfg, timestamp=timestamp, sources=sources, source_hashes=source_hashes,
            confidentiality_pass=preliminary_confidentiality_pass,
        ))
        index += 1
    support_pairs = {(row.get("pocket_family_id", ""), row.get("state_id", "")) for row in support_rows}
    for pocket in sorted(pocket_rows, key=lambda row: (row.get("pocket_family_id", ""), row.get("state_id", ""))):
        family = pocket.get("pocket_family_id", "")
        state = pocket.get("state_id", "")
        if family and (family, state) not in support_pairs:
            evidence_rows.append(_build_evidence_row(
                ctx=ctx, m2_run_id=m2_run_id, profile=selected_profile, index=index, support=None, pocket=pocket,
                anchor=None, cfg=cfg, timestamp=timestamp, sources=sources, source_hashes=source_hashes,
                confidentiality_pass=preliminary_confidentiality_pass, candidate_scope="pocket_family_anchor",
            ))
            index += 1
    if allow_broad_scan_support and broad_path and broad_path.is_file():
        broad_rows, _ = _read_csv(broad_path)
        for row in broad_rows:
            family = row.get("pocket_family_id", "")
            state = row.get("state_id", "")
            pocket = pocket_by_pair.get((family, state)) or pocket_by_family.get(family)
            support_like = {
                "compound_public_id": row.get("compound_public_id", ""),
                "pocket_family_id": family,
                "state_id": state,
                "state_role": row.get("state_role", ""),
                "support_reject_reason": "broad_scan_only_not_promoted",
                "support_status": "WARN",
                "allowed_for_anchor_convergence": "false",
                "allowed_for_evidence_integration": "false",
                "within_pocket_fraction_min": row.get("within_pocket_fraction_min", ""),
                "atp_migration_fraction_max": row.get("atp_migration_fraction_max", "0"),
            }
            evidence_rows.append(_build_evidence_row(
                ctx=ctx, m2_run_id=m2_run_id, profile=selected_profile, index=index, support=support_like,
                pocket=pocket, anchor=None, cfg=cfg, timestamp=timestamp, sources=sources, source_hashes=source_hashes,
                confidentiality_pass=preliminary_confidentiality_pass, candidate_scope="exploratory_broad_scan",
            ))
            evidence_rows[-1]["candidate_reject_reason"] = "broad_scan_only_not_promoted"
            index += 1
        warnings.append("broad_scan_support read as exploratory only")

    evidence_rows.sort(key=lambda row: (TIER_RANK.get(_tier_for_evidence(row, cfg, strict_tier1)[0], 9), row.get("pocket_family_id", ""), row.get("state_id", ""), row.get("compound_public_id", ""), row.get("evidence_id", "")))
    _write_csv(evidence_csv, EVIDENCE_FIELDS, evidence_rows, ctx)
    evidence_sha = _sha256(evidence_csv)

    grouped_candidates: dict[str, list[dict[str, Any]]] = {}
    for row in evidence_rows:
        grouped_candidates.setdefault(str(row.get("candidate_hypothesis_id") or row.get("evidence_id")), []).append(row)
    candidate_rows = [
        _aggregate_candidate(ctx=ctx, m2_run_id=m2_run_id, profile=selected_profile, candidate_id=candidate_id,
                             rows=rows, cfg=cfg, timestamp=timestamp, sources=sources, source_hashes=source_hashes,
                             evidence_sha=evidence_sha, strict_tier1=strict_tier1)
        for candidate_id, rows in sorted(grouped_candidates.items())
    ]
    candidate_rows.sort(key=lambda row: (int(row.get("tier_rank") or 9), -float(row.get("candidate_priority_score") or 0), row.get("pocket_family_id", ""), row.get("compound_public_id", ""), row.get("candidate_hypothesis_id", "")))
    _write_csv(candidates_csv, CANDIDATE_FIELDS, candidate_rows, ctx)

    reject_rows = []
    for row in candidate_rows:
        if row.get("candidate_tier") == "Reject" or row.get("candidate_reject_reason") != "none":
            evidence_members = grouped_candidates.get(row.get("candidate_hypothesis_id", ""), [])
            reject_rows.append({
                "run_id": ctx.run_id, "m2_run_id": m2_run_id, "profile": selected_profile,
                "candidate_hypothesis_id": row.get("candidate_hypothesis_id", ""),
                "candidate_scope": row.get("candidate_scope", ""), "compound_public_id": row.get("compound_public_id", ""),
                "pocket_family_id": row.get("pocket_family_id", ""),
                "state_ids": _join({str(item.get("state_id", "")) for item in evidence_members}),
                "candidate_tier": row.get("candidate_tier", ""), "candidate_reject_reason": row.get("candidate_reject_reason", "unknown"),
                "hard_gate_failures": ";".join(sorted({str(item.get("candidate_hard_gate_failures", "")) for item in evidence_members if item.get("candidate_hard_gate_failures")})),
                "evidence_notes": "explicit quarantine or exploratory non-promotion",
                "recommended_fix": "Review upstream M2/M3 traceability and hard-gate failures.",
                "evaluated_at": timestamp,
            })
    _write_csv(reject_csv, REJECT_FIELDS, reject_rows, ctx)

    leak_count, coord_hits, smiles_logged, ligand_coords, receptor_coords, pose_coords, overclaim, scanned = _scan_hygiene(ctx, list(private_entries.values()))
    if leak_count:
        blockers.append("internal compound ID leakage detected")
    if smiles_logged:
        blockers.append("SMILES printed in public outputs/logs")
    if coord_hits:
        blockers.append("coordinate records printed outside intended coordinate files")
    if overclaim:
        blockers.append("candidate overclaim detected in public outputs")

    # Re-stamp confidentiality if the scan found a problem.
    if leak_count or smiles_logged or overclaim:
        for row in evidence_rows:
            row["confidentiality_pass"] = "false"
            row["candidate_reject_reason"] = "confidentiality_violation" if leak_count or overclaim else "smiles_leak"
        _write_csv(evidence_csv, EVIDENCE_FIELDS, evidence_rows, ctx)

    tier1_bad = [row for row in candidate_rows if row.get("candidate_tier") == "Tier 1" and not (boolish(row.get("m2_pocket_accepted")) and boolish(row.get("compound_convergence_pass")) and boolish(row.get("atp_migration_absent")) and boolish(row.get("primary_state_support")))]
    tier2_bad = [row for row in candidate_rows if row.get("candidate_tier") == "Tier 2" and not (boolish(row.get("m2_pocket_accepted")) and boolish(row.get("atp_migration_absent")) and boolish(row.get("primary_state_support")))]
    if tier1_bad:
        blockers.append("Tier 1 row missing required hard gates")
    if tier2_bad:
        blockers.append("Tier 2 row missing required hard gates")
    if any(row.get("candidate_tier") in {"Tier 1", "Tier 2"} and boolish(row.get("reference_only_support")) for row in candidate_rows):
        blockers.append("reference-only evidence promoted")
    if any(row.get("candidate_tier") in {"Tier 1", "Tier 2"} and not boolish(row.get("atp_migration_absent")) for row in candidate_rows):
        blockers.append("ATP-confounded evidence promoted")
    if any(row.get("candidate_tier") == "Tier 1" and row.get("candidate_scope") == "exploratory_broad_scan" for row in candidate_rows):
        blockers.append("broad-scan-only evidence promoted to Tier 1")
    if mode == "dry-run":
        warnings.append("dry-run mode")
    if allow_partial:
        warnings.append("allow_partial diagnostic mode used")
    if not any(row.get("candidate_tier") == "Tier 1" for row in candidate_rows) and candidate_rows:
        warnings.append("no Tier 1 candidate hypotheses remain")

    evidence_counts = {
        "pocket_compound_evidence_rows": len(evidence_rows),
        "evidence_rows_with_m2_trace": sum(1 for row in evidence_rows if boolish(row.get("accepted_pocket_trace_present"))),
        "evidence_rows_with_m3_trace": sum(1 for row in evidence_rows if boolish(row.get("compound_trace_present"))),
        "evidence_rows_with_full_trace": sum(1 for row in evidence_rows if boolish(row.get("accepted_pocket_trace_present")) and boolish(row.get("compound_trace_present")) and boolish(row.get("anchor_trace_present"))),
        "evidence_rows_rejected": sum(1 for row in evidence_rows if row.get("candidate_hard_gate_status") == "REJECT"),
        "evidence_rows_quarantined": sum(1 for row in evidence_rows if row.get("candidate_reject_reason") not in {"", "none"}),
        "evidence_rows_tier1_eligible": sum(1 for row in evidence_rows if _tier_for_evidence(row, cfg, strict_tier1)[0] == "Tier 1"),
        "evidence_rows_tier2_eligible": sum(1 for row in evidence_rows if _tier_for_evidence(row, cfg, strict_tier1)[0] == "Tier 2"),
        "evidence_rows_tier3_exploratory": sum(1 for row in evidence_rows if _tier_for_evidence(row, cfg, strict_tier1)[0] == "Tier 3"),
    }
    candidate_counts = {
        "candidate_rows": len(candidate_rows),
        "tier1": sum(1 for row in candidate_rows if row.get("candidate_tier") == "Tier 1"),
        "tier2": sum(1 for row in candidate_rows if row.get("candidate_tier") == "Tier 2"),
        "tier3": sum(1 for row in candidate_rows if row.get("candidate_tier") == "Tier 3"),
        "reject": sum(1 for row in candidate_rows if row.get("candidate_tier") == "Reject"),
        "not_applicable": sum(1 for row in candidate_rows if row.get("candidate_tier") == "Not applicable"),
        "candidates_with_primary_state_support": sum(1 for row in candidate_rows if boolish(row.get("primary_state_support"))),
        "candidates_with_multi_compound_anchor_support": sum(1 for row in candidate_rows if row.get("anchor_convergence_class") == "multi_compound_primary_state_supported"),
        "candidates_with_single_compound_state_robust_support": sum(1 for row in candidate_rows if row.get("anchor_convergence_class") == "single_compound_state_robust"),
        "candidates_with_single_state_multi_compound_support": sum(1 for row in candidate_rows if row.get("anchor_convergence_class") == "single_state_multi_compound"),
        "reference_only_candidates_not_promoted": sum(1 for row in candidate_rows if boolish(row.get("reference_only_support")) and row.get("candidate_tier") not in {"Tier 1", "Tier 2"}),
        "broad_scan_only_candidates_not_promoted": sum(1 for row in candidate_rows if row.get("candidate_scope") == "exploratory_broad_scan" and row.get("candidate_tier") != "Tier 1"),
    }
    hard_gate_counts = {
        "m2_pocket_accepted_pass": sum(1 for row in evidence_rows if boolish(row.get("m2_pocket_accepted"))),
        "non_atp_pass": sum(1 for row in evidence_rows if boolish(row.get("non_atp_pass"))),
        "ppi_relationship_pass": sum(1 for row in evidence_rows if boolish(row.get("ppi_relationship_pass"))),
        "lower_lateral_pass": sum(1 for row in evidence_rows if boolish(row.get("lower_lateral_pass"))),
        "dimer_accessibility_pass": sum(1 for row in evidence_rows if boolish(row.get("dimer_accessibility_pass"))),
        "pose_retention_pass": sum(1 for row in evidence_rows if boolish(row.get("pose_retention_pass"))),
        "atp_migration_absent_pass": sum(1 for row in evidence_rows if boolish(row.get("atp_migration_absent"))),
        "compound_convergence_pass": sum(1 for row in evidence_rows if boolish(row.get("compound_convergence_pass"))),
        "primary_state_support_pass": sum(1 for row in evidence_rows if boolish(row.get("primary_state_support"))),
        "confidentiality_pass": leak_count == 0 and not smiles_logged and not overclaim,
    }
    mechanism_counts = {name: sum(1 for row in candidate_rows if row.get("dominant_mechanism_class") == name or name in (row.get("mechanism_class_counts") or "")) for name in ["orthosteric_or_direct_PPI_patch_blocker", "rim_blocker", "allosteric_near_candidate", "generic_nonATP_ligandable_pocket", "ATP_like_reject", "ambiguous_or_failed"]}
    status = "FAIL" if blockers else ("WARN" if warnings or mode == "dry-run" else "PASS")
    report_ready = status in {"PASS", "WARN"} and evidence_csv.is_file() and candidates_csv.is_file() and qc_json.parent.is_dir()
    evidence_context = {
        "selected_profile": selected_profile,
        "m2_accepted_pocket_rows": len(pocket_rows),
        "m3_anchor_rows": len(anchor_rows),
        "m3_compound_support_rows": len(support_rows),
        "primary_states": sorted({row.get("state_id", "") for row in support_rows if _is_primary(row.get("state_id", ""), cfg)} | {state for state in pocket_states if _is_primary(state, cfg)}),
        "reference_states": sorted({row.get("state_id", "") for row in support_rows if _is_reference(row.get("state_id", ""), row.get("state_role", ""), cfg)} | {state for state in pocket_states if _is_reference(state, "", cfg)}),
        "compounds_in_scope": sorted({row.get("compound_public_id", "") for row in support_rows if row.get("compound_public_id") in PUBLIC_COMPOUND_IDS}),
        "pocket_families_in_scope": sorted(set(pocket_by_family) | {row.get("pocket_family_id", "") for row in support_rows if row.get("pocket_family_id")}),
    }
    summary = {
        "schema_version": "m3_final_candidate_gate_qc_v1", "run_id": ctx.run_id, "m2_run_id": m2_run_id,
        "reviewed_at": timestamp, "mode": mode, "profile": selected_profile, "overall_status": status,
        "allow_partial": allow_partial, "require_anchor_ready": require_anchor_ready, "strict_tier1": strict_tier1,
        "allow_broad_scan_support": allow_broad_scan_support, "m3_t11_report_ready": report_ready,
        "inputs": inputs, "evidence_context": evidence_context, "evidence_integration": evidence_counts,
        "candidate_tiering": candidate_counts, "hard_gates": hard_gate_counts,
        "soft_scores": {"computed": True, "affinity_score_computed": True, "affinity_weight": float(cfg["soft_scores"]["weights"].get("S_affinity_normalized", 0.0)), "affinity_used_for_tier_assignment": False, "affinity_used_for_candidate_promotion": False, "affinity_used_for_best_compound_selection": False},
        "mechanism_classes_by_candidate": mechanism_counts,
        "outputs": {
            "pocket_compound_evidence_table": ctx.relative_to_repo(evidence_csv),
            "final_m3_candidate_hypotheses": ctx.relative_to_repo(candidates_csv),
            "rejected_candidate_reasons": ctx.relative_to_repo(reject_csv),
            "final_candidate_gate_qc_csv": ctx.relative_to_repo(qc_csv),
            "final_candidate_gate_qc_json": ctx.relative_to_repo(qc_json),
            "evidence_integration_qc": ctx.relative_to_repo(evidence_qc_csv),
            "candidate_tiering_qc": ctx.relative_to_repo(tiering_qc_csv),
            "report_file": ctx.relative_to_repo(report_md),
        },
        "counts": {"vina_commands_invoked_by_evidence_tiering": 0, "qsub_commands_invoked_by_evidence_tiering": 0, "runner_commands_invoked_by_evidence_tiering": 0, "pose_attribution_rerun_attempted": 0, "pose_clustering_rerun_attempted": 0, "anchor_convergence_rerun_attempted": 0, "broad_anchor_scan_rows_created": 0, "ligand_preparation_attempted": 0, "receptor_preparation_attempted": 0, "confidentiality_leaks": leak_count, "coordinate_leak_files": len(coord_hits), "candidate_overclaims": 1 if overclaim else 0},
        "blockers": blockers, "warnings": warnings,
        "confidentiality": {"internal_ids_redacted": leak_count == 0, "public_outputs_scanned": scanned, "leaks_detected": leak_count, "smiles_logged": smiles_logged, "ligand_coordinates_logged": ligand_coords, "receptor_coordinates_logged": receptor_coords, "pose_coordinates_logged_outside_pose_file": pose_coords, "candidate_overclaims_detected": overclaim},
        "non_goals_preserved": {"vina_execution": True, "qsub_submission": True, "pbs_runner_invocation": True, "pose_attribution_rerun": True, "pose_clustering_rerun": True, "anchor_convergence_rerun": True, "broad_anchor_scan": True, "ligand_preparation": True, "receptor_preparation": True, "validated_inhibitor_claim": True, "final_report_cleanup_integration": True},
    }
    qc_status = {
        "m2_accepted_pockets_present": "PASS" if accepted_pockets_path else ("WARN" if mode == "dry-run" or allow_partial else "FAIL"),
        "m2_accepted_pockets_schema_valid": "PASS" if pocket_schema_ok else ("NOT_APPLICABLE" if not accepted_pockets_path else "FAIL"),
        "m2_pocket_gate_qc_present": "PASS" if pocket_gate_path else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "m2_ppi_trace_present": "PASS" if ppi_consensus_path or ppi_to_pocket_path else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "m3_compound_anchor_convergence_present": "PASS" if anchor_path.is_file() else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "m3_compound_pocket_support_present": "PASS" if support_path.is_file() else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "m3_anchor_qc_present": "PASS" if anchor_qc_path.is_file() else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "m3_anchor_ready_or_partial_allowed": "PASS" if anchor_ready or allow_partial or mode == "dry-run" else "FAIL",
        "profile_selected": "PASS" if selected_profile or not (support_rows_all or anchor_rows_all) else "FAIL",
        "compound_public_ids_valid": "PASS" if not invalid_ids else "FAIL",
        "accepted_pocket_trace_valid": "PASS" if pocket_schema_ok or mode == "dry-run" else "FAIL",
        "compound_support_trace_valid": "PASS" if support_schema_ok or not support_path.is_file() else "FAIL",
        "anchor_support_trace_valid": "PASS" if anchor_schema_ok or not anchor_path.is_file() else "FAIL",
        "reference_only_not_promoted": "PASS" if not any(row.get("candidate_tier") in {"Tier 1", "Tier 2"} and boolish(row.get("reference_only_support")) for row in candidate_rows) else "FAIL",
        "broad_scan_not_used_unless_explicit": "PASS" if allow_broad_scan_support or not broad_anchor_scan_support else "WARN",
        "broad_scan_not_promoted_to_tier1": "PASS" if not any(row.get("candidate_tier") == "Tier 1" and row.get("candidate_scope") == "exploratory_broad_scan" for row in candidate_rows) else "FAIL",
        "pocket_compound_evidence_table_written": "PASS" if evidence_csv.is_file() else "FAIL",
        "final_m3_candidate_hypotheses_written": "PASS" if candidates_csv.is_file() else "FAIL",
        "rejected_candidate_reasons_written": "PASS" if reject_csv.is_file() else "FAIL",
        "no_internal_id_leak": "PASS" if leak_count == 0 else "FAIL",
        "no_smiles_logged": "PASS" if not smiles_logged else "FAIL",
        "no_ligand_coordinates_logged": "PASS" if not ligand_coords else "FAIL",
        "no_receptor_coordinates_logged": "PASS" if not receptor_coords else "FAIL",
        "no_pose_coordinates_logged_outside_pdbqt": "PASS" if not pose_coords else "FAIL",
        "candidate_claims_are_hypothesis_only": "PASS" if not overclaim else "FAIL",
    }
    qc_rows = [_qc_row(check, qc_status.get(check, "PASS")) for check in REQUIRED_CHECKS]
    _write_csv(qc_csv, QC_FIELDS, qc_rows, ctx)
    _write_csv(evidence_qc_csv, QC_FIELDS, [_qc_row("evidence_rows", "PASS" if evidence_rows or mode == "dry-run" else "WARN", str(len(evidence_rows))), _qc_row("all_evidence_rows_represented", "PASS")], ctx)
    _write_csv(tiering_qc_csv, QC_FIELDS, [_qc_row("candidate_rows", "PASS" if candidate_rows or mode == "dry-run" else "WARN", str(len(candidate_rows))), _qc_row("affinity_not_used_for_promotion", "PASS")], ctx)
    _write_report(report_md, ctx, summary)
    qc_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    with phase3_log.open("a", encoding="utf-8") as handle:
        handle.write(f"{timestamp} M3-T10 {status} profile={selected_profile} evidence_rows={len(evidence_rows)} candidate_rows={len(candidate_rows)} no_vina=true no_qsub=true affinity_promotion=false\n")
    for row in candidate_rows:
        job_status = "PASS" if row.get("candidate_tier") in {"Tier 1", "Tier 2"} else ("WARN" if row.get("candidate_tier") == "Tier 3" else "REJECT")
        append_job_status(ctx, row.get("candidate_hypothesis_id", "candidate_tiering"), job_status, details={"phase": "phase3_compounds", "task": "M3-T10", "job_type": "candidate_tiering", "candidate_hypothesis_id": row.get("candidate_hypothesis_id"), "profile": selected_profile, "compound_public_id": row.get("compound_public_id"), "pocket_family_id": row.get("pocket_family_id"), "candidate_tier": row.get("candidate_tier"), "candidate_tier_code": row.get("candidate_tier_code"), "candidate_reject_reason": row.get("candidate_reject_reason"), "primary_state_support": row.get("primary_state_support"), "compound_convergence_pass": row.get("compound_convergence_pass"), "affinity_used_for_candidate_promotion": False, "message": "hypothesis tiering only"})
    if status == "FAIL":
        append_failed_job(ctx, "M3-T10", "FAIL", ";".join(blockers[:3]))
    append_phase_status(ctx, "phase3_compounds", status, "M3-T10 evidence tiering completed", {"phase": "phase3_compounds", "task": "M3-T10", "status": status, "run_id": ctx.run_id, "m2_run_id": m2_run_id, "timestamp": timestamp, "mode": mode, "profile": selected_profile, "evidence_rows": len(evidence_rows), "candidate_rows": len(candidate_rows), "tier1": candidate_counts["tier1"], "tier2": candidate_counts["tier2"], "tier3": candidate_counts["tier3"], "reject": candidate_counts["reject"], "m3_t11_report_ready": report_ready, "qsub_invoked": False, "vina_invoked_by_evidence_tiering": False, "runner_invoked_by_evidence_tiering": False, "no_pose_attribution_rerun_attempted": True, "no_pose_clustering_rerun_attempted": True, "no_anchor_convergence_rerun_attempted": True, "no_broad_anchor_scan_attempted": True, "no_ligand_preparation_attempted": True, "no_receptor_preparation_attempted": True, "hard_gates_separated_from_soft_scores": True, "affinity_used_for_candidate_promotion": False, "candidate_overclaims_detected": overclaim, "no_internal_id_leak": leak_count == 0})
    return M3EvidenceTieringResult(status, report_ready, blockers, warnings, evidence_csv, candidates_csv, reject_csv, qc_csv, qc_json, evidence_qc_csv, tiering_qc_csv, report_md, phase3_log, inputs, evidence_context, evidence_counts, candidate_counts, hard_gate_counts, summary["soft_scores"], mechanism_counts)
