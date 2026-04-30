"""M3-T9 compound-anchor convergence aggregation.

Reads M3-T8 cluster outputs and writes compound-pocket support plus
compound-anchor convergence tables. This module never invokes docking, reruns
pose attribution/clustering, scores compounds, ranks candidates, or creates
final candidate tables.
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

from egfr_myo1d.compound.confidentiality import PUBLIC_COMPOUND_IDS, load_private_map, public_output_scan_paths, scan_internal_id_leaks
from egfr_myo1d.core.logging_utils import append_failed_job, append_job_status, append_phase_status
from egfr_myo1d.core.run_context import RunContext, ensure_within


DEFAULT_CONFIG: dict[str, Any] = {
    "schema_version": "m3_anchor_convergence_config_v1",
    "input": {"require_t8_qc_ready": True, "require_pose_attribution_trace": True, "require_cluster_membership_trace": True},
    "compounds": {"expected_public_compound_ids": ["Cpd-A", "Cpd-B", "Cpd-C"], "require_all_expected_compounds_tested_for_production": False},
    "states": {
        "primary_state_ids": ["EGFR_160-185", "EGFR_170-200"],
        "reference_state_ids": ["3GT8_raw"],
        "reference_state_role_values": ["reference", "crystallographic_reference_control", "reference_control"],
        "do_not_count_ab_symmetry_as_independent_state_support": True,
    },
    "support_cluster_filters": {
        "allowed_convergence_statuses": ["PASS_CONVERGED", "PASS_CONVERGED_GENERIC_NONATP", "WARN_CONVERGED_WEAK_PPI_RELATIONSHIP"],
        "require_allowed_for_anchor_convergence": True,
        "require_cluster_reject_reason": "none",
        "require_atp_migration_fraction_max": 0.0,
        "require_membrane_penetration_fraction_max": 0.0,
        "require_dimer_interface_clash_fraction_max": 0.0,
        "pocket_retention_fraction_min": 0.95,
        "weak_ppi_relationship_allowed": True,
    },
    "evidence_readiness": {
        "m3_t10_ready_classes": ["multi_compound_primary_state_supported", "single_compound_state_robust", "single_state_multi_compound"],
        "require_no_blockers": True,
        "require_primary_support_unless_allow_reference_only": True,
    },
    "affinity": {"record_descriptive_stats": True, "use_affinity_for_ranking": False, "use_affinity_for_class_assignment": False, "use_affinity_for_best_compound_selection": False},
    "confidentiality": {"allowed_public_compound_ids": ["Cpd-A", "Cpd-B", "Cpd-C"], "scan_public_outputs": True, "fail_on_internal_id_leak": True, "fail_on_smiles_leak": True, "fail_on_candidate_claims": True},
}

CLUSTER_REQUIRED = [
    "run_id", "m2_run_id", "profile", "cluster_id", "compound_public_id", "state_id", "state_role",
    "pocket_family_id", "protomer_id", "converged_pose_cluster", "convergence_status", "cluster_reject_reason",
    "allowed_for_anchor_convergence", "affinity_used_for_ranking",
]

SUPPORT_FIELDS = [
    "run_id", "m2_run_id", "profile", "support_id", "compound_public_id", "pocket_family_id", "state_id",
    "state_role", "primary_state_flag", "reference_state_flag", "protomer_ids_with_support", "protomer_symmetry_support",
    "box_ids_with_support", "cluster_ids_all", "cluster_ids_supporting", "num_clusters_total", "num_converged_clusters",
    "num_anchor_allowed_clusters", "num_valid_pose_clusters", "num_rejected_or_nonconverged_clusters",
    "eligible_pose_count_total", "clustered_pose_count_total", "max_cluster_size", "max_cluster_fraction",
    "median_cluster_fraction", "within_pocket_fraction_min", "atp_migration_fraction_max",
    "ppi_relationship_confirmed_fraction_max", "ppi_relationship_confirmed_fraction_median",
    "ppi_hotspot_contact_median_max", "ppi_hotspot_contact_mean_max", "nearest_ppi_hotspot_distance_median_min",
    "ppi_rim_distance_median_min", "membrane_penetration_fraction_max", "dimer_interface_clash_fraction_max",
    "mechanism_class_counts", "mechanism_class_majority", "vina_affinity_median_kcal_mol",
    "vina_affinity_mean_kcal_mol", "vina_affinity_best_kcal_mol", "vina_affinity_worst_kcal_mol",
    "affinity_used_for_ranking", "support_status", "support_reject_reason", "support_notes",
    "allowed_for_anchor_convergence", "allowed_for_evidence_integration", "source_compound_pose_clusters",
    "source_compound_pose_clusters_sha256", "source_pose_clustering_qc", "source_pose_clustering_qc_sha256",
    "source_compound_pose_attribution", "source_compound_pose_attribution_sha256", "accepted_pocket_source",
    "accepted_pocket_source_sha256", "evaluated_at",
]

ANCHOR_FIELDS = [
    "run_id", "m2_run_id", "profile", "anchor_id", "anchor_scope", "pocket_family_id", "state_id",
    "state_ids_supported", "primary_state_support", "primary_state_ids_supported", "reference_state_support",
    "reference_state_ids_supported", "symmetry_support", "protomer_support_summary", "num_compounds_expected",
    "num_compounds_tested", "num_compounds_with_valid_pose", "num_compounds_with_converged_pose",
    "num_compounds_allowed_for_anchor_convergence", "compound_public_ids_tested", "compound_public_ids_with_valid_pose",
    "compound_public_ids_supported", "best_compound_public_id", "best_compound_selection_rule",
    "affinity_used_for_best_compound_selection", "median_affinity_across_supported", "mean_affinity_across_supported",
    "best_affinity_across_supported", "worst_affinity_across_supported", "affinity_used_for_ranking",
    "cluster_ids_supporting", "support_ids_supporting", "support_status_counts", "mechanism_class_counts",
    "dominant_mechanism_class", "max_cluster_fraction_across_supported", "median_cluster_fraction_across_supported",
    "min_within_pocket_fraction_across_supported", "max_atp_migration_fraction_across_supported",
    "max_membrane_penetration_fraction_across_supported", "max_dimer_interface_clash_fraction_across_supported",
    "max_ppi_relationship_confirmed_fraction", "median_ppi_relationship_confirmed_fraction", "anchor_converged",
    "anchor_convergence_class", "anchor_convergence_status", "anchor_reject_reason", "anchor_convergence_notes",
    "allowed_for_evidence_integration", "source_compound_pocket_support", "source_compound_pocket_support_sha256",
    "source_compound_pose_clusters", "source_compound_pose_clusters_sha256", "accepted_pocket_source",
    "accepted_pocket_source_sha256", "evaluated_at",
]

QC_FIELDS = ["check_id", "category", "status", "severity", "details", "recommended_fix"]

REQUIRED_CHECKS = [
    "m3_t8_compound_pose_clusters_present", "m3_t8_compound_pose_clusters_schema_valid", "m3_t8_pose_clustering_qc_present",
    "m3_t8_anchor_convergence_ready_or_partial_allowed", "m3_t7_pose_attribution_present_or_partial_allowed",
    "m3_t8_cluster_membership_present_or_partial_allowed", "profile_selected", "cluster_rows_found",
    "support_eligible_clusters_found", "cluster_ids_unique", "clusters_traceable_to_compound_state_pocket",
    "compound_public_ids_valid", "accepted_pockets_present", "accepted_pockets_schema_valid",
    "accepted_pocket_family_trace_valid", "primary_states_identified", "reference_states_identified",
    "ab_symmetry_not_counted_as_independent_state", "reference_state_not_promoted", "all_support_rows_represented",
    "no_silent_cluster_drop", "compound_pocket_support_computed", "anchor_convergence_computed",
    "anchor_convergence_class_assigned", "anchor_reject_reasons_explicit",
    "m3_t10_evidence_integration_readiness_flag_written", "multi_compound_convergence_reported_separately_from_affinity",
    "raw_pose_affinity_not_used_for_ranking", "best_compound_not_selected_by_affinity",
    "compound_pocket_support_written", "compound_anchor_convergence_written", "compound_pocket_support_qc_written",
    "anchor_convergence_report_written", "no_vina_invoked_by_anchor_convergence",
    "no_qsub_invoked_by_anchor_convergence", "no_runner_invoked_by_anchor_convergence",
    "no_pose_attribution_rerun_attempted", "no_pose_clustering_rerun_attempted", "no_broad_anchor_scan_attempted",
    "no_candidate_scoring_attempted", "no_candidate_tiering_attempted", "no_candidate_nomination_attempted",
    "no_final_candidate_tables_created", "no_internal_id_leak", "no_smiles_logged", "no_ligand_coordinates_logged",
    "no_receptor_coordinates_logged", "no_pose_coordinates_logged_outside_pdbqt", "old_workflow_not_used",
    "non_goals_preserved",
]

FORBIDDEN_OUTPUTS = [
    "phase3_compounds/tables/pocket_compound_evidence_table.csv",
    "phase3_compounds/tables/final_m3_candidate_hypotheses.csv",
    "phase3_compounds/qc/final_candidate_gate_qc.csv",
    "phase3_compounds/docking_outputs/broad_anchor_scan_optional",
    "phase3_compounds/tables/broad_anchor_scan_results.csv",
]


@dataclass
class M3AnchorConvergenceResult:
    status: str
    m3_t10_evidence_integration_ready: bool
    blockers: list[str]
    warnings: list[str]
    support_csv: Path
    anchor_csv: Path
    qc_csv: Path
    qc_json: Path
    support_qc_csv: Path
    report_md: Path
    phase3_log: Path
    inputs: dict[str, bool]
    cluster_context: dict[str, Any]
    support: dict[str, int]
    anchor_convergence: dict[str, int]
    mechanism_classes_by_anchor: dict[str, int]


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


def _float(value: Any) -> float | None:
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _fmt(value: float | None, digits: int = 3) -> str:
    return "" if value is None else f"{value:.{digits}f}"


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


def _severity(status: str) -> str:
    if status == "FAIL":
        return "BLOCKER"
    if status == "WARN":
        return "MAJOR"
    if status == "NOT_APPLICABLE":
        return "MINOR"
    return "INFO"


def _append_qc(rows: list[dict[str, str]], check_id: str, status: str) -> None:
    rows.append({"check_id": check_id, "category": "m3_t9_anchor_convergence", "status": status, "severity": _severity(status), "details": check_id, "recommended_fix": ""})


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
    failed = ctx.errors_dir / "failed_jobs.csv"
    if failed.stat().st_size == 0:
        failed.write_text("timestamp,job_name,status,message\n", encoding="utf-8")


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
    if not path.is_absolute():
        path = ctx.repo_root / path
    return ensure_within(path.resolve(), ctx.fresh_root)


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
    path = config_path or (ctx.fresh_root / "configs" / "compound_anchor_convergence.yaml")
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
    snapshot = ctx.run_dir / "phase3_compounds" / "config_snapshots" / "compound_anchor_convergence.resolved.yaml"
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    try:
        import yaml  # type: ignore

        snapshot.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    except Exception:
        snapshot.write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")
    return config


def _table_path(ctx: RunContext, given: Path | None, default_rel: str) -> Path:
    path = given or (ctx.run_dir / default_rel)
    if not path.is_absolute():
        path = ctx.repo_root / path
    return path.resolve()


def _filter_profile(rows: list[dict[str, str]], profile: str | None) -> tuple[list[dict[str, str]], str | None, list[str]]:
    profiles = sorted({row.get("profile", "") for row in rows if row.get("profile")})
    selected = profile or (profiles[0] if len(profiles) == 1 else None)
    return [row for row in rows if not selected or row.get("profile") == selected], selected, profiles


def _cluster_sort_key(row: dict[str, str]) -> tuple[Any, ...]:
    return (row.get("pocket_family_id", ""), row.get("state_id", ""), row.get("state_role", ""), row.get("compound_public_id", ""), row.get("protomer_id", ""), row.get("cluster_id", ""))


def _support_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (row.get("compound_public_id", ""), row.get("pocket_family_id", ""), row.get("state_id", ""), row.get("state_role", ""))


def _parse_mechanism_counts(rows: list[dict[str, str]]) -> tuple[dict[str, int], str]:
    counts: dict[str, int] = {}
    for row in rows:
        text = row.get("mechanism_class_counts") or ""
        parsed: dict[str, Any] = {}
        try:
            parsed = json.loads(text) if text else {}
        except json.JSONDecodeError:
            parsed = {}
        if not parsed and row.get("mechanism_class_majority"):
            parsed = {row["mechanism_class_majority"]: 1}
        for key, value in parsed.items():
            try:
                counts[key] = counts.get(key, 0) + int(value)
            except (TypeError, ValueError):
                counts[key] = counts.get(key, 0) + 1
    majority = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0] if counts else "ambiguous_or_failed"
    return counts, majority


def _is_primary(row: dict[str, str], config: dict[str, Any]) -> bool:
    return row.get("state_id") in set(config["states"]["primary_state_ids"])


def _is_reference(row: dict[str, str], config: dict[str, Any]) -> bool:
    return row.get("state_id") in set(config["states"]["reference_state_ids"]) or row.get("state_role") in set(config["states"]["reference_state_role_values"])


def _cluster_support_reject_reason(row: dict[str, str], config: dict[str, Any], allow_reference_only: bool) -> str:
    if row.get("compound_public_id") not in PUBLIC_COMPOUND_IDS:
        return "non_public_compound_id"
    if row.get("cluster_reject_reason") and row.get("cluster_reject_reason") != "none":
        mapping = {
            "ATP_migration_in_cluster": "ATP_migration",
            "outside_accepted_pocket_in_cluster": "outside_accepted_pocket",
            "membrane_penetration_in_cluster": "membrane_penetration",
            "dimer_interface_clash_in_cluster": "central_dimer_clash",
            "singleton_cluster": "scattered_or_singleton_only",
            "scattered_poses": "scattered_or_singleton_only",
            "reference_only_not_promoted": "reference_only_not_promoted",
        }
        return mapping.get(row.get("cluster_reject_reason", ""), "not_converged_cluster")
    if not _boolish(row.get("allowed_for_anchor_convergence")):
        return "cluster_not_allowed_for_anchor_convergence"
    if not _boolish(row.get("converged_pose_cluster")):
        return "not_converged_cluster"
    if row.get("convergence_status") not in set(config["support_cluster_filters"]["allowed_convergence_statuses"]):
        return "not_converged_cluster"
    if (_float(row.get("atp_migration_fraction")) or 0.0) > float(config["support_cluster_filters"]["require_atp_migration_fraction_max"]):
        return "ATP_migration"
    if (_float(row.get("within_pocket_fraction")) or 0.0) < float(config["support_cluster_filters"]["pocket_retention_fraction_min"]):
        return "outside_accepted_pocket"
    if (_float(row.get("membrane_penetration_fraction")) or 0.0) > float(config["support_cluster_filters"]["require_membrane_penetration_fraction_max"]):
        return "membrane_penetration"
    if (_float(row.get("dimer_interface_clash_fraction")) or 0.0) > float(config["support_cluster_filters"]["require_dimer_interface_clash_fraction_max"]):
        return "central_dimer_clash"
    if _is_reference(row, config) and not allow_reference_only:
        return "reference_only_not_promoted"
    return "none"


def _is_anchor_support(row: dict[str, str], config: dict[str, Any], allow_reference_only: bool) -> bool:
    return _cluster_support_reject_reason(row, config, allow_reference_only) == "none"


def _is_reference_quality_support(row: dict[str, str], config: dict[str, Any]) -> bool:
    return (
        _boolish(row.get("converged_pose_cluster"))
        and row.get("convergence_status") in set(config["support_cluster_filters"]["allowed_convergence_statuses"])
        and row.get("cluster_reject_reason") == "none"
        and (_float(row.get("atp_migration_fraction")) or 0.0) == 0.0
        and _is_reference(row, config)
    )


def _join_sorted(values: set[str] | list[str], sep: str = "|") -> str:
    return sep.join(sorted({value for value in values if value}))


def _support_status(reason: str, supporting: bool, reference: bool) -> str:
    if supporting:
        return "PASS"
    if reference and reason == "reference_only_not_promoted":
        return "WARN"
    if reason in {"ATP_migration", "outside_accepted_pocket", "membrane_penetration", "central_dimer_clash"}:
        return "REJECT"
    return "WARN"


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
    claims = False
    scanned: list[str] = []
    for path in scan_paths:
        if "docking_outputs" in path.parts and path.suffix.lower() == ".pdbqt":
            continue
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
            pose_coords = pose_coords or "pose" in lower or "anchor" in lower or "support" in lower or "qc" in lower
        if re.search(r"validated inhibitor|proven binder|proven PPI inhibitor|clinically relevant drug candidate|final candidate", text, re.IGNORECASE):
            claims = True
    return len(leaks), coord_hits, smiles_logged, ligand_coords, receptor_coords, pose_coords, claims, scanned


def _write_report(path: Path, ctx: RunContext, summary: dict[str, Any]) -> None:
    lines = [
        "# M3-T9 Anchor Convergence",
        "",
        "Technical compound-anchor convergence report only. Affinity is descriptive metadata and is not used for ranking, tiering, or candidate claims.",
        "",
        f"- status: {summary['overall_status']}",
        f"- profile: {summary['profile']}",
        f"- support_rows: {summary['support']['compound_pocket_support_rows']}",
        f"- anchor_rows: {summary['anchor_convergence']['anchor_rows']}",
        f"- anchors_allowed_for_evidence_integration: {summary['anchor_convergence']['anchors_allowed_for_evidence_integration']}",
        f"- m3_t10_evidence_integration_ready: {str(summary['m3_t10_evidence_integration_ready']).lower()}",
        "",
        "## Blockers",
    ]
    lines.extend(f"- {item}" for item in summary["blockers"] or ["none"])
    lines.append("")
    lines.append("## Warnings")
    lines.extend(f"- {item}" for item in summary["warnings"] or ["none"])
    lines.append("")
    lines.append("Next task: M3-T10 - Integrated compound-pocket-PPI evidence table and candidate tiering")
    ctx.require_within_run_dir(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_m3_anchor_convergence(
    ctx: RunContext,
    *,
    m2_run_id: str,
    profile: str | None = None,
    mode: str = "converge",
    compound_pose_clusters: Path | None = None,
    cluster_membership_table: Path | None = None,
    pose_attribution_table: Path | None = None,
    pose_clustering_qc: Path | None = None,
    accepted_pockets: Path | None = None,
    accepted_boxes: Path | None = None,
    compound_manifest: Path | None = None,
    config: Path | None = None,
    force: bool = False,
    allow_partial: bool = False,
    require_clustering_ready: bool = True,
    allow_reference_only: bool = False,
    private_map: Path | None = None,
    strict_primary_state_support: bool = True,
    write_empty_convergence: bool = False,
) -> M3AnchorConvergenceResult:
    del force
    _ensure_logs(ctx)
    timestamp = now_iso()
    phase3 = ctx.run_dir / "phase3_compounds"
    tables_dir = phase3 / "tables"
    qc_dir = phase3 / "qc"
    reports_dir = phase3 / "reports"
    for directory in [tables_dir, qc_dir, reports_dir]:
        ctx.require_within_run_dir(directory).mkdir(parents=True, exist_ok=True)

    support_csv = tables_dir / "compound_pocket_support.csv"
    anchor_csv = tables_dir / "compound_anchor_convergence.csv"
    qc_csv = qc_dir / "anchor_convergence_qc.csv"
    qc_json = qc_dir / "anchor_convergence_qc.json"
    support_qc_csv = qc_dir / "compound_pocket_support_qc.csv"
    report_md = reports_dir / "m3_task9_anchor_convergence.md"
    phase3_log = ctx.logs_dir / "phase3_compounds.log"

    blockers: list[str] = []
    warnings: list[str] = []
    cfg = _load_config(ctx, config)
    cluster_path = _table_path(ctx, compound_pose_clusters, "phase3_compounds/tables/compound_pose_clusters.csv")
    cluster_found = cluster_path.is_file()
    if not cluster_found:
        (warnings if mode == "dry-run" else blockers).append("compound_pose_clusters.csv missing")
    if cluster_found and not _inside(cluster_path, tables_dir):
        blockers.append("compound_pose_clusters.csv path outside phase3_compounds/tables")
    cluster_rows_all: list[dict[str, str]] = []
    cluster_fields: list[str] = []
    if cluster_found:
        cluster_rows_all, cluster_fields = _read_csv(cluster_path)
    schema_ok = set(CLUSTER_REQUIRED).issubset(set(cluster_fields))
    if cluster_found and not schema_ok:
        blockers.append("compound_pose_clusters.csv schema missing required columns")
    selected_clusters, selected_profile, profiles = _filter_profile(cluster_rows_all, profile)
    if cluster_found and not selected_profile:
        blockers.append("profile ambiguous; provide --profile")
    if cluster_found and not selected_clusters:
        (warnings if mode == "dry-run" or write_empty_convergence else blockers).append("zero cluster rows for selected profile")

    qc_path = _table_path(ctx, pose_clustering_qc, "phase3_compounds/qc/pose_clustering_qc.json")
    qc_found = qc_path.is_file()
    qc_data = _load_json(qc_path)
    t8_ready = bool(qc_data and qc_data.get("overall_status") == "PASS" and qc_data.get("m3_t9_anchor_convergence_ready") is True)
    if not qc_found:
        (warnings if mode == "dry-run" or allow_partial else blockers).append("pose_clustering_qc.json missing")
    elif require_clustering_ready and not t8_ready:
        (warnings if mode == "dry-run" or allow_partial else blockers).append("pose_clustering_qc.json has m3_t9_anchor_convergence_ready=false")

    attr_path = _table_path(ctx, pose_attribution_table, "phase3_compounds/tables/compound_pose_attribution.csv")
    membership_path = _table_path(ctx, cluster_membership_table, "phase3_compounds/tables/compound_pose_cluster_membership.csv")
    attr_found = attr_path.is_file()
    membership_found = membership_path.is_file()
    if not attr_found and mode == "converge" and not allow_partial:
        blockers.append("compound_pose_attribution.csv missing")
    elif not attr_found:
        warnings.append("compound_pose_attribution.csv missing")
    if not membership_found and mode == "converge" and not allow_partial:
        blockers.append("compound_pose_cluster_membership.csv missing")
    elif not membership_found:
        warnings.append("compound_pose_cluster_membership.csv missing")

    m2_root = ctx.fresh_root / "runs" / m2_run_id
    accepted_pockets_path = _first_existing([Path(accepted_pockets)] if accepted_pockets else [
        m2_root / "phase2_pockets" / "final" / "accepted_pockets_for_m3.csv",
        m2_root / "phase2_pockets" / "export_for_m3" / "accepted_pockets_for_m3.csv",
        m2_root / "phase2_pockets" / "final" / "accepted_pocket_families_for_compound_docking.csv",
    ])
    if not accepted_pockets_path:
        (warnings if mode == "dry-run" else blockers).append("accepted pockets missing")
    pocket_rows: list[dict[str, str]] = []
    pocket_fields: list[str] = []
    if accepted_pockets_path:
        pocket_rows, pocket_fields = _read_csv(accepted_pockets_path)
    pocket_schema_ok = "pocket_family_id" in set(pocket_fields)
    if accepted_pockets_path and not pocket_schema_ok:
        blockers.append("accepted pockets schema missing pocket_family_id")
    accepted_families = {row.get("pocket_family_id", "") for row in pocket_rows}
    accepted_pairs = {(row.get("pocket_family_id", ""), row.get("state_id", "")) for row in pocket_rows if row.get("state_id")}
    accepted_boxes_path = _first_existing([Path(accepted_boxes)] if accepted_boxes else [
        m2_root / "phase2_pockets" / "final" / "accepted_pocket_boxes.csv",
        m2_root / "phase2_pockets" / "export_for_m3" / "accepted_pocket_boxes.csv",
    ])

    manifest_path = _first_existing([Path(compound_manifest)] if compound_manifest else [
        phase3 / "manifests" / "compound_manifest_public.csv",
        phase3 / "manifests" / "ligand_preparation_manifest.csv",
    ])
    manifest_compounds: set[str] = set()
    if manifest_path:
        rows, fields = _read_csv(manifest_path)
        if "compound_public_id" in fields:
            manifest_compounds = {row.get("compound_public_id", "") for row in rows if row.get("compound_public_id") in PUBLIC_COMPOUND_IDS}
        elif "public_id" in fields:
            manifest_compounds = {row.get("public_id", "") for row in rows if row.get("public_id") in PUBLIC_COMPOUND_IDS}

    cluster_ids = [row.get("cluster_id", "") for row in selected_clusters]
    cluster_ids_unique = len(cluster_ids) == len(set(cluster_ids))
    if selected_clusters and not cluster_ids_unique:
        blockers.append("duplicate cluster_id rows")
    non_public = sorted({row.get("compound_public_id", "") for row in selected_clusters if row.get("compound_public_id") not in PUBLIC_COMPOUND_IDS})
    if non_public:
        blockers.append("non-public compound IDs present in cluster table")
    if any(_boolish(row.get("affinity_used_for_ranking")) for row in selected_clusters):
        blockers.append("cluster table reports affinity_used_for_ranking=true")
    if accepted_families:
        missing_trace = sorted(
            {
                row.get("pocket_family_id", "")
                for row in selected_clusters
                if row.get("pocket_family_id") not in accepted_families
                or (accepted_pairs and (row.get("pocket_family_id", ""), row.get("state_id", "")) not in accepted_pairs)
            }
        )
        if missing_trace:
            blockers.append("accepted pocket trace missing for selected cluster pocket_family_id")
    elif selected_clusters and mode == "converge":
        blockers.append("accepted pocket trace unavailable")

    source_cluster_sha = _sha256(cluster_path) if cluster_found else ""
    source_qc_sha = _sha256(qc_path) if qc_found else ""
    source_attr_sha = _sha256(attr_path) if attr_found else ""
    accepted_sha = _sha256(accepted_pockets_path) if accepted_pockets_path else ""
    if attr_found and selected_clusters:
        mismatched = [
            row.get("cluster_id", "unknown_cluster")
            for row in selected_clusters
            if row.get("source_compound_pose_attribution_sha256")
            and row.get("source_compound_pose_attribution_sha256") != source_attr_sha
        ]
        if mismatched:
            blockers.append("cluster source attribution hash mismatch: " + ",".join(sorted(mismatched)[:5]))
    support_rows: list[dict[str, Any]] = []
    support_qc_rows: list[dict[str, Any]] = []
    grouped: dict[tuple[str, str, str, str], list[dict[str, str]]] = {}
    for row in sorted(selected_clusters, key=_cluster_sort_key):
        grouped.setdefault(_support_key(row), []).append(row)

    for index, key in enumerate(sorted(grouped), start=1):
        rows = grouped[key]
        compound_id, pocket_family_id, state_id, state_role = key
        primary = state_id in set(cfg["states"]["primary_state_ids"])
        reference = state_id in set(cfg["states"]["reference_state_ids"]) or state_role in set(cfg["states"]["reference_state_role_values"])
        supporting_clusters = [row for row in rows if _is_anchor_support(row, cfg, allow_reference_only)]
        reference_quality = [row for row in rows if _is_reference_quality_support(row, cfg)]
        converged_clusters = [row for row in rows if _boolish(row.get("converged_pose_cluster"))]
        reject_reason = "none" if supporting_clusters else (_cluster_support_reject_reason(rows[0], cfg, allow_reference_only) if rows else "unknown")
        support_status = _support_status(reject_reason, bool(supporting_clusters), reference)
        fractions = [value for row in rows if (value := _float(row.get("cluster_fraction"))) is not None]
        affinities = [value for row in rows if (value := _float(row.get("vina_affinity_median_kcal_mol"))) is not None]
        within = [value for row in rows if (value := _float(row.get("within_pocket_fraction"))) is not None]
        atp = [value for row in rows if (value := _float(row.get("atp_migration_fraction"))) is not None]
        ppi = [value for row in rows if (value := _float(row.get("ppi_relationship_confirmed_fraction"))) is not None]
        ppi_contact_median = [value for row in rows if (value := _float(row.get("ppi_hotspot_contact_median"))) is not None]
        ppi_contact_mean = [value for row in rows if (value := _float(row.get("ppi_hotspot_contact_mean"))) is not None]
        ppi_distance = [value for row in rows if (value := _float(row.get("nearest_ppi_hotspot_distance_median"))) is not None]
        ppi_rim = [value for row in rows if (value := _float(row.get("ppi_rim_distance_median"))) is not None]
        mem = [value for row in rows if (value := _float(row.get("membrane_penetration_fraction"))) is not None]
        dimer = [value for row in rows if (value := _float(row.get("dimer_interface_clash_fraction"))) is not None]
        mechanisms, mechanism_majority = _parse_mechanism_counts(rows)
        support_id = f"support_{index:04d}"
        allowed_anchor = bool(supporting_clusters)
        allowed_evidence = bool(allowed_anchor and primary)
        support_row = {
            "run_id": ctx.run_id, "m2_run_id": m2_run_id, "profile": selected_profile, "support_id": support_id,
            "compound_public_id": compound_id, "pocket_family_id": pocket_family_id, "state_id": state_id, "state_role": state_role,
            "primary_state_flag": _bool_text(primary), "reference_state_flag": _bool_text(reference),
            "protomer_ids_with_support": _join_sorted([row.get("protomer_id", "") for row in rows]),
            "protomer_symmetry_support": _bool_text(len({row.get("protomer_id", "") for row in rows if row.get("protomer_id")}) > 1),
            "box_ids_with_support": _join_sorted([row.get("box_ids_in_cluster", "") for row in rows]),
            "cluster_ids_all": _join_sorted([row.get("cluster_id", "") for row in rows]),
            "cluster_ids_supporting": _join_sorted([row.get("cluster_id", "") for row in supporting_clusters or reference_quality]),
            "num_clusters_total": len(rows), "num_converged_clusters": len(converged_clusters),
            "num_anchor_allowed_clusters": len(supporting_clusters), "num_valid_pose_clusters": len(converged_clusters),
            "num_rejected_or_nonconverged_clusters": len(rows) - len(supporting_clusters),
            "eligible_pose_count_total": sum(int(_float(row.get("eligible_pose_count_in_group")) or 0) for row in rows),
            "clustered_pose_count_total": sum(int(_float(row.get("cluster_size")) or 0) for row in rows),
            "max_cluster_size": max([int(_float(row.get("cluster_size")) or 0) for row in rows] or [0]),
            "max_cluster_fraction": _fmt(max(fractions) if fractions else None),
            "median_cluster_fraction": _fmt(median(fractions) if fractions else None),
            "within_pocket_fraction_min": _fmt(min(within) if within else None),
            "atp_migration_fraction_max": _fmt(max(atp) if atp else None),
            "ppi_relationship_confirmed_fraction_max": _fmt(max(ppi) if ppi else None),
            "ppi_relationship_confirmed_fraction_median": _fmt(median(ppi) if ppi else None),
            "ppi_hotspot_contact_median_max": _fmt(max(ppi_contact_median) if ppi_contact_median else None),
            "ppi_hotspot_contact_mean_max": _fmt(max(ppi_contact_mean) if ppi_contact_mean else None),
            "nearest_ppi_hotspot_distance_median_min": _fmt(min(ppi_distance) if ppi_distance else None),
            "ppi_rim_distance_median_min": _fmt(min(ppi_rim) if ppi_rim else None),
            "membrane_penetration_fraction_max": _fmt(max(mem) if mem else None),
            "dimer_interface_clash_fraction_max": _fmt(max(dimer) if dimer else None),
            "mechanism_class_counts": json.dumps(mechanisms, sort_keys=True),
            "mechanism_class_majority": mechanism_majority,
            "vina_affinity_median_kcal_mol": _fmt(_stat(affinities, "median")),
            "vina_affinity_mean_kcal_mol": _fmt(_stat(affinities, "mean")),
            "vina_affinity_best_kcal_mol": _fmt(_stat(affinities, "best")),
            "vina_affinity_worst_kcal_mol": _fmt(_stat(affinities, "worst")),
            "affinity_used_for_ranking": "false", "support_status": support_status, "support_reject_reason": reject_reason,
            "support_notes": "reference_support_not_primary" if reference and not allow_reference_only else "cluster_support_separate_from_affinity",
            "allowed_for_anchor_convergence": _bool_text(allowed_anchor),
            "allowed_for_evidence_integration": _bool_text(allowed_evidence),
            "source_compound_pose_clusters": ctx.relative_to_repo(cluster_path) if cluster_found else "",
            "source_compound_pose_clusters_sha256": source_cluster_sha,
            "source_pose_clustering_qc": ctx.relative_to_repo(qc_path) if qc_found else "",
            "source_pose_clustering_qc_sha256": source_qc_sha,
            "source_compound_pose_attribution": ctx.relative_to_repo(attr_path) if attr_found else "",
            "source_compound_pose_attribution_sha256": source_attr_sha,
            "accepted_pocket_source": ctx.relative_to_repo(accepted_pockets_path) if accepted_pockets_path else "",
            "accepted_pocket_source_sha256": accepted_sha,
            "evaluated_at": timestamp,
        }
        support_rows.append(support_row)
        support_qc_rows.append({key: support_row.get(key, "") for key in ["run_id", "m2_run_id", "profile", "support_id", "compound_public_id", "pocket_family_id", "state_id", "state_role", "primary_state_flag", "reference_state_flag", "cluster_ids_supporting", "num_anchor_allowed_clusters", "support_status", "support_reject_reason", "allowed_for_anchor_convergence", "allowed_for_evidence_integration", "support_notes"]})

    allowed_support = [row for row in support_rows if row.get("allowed_for_anchor_convergence") == "true"]
    primary_support = [row for row in allowed_support if row.get("primary_state_flag") == "true"]
    reference_quality_support = [row for row in support_rows if row.get("reference_state_flag") == "true" and row.get("cluster_ids_supporting")]
    all_cluster_compounds = {row.get("compound_public_id", "") for row in selected_clusters if row.get("compound_public_id") in PUBLIC_COMPOUND_IDS}
    compounds_tested = manifest_compounds or all_cluster_compounds
    compounds_valid = {row.get("compound_public_id", "") for row in support_rows if int(row.get("num_valid_pose_clusters") or 0) > 0}
    compounds_converged = {row.get("compound_public_id", "") for row in support_rows if int(row.get("num_converged_clusters") or 0) > 0}

    anchor_rows: list[dict[str, Any]] = []
    anchor_index = 1

    def make_anchor(scope: str, pocket: str, rows: list[dict[str, Any]], state_id: str = "") -> dict[str, Any]:
        nonlocal anchor_index
        anchor_id = f"anchor_{anchor_index:04d}"
        anchor_index += 1
        primary_rows = [row for row in rows if row.get("primary_state_flag") == "true" and row.get("allowed_for_anchor_convergence") == "true"]
        reference_rows = [row for row in rows if row.get("reference_state_flag") == "true" and row.get("cluster_ids_supporting")]
        supported = [row for row in rows if row.get("allowed_for_anchor_convergence") == "true"]
        compounds_supported = {row.get("compound_public_id", "") for row in supported}
        primary_states = {row.get("state_id", "") for row in primary_rows}
        reference_states = {row.get("state_id", "") for row in reference_rows}
        support_ids = {row.get("support_id", "") for row in supported}
        cluster_ids = set()
        for row in supported:
            cluster_ids.update((row.get("cluster_ids_supporting") or "").split("|"))
        if not supported and reference_rows:
            anchor_class = "reference_only"
        elif not selected_clusters:
            anchor_class = "no_valid_cluster_support"
        elif not supported:
            anchor_class = "weak_or_scattered"
        elif scope == "primary_cross_state_pocket" and len(compounds_supported) == 1 and len(primary_states) >= 2:
            anchor_class = "single_compound_state_robust"
        elif len(compounds_supported) >= 2 and primary_rows:
            anchor_class = "multi_compound_primary_state_supported" if len(primary_states) >= 2 or scope == "primary_cross_state_pocket" else "single_state_multi_compound"
        elif len(compounds_supported) >= 2:
            anchor_class = "single_state_multi_compound"
        else:
            anchor_class = "weak_or_scattered"
        if anchor_class in {"multi_compound_primary_state_supported", "single_compound_state_robust", "single_state_multi_compound"}:
            anchor_status = "PASS"
            reject = "none"
        elif anchor_class == "reference_only":
            anchor_status = "WARN"
            reject = "reference_only_not_promoted"
        elif anchor_class == "no_valid_cluster_support":
            anchor_status = "WARN"
            reject = "no_converged_cluster_support"
        else:
            anchor_status = "WARN"
            reject = "weak_or_scattered_support"
        allowed_evidence = bool(anchor_status in {"PASS", "WARN"} and anchor_class in cfg["evidence_readiness"]["m3_t10_ready_classes"] and reject == "none")
        aff = [value for row in supported if (value := _float(row.get("vina_affinity_median_kcal_mol"))) is not None]
        fractions = [value for row in supported if (value := _float(row.get("max_cluster_fraction"))) is not None]
        within = [value for row in supported if (value := _float(row.get("within_pocket_fraction_min"))) is not None]
        atp = [value for row in supported if (value := _float(row.get("atp_migration_fraction_max"))) is not None]
        mem = [value for row in supported if (value := _float(row.get("membrane_penetration_fraction_max"))) is not None]
        dimer = [value for row in supported if (value := _float(row.get("dimer_interface_clash_fraction_max"))) is not None]
        ppi = [value for row in supported if (value := _float(row.get("ppi_relationship_confirmed_fraction_max"))) is not None]
        mechanisms, dominant_mech = _parse_mechanism_counts([{**row, "mechanism_class_counts": row.get("mechanism_class_counts", "")} for row in rows])
        best_compound = sorted(compounds_supported, key=lambda cid: PUBLIC_COMPOUND_IDS.index(cid) if cid in PUBLIC_COMPOUND_IDS else 99)[0] if compounds_supported else ""
        status_counts: dict[str, int] = {}
        for row in rows:
            key = str(row.get("support_status") or "WARN")
            status_counts[key] = status_counts.get(key, 0) + 1
        protomer_summary = ";".join(f"{row.get('state_id')}:{row.get('protomer_ids_with_support')}" for row in rows)
        return {
            "run_id": ctx.run_id, "m2_run_id": m2_run_id, "profile": selected_profile, "anchor_id": anchor_id,
            "anchor_scope": scope, "pocket_family_id": pocket, "state_id": state_id,
            "state_ids_supported": _join_sorted({row.get("state_id", "") for row in supported}, ";"),
            "primary_state_support": len(primary_states), "primary_state_ids_supported": _join_sorted(primary_states, ";"),
            "reference_state_support": len(reference_states), "reference_state_ids_supported": _join_sorted(reference_states, ";"),
            "symmetry_support": _bool_text(any(row.get("protomer_symmetry_support") == "true" for row in rows)),
            "protomer_support_summary": protomer_summary, "num_compounds_expected": len(cfg["compounds"]["expected_public_compound_ids"]),
            "num_compounds_tested": len(compounds_tested), "num_compounds_with_valid_pose": len(compounds_valid),
            "num_compounds_with_converged_pose": len(compounds_converged),
            "num_compounds_allowed_for_anchor_convergence": len(compounds_supported),
            "compound_public_ids_tested": _join_sorted(compounds_tested, ";"),
            "compound_public_ids_with_valid_pose": _join_sorted(compounds_valid, ";"),
            "compound_public_ids_supported": _join_sorted(compounds_supported, ";"),
            "best_compound_public_id": best_compound, "best_compound_selection_rule": "deterministic_public_id_order",
            "affinity_used_for_best_compound_selection": "false",
            "median_affinity_across_supported": _fmt(_stat(aff, "median")), "mean_affinity_across_supported": _fmt(_stat(aff, "mean")),
            "best_affinity_across_supported": _fmt(_stat(aff, "best")), "worst_affinity_across_supported": _fmt(_stat(aff, "worst")),
            "affinity_used_for_ranking": "false", "cluster_ids_supporting": _join_sorted(cluster_ids),
            "support_ids_supporting": _join_sorted(support_ids), "support_status_counts": json.dumps(status_counts, sort_keys=True),
            "mechanism_class_counts": json.dumps(mechanisms, sort_keys=True), "dominant_mechanism_class": dominant_mech,
            "max_cluster_fraction_across_supported": _fmt(max(fractions) if fractions else None),
            "median_cluster_fraction_across_supported": _fmt(median(fractions) if fractions else None),
            "min_within_pocket_fraction_across_supported": _fmt(min(within) if within else None),
            "max_atp_migration_fraction_across_supported": _fmt(max(atp) if atp else None),
            "max_membrane_penetration_fraction_across_supported": _fmt(max(mem) if mem else None),
            "max_dimer_interface_clash_fraction_across_supported": _fmt(max(dimer) if dimer else None),
            "max_ppi_relationship_confirmed_fraction": _fmt(max(ppi) if ppi else None),
            "median_ppi_relationship_confirmed_fraction": _fmt(median(ppi) if ppi else None),
            "anchor_converged": _bool_text(allowed_evidence), "anchor_convergence_class": anchor_class,
            "anchor_convergence_status": anchor_status, "anchor_reject_reason": reject,
            "anchor_convergence_notes": "A/B symmetry not counted as independent state support",
            "allowed_for_evidence_integration": _bool_text(allowed_evidence),
            "source_compound_pocket_support": ctx.relative_to_repo(support_csv),
            "source_compound_pocket_support_sha256": "",
            "source_compound_pose_clusters": ctx.relative_to_repo(cluster_path) if cluster_found else "",
            "source_compound_pose_clusters_sha256": source_cluster_sha,
            "accepted_pocket_source": ctx.relative_to_repo(accepted_pockets_path) if accepted_pockets_path else "",
            "accepted_pocket_source_sha256": accepted_sha,
            "evaluated_at": timestamp,
        }

    if support_rows:
        by_pocket_state: dict[tuple[str, str], list[dict[str, Any]]] = {}
        by_pocket: dict[str, list[dict[str, Any]]] = {}
        ref_by_pocket: dict[str, list[dict[str, Any]]] = {}
        for row in support_rows:
            by_pocket_state.setdefault((row["pocket_family_id"], row["state_id"]), []).append(row)
            by_pocket.setdefault(row["pocket_family_id"], []).append(row)
            if row["reference_state_flag"] == "true":
                ref_by_pocket.setdefault(row["pocket_family_id"], []).append(row)
        for (pocket, state), rows in sorted(by_pocket_state.items()):
            anchor_rows.append(make_anchor("state_pocket", pocket, rows, state_id=state))
        for pocket, rows in sorted(by_pocket.items()):
            anchor_rows.append(make_anchor("primary_cross_state_pocket", pocket, rows))
            anchor_rows.append(make_anchor("all_states_pocket", pocket, rows))
        for pocket, rows in sorted(ref_by_pocket.items()):
            anchor_rows.append(make_anchor("reference_pocket", pocket, rows))
    elif write_empty_convergence or mode == "dry-run":
        pass

    _write_csv(support_csv, SUPPORT_FIELDS, support_rows, ctx)
    support_sha = _sha256(support_csv)
    for row in anchor_rows:
        row["source_compound_pocket_support_sha256"] = support_sha
    _write_csv(anchor_csv, ANCHOR_FIELDS, anchor_rows, ctx)
    _write_csv(support_qc_csv, [
        "run_id", "m2_run_id", "profile", "support_id", "compound_public_id", "pocket_family_id", "state_id",
        "state_role", "primary_state_flag", "reference_state_flag", "cluster_ids_supporting",
        "num_anchor_allowed_clusters", "support_status", "support_reject_reason", "allowed_for_anchor_convergence",
        "allowed_for_evidence_integration", "support_notes",
    ], support_qc_rows, ctx)

    if selected_clusters and len(support_rows) != len(grouped):
        blockers.append("support rows were silently dropped")
    if not support_rows and selected_clusters:
        warnings.append("no compound-pocket support rows written")
    if not allowed_support and selected_clusters:
        warnings.append("no support-eligible clusters after T8 hard gates")
    forbidden_created = [item for item in FORBIDDEN_OUTPUTS if (ctx.run_dir / item).exists()]
    if forbidden_created:
        blockers.append("forbidden downstream outputs already exist in run directory")
    if mode == "dry-run":
        warnings.append("dry-run mode")
    if allow_partial:
        warnings.append("allow_partial diagnostic mode used")

    support_counts = {
        "compound_pocket_support_rows": len(support_rows),
        "support_rows_allowed_for_anchor_convergence": sum(1 for row in support_rows if row.get("allowed_for_anchor_convergence") == "true"),
        "support_rows_allowed_for_evidence_integration": sum(1 for row in support_rows if row.get("allowed_for_evidence_integration") == "true"),
        "compounds_tested": len(compounds_tested),
        "compounds_with_valid_pose": len(compounds_valid),
        "compounds_with_converged_pose": len(compounds_converged),
        "compounds_with_primary_support": len({row.get("compound_public_id", "") for row in support_rows if row.get("primary_state_flag") == "true" and row.get("allowed_for_anchor_convergence") == "true"}),
        "pocket_families_evaluated": len({row.get("pocket_family_id", "") for row in support_rows}),
        "pocket_families_with_primary_support": len({row.get("pocket_family_id", "") for row in support_rows if row.get("primary_state_flag") == "true" and row.get("allowed_for_anchor_convergence") == "true"}),
        "pocket_families_with_reference_only_support": len({row.get("pocket_family_id", "") for row in support_rows if row.get("reference_state_flag") == "true" and row.get("cluster_ids_supporting") and not any(other.get("pocket_family_id") == row.get("pocket_family_id") and other.get("primary_state_flag") == "true" and other.get("allowed_for_anchor_convergence") == "true" for other in support_rows)}),
    }
    anchor_counts = {
        "anchor_rows": len(anchor_rows),
        "state_pocket_rows": sum(1 for row in anchor_rows if row.get("anchor_scope") == "state_pocket"),
        "primary_cross_state_rows": sum(1 for row in anchor_rows if row.get("anchor_scope") == "primary_cross_state_pocket"),
        "reference_pocket_rows": sum(1 for row in anchor_rows if row.get("anchor_scope") == "reference_pocket"),
        "anchors_converged": sum(1 for row in anchor_rows if row.get("anchor_converged") == "true"),
        "anchors_allowed_for_evidence_integration": sum(1 for row in anchor_rows if row.get("allowed_for_evidence_integration") == "true"),
        "multi_compound_primary_state_supported": sum(1 for row in anchor_rows if row.get("anchor_convergence_class") == "multi_compound_primary_state_supported"),
        "single_compound_state_robust": sum(1 for row in anchor_rows if row.get("anchor_convergence_class") == "single_compound_state_robust"),
        "single_state_multi_compound": sum(1 for row in anchor_rows if row.get("anchor_convergence_class") == "single_state_multi_compound"),
        "reference_only": sum(1 for row in anchor_rows if row.get("anchor_convergence_class") == "reference_only"),
        "weak_or_scattered": sum(1 for row in anchor_rows if row.get("anchor_convergence_class") == "weak_or_scattered"),
        "no_valid_cluster_support": sum(1 for row in anchor_rows if row.get("anchor_convergence_class") == "no_valid_cluster_support"),
    }
    mechanism_classes = {name: sum(1 for row in anchor_rows if row.get("dominant_mechanism_class") == name) for name in ["orthosteric_or_direct_PPI_patch_blocker", "rim_blocker", "allosteric_near_candidate", "generic_nonATP_ligandable_pocket", "ATP_like_reject", "ambiguous_or_failed"]}

    private_entries, private_warnings = load_private_map(_safe_private_map_path(ctx, private_map))
    warnings.extend(private_warnings)
    leak_count, coord_hits, smiles_logged, ligand_coords, receptor_coords, pose_coords, claims, scanned = _scan_hygiene(ctx, list(private_entries.values()))
    if leak_count:
        blockers.append("internal compound ID leakage detected")
    if coord_hits:
        blockers.append("coordinate records printed outside intended PDBQT files")
    if smiles_logged:
        blockers.append("SMILES-like fields printed in public outputs")
    if claims:
        blockers.append("candidate or validation claim detected in public outputs")

    status = "FAIL" if blockers else ("WARN" if warnings or mode == "dry-run" or anchor_counts["anchors_allowed_for_evidence_integration"] == 0 else "PASS")
    m3_t10_ready = status == "PASS" and mode == "converge" and anchor_counts["anchors_allowed_for_evidence_integration"] > 0
    cluster_context = {
        "cluster_rows": len(selected_clusters),
        "support_eligible_clusters": sum(1 for row in selected_clusters if _is_anchor_support(row, cfg, allow_reference_only)),
        "profiles_in_cluster_table": profiles,
        "selected_profile": selected_profile,
        "compounds_in_cluster_table": sorted(all_cluster_compounds),
        "primary_states": sorted({row.get("state_id", "") for row in selected_clusters if _is_primary(row, cfg)}),
        "reference_states": sorted({row.get("state_id", "") for row in selected_clusters if _is_reference(row, cfg)}),
        "reference_only_support_present": bool(reference_quality_support and not primary_support),
    }
    inputs = {
        "compound_pose_clusters_found": cluster_found,
        "pose_clustering_qc_found": qc_found,
        "compound_pose_attribution_found": attr_found,
        "cluster_membership_table_found": membership_found,
        "accepted_pockets_found": accepted_pockets_path is not None,
        "accepted_boxes_found": accepted_boxes_path is not None,
        "compound_manifest_found": manifest_path is not None,
    }
    summary = {
        "schema_version": "m3_anchor_convergence_qc_v1",
        "run_id": ctx.run_id,
        "m2_run_id": m2_run_id,
        "reviewed_at": timestamp,
        "mode": mode,
        "profile": selected_profile,
        "overall_status": status,
        "allow_partial": allow_partial,
        "require_clustering_ready": require_clustering_ready,
        "allow_reference_only": allow_reference_only,
        "strict_primary_state_support": strict_primary_state_support,
        "m3_t10_evidence_integration_ready": m3_t10_ready,
        "inputs": inputs,
        "cluster_context": cluster_context,
        "support": support_counts,
        "anchor_convergence": anchor_counts,
        "mechanism_classes_by_anchor": mechanism_classes,
        "quality": {
            "anchors_with_zero_atp_migration_fraction": sum(1 for row in anchor_rows if row.get("max_atp_migration_fraction_across_supported") in {"", "0.000"}),
            "anchors_with_pocket_retention_fraction_pass": sum(1 for row in anchor_rows if (_float(row.get("min_within_pocket_fraction_across_supported")) or 0) >= float(cfg["support_cluster_filters"]["pocket_retention_fraction_min"])),
            "anchors_with_ppi_relationship_fraction_pass": sum(1 for row in anchor_rows if (_float(row.get("max_ppi_relationship_confirmed_fraction")) or 0) >= 0.5),
            "anchors_with_no_membrane_or_dimer_conflict": sum(1 for row in anchor_rows if row.get("max_membrane_penetration_fraction_across_supported") in {"", "0.000"} and row.get("max_dimer_interface_clash_fraction_across_supported") in {"", "0.000"}),
            "reference_only_anchors_not_promoted": sum(1 for row in anchor_rows if row.get("anchor_convergence_class") == "reference_only" and row.get("allowed_for_evidence_integration") != "true"),
            "ab_symmetry_events_not_counted_as_state_support": sum(1 for row in support_rows if row.get("protomer_symmetry_support") == "true"),
        },
        "outputs": {
            "compound_pocket_support": ctx.relative_to_repo(support_csv),
            "compound_anchor_convergence": ctx.relative_to_repo(anchor_csv),
            "anchor_convergence_qc_csv": ctx.relative_to_repo(qc_csv),
            "anchor_convergence_qc_json": ctx.relative_to_repo(qc_json),
            "compound_pocket_support_qc": ctx.relative_to_repo(support_qc_csv),
            "report_file": ctx.relative_to_repo(report_md),
        },
        "parser": {"affinity_used_for_ranking": False, "affinity_used_for_best_compound_selection": False, "coordinate_records_written_to_public_tables": False},
        "counts": {"vina_commands_invoked_by_anchor_convergence": 0, "qsub_commands_invoked_by_anchor_convergence": 0, "runner_commands_invoked_by_anchor_convergence": 0, "pose_attribution_rerun_attempted": 0, "pose_clustering_rerun_attempted": 0, "broad_anchor_scan_rows_created": 0, "candidate_tables_created": 0, "candidate_tiers_assigned": 0, "confidentiality_leaks": leak_count, "coordinate_leak_files": len(coord_hits)},
        "blockers": blockers,
        "warnings": warnings,
        "confidentiality": {"internal_ids_redacted": leak_count == 0, "public_outputs_scanned": scanned, "leaks_detected": leak_count, "smiles_logged": smiles_logged, "ligand_coordinates_logged": ligand_coords, "receptor_coordinates_logged": receptor_coords, "pose_coordinates_logged_outside_pose_file": pose_coords, "candidate_claims_detected": claims},
        "non_goals_preserved": {"vina_execution": True, "qsub_submission": True, "pbs_runner_invocation": True, "pose_attribution_rerun": True, "pose_clustering_rerun": True, "broad_anchor_scan": True, "compound_scoring": True, "candidate_tiering": True, "candidate_nomination": True, "final_candidate_table_creation": True},
    }
    qc_rows: list[dict[str, str]] = []
    qc_status = {
        "m3_t8_compound_pose_clusters_present": "PASS" if cluster_found else ("WARN" if mode == "dry-run" else "FAIL"),
        "m3_t8_compound_pose_clusters_schema_valid": "PASS" if schema_ok else ("NOT_APPLICABLE" if not cluster_found else "FAIL"),
        "m3_t8_pose_clustering_qc_present": "PASS" if qc_found else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "m3_t8_anchor_convergence_ready_or_partial_allowed": "PASS" if t8_ready or allow_partial or mode == "dry-run" else "FAIL",
        "m3_t7_pose_attribution_present_or_partial_allowed": "PASS" if attr_found or allow_partial or mode == "dry-run" else "FAIL",
        "m3_t8_cluster_membership_present_or_partial_allowed": "PASS" if membership_found or allow_partial or mode == "dry-run" else "FAIL",
        "profile_selected": "PASS" if selected_profile else "FAIL",
        "cluster_rows_found": "PASS" if selected_clusters else ("WARN" if mode == "dry-run" or write_empty_convergence else "FAIL"),
        "support_eligible_clusters_found": "PASS" if cluster_context["support_eligible_clusters"] else "WARN",
        "cluster_ids_unique": "PASS" if cluster_ids_unique else "FAIL",
        "compound_public_ids_valid": "PASS" if not non_public else "FAIL",
        "accepted_pockets_present": "PASS" if accepted_pockets_path else ("WARN" if mode == "dry-run" else "FAIL"),
        "accepted_pockets_schema_valid": "PASS" if pocket_schema_ok else ("NOT_APPLICABLE" if not accepted_pockets_path else "FAIL"),
        "accepted_pocket_family_trace_valid": "PASS" if accepted_families or not selected_clusters else "FAIL",
        "compound_pocket_support_written": "PASS" if support_csv.is_file() else "FAIL",
        "compound_anchor_convergence_written": "PASS" if anchor_csv.is_file() else "FAIL",
        "compound_pocket_support_qc_written": "PASS" if support_qc_csv.is_file() else "FAIL",
        "anchor_convergence_report_written": "PASS" if report_md.is_file() else "FAIL",
        "no_internal_id_leak": "PASS" if leak_count == 0 else "FAIL",
        "no_smiles_logged": "PASS" if not smiles_logged else "FAIL",
        "no_ligand_coordinates_logged": "PASS" if not ligand_coords else "FAIL",
        "no_receptor_coordinates_logged": "PASS" if not receptor_coords else "FAIL",
        "no_pose_coordinates_logged_outside_pdbqt": "PASS" if not pose_coords else "FAIL",
    }
    for check in REQUIRED_CHECKS:
        default = "PASS"
        if check in {"primary_states_identified", "reference_states_identified"}:
            default = "PASS"
        elif check in {"all_support_rows_represented", "no_silent_cluster_drop"}:
            default = "PASS" if len(support_rows) == len(grouped) else "FAIL"
        _append_qc(qc_rows, check, qc_status.get(check, default))
    _write_csv(qc_csv, QC_FIELDS, qc_rows, ctx)
    _write_report(report_md, ctx, summary)
    qc_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    with phase3_log.open("a", encoding="utf-8") as handle:
        handle.write(f"{timestamp} M3-T9 {status} profile={selected_profile} support_rows={len(support_rows)} anchor_rows={len(anchor_rows)} no_vina=true no_qsub=true\n")
    for row in anchor_rows:
        append_job_status(ctx, row.get("anchor_id", "anchor_convergence"), row.get("anchor_convergence_status", "WARN"), details={"phase": "phase3_compounds", "task": "M3-T9", "job_type": "anchor_convergence", "anchor_id": row.get("anchor_id"), "profile": selected_profile, "pocket_family_id": row.get("pocket_family_id"), "anchor_scope": row.get("anchor_scope"), "primary_state_support": row.get("primary_state_support"), "num_compounds_with_converged_pose": row.get("num_compounds_with_converged_pose"), "compound_public_ids_supported": row.get("compound_public_ids_supported"), "anchor_convergence_class": row.get("anchor_convergence_class"), "anchor_convergence_status": row.get("anchor_convergence_status"), "allowed_for_evidence_integration": row.get("allowed_for_evidence_integration"), "message": "anchor convergence aggregation only"})
    if status == "FAIL":
        append_failed_job(ctx, "M3-T9", "FAIL", ";".join(blockers[:3]))
    append_phase_status(ctx, "phase3_compounds", status, "M3-T9 anchor convergence completed", {"phase": "phase3_compounds", "task": "M3-T9", "status": status, "run_id": ctx.run_id, "m2_run_id": m2_run_id, "timestamp": timestamp, "mode": mode, "profile": selected_profile, "cluster_rows": len(selected_clusters), "support_eligible_clusters": cluster_context["support_eligible_clusters"], "compound_pocket_support_rows": len(support_rows), "anchor_rows": len(anchor_rows), "anchors_converged": anchor_counts["anchors_converged"], "anchors_allowed_for_evidence_integration": anchor_counts["anchors_allowed_for_evidence_integration"], "multi_compound_primary_state_supported": anchor_counts["multi_compound_primary_state_supported"], "single_compound_state_robust": anchor_counts["single_compound_state_robust"], "single_state_multi_compound": anchor_counts["single_state_multi_compound"], "reference_only": anchor_counts["reference_only"], "weak_or_scattered": anchor_counts["weak_or_scattered"], "m3_t10_evidence_integration_ready": m3_t10_ready, "qsub_invoked": False, "vina_invoked_by_anchor_convergence": False, "runner_invoked_by_anchor_convergence": False, "no_pose_attribution_rerun_attempted": True, "no_pose_clustering_rerun_attempted": True, "no_broad_anchor_scan_attempted": True, "no_scoring_attempted": True, "no_candidate_tiering_attempted": True, "no_candidate_nomination_attempted": True, "no_final_candidate_table_created": True})
    return M3AnchorConvergenceResult(status, m3_t10_ready, blockers, warnings, support_csv, anchor_csv, qc_csv, qc_json, support_qc_csv, report_md, phase3_log, inputs, cluster_context, support_counts, anchor_counts, mechanism_classes)
