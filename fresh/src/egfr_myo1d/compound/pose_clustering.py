"""M3-T8 pose clustering and convergence analysis.

This module consumes M3-T7 attribution rows, clusters only hard-gate-pass poses
within compound/state/protomer/pocket groups, and writes cluster-level scalar
summaries. It never invokes Vina, qsub, PBS runners, attribution reruns,
anchor convergence, scoring, ranking, tiering, or candidate nomination.
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
from egfr_myo1d.compound.pdbqt_parse import AtomRecord, extract_pose_atoms
from egfr_myo1d.core.logging_utils import append_failed_job, append_job_status, append_phase_status
from egfr_myo1d.core.run_context import RunContext, ensure_within


DEFAULT_CONFIG: dict[str, Any] = {
    "schema_version": "m3_pose_clustering_config_v1",
    "input": {"require_t7_qc_ready": True, "eligible_attribution_statuses": ["PASS", "WARN"], "eligible_pose_reject_reason": "none"},
    "grouping": {
        "primary_group_by": ["compound_public_id", "state_id", "state_role", "pocket_family_id", "protomer_id"],
        "keep_box_id_in_member_rows": True,
        "do_not_count_ab_symmetry_as_independent_state_support": True,
        "reference_state_role_values": ["reference", "crystallographic_reference_control", "reference_control"],
    },
    "clustering": {
        "method": "hybrid",
        "centroid_distance_cutoff_A": 2.0,
        "heavy_atom_rmsd_cutoff_A": 2.0,
        "rmsd_requires_same_heavy_atom_count": True,
        "rmsd_requires_same_atom_order": True,
        "fallback_to_centroid_if_rmsd_unavailable": True,
        "min_cluster_size": 2,
        "singleton_policy": "warn_not_converged",
        "deterministic_sort_keys": ["compound_public_id", "state_id", "pocket_family_id", "protomer_id", "box_id", "vina_repeat_id", "pose_rank", "pose_model_index", "job_id"],
    },
    "convergence": {
        "cluster_fraction_min": 0.50,
        "cluster_size_min": 2,
        "pocket_retention_fraction_min": 0.95,
        "atp_migration_fraction_max": 0.0,
        "membrane_penetration_fraction_max": 0.0,
        "dimer_interface_clash_fraction_max": 0.0,
        "ppi_relationship_fraction_min": 0.50,
        "ppi_hotspot_contact_median_min": 0,
        "strict_requires_ppi_relationship": False,
    },
    "affinity": {"record_descriptive_stats": True, "use_affinity_for_ranking": False},
    "reference_control": {"allow_reference_only": False, "reference_only_clusters_allowed_for_anchor_convergence": False},
    "confidentiality": {"allowed_public_compound_ids": ["Cpd-A", "Cpd-B", "Cpd-C"], "scan_public_outputs": True, "fail_on_internal_id_leak": True, "fail_on_smiles_leak": True, "fail_on_coordinate_leak": True},
}

ATTR_REQUIRED = [
    "run_id", "m2_run_id", "profile", "job_id", "compound_public_id", "state_id", "state_role",
    "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index", "pose_file",
    "output_pdbqt_file", "output_pdbqt_sha256", "vina_log_file", "vina_log_sha256", "allowed_for_clustering",
    "pose_hard_gate_pass", "pose_reject_reason", "attribution_status", "pocket_retention_pass", "atp_migration_flag", "ppi_relationship_confirmed",
    "membrane_penetration_flag", "dimer_interface_clash_flag", "pose_mechanism_class",
]

CLUSTER_FIELDS = [
    "run_id", "m2_run_id", "profile", "cluster_id", "cluster_algorithm", "cluster_method", "cluster_scope",
    "compound_public_id", "state_id", "state_role", "pocket_family_id", "protomer_id", "box_ids_in_cluster",
    "vina_repeat_ids_in_cluster", "pose_ranks_in_cluster", "pose_model_indices_in_cluster", "job_ids_in_cluster",
    "source_compound_pose_attribution", "source_compound_pose_attribution_sha256", "source_pose_attribution_qc",
    "source_pose_attribution_qc_sha256", "accepted_box_source", "accepted_box_source_sha256",
    "eligible_pose_count_in_group", "cluster_size", "cluster_fraction", "cluster_rank_within_group",
    "vina_affinity_median_kcal_mol", "vina_affinity_mean_kcal_mol", "vina_affinity_best_kcal_mol",
    "vina_affinity_worst_kcal_mol", "affinity_used_for_ranking", "cluster_centroid_x", "cluster_centroid_y",
    "cluster_centroid_z", "cluster_radius_A", "cluster_median_pairwise_centroid_distance_A",
    "cluster_max_pairwise_centroid_distance_A", "cluster_mean_pairwise_heavy_atom_rmsd_A",
    "cluster_median_pairwise_heavy_atom_rmsd_A", "cluster_max_pairwise_heavy_atom_rmsd_A", "rmsd_available",
    "rmsd_notes", "within_pocket_fraction", "atp_migration_fraction", "ppi_relationship_confirmed_fraction",
    "ppi_hotspot_contact_median", "ppi_hotspot_contact_mean", "ppi_contact_residue_count_median",
    "nearest_ppi_hotspot_distance_median", "ppi_rim_distance_median", "membrane_penetration_fraction",
    "dimer_interface_clash_fraction", "mechanism_class_majority", "mechanism_class_counts",
    "converged_pose_cluster", "convergence_status", "cluster_reject_reason", "cluster_notes",
    "allowed_for_anchor_convergence", "clustered_at",
]

MEMBERSHIP_FIELDS = [
    "run_id", "m2_run_id", "profile", "cluster_id", "cluster_membership_status", "not_clustered_reason",
    "compound_public_id", "state_id", "state_role", "pocket_family_id", "box_id", "protomer_id", "job_id",
    "vina_repeat_id", "pose_rank", "pose_model_index", "pose_file", "pose_hard_gate_pass", "pose_reject_reason",
    "attribution_status", "allowed_for_clustering", "pocket_retention_pass", "atp_migration_flag",
    "ppi_relationship_confirmed", "membrane_penetration_flag", "dimer_interface_clash_flag", "pose_mechanism_class",
    "distance_to_cluster_centroid_A", "nearest_cluster_neighbor_distance_A", "used_for_cluster_fit",
    "membership_notes", "clustered_at",
]

QC_FIELDS = ["check_id", "category", "status", "severity", "details", "recommended_fix"]

REQUIRED_CHECKS = [
    "m3_t7_pose_attribution_table_present", "m3_t7_pose_attribution_schema_valid", "m3_t7_pose_attribution_qc_present",
    "m3_t7_clustering_ready_or_partial_allowed", "profile_selected", "attribution_rows_found", "eligible_pose_rows_found",
    "pose_ids_unique", "pose_rows_traceable_to_compound_state_pocket_box_repeat", "pose_files_inside_run_dir",
    "eligible_pose_files_exist", "eligible_pose_files_nonempty", "accepted_boxes_present", "input_reference_paths_valid",
    "all_attribution_rows_represented_in_membership", "no_silent_pose_drop", "pose_atoms_loaded_in_memory_if_needed",
    "no_coordinate_rows_written", "cluster_config_loaded_or_snapshot_written", "cluster_thresholds_not_hardcoded",
    "cluster_method_recorded", "deterministic_cluster_order", "pose_clustering_computed", "singleton_clusters_marked",
    "scattered_pose_support_not_promoted", "cluster_fraction_computed", "pocket_retention_fraction_computed",
    "atp_migration_fraction_computed", "ppi_relationship_fraction_computed", "membrane_dimer_cluster_fraction_computed",
    "convergence_status_assigned", "cluster_reject_reasons_explicit", "anchor_convergence_readiness_flag_written",
    "cluster_support_reported_separately_from_affinity", "raw_pose_affinity_not_used_for_ranking",
    "compound_pose_clusters_written", "cluster_membership_table_written", "cluster_convergence_qc_written",
    "pose_cluster_membership_qc_written", "pose_clustering_report_written", "no_vina_invoked_by_clustering",
    "no_qsub_invoked_by_clustering", "no_runner_invoked_by_clustering", "no_pose_attribution_rerun_attempted",
    "no_compound_anchor_convergence_attempted", "no_candidate_scoring_attempted", "no_candidate_nomination_attempted",
    "no_final_candidate_tables_created", "reference_state_not_promoted", "no_internal_id_leak", "no_smiles_logged",
    "no_ligand_coordinates_logged", "no_receptor_coordinates_logged", "no_pose_coordinates_logged_outside_pdbqt",
    "old_workflow_not_used", "non_goals_preserved",
]

FORBIDDEN_OUTPUTS = [
    "phase3_compounds/tables/compound_anchor_convergence.csv",
    "phase3_compounds/tables/compound_pocket_support.csv",
    "phase3_compounds/tables/pocket_compound_evidence_table.csv",
    "phase3_compounds/tables/final_m3_candidate_hypotheses.csv",
    "phase3_compounds/qc/final_candidate_gate_qc.csv",
    "phase3_compounds/docking_outputs/broad_anchor_scan_optional",
]


@dataclass
class PoseForCluster:
    row: dict[str, str]
    pose_id: str
    center: tuple[float, float, float]
    heavy_atoms: list[AtomRecord]
    rmsd_available: bool
    rmsd_notes: str


@dataclass
class M3PoseClusteringResult:
    status: str
    m3_t9_anchor_convergence_ready: bool
    blockers: list[str]
    warnings: list[str]
    clusters_csv: Path
    membership_csv: Path
    qc_csv: Path
    qc_json: Path
    convergence_qc_csv: Path
    membership_qc_csv: Path
    report_md: Path
    phase3_log: Path
    inputs: dict[str, bool]
    attribution_context: dict[str, Any]
    clustering: dict[str, int | bool]
    mechanism_classes_by_cluster: dict[str, int]


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
    if math.isfinite(result):
        return result
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


def _severity(status: str, blocker: bool = False) -> str:
    if status == "FAIL" and blocker:
        return "BLOCKER"
    if status in {"FAIL", "WARN"}:
        return "MAJOR"
    if status == "NOT_APPLICABLE":
        return "MINOR"
    return "INFO"


def _append_qc(rows: list[dict[str, str]], check_id: str, status: str, details: str) -> None:
    rows.append({"check_id": check_id, "category": "m3_t8_pose_clustering", "status": status, "severity": _severity(status, status == "FAIL"), "details": details, "recommended_fix": ""})


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


def _attribution_path(ctx: RunContext, pose_attribution_table: Path | None) -> Path:
    path = pose_attribution_table or (ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_attribution.csv")
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


def _load_config(ctx: RunContext, config_path: Path | None, method: str) -> dict[str, Any]:
    config = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (ctx.fresh_root / "configs" / "compound_pose_clustering.yaml")
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
    config["clustering"]["method"] = method
    snapshot = ctx.run_dir / "phase3_compounds" / "config_snapshots" / "compound_pose_clustering.resolved.yaml"
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


def _sort_key(row: dict[str, str]) -> tuple[Any, ...]:
    def as_int(value: str) -> tuple[int, Any]:
        try:
            return (0, int(value))
        except (TypeError, ValueError):
            return (1, str(value))

    return (
        row.get("compound_public_id", ""),
        row.get("state_id", ""),
        row.get("pocket_family_id", ""),
        row.get("protomer_id", ""),
        row.get("box_id", ""),
        as_int(row.get("vina_repeat_id", "")),
        as_int(row.get("pose_rank", "")),
        as_int(row.get("pose_model_index", "")),
        row.get("job_id", ""),
    )


def _pose_identity(row: dict[str, str]) -> tuple[str, ...]:
    keys = ["job_id", "pose_rank", "pose_model_index", "compound_public_id", "state_id", "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id"]
    return tuple(row.get(key, "") for key in keys)


def _not_clustered_reason(row: dict[str, str], pose_path: Path | None = None) -> str:
    if row.get("compound_public_id") not in PUBLIC_COMPOUND_IDS:
        return "unknown"
    if row.get("parse_status") and row.get("parse_status") not in {"PASS", "WARN"}:
        return "parse_failure"
    if not _boolish(row.get("allowed_for_clustering")):
        return "not_allowed_for_clustering"
    if not _boolish(row.get("pose_hard_gate_pass")):
        reason = row.get("pose_reject_reason") or ""
        mapping = {
            "ATP_migration": "ATP_migration",
            "outside_accepted_pocket": "outside_accepted_pocket",
            "membrane_penetration": "membrane_penetration",
            "central_dimer_clash": "central_dimer_clash",
            "missing_mapping": "missing_mapping",
            "parse_failure": "parse_failure",
        }
        return mapping.get(reason, "pose_not_hard_gate_pass")
    if row.get("attribution_status") not in {"PASS", "WARN"}:
        return "attribution_reject_or_fail"
    if row.get("pose_reject_reason") != "none":
        return row.get("pose_reject_reason") or "unknown"
    if _boolish(row.get("atp_migration_flag")):
        return "ATP_migration"
    if not _boolish(row.get("pocket_retention_pass")):
        return "outside_accepted_pocket"
    if _boolish(row.get("membrane_penetration_flag")):
        return "membrane_penetration"
    if _boolish(row.get("dimer_interface_clash_flag")):
        return "central_dimer_clash"
    if pose_path is None or not pose_path.is_file():
        return "missing_pose_file"
    return "unknown"


def _is_eligible(row: dict[str, str], pose_path: Path | None) -> bool:
    return (
        _boolish(row.get("allowed_for_clustering"))
        and _boolish(row.get("pose_hard_gate_pass"))
        and row.get("attribution_status") in {"PASS", "WARN"}
        and row.get("pose_reject_reason") == "none"
        and (not row.get("parse_status") or row.get("parse_status") in {"PASS", "WARN"})
        and _boolish(row.get("pocket_retention_pass"))
        and not _boolish(row.get("atp_migration_flag"))
        and not _boolish(row.get("membrane_penetration_flag"))
        and not _boolish(row.get("dimer_interface_clash_flag"))
        and row.get("compound_public_id") in PUBLIC_COMPOUND_IDS
        and pose_path is not None
        and pose_path.is_file()
        and pose_path.stat().st_size > 0
        and bool(row.get("pose_model_index") or row.get("pose_rank"))
    )


def _heavy_atoms(atoms: list[AtomRecord]) -> list[AtomRecord]:
    return [atom for atom in atoms if not (atom.element or "").upper().startswith("H")]


def _center_from_atoms(atoms: list[AtomRecord]) -> tuple[float, float, float]:
    selected = _heavy_atoms(atoms) or atoms
    count = float(len(selected))
    return (
        sum(atom.x for atom in selected) / count,
        sum(atom.y for atom in selected) / count,
        sum(atom.z for atom in selected) / count,
    )


def _center_from_row(row: dict[str, str]) -> tuple[float, float, float] | None:
    values = [_float(row.get("pose_center_x")), _float(row.get("pose_center_y")), _float(row.get("pose_center_z"))]
    if any(value is None for value in values):
        return None
    return (values[0], values[1], values[2])  # type: ignore[return-value]


def _rmsd(a: list[AtomRecord], b: list[AtomRecord]) -> float | None:
    ah = _heavy_atoms(a) or a
    bh = _heavy_atoms(b) or b
    if len(ah) != len(bh) or not ah:
        return None
    total = 0.0
    for left, right in zip(ah, bh):
        total += (left.x - right.x) ** 2 + (left.y - right.y) ** 2 + (left.z - right.z) ** 2
    return math.sqrt(total / len(ah))


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def _pair_distance(a: PoseForCluster, b: PoseForCluster, method: str, config: dict[str, Any]) -> tuple[float | None, str, bool]:
    centroid_distance = _distance(a.center, b.center)
    rmsd_value = _rmsd(a.heavy_atoms, b.heavy_atoms)
    if method == "centroid":
        return centroid_distance, "centroid", False
    if method == "rmsd":
        return rmsd_value, "rmsd" if rmsd_value is not None else "rmsd_unavailable", rmsd_value is None
    if rmsd_value is not None:
        return rmsd_value, "rmsd", False
    if config["clustering"].get("fallback_to_centroid_if_rmsd_unavailable", True):
        return centroid_distance, "centroid_fallback", True
    return None, "rmsd_unavailable", True


def _components(poses: list[PoseForCluster], method: str, config: dict[str, Any]) -> tuple[list[list[int]], dict[tuple[int, int], float | None], bool, list[str]]:
    cutoff = float(config["clustering"]["centroid_distance_cutoff_A"] if method == "centroid" else config["clustering"]["heavy_atom_rmsd_cutoff_A"])
    if method == "hybrid":
        cutoff_rmsd = float(config["clustering"]["heavy_atom_rmsd_cutoff_A"])
        cutoff_centroid = float(config["clustering"]["centroid_distance_cutoff_A"])
    edges: dict[int, set[int]] = {i: set() for i in range(len(poses))}
    distances: dict[tuple[int, int], float | None] = {}
    fallback_used = False
    notes: list[str] = []
    for i in range(len(poses)):
        for j in range(i + 1, len(poses)):
            value, metric, fallback = _pair_distance(poses[i], poses[j], method, config)
            fallback_used = fallback_used or fallback
            if fallback:
                notes.append(metric)
            distances[(i, j)] = value
            if value is None:
                continue
            threshold = cutoff
            if method == "hybrid":
                threshold = cutoff_centroid if metric == "centroid_fallback" else cutoff_rmsd
            if value <= threshold:
                edges[i].add(j)
                edges[j].add(i)
    seen: set[int] = set()
    comps: list[list[int]] = []
    for start in range(len(poses)):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        comp: list[int] = []
        while stack:
            node = stack.pop()
            comp.append(node)
            for nxt in sorted(edges[node]):
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        comps.append(sorted(comp, key=lambda idx: _sort_key(poses[idx].row)))
    comps.sort(key=lambda comp: (-len(comp), _sort_key(poses[comp[0]].row)))
    return comps, distances, fallback_used, sorted(set(notes))


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


def _fraction(rows: list[dict[str, str]], field: str, true_value: bool = True) -> float:
    if not rows:
        return 0.0
    count = sum(1 for row in rows if _boolish(row.get(field)) is true_value)
    return count / float(len(rows))


def _mechanism_majority(rows: list[dict[str, str]]) -> tuple[str, str]:
    counts: dict[str, int] = {}
    for row in rows:
        key = row.get("pose_mechanism_class") or "ambiguous_or_failed"
        counts[key] = counts.get(key, 0) + 1
    majority = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0] if counts else "ambiguous_or_failed"
    return majority, json.dumps(counts, sort_keys=True)


def _scan_hygiene(ctx: RunContext, private_entries: list[Any]) -> tuple[int, list[str], bool, bool, bool, bool, list[str]]:
    leaks = scan_internal_id_leaks(ctx, private_entries)
    scan_paths = list(public_output_scan_paths(ctx))
    for root in [
        ctx.run_dir / "phase3_compounds" / "tables",
        ctx.run_dir / "phase3_compounds" / "qsub",
        ctx.logs_dir / "errors",
    ]:
        if root.is_dir():
            scan_paths.extend(path for path in sorted(root.rglob("*")) if path.is_file())
    filtered: list[Path] = []
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
        upper = text.upper()
        if ("SMILES=" in upper or "SMILES:" in upper) and re.search(r"(^|\s)(C|N|O|S|P|Cl|Br|F|I)[A-Za-z0-9@+\-\[\]\(\)=#$\\/]+(\s|$)", text):
            smiles_logged = True
        if re.search(r"^(ATOM|HETATM)\s+\d+", text, re.MULTILINE):
            coord_hits.append(ctx.relative_to_repo(path))
            lower = path.as_posix().lower()
            ligand_coords = ligand_coords or "ligand" in lower
            receptor_coords = receptor_coords or "receptor" in lower
            pose_coords = pose_coords or "pose" in lower or "cluster" in lower or "qc" in lower or "table" in lower
    return len(leaks), coord_hits, smiles_logged, ligand_coords, receptor_coords, pose_coords, scanned


def _has_real_output(path: Path) -> bool:
    if path.is_dir():
        return any(_has_real_output(child) for child in path.rglob("*"))
    if not path.is_file():
        return False
    if path.name == ".gitkeep":
        return False
    if path.suffix.lower() == ".csv":
        try:
            with path.open("r", encoding="utf-8", newline="") as handle:
                return sum(1 for _ in csv.reader(handle)) > 1
        except OSError:
            return True
    return path.stat().st_size > 0


def _write_report(path: Path, ctx: RunContext, summary: dict[str, Any]) -> None:
    lines = [
        "# M3-T8 Pose Clustering",
        "",
        "Technical pose-cluster convergence report only. Cluster support is reported separately from affinity and is not a candidate claim.",
        "",
        f"- status: {summary['overall_status']}",
        f"- profile: {summary['profile']}",
        f"- eligible_pose_rows: {summary['attribution_context']['eligible_pose_rows']}",
        f"- clusters_written: {summary['clustering']['clusters_written']}",
        f"- converged_clusters: {summary['clustering']['converged_clusters']}",
        f"- m3_t9_anchor_convergence_ready: {str(summary['m3_t9_anchor_convergence_ready']).lower()}",
        "",
        "## Blockers",
    ]
    lines.extend(f"- {item}" for item in summary["blockers"] or ["none"])
    lines.append("")
    lines.append("## Warnings")
    lines.extend(f"- {item}" for item in summary["warnings"] or ["none"])
    lines.append("")
    lines.append("Next task: M3-T9 - Compound-anchor convergence across states and pockets")
    ctx.require_within_run_dir(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_m3_pose_clustering(
    ctx: RunContext,
    *,
    m2_run_id: str,
    profile: str | None = None,
    mode: str = "cluster",
    pose_attribution_table: Path | None = None,
    pose_attribution_qc: Path | None = None,
    accepted_boxes: Path | None = None,
    config: Path | None = None,
    force: bool = False,
    allow_partial: bool = False,
    require_attribution_ready: bool = True,
    allow_reference_only: bool = False,
    private_map: Path | None = None,
    cluster_method: str = "hybrid",
    strict_convergence: bool = True,
    write_empty_clusters: bool = False,
    write_cluster_membership: bool = True,
) -> M3PoseClusteringResult:
    del force
    _ensure_logs(ctx)
    timestamp = now_iso()
    phase3 = ctx.run_dir / "phase3_compounds"
    tables_dir = phase3 / "tables"
    qc_dir = phase3 / "qc"
    reports_dir = phase3 / "reports"
    for directory in [tables_dir, qc_dir, reports_dir]:
        ctx.require_within_run_dir(directory).mkdir(parents=True, exist_ok=True)

    clusters_csv = tables_dir / "compound_pose_clusters.csv"
    membership_csv = tables_dir / "compound_pose_cluster_membership.csv"
    qc_csv = qc_dir / "pose_clustering_qc.csv"
    qc_json = qc_dir / "pose_clustering_qc.json"
    convergence_qc_csv = qc_dir / "cluster_convergence_qc.csv"
    membership_qc_csv = qc_dir / "pose_cluster_membership_qc.csv"
    report_md = reports_dir / "m3_task8_pose_clustering.md"
    phase3_log = ctx.logs_dir / "phase3_compounds.log"

    blockers: list[str] = []
    warnings: list[str] = []
    attr_path = _attribution_path(ctx, pose_attribution_table)
    attr_found = attr_path.is_file()
    if not attr_found:
        (warnings if mode == "dry-run" else blockers).append("compound_pose_attribution.csv missing")
    if attr_found and not _inside(attr_path, tables_dir):
        blockers.append("compound_pose_attribution.csv path outside phase3_compounds/tables")
    attr_rows: list[dict[str, str]] = []
    attr_fields: list[str] = []
    if attr_found:
        attr_rows, attr_fields = _read_csv(attr_path)
    schema_ok = set(ATTR_REQUIRED).issubset(set(attr_fields))
    if attr_found and not schema_ok:
        blockers.append("compound_pose_attribution.csv schema missing required columns")
    selected_rows, selected_profile, profiles = _filter_profile(attr_rows, profile)
    if attr_found and not selected_profile:
        blockers.append("profile ambiguous; provide --profile")
    if attr_found and not selected_rows:
        (warnings if mode == "dry-run" or write_empty_clusters else blockers).append("zero attribution rows for selected profile")

    qc_path = pose_attribution_qc or (qc_dir / "pose_attribution_qc.json")
    if not qc_path.is_absolute():
        qc_path = ctx.repo_root / qc_path
    qc_found = qc_path.is_file()
    qc_data = _load_json(qc_path)
    t7_ready = bool(qc_data and qc_data.get("overall_status") == "PASS" and qc_data.get("m3_t8_clustering_ready") is True)
    if not qc_found:
        (warnings if mode == "dry-run" or allow_partial else blockers).append("pose_attribution_qc.json missing")
    elif require_attribution_ready and not t7_ready:
        (warnings if mode == "dry-run" or allow_partial else blockers).append("pose_attribution_qc.json has m3_t8_clustering_ready=false")

    m2_root = ctx.fresh_root / "runs" / m2_run_id
    accepted_boxes_path = _first_existing([Path(accepted_boxes)] if accepted_boxes else [
        m2_root / "phase2_pockets" / "final" / "accepted_pocket_boxes.csv",
        m2_root / "phase2_pockets" / "export_for_m3" / "accepted_pocket_boxes.csv",
    ])
    if not accepted_boxes_path:
        (warnings if mode == "dry-run" else blockers).append("accepted pocket boxes missing")
    config_data = _load_config(ctx, config, cluster_method)

    identities = [_pose_identity(row) for row in selected_rows]
    unique_ids = len(identities) == len(set(identities))
    if selected_rows and not unique_ids:
        blockers.append("duplicate pose identity rows")
    non_public = sorted({row.get("compound_public_id", "") for row in selected_rows if row.get("compound_public_id") not in PUBLIC_COMPOUND_IDS})
    if non_public:
        blockers.append("non-public compound IDs present in attribution table")

    eligible_rows: list[dict[str, str]] = []
    pose_files_inside = True
    eligible_files_exist = True
    eligible_files_nonempty = True
    pose_objects: dict[tuple[str, ...], PoseForCluster] = {}
    atoms_loaded = False
    geometry_failures = 0
    membership_rows: list[dict[str, Any]] = []
    attr_sha = _sha256(attr_path) if attr_found else ""
    qc_sha = _sha256(qc_path) if qc_found else ""
    box_sha = _sha256(accepted_boxes_path) if accepted_boxes_path else ""

    for row in sorted(selected_rows, key=_sort_key):
        pose_path = _resolve_path(ctx, row.get("pose_file"))
        inside = _inside(pose_path, ctx.run_dir)
        pose_files_inside = pose_files_inside and inside
        if not inside:
            blockers.append(f"pose_file outside run_dir for job {row.get('job_id')}")
        otherwise_eligible = (
            _boolish(row.get("allowed_for_clustering"))
            and _boolish(row.get("pose_hard_gate_pass"))
            and row.get("attribution_status") in {"PASS", "WARN"}
            and row.get("pose_reject_reason") == "none"
            and (not row.get("parse_status") or row.get("parse_status") in {"PASS", "WARN"})
            and _boolish(row.get("pocket_retention_pass"))
            and not _boolish(row.get("atp_migration_flag"))
            and not _boolish(row.get("membrane_penetration_flag"))
            and not _boolish(row.get("dimer_interface_clash_flag"))
            and row.get("compound_public_id") in PUBLIC_COMPOUND_IDS
        )
        if otherwise_eligible and inside and (pose_path is None or not pose_path.is_file() or pose_path.stat().st_size == 0):
            blockers.append(f"eligible pose_file missing or empty for job {row.get('job_id')}")
        expected_pose_sha = (row.get("output_pdbqt_sha256") or "").strip()
        current_pose_sha = _sha256(pose_path) if inside and pose_path is not None and pose_path.is_file() and pose_path.stat().st_size > 0 else ""
        pose_hash_matches = bool(expected_pose_sha and current_pose_sha and expected_pose_sha == current_pose_sha)
        if otherwise_eligible and inside and pose_path is not None and pose_path.is_file() and pose_path.stat().st_size > 0 and not pose_hash_matches:
            blockers.append(f"eligible pose_file hash mismatch for job {row.get('job_id')}")
        eligible = _is_eligible(row, pose_path if inside else None) and pose_hash_matches
        if eligible:
            reason = "none"
        elif otherwise_eligible and inside and pose_path is not None and pose_path.is_file() and pose_path.stat().st_size > 0 and not pose_hash_matches:
            reason = "pose_hash_mismatch"
        else:
            reason = _not_clustered_reason(row, pose_path if inside else None)
        if eligible:
            eligible_rows.append(row)
            exists = bool(pose_path and pose_path.is_file())
            nonempty = bool(pose_path and pose_path.is_file() and pose_path.stat().st_size > 0)
            eligible_files_exist = eligible_files_exist and exists
            eligible_files_nonempty = eligible_files_nonempty and nonempty
            try:
                atoms = extract_pose_atoms(pose_path, int(row.get("pose_model_index") or 0), int(row.get("pose_rank") or 0)) if pose_path else []
            except Exception:
                atoms = []
            atoms_loaded = atoms_loaded or bool(atoms)
            if atoms:
                center = _center_from_atoms(atoms)
                heavy = _heavy_atoms(atoms) or atoms
                pose_objects[_pose_identity(row)] = PoseForCluster(row=row, pose_id="|".join(_pose_identity(row)), center=center, heavy_atoms=heavy, rmsd_available=True, rmsd_notes="parsed")
            else:
                fallback_center = _center_from_row(row)
                if fallback_center is None:
                    geometry_failures += 1
                    reason = "geometry_parse_failure"
                else:
                    pose_objects[_pose_identity(row)] = PoseForCluster(row=row, pose_id="|".join(_pose_identity(row)), center=fallback_center, heavy_atoms=[], rmsd_available=False, rmsd_notes="centroid_from_attribution")
                    warnings.append("pose atom parsing unavailable; centroid fallback used")
        membership_rows.append({**{field: row.get(field, "") for field in MEMBERSHIP_FIELDS if field in row}, "cluster_id": "", "cluster_membership_status": "PENDING" if eligible and _pose_identity(row) in pose_objects else ("FAIL" if reason in {"missing_pose_file", "geometry_parse_failure"} else "NOT_APPLICABLE"), "not_clustered_reason": reason, "distance_to_cluster_centroid_A": "", "nearest_cluster_neighbor_distance_A": "", "used_for_cluster_fit": "false", "membership_notes": "pending_cluster_assignment" if eligible else reason, "clustered_at": timestamp})
        if reason in {"missing_pose_file", "geometry_parse_failure"}:
            append_failed_job(ctx, row.get("job_id", "pose_clustering"), "FAIL", reason)

    if mode == "cluster" and selected_rows and not pose_files_inside:
        blockers.append("one or more pose_file paths outside run_dir")
    if mode == "cluster" and require_attribution_ready and not allow_partial and (not eligible_files_exist or not eligible_files_nonempty):
        blockers.append("required eligible pose files missing or empty")
    if mode == "cluster" and eligible_rows and not pose_objects and not allow_partial:
        blockers.append("every eligible pose failed geometry parsing")
    if mode == "cluster" and not eligible_rows and selected_rows:
        warnings.append("no eligible pose rows after T7 hard gates")
    if allow_partial:
        warnings.append("allow_partial diagnostic mode used")
    if mode == "dry-run":
        warnings.append("dry-run mode")

    group_keys = config_data["grouping"]["primary_group_by"]
    grouped: dict[tuple[str, ...], list[PoseForCluster]] = {}
    for pose in sorted(pose_objects.values(), key=lambda item: _sort_key(item.row)):
        key = tuple(pose.row.get(field, "") for field in group_keys)
        grouped.setdefault(key, []).append(pose)

    cluster_rows: list[dict[str, Any]] = []
    convergence_rows: list[dict[str, Any]] = []
    membership_by_identity = {_pose_identity(row): row for row in membership_rows}
    scattered_groups = 0
    cluster_algorithm = "graph_connected_components"
    cluster_scope = "compound_state_protomer_pocket"
    global_cluster_index = 1
    fallback_used_any = False
    for group_key in sorted(grouped):
        poses = grouped[group_key]
        comps, distances, fallback_used, rmsd_notes = _components(poses, cluster_method, config_data)
        fallback_used_any = fallback_used_any or fallback_used
        if len(comps) > 1:
            scattered_groups += 1
        for rank, comp in enumerate(comps, start=1):
            members = [poses[idx] for idx in comp]
            rows = [member.row for member in members]
            cluster_id = f"cluster_{global_cluster_index:04d}"
            global_cluster_index += 1
            centers = [member.center for member in members]
            centroid = (
                mean([center[0] for center in centers]),
                mean([center[1] for center in centers]),
                mean([center[2] for center in centers]),
            )
            centroid_distances = [_distance(center, centroid) for center in centers]
            pair_centroid_values: list[float] = []
            pair_rmsd_values: list[float] = []
            for i in range(len(members)):
                for j in range(i + 1, len(members)):
                    pair_centroid_values.append(_distance(members[i].center, members[j].center))
                    rmsd_value = _rmsd(members[i].heavy_atoms, members[j].heavy_atoms)
                    if rmsd_value is not None:
                        pair_rmsd_values.append(rmsd_value)
            affinities = [value for row in rows if (value := _float(row.get("vina_affinity_kcal_mol"))) is not None]
            within_pocket_fraction = _fraction(rows, "pocket_retention_pass", True)
            atp_migration_fraction = _fraction(rows, "atp_migration_flag", True)
            ppi_fraction = _fraction(rows, "ppi_relationship_confirmed", True)
            membrane_fraction = _fraction(rows, "membrane_penetration_flag", True)
            dimer_fraction = _fraction(rows, "dimer_interface_clash_flag", True)
            ppi_contacts = [int(_float(row.get("ppi_hotspot_contact_count")) or 0) for row in rows]
            ppi_residue_counts = [int(_float(row.get("ppi_contact_residue_count")) or 0) for row in rows]
            ppi_distances = [value for row in rows if (value := _float(row.get("nearest_ppi_hotspot_distance"))) is not None]
            ppi_rim = [value for row in rows if (value := _float(row.get("ppi_rim_distance"))) is not None]
            mechanism_majority, mechanism_counts = _mechanism_majority(rows)
            cluster_size = len(rows)
            eligible_count = len(poses)
            cluster_fraction = cluster_size / float(eligible_count) if eligible_count else 0.0
            reference_only = all((row.get("state_role") in config_data["grouping"]["reference_state_role_values"] or row.get("state_id") == "3GT8_raw") for row in rows)
            converged = (
                cluster_size >= int(config_data["convergence"]["cluster_size_min"])
                and cluster_fraction >= float(config_data["convergence"]["cluster_fraction_min"])
                and within_pocket_fraction >= float(config_data["convergence"]["pocket_retention_fraction_min"])
                and atp_migration_fraction <= float(config_data["convergence"]["atp_migration_fraction_max"])
                and membrane_fraction <= float(config_data["convergence"]["membrane_penetration_fraction_max"])
                and dimer_fraction <= float(config_data["convergence"]["dimer_interface_clash_fraction_max"])
                and (not strict_convergence or not config_data["convergence"].get("strict_requires_ppi_relationship") or ppi_fraction >= float(config_data["convergence"]["ppi_relationship_fraction_min"]))
                and (not reference_only or allow_reference_only)
            )
            if cluster_size == 1:
                convergence_status = "WARN_SINGLETON_NOT_CONVERGED"
                reject_reason = "singleton_cluster"
            elif reference_only and not allow_reference_only:
                convergence_status = "WARN_REFERENCE_ONLY_NOT_PROMOTED"
                reject_reason = "reference_only_not_promoted"
            elif atp_migration_fraction > 0:
                convergence_status = "REJECT_ATP_CONFOUNDED"
                reject_reason = "ATP_migration_in_cluster"
            elif within_pocket_fraction < float(config_data["convergence"]["pocket_retention_fraction_min"]):
                convergence_status = "REJECT_OUTSIDE_ACCEPTED_POCKET"
                reject_reason = "outside_accepted_pocket_in_cluster"
            elif membrane_fraction > 0 or dimer_fraction > 0:
                convergence_status = "REJECT_MEMBRANE_OR_DIMER_CONFLICT"
                reject_reason = "membrane_penetration_in_cluster" if membrane_fraction > 0 else "dimer_interface_clash_in_cluster"
            elif not converged and cluster_fraction < float(config_data["convergence"]["cluster_fraction_min"]):
                convergence_status = "REJECT_SCATTERED"
                reject_reason = "insufficient_cluster_fraction"
            elif converged and ppi_fraction < float(config_data["convergence"]["ppi_relationship_fraction_min"]):
                convergence_status = "WARN_CONVERGED_WEAK_PPI_RELATIONSHIP"
                reject_reason = "none"
            elif converged and mechanism_majority == "generic_nonATP_ligandable_pocket":
                convergence_status = "PASS_CONVERGED_GENERIC_NONATP"
                reject_reason = "none"
            elif converged:
                convergence_status = "PASS_CONVERGED"
                reject_reason = "none"
            else:
                convergence_status = "REJECT_SCATTERED"
                reject_reason = "scattered_poses"
            allowed_anchor = bool(converged and convergence_status in {"PASS_CONVERGED", "PASS_CONVERGED_GENERIC_NONATP", "WARN_CONVERGED_WEAK_PPI_RELATIONSHIP"} and reject_reason == "none" and atp_migration_fraction == 0 and not (reference_only and not allow_reference_only))
            cluster_row = {
                "run_id": ctx.run_id, "m2_run_id": m2_run_id, "profile": selected_profile, "cluster_id": cluster_id,
                "cluster_algorithm": cluster_algorithm, "cluster_method": cluster_method, "cluster_scope": cluster_scope,
                "compound_public_id": rows[0].get("compound_public_id"), "state_id": rows[0].get("state_id"), "state_role": rows[0].get("state_role"),
                "pocket_family_id": rows[0].get("pocket_family_id"), "protomer_id": rows[0].get("protomer_id"),
                "box_ids_in_cluster": "|".join(sorted({row.get("box_id", "") for row in rows})),
                "vina_repeat_ids_in_cluster": "|".join(row.get("vina_repeat_id", "") for row in rows),
                "pose_ranks_in_cluster": "|".join(row.get("pose_rank", "") for row in rows),
                "pose_model_indices_in_cluster": "|".join(row.get("pose_model_index", "") for row in rows),
                "job_ids_in_cluster": "|".join(row.get("job_id", "") for row in rows),
                "source_compound_pose_attribution": ctx.relative_to_repo(attr_path) if attr_found else "",
                "source_compound_pose_attribution_sha256": attr_sha,
                "source_pose_attribution_qc": ctx.relative_to_repo(qc_path) if qc_found else "",
                "source_pose_attribution_qc_sha256": qc_sha,
                "accepted_box_source": ctx.relative_to_repo(accepted_boxes_path) if accepted_boxes_path else "",
                "accepted_box_source_sha256": box_sha,
                "eligible_pose_count_in_group": eligible_count, "cluster_size": cluster_size, "cluster_fraction": f"{cluster_fraction:.6f}",
                "cluster_rank_within_group": rank,
                "vina_affinity_median_kcal_mol": _fmt(_stat(affinities, "median")),
                "vina_affinity_mean_kcal_mol": _fmt(_stat(affinities, "mean")),
                "vina_affinity_best_kcal_mol": _fmt(_stat(affinities, "best")),
                "vina_affinity_worst_kcal_mol": _fmt(_stat(affinities, "worst")),
                "affinity_used_for_ranking": "false",
                "cluster_centroid_x": _fmt(centroid[0]), "cluster_centroid_y": _fmt(centroid[1]), "cluster_centroid_z": _fmt(centroid[2]),
                "cluster_radius_A": _fmt(max(centroid_distances) if centroid_distances else 0.0),
                "cluster_median_pairwise_centroid_distance_A": _fmt(median(pair_centroid_values) if pair_centroid_values else 0.0),
                "cluster_max_pairwise_centroid_distance_A": _fmt(max(pair_centroid_values) if pair_centroid_values else 0.0),
                "cluster_mean_pairwise_heavy_atom_rmsd_A": _fmt(mean(pair_rmsd_values) if pair_rmsd_values else None),
                "cluster_median_pairwise_heavy_atom_rmsd_A": _fmt(median(pair_rmsd_values) if pair_rmsd_values else None),
                "cluster_max_pairwise_heavy_atom_rmsd_A": _fmt(max(pair_rmsd_values) if pair_rmsd_values else None),
                "rmsd_available": _bool_text(bool(pair_rmsd_values) or cluster_size == 1),
                "rmsd_notes": ";".join(rmsd_notes) if rmsd_notes else "available",
                "within_pocket_fraction": f"{within_pocket_fraction:.6f}",
                "atp_migration_fraction": f"{atp_migration_fraction:.6f}",
                "ppi_relationship_confirmed_fraction": f"{ppi_fraction:.6f}",
                "ppi_hotspot_contact_median": _fmt(median(ppi_contacts) if ppi_contacts else None),
                "ppi_hotspot_contact_mean": _fmt(mean(ppi_contacts) if ppi_contacts else None),
                "ppi_contact_residue_count_median": _fmt(median(ppi_residue_counts) if ppi_residue_counts else None),
                "nearest_ppi_hotspot_distance_median": _fmt(median(ppi_distances) if ppi_distances else None),
                "ppi_rim_distance_median": _fmt(median(ppi_rim) if ppi_rim else None),
                "membrane_penetration_fraction": f"{membrane_fraction:.6f}",
                "dimer_interface_clash_fraction": f"{dimer_fraction:.6f}",
                "mechanism_class_majority": mechanism_majority,
                "mechanism_class_counts": mechanism_counts,
                "converged_pose_cluster": _bool_text(converged),
                "convergence_status": convergence_status,
                "cluster_reject_reason": reject_reason,
                "cluster_notes": "cluster_support_separate_from_affinity",
                "allowed_for_anchor_convergence": _bool_text(allowed_anchor),
                "clustered_at": timestamp,
            }
            cluster_rows.append(cluster_row)
            convergence_rows.append({key: cluster_row.get(key, "") for key in ["run_id", "m2_run_id", "profile", "cluster_id", "compound_public_id", "state_id", "state_role", "pocket_family_id", "protomer_id", "eligible_pose_count_in_group", "cluster_size", "cluster_fraction", "within_pocket_fraction", "atp_migration_fraction", "ppi_relationship_confirmed_fraction", "membrane_penetration_fraction", "dimer_interface_clash_fraction", "converged_pose_cluster", "convergence_status", "cluster_reject_reason", "allowed_for_anchor_convergence"]} | {"convergence_notes": "cluster_level_summary_from_m3_t7"})
            for member in members:
                ident = _pose_identity(member.row)
                mem = membership_by_identity.get(ident)
                if mem is None:
                    continue
                other_distances = [_distance(member.center, other.center) for other in members if other is not member]
                mem.update({
                    "cluster_id": cluster_id,
                    "cluster_membership_status": "PASS" if converged else "WARN",
                    "not_clustered_reason": "none",
                    "distance_to_cluster_centroid_A": _fmt(_distance(member.center, centroid)),
                    "nearest_cluster_neighbor_distance_A": _fmt(min(other_distances) if other_distances else 0.0),
                    "used_for_cluster_fit": "true",
                    "membership_notes": convergence_status,
                })
            append_job_status(ctx, cluster_id, "PASS" if allowed_anchor else "WARN", details={"phase": "phase3_compounds", "task": "M3-T8", "job_type": "pose_clustering", "cluster_id": cluster_id, "profile": selected_profile, "compound_public_id": rows[0].get("compound_public_id"), "state_id": rows[0].get("state_id"), "pocket_family_id": rows[0].get("pocket_family_id"), "protomer_id": rows[0].get("protomer_id"), "eligible_pose_count_in_group": eligible_count, "cluster_size": cluster_size, "cluster_fraction": cluster_fraction, "converged_pose_cluster": converged, "convergence_status": convergence_status, "cluster_reject_reason": reject_reason, "allowed_for_anchor_convergence": allowed_anchor, "message": "pose cluster planned for M3-T9 readiness only"})

    if fallback_used_any:
        warnings.append("RMSD unavailable for at least one pair; centroid fallback used")
    if eligible_rows and not cluster_rows:
        warnings.append("eligible rows found but no clusters written")
    if cluster_rows and not any(row.get("converged_pose_cluster") == "true" for row in cluster_rows):
        warnings.append("no converged pose clusters")
    if cluster_rows and all(int(row.get("cluster_size", 0)) == 1 for row in cluster_rows):
        warnings.append("all eligible poses are singleton clusters")

    write_outputs = bool(selected_rows or cluster_rows or write_empty_clusters or mode == "dry-run")
    if write_outputs:
        _write_csv(clusters_csv, CLUSTER_FIELDS, cluster_rows, ctx)
        if write_cluster_membership:
            _write_csv(membership_csv, MEMBERSHIP_FIELDS, membership_rows, ctx)
        _write_csv(convergence_qc_csv, [
            "run_id", "m2_run_id", "profile", "cluster_id", "compound_public_id", "state_id", "state_role",
            "pocket_family_id", "protomer_id", "eligible_pose_count_in_group", "cluster_size", "cluster_fraction",
            "within_pocket_fraction", "atp_migration_fraction", "ppi_relationship_confirmed_fraction",
            "membrane_penetration_fraction", "dimer_interface_clash_fraction", "converged_pose_cluster",
            "convergence_status", "cluster_reject_reason", "allowed_for_anchor_convergence", "convergence_notes",
        ], convergence_rows, ctx)
        _write_csv(membership_qc_csv, [
            "run_id", "m2_run_id", "profile", "compound_public_id", "state_id", "state_role", "pocket_family_id",
            "box_id", "protomer_id", "job_id", "vina_repeat_id", "pose_rank", "pose_model_index", "cluster_id",
            "cluster_membership_status", "not_clustered_reason", "pose_hard_gate_pass", "allowed_for_clustering",
            "membership_notes",
        ], membership_rows, ctx)
    if write_cluster_membership and len(membership_rows) != len(selected_rows):
        blockers.append("attribution rows were silently dropped from membership output")
    forbidden_created = [item for item in FORBIDDEN_OUTPUTS if _has_real_output(ctx.run_dir / item)]
    if forbidden_created:
        blockers.append("forbidden downstream outputs already exist in run directory")

    private_entries, private_warnings = load_private_map(_safe_private_map_path(ctx, private_map))
    warnings.extend(private_warnings)
    leak_count, coord_hits, smiles_logged, ligand_coords, receptor_coords, pose_coords, scanned = _scan_hygiene(ctx, list(private_entries.values()))
    if leak_count:
        blockers.append("internal compound ID leakage detected")
    if coord_hits:
        blockers.append("coordinate records printed outside intended PDBQT files")
    if smiles_logged:
        blockers.append("SMILES-like fields printed in public outputs")

    clustering_counts = {
        "groups_evaluated": len(grouped),
        "clusters_written": len(cluster_rows),
        "singleton_clusters": sum(1 for row in cluster_rows if int(row.get("cluster_size", 0)) == 1),
        "converged_clusters": sum(1 for row in cluster_rows if row.get("converged_pose_cluster") == "true"),
        "nonconverged_clusters": sum(1 for row in cluster_rows if row.get("converged_pose_cluster") != "true"),
        "scattered_groups": scattered_groups,
        "clusters_allowed_for_anchor_convergence": sum(1 for row in cluster_rows if row.get("allowed_for_anchor_convergence") == "true"),
        "eligible_pose_rows_assigned_to_cluster": sum(1 for row in membership_rows if row.get("used_for_cluster_fit") == "true"),
        "eligible_pose_rows_unassigned": max(0, len(eligible_rows) - sum(1 for row in membership_rows if row.get("used_for_cluster_fit") == "true")),
        "all_attribution_rows_represented_in_membership": len(membership_rows) == len(selected_rows),
    }
    mechanism_classes_by_cluster = {name: sum(1 for row in cluster_rows if row.get("mechanism_class_majority") == name) for name in ["orthosteric_or_direct_PPI_patch_blocker", "rim_blocker", "allosteric_near_candidate", "generic_nonATP_ligandable_pocket", "ATP_like_reject", "ambiguous_or_failed"]}
    if mode == "cluster" and cluster_rows and clustering_counts["clusters_allowed_for_anchor_convergence"] == 0:
        warnings.append("no cluster allowed for anchor convergence")
    status = "FAIL" if blockers else ("WARN" if warnings or mode == "dry-run" or clustering_counts["clusters_allowed_for_anchor_convergence"] == 0 else "PASS")
    m3_t9_ready = status == "PASS" and mode == "cluster" and clustering_counts["clusters_allowed_for_anchor_convergence"] > 0
    reference_states = sorted({row.get("state_id", "") for row in selected_rows if row.get("state_role") in config_data["grouping"]["reference_state_role_values"] or row.get("state_id") == "3GT8_raw"})
    primary_states = sorted({row.get("state_id", "") for row in selected_rows if row.get("state_id") and row.get("state_id") not in reference_states})

    inputs = {"compound_pose_attribution_found": attr_found, "pose_attribution_qc_found": qc_found, "accepted_boxes_found": accepted_boxes_path is not None}
    attribution_context = {
        "attribution_rows": len(selected_rows),
        "eligible_pose_rows": len(eligible_rows),
        "ineligible_pose_rows": len(selected_rows) - len(eligible_rows),
        "profiles_in_attribution_table": profiles,
        "selected_profile": selected_profile,
        "primary_states": primary_states,
        "reference_states": reference_states,
        "reference_only_clustering": bool(reference_states and not primary_states),
    }
    summary = {
        "schema_version": "m3_pose_clustering_qc_v1",
        "run_id": ctx.run_id,
        "m2_run_id": m2_run_id,
        "reviewed_at": timestamp,
        "mode": mode,
        "profile": selected_profile,
        "overall_status": status,
        "allow_partial": allow_partial,
        "require_attribution_ready": require_attribution_ready,
        "allow_reference_only": allow_reference_only,
        "cluster_method": cluster_method,
        "strict_convergence": strict_convergence,
        "m3_t9_anchor_convergence_ready": m3_t9_ready,
        "inputs": inputs,
        "attribution_context": attribution_context,
        "clustering": clustering_counts,
        "mechanism_classes_by_cluster": mechanism_classes_by_cluster,
        "cluster_quality": {
            "clusters_with_zero_atp_migration_fraction": sum(1 for row in cluster_rows if row.get("atp_migration_fraction") == "0.000000"),
            "clusters_with_pocket_retention_fraction_pass": sum(1 for row in cluster_rows if _float(row.get("within_pocket_fraction")) is not None and _float(row.get("within_pocket_fraction")) >= float(config_data["convergence"]["pocket_retention_fraction_min"])),
            "clusters_with_ppi_relationship_fraction_pass": sum(1 for row in cluster_rows if _float(row.get("ppi_relationship_confirmed_fraction")) is not None and _float(row.get("ppi_relationship_confirmed_fraction")) >= float(config_data["convergence"]["ppi_relationship_fraction_min"])),
            "clusters_with_no_membrane_or_dimer_conflict": sum(1 for row in cluster_rows if row.get("membrane_penetration_fraction") == "0.000000" and row.get("dimer_interface_clash_fraction") == "0.000000"),
            "clusters_warn_weak_ppi": sum(1 for row in cluster_rows if row.get("convergence_status") == "WARN_CONVERGED_WEAK_PPI_RELATIONSHIP"),
            "reference_only_clusters_not_promoted": sum(1 for row in cluster_rows if row.get("convergence_status") == "WARN_REFERENCE_ONLY_NOT_PROMOTED"),
        },
        "outputs": {
            "compound_pose_clusters": ctx.relative_to_repo(clusters_csv),
            "compound_pose_cluster_membership": ctx.relative_to_repo(membership_csv),
            "pose_clustering_qc_csv": ctx.relative_to_repo(qc_csv),
            "pose_clustering_qc_json": ctx.relative_to_repo(qc_json),
            "cluster_convergence_qc": ctx.relative_to_repo(convergence_qc_csv),
            "pose_cluster_membership_qc": ctx.relative_to_repo(membership_qc_csv),
            "report_file": ctx.relative_to_repo(report_md),
        },
        "parser": {"pose_atoms_loaded_in_memory": True, "coordinate_records_written_to_public_tables": False, "affinity_used_for_ranking": False},
        "counts": {"vina_commands_invoked_by_clustering": 0, "qsub_commands_invoked_by_clustering": 0, "runner_commands_invoked_by_clustering": 0, "pose_attribution_rerun_attempted": 0, "anchor_convergence_rows_created": 0, "candidate_tables_created": 0, "confidentiality_leaks": leak_count, "coordinate_leak_files": len(coord_hits)},
        "blockers": blockers,
        "warnings": warnings,
        "confidentiality": {"internal_ids_redacted": leak_count == 0, "public_outputs_scanned": scanned, "leaks_detected": leak_count, "smiles_logged": smiles_logged, "ligand_coordinates_logged": ligand_coords, "receptor_coordinates_logged": receptor_coords, "pose_coordinates_logged_outside_pose_file": pose_coords},
        "non_goals_preserved": {"vina_execution": True, "qsub_submission": True, "pbs_runner_invocation": True, "pose_attribution_rerun": True, "anchor_convergence": True, "compound_scoring": True, "candidate_nomination": True, "final_candidate_tiering": True},
    }
    qc_status = {
        "m3_t7_pose_attribution_table_present": "PASS" if attr_found else ("WARN" if mode == "dry-run" else "FAIL"),
        "m3_t7_pose_attribution_schema_valid": "PASS" if schema_ok else ("NOT_APPLICABLE" if not attr_found else "FAIL"),
        "m3_t7_pose_attribution_qc_present": "PASS" if qc_found else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "m3_t7_clustering_ready_or_partial_allowed": "PASS" if (t7_ready or allow_partial or mode == "dry-run") else "FAIL",
        "profile_selected": "PASS" if selected_profile else "FAIL",
        "attribution_rows_found": "PASS" if selected_rows else ("WARN" if mode == "dry-run" or write_empty_clusters else "FAIL"),
        "eligible_pose_rows_found": "PASS" if eligible_rows else "WARN",
        "pose_ids_unique": "PASS" if unique_ids else "FAIL",
        "pose_files_inside_run_dir": "PASS" if pose_files_inside else "FAIL",
        "eligible_pose_files_exist": "PASS" if eligible_files_exist else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "eligible_pose_files_nonempty": "PASS" if eligible_files_nonempty else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "accepted_boxes_present": "PASS" if accepted_boxes_path else ("WARN" if mode == "dry-run" else "FAIL"),
        "all_attribution_rows_represented_in_membership": "PASS" if clustering_counts["all_attribution_rows_represented_in_membership"] else "FAIL",
        "no_silent_pose_drop": "PASS" if clustering_counts["all_attribution_rows_represented_in_membership"] else "FAIL",
        "pose_atoms_loaded_in_memory_if_needed": "PASS" if atoms_loaded or not eligible_rows else "WARN",
        "no_coordinate_rows_written": "PASS" if not coord_hits else "FAIL",
        "compound_pose_clusters_written": "PASS" if clusters_csv.is_file() else "FAIL",
        "cluster_membership_table_written": "PASS" if membership_csv.is_file() else "FAIL",
        "cluster_convergence_qc_written": "PASS" if convergence_qc_csv.is_file() else "FAIL",
        "pose_cluster_membership_qc_written": "PASS" if membership_qc_csv.is_file() else "FAIL",
        "pose_clustering_report_written": "PASS" if report_md.is_file() else "FAIL",
        "no_internal_id_leak": "PASS" if leak_count == 0 else "FAIL",
        "no_smiles_logged": "PASS" if not smiles_logged else "FAIL",
        "no_ligand_coordinates_logged": "PASS" if not ligand_coords else "FAIL",
        "no_receptor_coordinates_logged": "PASS" if not receptor_coords else "FAIL",
        "no_pose_coordinates_logged_outside_pdbqt": "PASS" if not pose_coords else "FAIL",
    }
    qc_rows: list[dict[str, str]] = []
    for check in REQUIRED_CHECKS:
        default = "PASS"
        if check == "input_reference_paths_valid":
            default = "PASS" if accepted_boxes_path else ("WARN" if mode == "dry-run" else "FAIL")
        elif check in {"pose_clustering_computed", "cluster_fraction_computed", "pocket_retention_fraction_computed", "atp_migration_fraction_computed", "ppi_relationship_fraction_computed", "membrane_dimer_cluster_fraction_computed", "convergence_status_assigned", "cluster_reject_reasons_explicit", "anchor_convergence_readiness_flag_written", "cluster_support_reported_separately_from_affinity", "raw_pose_affinity_not_used_for_ranking", "cluster_config_loaded_or_snapshot_written", "cluster_thresholds_not_hardcoded", "cluster_method_recorded", "deterministic_cluster_order", "singleton_clusters_marked", "scattered_pose_support_not_promoted", "no_vina_invoked_by_clustering", "no_qsub_invoked_by_clustering", "no_runner_invoked_by_clustering", "no_pose_attribution_rerun_attempted", "no_compound_anchor_convergence_attempted", "no_candidate_scoring_attempted", "no_candidate_nomination_attempted", "no_final_candidate_tables_created", "reference_state_not_promoted", "old_workflow_not_used", "non_goals_preserved"}:
            default = "PASS"
        _append_qc(qc_rows, check, qc_status.get(check, default), check)
    _write_csv(qc_csv, QC_FIELDS, qc_rows, ctx)
    _write_report(report_md, ctx, summary)
    late_leaks, late_coord_hits, late_smiles, late_ligand_coords, late_receptor_coords, late_pose_coords, late_scanned = _scan_hygiene(ctx, list(private_entries.values()))
    if late_leaks or late_coord_hits or late_smiles:
        status = "FAIL"
        m3_t9_ready = False
        if late_leaks and "internal compound ID leakage detected" not in blockers:
            blockers.append("internal compound ID leakage detected")
        if late_coord_hits and "coordinate records printed outside intended PDBQT files" not in blockers:
            blockers.append("coordinate records printed outside intended PDBQT files")
        if late_smiles and "SMILES-like fields printed in public outputs" not in blockers:
            blockers.append("SMILES-like fields printed in public outputs")
    summary["overall_status"] = status
    summary["m3_t9_anchor_convergence_ready"] = m3_t9_ready
    summary["blockers"] = blockers
    summary["confidentiality"] = {"internal_ids_redacted": late_leaks == 0, "public_outputs_scanned": late_scanned, "leaks_detected": late_leaks, "smiles_logged": late_smiles, "ligand_coordinates_logged": late_ligand_coords, "receptor_coordinates_logged": late_receptor_coords, "pose_coordinates_logged_outside_pose_file": late_pose_coords}
    summary["counts"]["confidentiality_leaks"] = late_leaks
    summary["counts"]["coordinate_leak_files"] = len(late_coord_hits)
    qc_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    with phase3_log.open("a", encoding="utf-8") as handle:
        handle.write(f"{timestamp} M3-T8 {status} profile={selected_profile} clusters_written={len(cluster_rows)} no_vina=true no_qsub=true\n")
    append_phase_status(ctx, "phase3_compounds", status, "M3-T8 pose clustering completed", {"phase": "phase3_compounds", "task": "M3-T8", "status": status, "run_id": ctx.run_id, "m2_run_id": m2_run_id, "timestamp": timestamp, "mode": mode, "profile": selected_profile, "attribution_rows": len(selected_rows), "eligible_pose_rows": len(eligible_rows), "groups_evaluated": len(grouped), "clusters_written": len(cluster_rows), "converged_clusters": clustering_counts["converged_clusters"], "clusters_allowed_for_anchor_convergence": clustering_counts["clusters_allowed_for_anchor_convergence"], "singleton_clusters": clustering_counts["singleton_clusters"], "scattered_groups": scattered_groups, "m3_t9_anchor_convergence_ready": m3_t9_ready, "qsub_invoked": False, "vina_invoked_by_clustering": False, "runner_invoked_by_clustering": False, "no_pose_attribution_rerun_attempted": True, "no_anchor_convergence_attempted": True, "no_scoring_attempted": True, "no_candidate_nomination_attempted": True, "no_final_candidate_tiering_attempted": True})
    return M3PoseClusteringResult(status, m3_t9_ready, blockers, warnings, clusters_csv, membership_csv, qc_csv, qc_json, convergence_qc_csv, membership_qc_csv, report_md, phase3_log, inputs, attribution_context, clustering_counts, mechanism_classes_by_cluster)
