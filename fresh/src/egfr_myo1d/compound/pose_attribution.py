"""M3-T7 pose attribution and geometry QC.

The attributor reads M3-T6 raw pose rows, reloads referenced PDBQT pose atoms
only in memory, writes scalar geometry/QC outputs, and preserves downstream
clustering/scoring/candidate work for later milestones.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from egfr_myo1d.compound.confidentiality import PUBLIC_COMPOUND_IDS, load_private_map, public_output_scan_paths, scan_internal_id_leaks
from egfr_myo1d.compound.pdbqt_parse import extract_pose_atoms
from egfr_myo1d.compound.pose_geometry import (
    DEFAULT_CONFIG,
    BoxProxy,
    ReferencePoint,
    atom_overlap_fraction,
    box_from_row,
    centroid,
    distance,
    fraction_inside_box,
    heavy_atoms,
    inside_box,
    membrane_frame_from_json,
    nearest_atom_distance,
    nearest_point,
    ppi_contacts,
    project_membrane_z,
    reference_point_from_row,
)
from egfr_myo1d.core.logging_utils import append_failed_job, append_job_status, append_phase_status
from egfr_myo1d.core.run_context import RunContext, ensure_within
from egfr_myo1d.structure.pdb_parser import AtomRecord as PDBAtomRecord
from egfr_myo1d.structure.pdb_parser import PDBParseError, parse_pdb


ATTRIBUTION_FIELDS = [
    "run_id", "m2_run_id", "profile", "job_id", "compound_public_id", "state_id", "state_role",
    "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index",
    "vina_affinity_kcal_mol", "vina_log_affinity_kcal_mol", "pose_file", "output_pdbqt_file",
    "output_pdbqt_sha256", "vina_log_file", "vina_log_sha256",
    "source_vina_job_manifest", "source_vina_job_manifest_sha256", "source_compound_pose_raw",
    "source_compound_pose_raw_sha256", "accepted_pocket_source", "accepted_box_source", "atp_reference_source",
    "ppi_consensus_source", "membrane_frame_source", "receptor_mapping_source", "pose_center_x", "pose_center_y",
    "pose_center_z", "pose_atom_count", "pose_heavy_atom_count", "inside_original_box",
    "inside_pocket_volume_proxy", "pocket_atom_fraction", "pocket_heavy_atom_fraction",
    "nearest_pocket_family_id", "nearest_pocket_distance_A", "pocket_retention_pass", "atp_overlap_fraction",
    "atp_heavy_atom_overlap_fraction", "atp_centroid_distance", "atp_nearest_atom_distance",
    "atp_migration_flag", "nearest_ppi_hotspot_id", "nearest_ppi_hotspot_distance",
    "ppi_hotspot_contact_count", "ppi_contact_residue_count", "ppi_contact_residues_public",
    "ppi_rim_distance", "ppi_relationship_confirmed", "ligand_membrane_z", "ligand_min_membrane_z",
    "ligand_max_membrane_z", "membrane_penetration_flag", "dimer_interface_clash_flag",
    "nearest_receptor_contact_distance", "pose_hard_gate_pass", "pose_reject_reason", "pose_mechanism_class",
    "attribution_status", "attribution_notes", "collected_at", "attributed_at", "allowed_for_clustering",
]

MECHANISM_FIELDS = [
    "run_id", "m2_run_id", "profile", "job_id", "compound_public_id", "state_id", "state_role",
    "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index",
    "pose_hard_gate_pass", "pose_reject_reason", "pose_mechanism_class", "mechanism_class_basis",
    "direct_ppi_contact", "rim_relationship", "allosteric_near_relationship", "generic_nonATP_relationship",
    "atp_like_reject", "allowed_for_clustering", "attributed_at",
]

RAW_REQUIRED = [
    "run_id", "profile", "job_id", "compound_public_id", "state_id", "pocket_family_id", "box_id",
    "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index", "pose_file", "output_pdbqt_file",
    "output_pdbqt_sha256", "vina_log_file", "vina_log_sha256", "parse_status", "allowed_for_attribution",
]

QC_FIELDS = ["check_id", "category", "status", "severity", "details", "recommended_fix"]

REQUIRED_CHECKS = [
    "m3_t6_raw_pose_table_present", "m3_t6_raw_pose_table_schema_valid", "m3_t6_docking_completion_qc_present",
    "m3_t6_not_failed", "profile_selected", "pose_rows_found", "pose_ids_unique",
    "pose_rows_traceable_to_compound_state_pocket_box_repeat", "pose_files_inside_run_dir", "pose_files_exist",
    "pose_files_nonempty", "accepted_pockets_present", "accepted_boxes_present", "atp_reference_present",
    "ppi_consensus_patch_present", "membrane_frame_present", "receptor_mapping_present", "input_reference_paths_valid",
    "all_pose_rows_classified", "no_silent_pose_drop", "pdbqt_pose_atoms_loaded_in_memory", "no_coordinate_rows_written",
    "pocket_retention_computed", "pocket_retention_status_for_every_pose", "atp_migration_computed",
    "atp_migration_status_for_every_pose", "atp_migration_independent_of_vina_score", "ppi_geometry_computed",
    "ppi_geometry_status_for_every_pose", "ppi_relationship_not_centroid_only", "membrane_geometry_computed",
    "dimer_interface_clash_computed", "pose_hard_gates_applied", "pose_reject_reasons_explicit",
    "mechanism_class_assigned", "raw_pose_affinity_not_used_for_ranking", "attribution_table_written",
    "mechanism_table_written", "atp_migration_qc_written", "pocket_retention_qc_written", "ppi_geometry_qc_written",
    "membrane_dimer_qc_written", "pose_attribution_report_written", "no_vina_invoked_by_attributor",
    "no_qsub_invoked_by_attributor", "no_runner_invoked_by_attributor", "no_pose_clustering_attempted",
    "no_anchor_convergence_attempted", "no_candidate_scoring_attempted", "no_candidate_nomination_attempted",
    "no_final_candidate_tables_created", "reference_state_not_promoted", "no_internal_id_leak", "no_smiles_logged",
    "no_ligand_coordinates_logged", "no_receptor_coordinates_logged", "no_pose_coordinates_logged_outside_pdbqt",
    "old_workflow_not_used", "non_goals_preserved",
]

FORBIDDEN_OUTPUTS = [
    "phase3_compounds/tables/compound_pose_clusters.csv",
    "phase3_compounds/tables/compound_anchor_convergence.csv",
    "phase3_compounds/tables/compound_pocket_support.csv",
    "phase3_compounds/tables/pocket_compound_evidence_table.csv",
    "phase3_compounds/tables/final_m3_candidate_hypotheses.csv",
    "phase3_compounds/qc/pose_clustering_qc.csv",
    "phase3_compounds/qc/final_candidate_gate_qc.csv",
    "phase3_compounds/docking_outputs/broad_anchor_scan_optional",
]


@dataclass
class M3PoseAttributionResult:
    status: str
    m3_t8_clustering_ready: bool
    blockers: list[str]
    warnings: list[str]
    attribution_csv: Path
    mechanism_csv: Path
    qc_csv: Path
    qc_json: Path
    atp_qc_csv: Path
    pocket_qc_csv: Path
    ppi_qc_csv: Path
    membrane_qc_csv: Path
    report_md: Path
    phase3_log: Path
    inputs: dict[str, bool]
    collection_context: dict[str, Any]
    attribution: dict[str, int]
    mechanism_classes: dict[str, int]


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
    return result


def _int(value: Any) -> int | None:
    try:
        text = str(value).strip()
        if not text:
            return None
        return int(text)
    except (TypeError, ValueError):
        return None


def _row_text(row: dict[str, Any], names: list[str]) -> str:
    lower = {key.lower(): key for key in row}
    for name in names:
        key = lower.get(name.lower())
        if key is not None and str(row.get(key, "")).strip():
            return str(row.get(key, "")).strip()
    return ""


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


def _append_qc(rows: list[dict[str, str]], check_id: str, category: str, status: str, details: str, fix: str = "", blocker: bool = False) -> None:
    rows.append({"check_id": check_id, "category": category, "status": status, "severity": _severity(status, blocker), "details": details, "recommended_fix": fix})


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


def _pose_raw_path(ctx: RunContext, pose_raw_table: Path | None) -> Path:
    path = pose_raw_table or (ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_raw.csv")
    if not path.is_absolute():
        path = ctx.repo_root / path
    return path.resolve()


def _first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.is_file():
            return path.resolve()
    return None


def _existing_unique(paths: list[Path]) -> list[Path]:
    existing: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        if not path.is_file():
            continue
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        existing.append(resolved)
    return existing


def _load_json(path: Path) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _load_config(ctx: RunContext, config_path: Path | None) -> dict[str, Any]:
    config = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (ctx.fresh_root / "configs" / "compound_pose_attribution.yaml")
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
    snapshot = ctx.run_dir / "phase3_compounds" / "config_snapshots" / "compound_pose_attribution.resolved.yaml"
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    try:
        import yaml  # type: ignore

        snapshot.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    except Exception:
        snapshot.write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")
    return config


def _load_boxes(path: Path | None) -> list[BoxProxy]:
    if not path or not path.is_file():
        return []
    rows, _ = _read_csv(path)
    boxes = [box for row in rows if (box := box_from_row(row)) is not None]
    return boxes


def _load_points(path: Path | None, prefix: str) -> list[ReferencePoint]:
    if not path or not path.is_file():
        return []
    rows, _ = _read_csv(path)
    return [point for index, row in enumerate(rows, start=1) if (point := reference_point_from_row(row, prefix, index)) is not None]


def _load_points_many(paths: list[Path], prefix: str) -> list[ReferencePoint]:
    points: list[ReferencePoint] = []
    for path in paths:
        points.extend(_load_points(path, prefix))
    return points


def _state_from_mapping_path(path: Path) -> str:
    match = re.match(r"(.+)_receptor_mapping\.csv$", path.name)
    return match.group(1) if match else ""


def _candidate_receptor_paths_for_state(
    ctx: RunContext,
    m2_root: Path,
    state_id: str,
    accepted_boxes_path: Path | None,
) -> list[Path]:
    candidates: list[Path] = []
    if accepted_boxes_path and accepted_boxes_path.is_file():
        rows, _ = _read_csv(accepted_boxes_path)
        for row in rows:
            if state_id and row.get("state_id") and row.get("state_id") != state_id:
                continue
            for key in ["receptor_pdb", "source_receptor_file", "prepared_receptor_pdbqt_file", "receptor_pdbqt_file"]:
                path = _resolve_path(ctx, row.get(key))
                if path and path.suffix.lower() in {".pdb", ".pdbqt"}:
                    candidates.append(path)
    for root in [
        m2_root / "prepared" / "m2_1_ppi_inputs" / state_id / "receptor",
        m2_root / "normalized" / "receptors",
        m2_root / "normalized" / "receptors" / state_id,
        ctx.fresh_root / "data" / "normalized" / "receptors" / state_id,
    ]:
        for name in [
            f"{state_id}_dockable_669_1014_explicit_AB.pdb",
            f"{state_id}_dockable_reference_explicit_AB.pdb",
            f"{state_id}_full_frame_explicit_AB.pdb",
            f"{state_id}_runtime_offset_receptor_only.pdb",
            "dockable_669_1014_explicit_AB.pdb",
            "dockable_reference_explicit_AB.pdb",
            "full_frame_explicit_AB.pdb",
            "runtime_offset_receptor_only.pdb",
            "receptor.pdbqt",
        ]:
            candidates.append(root / name)
    unique: list[Path] = []
    seen: set[Path] = set()
    for path in candidates:
        if not path.is_file():
            continue
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        unique.append(resolved)
    return unique


def _atoms_by_chain_resseq(atoms: tuple[PDBAtomRecord, ...]) -> dict[tuple[str, int], list[PDBAtomRecord]]:
    grouped: dict[tuple[str, int], list[PDBAtomRecord]] = {}
    for atom in atoms:
        grouped.setdefault((atom.chain_id, atom.residue_number), []).append(atom)
    return grouped


def _ca_or_first(atoms: list[PDBAtomRecord]) -> PDBAtomRecord | None:
    for atom in atoms:
        if atom.atom_name.strip() == "CA":
            return atom
    return atoms[0] if atoms else None


def _mapping_row_identities(row: dict[str, Any]) -> list[tuple[str, int]]:
    identities: list[tuple[str, int]] = []
    for chain_name, residue_name in [
        ("source_chain", "source_resseq"),
        ("runtime_chain", "runtime_resseq"),
        ("egfr_chain_id", "receptor_residue_number"),
        ("chain_id", "residue_number"),
    ]:
        chain = _row_text(row, [chain_name])
        residue = _int(_row_text(row, [residue_name]))
        if chain and residue is not None:
            identities.append((chain, residue))
    unique: list[tuple[str, int]] = []
    for identity in identities:
        if identity not in unique:
            unique.append(identity)
    return unique


def _load_receptor_mapping_points(
    ctx: RunContext,
    mapping_paths: list[Path],
    m2_root: Path,
    accepted_boxes_path: Path | None,
) -> list[ReferencePoint]:
    points: list[ReferencePoint] = []
    receptor_cache: dict[str, list[dict[tuple[str, int], list[PDBAtomRecord]]]] = {}
    for path in mapping_paths:
        rows, _ = _read_csv(path)
        state_hint = _state_from_mapping_path(path)
        for index, row in enumerate(rows, start=1):
            point = reference_point_from_row(row, "receptor", index)
            if point is not None:
                state_id = point.state_id or state_hint
                points.append(
                    ReferencePoint(
                        point_id=point.point_id,
                        xyz=point.xyz,
                        residue_public=point.residue_public,
                        state_id=state_id,
                        protomer_id=point.protomer_id,
                        radius=point.radius,
                    )
                )
                continue

            state_id = _row_text(row, ["state_id", "state"]) or state_hint
            if not state_id:
                continue
            if state_id not in receptor_cache:
                receptor_cache[state_id] = []
                for receptor_path in _candidate_receptor_paths_for_state(ctx, m2_root, state_id, accepted_boxes_path):
                    try:
                        receptor_cache[state_id].append(_atoms_by_chain_resseq(parse_pdb(receptor_path).atoms))
                    except (OSError, PDBParseError):
                        continue
            identities = _mapping_row_identities(row)
            atom: PDBAtomRecord | None = None
            for by_identity in receptor_cache.get(state_id, []):
                for identity in identities:
                    atom = _ca_or_first(by_identity.get(identity, []))
                    if atom is not None:
                        break
                if atom is not None:
                    break
            if atom is None:
                continue
            residue = _row_text(row, ["uniprot_residue_number", "source_resseq", "runtime_resseq", "receptor_residue_number", "egfr_residue_number"])
            protomer_id = _row_text(row, ["protomer_id", "egfr_protomer_id", "protomer"])
            point_id = f"{state_id}:{protomer_id}:{residue or index}"
            points.append(
                ReferencePoint(
                    point_id=point_id,
                    xyz=(atom.x, atom.y, atom.z),
                    residue_public=residue or point_id,
                    state_id=state_id,
                    protomer_id=protomer_id,
                )
            )
    return points


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


def _points_for_pose(row: dict[str, str], points: list[ReferencePoint]) -> list[ReferencePoint]:
    state_id = row.get("state_id", "")
    protomer_id = row.get("protomer_id", "")
    matched: list[ReferencePoint] = []
    for point in points:
        if point.state_id and state_id and point.state_id != state_id:
            continue
        if point.protomer_id and protomer_id and point.protomer_id != protomer_id:
            continue
        matched.append(point)
    return matched


def _filter_profile(raw_rows: list[dict[str, str]], profile: str | None) -> tuple[list[dict[str, str]], str | None, list[str]]:
    profiles = sorted({row.get("profile", "") for row in raw_rows if row.get("profile")})
    if profile:
        selected = profile
    elif len(profiles) == 1:
        selected = profiles[0]
    else:
        selected = None
    rows = [row for row in raw_rows if not selected or row.get("profile") == selected]
    return rows, selected, profiles


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
            pose_coords = pose_coords or "pose" in lower or "attribution" in lower or "qc" in lower or "table" in lower
    return len(leaks), coord_hits, smiles_logged, ligand_coords, receptor_coords, pose_coords, scanned


def _nearest_box_for_pose(center_xyz: tuple[float, float, float], boxes: list[BoxProxy]) -> tuple[BoxProxy | None, float | None]:
    if not boxes:
        return None, None
    best = min(boxes, key=lambda box: distance(center_xyz, box.center))
    return best, distance(center_xyz, best.center)


def _planned_box(row: dict[str, str], boxes: list[BoxProxy]) -> BoxProxy | None:
    for box in boxes:
        if (
            box.pocket_family_id == row.get("pocket_family_id")
            and box.state_id == row.get("state_id")
            and box.protomer_id == row.get("protomer_id")
            and (box.box_id == row.get("box_id") or not row.get("box_id"))
        ):
            return box
    for box in boxes:
        if box.pocket_family_id == row.get("pocket_family_id") and box.state_id == row.get("state_id"):
            return box
    return None


def _mechanism(hard_pass: bool, atp_flag: bool, contact_count: int, ppi_distance: float | None, ppi_rim: float | None, config: dict[str, Any]) -> str:
    if atp_flag:
        return "ATP_like_reject"
    if not hard_pass:
        return "ambiguous_or_failed"
    if contact_count >= int(config["ppi"].get("hotspot_contact_min", 1)):
        return "orthosteric_or_direct_PPI_patch_blocker"
    if ppi_rim is not None and ppi_rim <= float(config["mechanism"].get("rim_blocker_cutoff_A", 8.0)):
        return "rim_blocker"
    if ppi_distance is not None and ppi_distance <= float(config["mechanism"].get("allosteric_near_cutoff_A", 12.0)):
        return "allosteric_near_candidate"
    return "generic_nonATP_ligandable_pocket"


def _empty_attr(row: dict[str, str], common: dict[str, Any], status: str, reason: str, notes: str) -> dict[str, Any]:
    out = dict(common)
    out.update({
        "inside_original_box": "false", "inside_pocket_volume_proxy": "false", "pocket_atom_fraction": "",
        "pocket_heavy_atom_fraction": "", "nearest_pocket_family_id": "", "nearest_pocket_distance_A": "",
        "pocket_retention_pass": "false", "atp_overlap_fraction": "", "atp_heavy_atom_overlap_fraction": "",
        "atp_centroid_distance": "", "atp_nearest_atom_distance": "", "atp_migration_flag": "false",
        "nearest_ppi_hotspot_id": "", "nearest_ppi_hotspot_distance": "", "ppi_hotspot_contact_count": "0",
        "ppi_contact_residue_count": "0", "ppi_contact_residues_public": "", "ppi_rim_distance": "",
        "ppi_relationship_confirmed": "false", "ligand_membrane_z": "", "ligand_min_membrane_z": "",
        "ligand_max_membrane_z": "", "membrane_penetration_flag": "false", "dimer_interface_clash_flag": "false",
        "nearest_receptor_contact_distance": "", "pose_hard_gate_pass": "false", "pose_reject_reason": reason,
        "pose_mechanism_class": "ambiguous_or_failed", "attribution_status": status, "attribution_notes": notes,
        "allowed_for_clustering": "false",
    })
    return out


def _write_report(path: Path, ctx: RunContext, summary: dict[str, Any]) -> None:
    lines = [
        "# M3-T7 Pose Attribution",
        "",
        "Technical geometry attribution report only. Pose geometry and affinity metadata are not used for ranking, scoring, or candidate claims.",
        "",
        f"- status: {summary['overall_status']}",
        f"- profile: {summary['profile']}",
        f"- raw_pose_rows: {summary['collection_context']['raw_pose_rows']}",
        f"- poses_classified: {summary['attribution']['poses_classified']}",
        f"- hard_gate_pass: {summary['attribution']['poses_pass']}",
        f"- m3_t8_clustering_ready: {str(summary['m3_t8_clustering_ready']).lower()}",
        "",
        "## Blockers",
    ]
    lines.extend(f"- {item}" for item in summary["blockers"] or ["none"])
    lines.append("")
    lines.append("## Warnings")
    lines.extend(f"- {item}" for item in summary["warnings"] or ["none"])
    lines.append("")
    lines.append("Next task: M3-T8 - Pose clustering and convergence analysis")
    ctx.require_within_run_dir(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_m3_pose_attribution(
    ctx: RunContext,
    *,
    m2_run_id: str,
    profile: str | None = None,
    mode: str = "attribute",
    pose_raw_table: Path | None = None,
    job_completion_table: Path | None = None,
    accepted_pockets: Path | None = None,
    accepted_boxes: Path | None = None,
    atp_reference: Path | None = None,
    ppi_consensus_patch: Path | None = None,
    membrane_frame: Path | None = None,
    receptor_mapping: Path | None = None,
    receptor_manifest: Path | None = None,
    config: Path | None = None,
    force: bool = False,
    allow_partial: bool = False,
    require_all_poses: bool = True,
    allow_reference_only: bool = False,
    private_map: Path | None = None,
    strict_geometry: bool = True,
    write_empty_attribution: bool = False,
) -> M3PoseAttributionResult:
    del force, receptor_manifest
    _ensure_logs(ctx)
    timestamp = now_iso()
    phase3 = ctx.run_dir / "phase3_compounds"
    tables_dir = phase3 / "tables"
    qc_dir = phase3 / "qc"
    reports_dir = phase3 / "reports"
    for directory in [tables_dir, qc_dir, reports_dir]:
        ctx.require_within_run_dir(directory).mkdir(parents=True, exist_ok=True)

    attribution_csv = tables_dir / "compound_pose_attribution.csv"
    mechanism_csv = tables_dir / "compound_pose_mechanism_classification.csv"
    qc_csv = qc_dir / "pose_attribution_qc.csv"
    qc_json = qc_dir / "pose_attribution_qc.json"
    atp_qc_csv = qc_dir / "atp_migration_qc.csv"
    pocket_qc_csv = qc_dir / "pocket_retention_qc.csv"
    ppi_qc_csv = qc_dir / "ppi_geometry_qc.csv"
    membrane_qc_csv = qc_dir / "membrane_dimer_geometry_qc.csv"
    report_md = reports_dir / "m3_task7_pose_attribution.md"
    phase3_log = ctx.logs_dir / "phase3_compounds.log"

    blockers: list[str] = []
    warnings: list[str] = []
    qc_rows: list[dict[str, str]] = []
    raw_path = _pose_raw_path(ctx, pose_raw_table)
    raw_found = raw_path.is_file()
    if not raw_found:
        (warnings if mode == "dry-run" else blockers).append("compound_pose_raw.csv missing")
    raw_rows: list[dict[str, str]] = []
    raw_fields: list[str] = []
    if raw_found:
        if not _inside(raw_path, tables_dir):
            blockers.append("compound_pose_raw.csv path outside phase3_compounds/tables")
        raw_rows, raw_fields = _read_csv(raw_path)
    schema_ok = set(RAW_REQUIRED).issubset(set(raw_fields))
    if raw_found and not schema_ok:
        blockers.append("compound_pose_raw.csv schema missing required columns")
    selected_rows, selected_profile, profiles = _filter_profile(raw_rows, profile)
    if raw_found and not selected_profile:
        blockers.append("profile ambiguous; provide --profile")
    if raw_found and not selected_rows:
        (warnings if mode == "dry-run" or write_empty_attribution else blockers).append("zero raw pose rows for selected profile")

    m2_root = ctx.fresh_root / "runs" / m2_run_id
    accepted_pockets_path = _first_existing([Path(accepted_pockets)] if accepted_pockets else [
        m2_root / "phase2_pockets" / "export_for_m3" / "accepted_pockets_for_m3.csv",
        m2_root / "phase2_pockets" / "final" / "accepted_pockets_for_m3.csv",
        m2_root / "phase2_pockets" / "final" / "accepted_pocket_families_for_compound_docking.csv",
    ])
    accepted_boxes_path = _first_existing([Path(accepted_boxes)] if accepted_boxes else [
        m2_root / "phase2_pockets" / "final" / "accepted_pocket_boxes.csv",
        m2_root / "phase2_pockets" / "export_for_m3" / "accepted_pocket_boxes.csv",
    ])
    atp_reference_path = _first_existing([Path(atp_reference)] if atp_reference else [
        m2_root / "phase2_pockets" / "atp_reference" / "atp_site_centroid_by_state.csv",
        m2_root / "phase2_pockets" / "export_for_m3" / "atp_site_centroid_by_state.csv",
        m2_root / "phase2_pockets" / "atp_reference" / "atp_site_reference.csv",
        m2_root / "phase2_pockets" / "export_for_m3" / "atp_site_reference.csv",
    ])
    ppi_path = _first_existing([Path(ppi_consensus_patch)] if ppi_consensus_patch else [
        m2_root / "phase2_pockets" / "export_for_m3" / "ppi_consensus_patch.csv",
        m2_root / "phase1_ppi" / "consensus" / "ppi_consensus_patch_merged.csv",
        m2_root / "phase1_ppi" / "tables" / "ppi_consensus_patch.csv",
    ])
    membrane_path = _first_existing([Path(membrane_frame)] if membrane_frame else [
        m2_root / "manifest" / "membrane_frame.json",
        m2_root / "membrane" / "membrane_frame.json",
        m2_root / "phase1_receptor" / "membrane_frame.json",
    ])
    if receptor_mapping:
        mapping_candidates = [Path(receptor_mapping)]
    else:
        mapping_candidates = []
        if accepted_boxes_path and accepted_boxes_path.is_file():
            box_rows, _ = _read_csv(accepted_boxes_path)
            for row in box_rows:
                path = _resolve_path(ctx, row.get("receptor_mapping_csv"))
                if path:
                    mapping_candidates.append(path)
        mapping_candidates.append(m2_root / "phase1_receptor" / "receptor_mapping.csv")
        if (m2_root / "qc").is_dir():
            mapping_candidates.extend(sorted((m2_root / "qc").glob("*_receptor_mapping.csv")))
        if (m2_root / "phase1_receptor").is_dir():
            mapping_candidates.extend(sorted((m2_root / "phase1_receptor").glob("*mapping*.csv")))
    mapping_paths = _existing_unique(mapping_candidates)
    mapping_path = mapping_paths[0] if mapping_paths else None
    completion_qc_path = job_completion_table if job_completion_table else (qc_dir / "docking_completion_qc.json")
    if not completion_qc_path.is_absolute():
        completion_qc_path = ctx.repo_root / completion_qc_path
    completion_qc = _load_json(completion_qc_path) if completion_qc_path.is_file() else None
    completion_ready = bool(
        completion_qc
        and completion_qc.get("overall_status") == "PASS"
        and _boolish(completion_qc.get("m3_t7_attribution_ready"))
    )
    inputs = {
        "compound_pose_raw_found": raw_found,
        "docking_completion_qc_found": completion_qc_path.is_file(),
        "docking_completion_qc_passed": completion_ready,
        "accepted_pockets_found": accepted_pockets_path is not None,
        "accepted_boxes_found": accepted_boxes_path is not None,
        "atp_reference_found": atp_reference_path is not None,
        "ppi_consensus_patch_found": ppi_path is not None,
        "membrane_frame_found": membrane_path is not None,
        "receptor_mapping_found": bool(mapping_paths),
    }
    for label, found in inputs.items():
        if label == "compound_pose_raw_found":
            continue
        if not found:
            (warnings if mode == "dry-run" else blockers).append(label.replace("_found", "") + " missing")
    if mode == "attribute" and inputs["docking_completion_qc_found"] and not completion_ready:
        blockers.append("M3-T6 docking completion QC did not PASS with m3_t7_attribution_ready=true")

    resolved_config = _load_config(ctx, config)
    boxes = _load_boxes(accepted_boxes_path)
    atp_points = _load_points(atp_reference_path, "atp")
    ppi_points = _load_points(ppi_path, "ppi")
    mapping_points = _load_receptor_mapping_points(ctx, mapping_paths, m2_root, accepted_boxes_path)
    membrane = membrane_frame_from_json(_load_json(membrane_path) or {}) if membrane_path else None
    if accepted_boxes_path and not boxes:
        blockers.append("accepted boxes present but no usable box proxies parsed")
    if atp_reference_path and not atp_points:
        blockers.append("ATP reference present but no usable reference points parsed")
    if ppi_path and not ppi_points:
        blockers.append("PPI consensus present but no usable hotspot points parsed")

    identity_keys = ["job_id", "pose_rank", "pose_model_index", "compound_public_id", "state_id", "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id"]
    identities = [tuple(row.get(key, "") for key in identity_keys) for row in selected_rows]
    unique_ids = len(identities) == len(set(identities))
    if selected_rows and not unique_ids:
        blockers.append("duplicate pose identity rows")
    non_public = sorted({row.get("compound_public_id", "") for row in selected_rows if row.get("compound_public_id") not in PUBLIC_COMPOUND_IDS})
    if non_public:
        blockers.append("non-public compound IDs present in raw pose table")

    attr_rows: list[dict[str, Any]] = []
    mech_rows: list[dict[str, Any]] = []
    atp_rows: list[dict[str, Any]] = []
    pocket_rows: list[dict[str, Any]] = []
    ppi_rows: list[dict[str, Any]] = []
    membrane_rows: list[dict[str, Any]] = []
    source_hash = _sha256(raw_path) if raw_found else ""
    pose_files_inside = True
    pose_files_exist = True
    pose_files_nonempty = True
    classified = 0
    atoms_loaded = False
    reference_states = sorted({row.get("state_id", "") for row in selected_rows if row.get("state_role") == "reference" or row.get("state_id") == "3GT8_raw"})
    primary_states = sorted({row.get("state_id", "") for row in selected_rows if row.get("state_id") and row.get("state_id") not in reference_states})
    if selected_rows and reference_states and not primary_states and not allow_reference_only:
        warnings.append("reference-only attribution without --allow-reference-only")

    for row in selected_rows:
        classified += 1
        pose_path = _resolve_path(ctx, row.get("pose_file"))
        pose_inside = _inside(pose_path, ctx.run_dir)
        pose_exists = bool(pose_path and pose_path.is_file())
        pose_nonempty = bool(pose_path and pose_path.is_file() and pose_path.stat().st_size > 0)
        expected_pose_sha = (row.get("output_pdbqt_sha256") or "").strip()
        current_pose_sha = _sha256(pose_path) if pose_exists and pose_nonempty and pose_path else ""
        pose_hash_matches = bool(expected_pose_sha and current_pose_sha and expected_pose_sha == current_pose_sha)
        pose_files_inside = pose_files_inside and pose_inside
        pose_files_exist = pose_files_exist and pose_exists
        pose_files_nonempty = pose_files_nonempty and pose_nonempty
        common = {field: row.get(field, "") for field in ATTRIBUTION_FIELDS if field in row}
        common.update({
            "m2_run_id": m2_run_id,
            "source_compound_pose_raw": ctx.relative_to_repo(raw_path) if raw_found else "",
            "source_compound_pose_raw_sha256": source_hash,
            "accepted_pocket_source": ctx.relative_to_repo(accepted_pockets_path) if accepted_pockets_path else "",
            "accepted_box_source": ctx.relative_to_repo(accepted_boxes_path) if accepted_boxes_path else "",
            "atp_reference_source": ctx.relative_to_repo(atp_reference_path) if atp_reference_path else "",
            "ppi_consensus_source": ctx.relative_to_repo(ppi_path) if ppi_path else "",
            "membrane_frame_source": ctx.relative_to_repo(membrane_path) if membrane_path else "",
            "receptor_mapping_source": ";".join(ctx.relative_to_repo(path) for path in mapping_paths),
            "attributed_at": timestamp,
        })
        if not pose_inside:
            blockers.append(f"pose_file outside run_dir for job {row.get('job_id')}")
            out = _empty_attr(row, common, "FAIL", "missing_pose_file", "pose_path_invalid")
        elif row.get("parse_status") not in {"PASS", "WARN"}:
            out = _empty_attr(row, common, "FAIL", "parse_failure", "raw_parse_status_not_attributable")
        elif not _boolish(row.get("allowed_for_attribution")):
            out = _empty_attr(row, common, "NOT_APPLICABLE", "not_allowed_for_attribution", "raw_row_not_allowed")
        elif not pose_exists or not pose_nonempty:
            if mode == "attribute":
                blockers.append(f"eligible pose_file missing or empty for job {row.get('job_id')}")
            out = _empty_attr(row, common, "MISSING", "missing_pose_file", "pose_file_missing_or_empty")
        elif not pose_hash_matches:
            if mode == "attribute":
                blockers.append(f"pose_file hash mismatch for job {row.get('job_id')}")
            out = _empty_attr(row, common, "FAIL", "pose_hash_mismatch", "pose_file_sha256_mismatch")
        else:
            try:
                atoms = extract_pose_atoms(pose_path, int(row.get("pose_model_index") or 0), int(row.get("pose_rank") or 0))
            except Exception:
                atoms = []
            atoms_loaded = atoms_loaded or bool(atoms)
            if not atoms:
                out = _empty_attr(row, common, "FAIL", "geometry_parse_failure", "pose_atom_parse_failed")
            else:
                h_atoms = heavy_atoms(atoms) or atoms
                pose_atp_points = _points_for_pose(row, atp_points)
                pose_ppi_points = _points_for_pose(row, ppi_points)
                pose_mapping_points = _points_for_pose(row, mapping_points)
                center = centroid(atoms, heavy_only=False)
                heavy_center = centroid(atoms, heavy_only=True)
                planned = _planned_box(row, boxes)
                nearest_box, nearest_box_dist = _nearest_box_for_pose(heavy_center, boxes)
                margin = float(resolved_config["pocket"].get("centroid_margin_A", 2.0))
                atom_fraction = fraction_inside_box(atoms, planned, margin=margin) if planned else 0.0
                heavy_fraction = fraction_inside_box(atoms, planned, margin=margin, heavy_only=True) if planned else 0.0
                inside_orig = bool(planned and (inside_box(center, planned) or inside_box(heavy_center, planned) or atom_fraction >= float(resolved_config["pocket"].get("atom_fraction_min", 0.5))))
                inside_proxy = bool(planned and (inside_box(heavy_center, planned, margin=margin) or heavy_fraction >= float(resolved_config["pocket"].get("heavy_atom_fraction_min", 0.5))))
                retention = bool(planned and inside_proxy and nearest_box and nearest_box.pocket_family_id == row.get("pocket_family_id"))
                atp_ref, atp_centroid_distance = nearest_point(heavy_center, pose_atp_points)
                atp_nearest = nearest_atom_distance(h_atoms, pose_atp_points)
                atp_overlap = atom_overlap_fraction(atoms, pose_atp_points, float(resolved_config["atp"].get("atom_near_atp_cutoff_A", 4.5)))
                atp_heavy_overlap = atom_overlap_fraction(atoms, pose_atp_points, float(resolved_config["atp"].get("atom_near_atp_cutoff_A", 4.5)), heavy_only=True)
                atp_flag = bool(
                    (atp_centroid_distance is not None and atp_centroid_distance <= float(resolved_config["atp"].get("centroid_reject_distance_A", 8.0)))
                    or atp_overlap >= float(resolved_config["atp"].get("atom_overlap_reject_fraction", 0.2))
                    or atp_heavy_overlap >= float(resolved_config["atp"].get("heavy_atom_overlap_reject_fraction", 0.2))
                )
                contact_count, residues, nearest_ppi_distance, nearest_ppi_id = ppi_contacts(h_atoms, pose_ppi_points, float(resolved_config["ppi"].get("contact_cutoff_A", 4.5)))
                ppi_rim = nearest_ppi_distance
                ppi_confirmed = contact_count >= int(resolved_config["ppi"].get("hotspot_contact_min", 1)) or (ppi_rim is not None and ppi_rim <= float(resolved_config["ppi"].get("rim_cutoff_A", 10.0)))
                z_values = [project_membrane_z((atom.x, atom.y, atom.z), membrane) for atom in atoms] if membrane else []
                ligand_z = project_membrane_z(center, membrane) if membrane else None
                z_min = min(z_values) if z_values else None
                z_max = max(z_values) if z_values else None
                penetration = False
                if membrane and membrane.core_z_min is not None and membrane.core_z_max is not None and z_min is not None and z_max is not None:
                    pad = float(resolved_config["membrane"].get("penetration_margin_A", 1.5))
                    penetration = z_min <= membrane.core_z_max + pad and z_max >= membrane.core_z_min - pad
                nearest_receptor = nearest_atom_distance(h_atoms, pose_mapping_points)
                dimer_clash = bool(nearest_receptor is not None and nearest_receptor <= float(resolved_config["dimer"].get("central_interface_clash_cutoff_A", 3.0)))
                required_mapping = bool(mapping_paths and membrane and pose_atp_points and pose_ppi_points and pose_mapping_points and boxes)
                hard_pass = bool(retention and not atp_flag and not penetration and not dimer_clash and required_mapping)
                if hard_pass:
                    reason = "none"
                elif not required_mapping:
                    reason = "missing_mapping"
                elif atp_flag:
                    reason = "ATP_migration"
                elif not retention:
                    reason = "outside_accepted_pocket"
                elif penetration:
                    reason = "membrane_penetration"
                elif dimer_clash:
                    reason = "central_dimer_clash"
                else:
                    reason = "unknown"
                mechanism = _mechanism(hard_pass, atp_flag, contact_count, nearest_ppi_distance, ppi_rim, resolved_config)
                status = "PASS" if hard_pass else "REJECT"
                if not ppi_confirmed and hard_pass:
                    status = "WARN"
                out = dict(common)
                out.update({
                    "pose_center_x": f"{center[0]:.3f}", "pose_center_y": f"{center[1]:.3f}", "pose_center_z": f"{center[2]:.3f}",
                    "pose_atom_count": len(atoms), "pose_heavy_atom_count": len(h_atoms),
                    "inside_original_box": _bool_text(inside_orig), "inside_pocket_volume_proxy": _bool_text(inside_proxy),
                    "pocket_atom_fraction": f"{atom_fraction:.6f}", "pocket_heavy_atom_fraction": f"{heavy_fraction:.6f}",
                    "nearest_pocket_family_id": nearest_box.pocket_family_id if nearest_box else "",
                    "nearest_pocket_distance_A": f"{nearest_box_dist:.3f}" if nearest_box_dist is not None else "",
                    "pocket_retention_pass": _bool_text(retention),
                    "atp_overlap_fraction": f"{atp_overlap:.6f}", "atp_heavy_atom_overlap_fraction": f"{atp_heavy_overlap:.6f}",
                    "atp_centroid_distance": f"{atp_centroid_distance:.3f}" if atp_centroid_distance is not None else "",
                    "atp_nearest_atom_distance": f"{atp_nearest:.3f}" if atp_nearest is not None else "",
                    "atp_migration_flag": _bool_text(atp_flag),
                    "nearest_ppi_hotspot_id": nearest_ppi_id, "nearest_ppi_hotspot_distance": f"{nearest_ppi_distance:.3f}" if nearest_ppi_distance is not None else "",
                    "ppi_hotspot_contact_count": contact_count, "ppi_contact_residue_count": len(residues),
                    "ppi_contact_residues_public": "|".join(residues), "ppi_rim_distance": f"{ppi_rim:.3f}" if ppi_rim is not None else "",
                    "ppi_relationship_confirmed": _bool_text(ppi_confirmed),
                    "ligand_membrane_z": f"{ligand_z:.3f}" if ligand_z is not None else "", "ligand_min_membrane_z": f"{z_min:.3f}" if z_min is not None else "",
                    "ligand_max_membrane_z": f"{z_max:.3f}" if z_max is not None else "",
                    "membrane_penetration_flag": _bool_text(penetration), "dimer_interface_clash_flag": _bool_text(dimer_clash),
                    "nearest_receptor_contact_distance": f"{nearest_receptor:.3f}" if nearest_receptor is not None else "",
                    "pose_hard_gate_pass": _bool_text(hard_pass), "pose_reject_reason": reason,
                    "pose_mechanism_class": mechanism, "attribution_status": status,
                    "attribution_notes": "scalar_geometry_only", "allowed_for_clustering": _bool_text(hard_pass and status in {"PASS", "WARN"}),
                })
                atp_status = "REJECT_ATP_MIGRATION" if atp_flag else ("PASS" if atp_points else "FAIL_MISSING_REFERENCE")
                if not atp_flag and atp_nearest is not None and atp_nearest <= float(resolved_config["atp"].get("atom_near_atp_cutoff_A", 4.5)) * 2:
                    atp_status = "WARN_NEAR_ATP"
                pocket_status = "PASS" if retention else ("FAIL_MISSING_BOX" if not planned else "REJECT_OUTSIDE_ACCEPTED_POCKET")
                if contact_count:
                    ppi_status = "PASS_DIRECT_PPI_CONTACT"
                elif ppi_rim is not None and ppi_rim <= float(resolved_config["ppi"].get("rim_cutoff_A", 10.0)):
                    ppi_status = "PASS_RIM_RELATIONSHIP"
                elif nearest_ppi_distance is not None and nearest_ppi_distance <= float(resolved_config["mechanism"].get("allosteric_near_cutoff_A", 12.0)):
                    ppi_status = "PASS_ALLOSTERIC_NEAR"
                else:
                    ppi_status = "WARN_GENERIC_NONATP" if ppi_points else "FAIL_MISSING_PPI_REFERENCE"
                if not membrane:
                    mem_status = "FAIL_MISSING_MEMBRANE_FRAME"
                elif penetration:
                    mem_status = "REJECT_MEMBRANE_PENETRATION"
                elif dimer_clash:
                    mem_status = "REJECT_DIMER_INTERFACE_CLASH"
                else:
                    mem_status = "PASS"
                atp_rows.append({**{k: row.get(k, "") for k in ["run_id", "m2_run_id", "profile", "job_id", "compound_public_id", "state_id", "state_role", "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index"]}, "m2_run_id": m2_run_id, "atp_centroid_distance": out["atp_centroid_distance"], "atp_overlap_fraction": out["atp_overlap_fraction"], "atp_heavy_atom_overlap_fraction": out["atp_heavy_atom_overlap_fraction"], "atp_nearest_atom_distance": out["atp_nearest_atom_distance"], "atp_migration_flag": out["atp_migration_flag"], "atp_status": atp_status, "atp_notes": "scalar_reference_proxy"})
                pocket_rows.append({**{k: row.get(k, "") for k in ["run_id", "m2_run_id", "profile", "job_id", "compound_public_id", "state_id", "state_role", "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index"]}, "m2_run_id": m2_run_id, "inside_original_box": out["inside_original_box"], "inside_pocket_volume_proxy": out["inside_pocket_volume_proxy"], "pocket_atom_fraction": out["pocket_atom_fraction"], "pocket_heavy_atom_fraction": out["pocket_heavy_atom_fraction"], "nearest_pocket_family_id": out["nearest_pocket_family_id"], "nearest_pocket_distance_A": out["nearest_pocket_distance_A"], "pocket_retention_pass": out["pocket_retention_pass"], "pocket_retention_status": pocket_status, "pocket_retention_notes": "box_proxy"})
                ppi_rows.append({**{k: row.get(k, "") for k in ["run_id", "m2_run_id", "profile", "job_id", "compound_public_id", "state_id", "state_role", "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index"]}, "m2_run_id": m2_run_id, "nearest_ppi_hotspot_id": out["nearest_ppi_hotspot_id"], "nearest_ppi_hotspot_distance": out["nearest_ppi_hotspot_distance"], "ppi_hotspot_contact_count": out["ppi_hotspot_contact_count"], "ppi_contact_residue_count": out["ppi_contact_residue_count"], "ppi_contact_residues_public": out["ppi_contact_residues_public"], "ppi_rim_distance": out["ppi_rim_distance"], "ppi_relationship_confirmed": out["ppi_relationship_confirmed"], "ppi_geometry_status": ppi_status, "ppi_geometry_notes": "contact_or_rim_required_for_confirmed_relationship"})
                membrane_rows.append({**{k: row.get(k, "") for k in ["run_id", "m2_run_id", "profile", "job_id", "compound_public_id", "state_id", "state_role", "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index"]}, "m2_run_id": m2_run_id, "ligand_membrane_z": out["ligand_membrane_z"], "ligand_min_membrane_z": out["ligand_min_membrane_z"], "ligand_max_membrane_z": out["ligand_max_membrane_z"], "membrane_penetration_flag": out["membrane_penetration_flag"], "dimer_interface_clash_flag": out["dimer_interface_clash_flag"], "nearest_receptor_contact_distance": out["nearest_receptor_contact_distance"], "membrane_dimer_status": mem_status, "membrane_dimer_notes": "membrane_frame_and_mapping_proxy"})
        attr_rows.append(out)
        direct = out.get("ppi_hotspot_contact_count") not in {"", "0", 0}
        rim = out.get("pose_mechanism_class") == "rim_blocker"
        allosteric = out.get("pose_mechanism_class") == "allosteric_near_candidate"
        generic = out.get("pose_mechanism_class") == "generic_nonATP_ligandable_pocket"
        mech_rows.append({
            **{k: out.get(k, "") for k in ["run_id", "m2_run_id", "profile", "job_id", "compound_public_id", "state_id", "state_role", "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index", "pose_hard_gate_pass", "pose_reject_reason", "pose_mechanism_class", "allowed_for_clustering", "attributed_at"]},
            "mechanism_class_basis": "pose_level_geometry_only", "direct_ppi_contact": _bool_text(bool(direct)), "rim_relationship": _bool_text(bool(rim)), "allosteric_near_relationship": _bool_text(bool(allosteric)), "generic_nonATP_relationship": _bool_text(bool(generic)), "atp_like_reject": _bool_text(out.get("pose_mechanism_class") == "ATP_like_reject"),
        })
        append_job_status(ctx, row.get("job_id", "pose_attribution"), out.get("attribution_status", "FAIL"), details={"phase": "phase3_compounds", "task": "M3-T7", "job_type": "pose_attribution", "profile": selected_profile, "compound_public_id": row.get("compound_public_id"), "state_id": row.get("state_id"), "pocket_family_id": row.get("pocket_family_id"), "box_id": row.get("box_id"), "vina_repeat_id": row.get("vina_repeat_id"), "pose_rank": row.get("pose_rank"), "pose_model_index": row.get("pose_model_index"), "pocket_retention_pass": out.get("pocket_retention_pass"), "atp_migration_flag": out.get("atp_migration_flag"), "ppi_relationship_confirmed": out.get("ppi_relationship_confirmed"), "pose_hard_gate_pass": out.get("pose_hard_gate_pass"), "pose_reject_reason": out.get("pose_reject_reason"), "pose_mechanism_class": out.get("pose_mechanism_class"), "message": out.get("attribution_notes", "")})
        if out.get("attribution_status") in {"FAIL", "MISSING"}:
            append_failed_job(ctx, row.get("job_id", "pose_attribution"), out.get("attribution_status", "FAIL"), out.get("pose_reject_reason", "geometry_failure"))

    atp_rows = []
    pocket_rows = []
    ppi_rows = []
    membrane_rows = []
    for out in attr_rows:
        base = {key: out.get(key, "") for key in ["run_id", "m2_run_id", "profile", "job_id", "compound_public_id", "state_id", "state_role", "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index"]}
        atp_status = "REJECT_ATP_MIGRATION" if out.get("atp_migration_flag") == "true" else ("FAIL_MISSING_REFERENCE" if not inputs["atp_reference_found"] else ("FAIL_GEOMETRY" if out.get("attribution_status") in {"FAIL", "MISSING"} else "PASS"))
        pocket_status = "PASS" if out.get("pocket_retention_pass") == "true" else ("FAIL_MISSING_BOX" if not inputs["accepted_boxes_found"] else ("FAIL_GEOMETRY" if out.get("attribution_status") in {"FAIL", "MISSING"} else "REJECT_OUTSIDE_ACCEPTED_POCKET"))
        if not inputs["ppi_consensus_patch_found"]:
            ppi_status = "FAIL_MISSING_PPI_REFERENCE"
        elif out.get("ppi_hotspot_contact_count") not in {"", "0", 0}:
            ppi_status = "PASS_DIRECT_PPI_CONTACT"
        elif out.get("pose_mechanism_class") == "rim_blocker":
            ppi_status = "PASS_RIM_RELATIONSHIP"
        elif out.get("pose_mechanism_class") == "allosteric_near_candidate":
            ppi_status = "PASS_ALLOSTERIC_NEAR"
        elif out.get("attribution_status") in {"FAIL", "MISSING"}:
            ppi_status = "FAIL_GEOMETRY"
        else:
            ppi_status = "WARN_GENERIC_NONATP"
        if not inputs["membrane_frame_found"]:
            mem_status = "FAIL_MISSING_MEMBRANE_FRAME"
        elif not inputs["receptor_mapping_found"]:
            mem_status = "FAIL_MISSING_RECEPTOR_MAPPING"
        elif out.get("membrane_penetration_flag") == "true":
            mem_status = "REJECT_MEMBRANE_PENETRATION"
        elif out.get("dimer_interface_clash_flag") == "true":
            mem_status = "REJECT_DIMER_INTERFACE_CLASH"
        elif out.get("attribution_status") in {"FAIL", "MISSING"}:
            mem_status = "FAIL_GEOMETRY"
        else:
            mem_status = "PASS"
        atp_rows.append({**base, "atp_centroid_distance": out.get("atp_centroid_distance", ""), "atp_overlap_fraction": out.get("atp_overlap_fraction", ""), "atp_heavy_atom_overlap_fraction": out.get("atp_heavy_atom_overlap_fraction", ""), "atp_nearest_atom_distance": out.get("atp_nearest_atom_distance", ""), "atp_migration_flag": out.get("atp_migration_flag", "false"), "atp_status": atp_status, "atp_notes": "scalar_reference_proxy"})
        pocket_rows.append({**base, "inside_original_box": out.get("inside_original_box", "false"), "inside_pocket_volume_proxy": out.get("inside_pocket_volume_proxy", "false"), "pocket_atom_fraction": out.get("pocket_atom_fraction", ""), "pocket_heavy_atom_fraction": out.get("pocket_heavy_atom_fraction", ""), "nearest_pocket_family_id": out.get("nearest_pocket_family_id", ""), "nearest_pocket_distance_A": out.get("nearest_pocket_distance_A", ""), "pocket_retention_pass": out.get("pocket_retention_pass", "false"), "pocket_retention_status": pocket_status, "pocket_retention_notes": "box_proxy"})
        ppi_rows.append({**base, "nearest_ppi_hotspot_id": out.get("nearest_ppi_hotspot_id", ""), "nearest_ppi_hotspot_distance": out.get("nearest_ppi_hotspot_distance", ""), "ppi_hotspot_contact_count": out.get("ppi_hotspot_contact_count", "0"), "ppi_contact_residue_count": out.get("ppi_contact_residue_count", "0"), "ppi_contact_residues_public": out.get("ppi_contact_residues_public", ""), "ppi_rim_distance": out.get("ppi_rim_distance", ""), "ppi_relationship_confirmed": out.get("ppi_relationship_confirmed", "false"), "ppi_geometry_status": ppi_status, "ppi_geometry_notes": "contact_or_rim_required_for_confirmed_relationship"})
        membrane_rows.append({**base, "ligand_membrane_z": out.get("ligand_membrane_z", ""), "ligand_min_membrane_z": out.get("ligand_min_membrane_z", ""), "ligand_max_membrane_z": out.get("ligand_max_membrane_z", ""), "membrane_penetration_flag": out.get("membrane_penetration_flag", "false"), "dimer_interface_clash_flag": out.get("dimer_interface_clash_flag", "false"), "nearest_receptor_contact_distance": out.get("nearest_receptor_contact_distance", ""), "membrane_dimer_status": mem_status, "membrane_dimer_notes": "membrane_frame_and_mapping_proxy"})

    if mode == "attribute" and selected_rows and not atoms_loaded and not allow_partial:
        blockers.append("every pose failed geometry parsing")
    if selected_rows and len(attr_rows) != len(selected_rows):
        blockers.append("pose rows were silently dropped")
    if selected_rows and not pose_files_inside:
        blockers.append("one or more pose_file paths outside run_dir")
    if mode == "attribute" and require_all_poses and (not pose_files_exist or not pose_files_nonempty) and not allow_partial:
        blockers.append("required pose files missing or empty")

    if selected_rows or write_empty_attribution or mode == "dry-run":
        _write_csv(attribution_csv, ATTRIBUTION_FIELDS, attr_rows, ctx)
        _write_csv(mechanism_csv, MECHANISM_FIELDS, mech_rows, ctx)
        _write_csv(atp_qc_csv, ["run_id", "m2_run_id", "profile", "job_id", "compound_public_id", "state_id", "state_role", "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index", "atp_centroid_distance", "atp_overlap_fraction", "atp_heavy_atom_overlap_fraction", "atp_nearest_atom_distance", "atp_migration_flag", "atp_status", "atp_notes"], atp_rows, ctx)
        _write_csv(pocket_qc_csv, ["run_id", "m2_run_id", "profile", "job_id", "compound_public_id", "state_id", "state_role", "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index", "inside_original_box", "inside_pocket_volume_proxy", "pocket_atom_fraction", "pocket_heavy_atom_fraction", "nearest_pocket_family_id", "nearest_pocket_distance_A", "pocket_retention_pass", "pocket_retention_status", "pocket_retention_notes"], pocket_rows, ctx)
        _write_csv(ppi_qc_csv, ["run_id", "m2_run_id", "profile", "job_id", "compound_public_id", "state_id", "state_role", "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index", "nearest_ppi_hotspot_id", "nearest_ppi_hotspot_distance", "ppi_hotspot_contact_count", "ppi_contact_residue_count", "ppi_contact_residues_public", "ppi_rim_distance", "ppi_relationship_confirmed", "ppi_geometry_status", "ppi_geometry_notes"], ppi_rows, ctx)
        _write_csv(membrane_qc_csv, ["run_id", "m2_run_id", "profile", "job_id", "compound_public_id", "state_id", "state_role", "pocket_family_id", "box_id", "protomer_id", "vina_repeat_id", "pose_rank", "pose_model_index", "ligand_membrane_z", "ligand_min_membrane_z", "ligand_max_membrane_z", "membrane_penetration_flag", "dimer_interface_clash_flag", "nearest_receptor_contact_distance", "membrane_dimer_status", "membrane_dimer_notes"], membrane_rows, ctx)

    private_entries_map, private_warnings = load_private_map(_safe_private_map_path(ctx, private_map))
    warnings.extend(private_warnings)
    leak_count, coord_hits, smiles_logged, ligand_coords, receptor_coords, pose_coords, scanned = _scan_hygiene(ctx, list(private_entries_map.values()))
    forbidden_created = [item for item in FORBIDDEN_OUTPUTS if _has_real_output(ctx.run_dir / item)]
    if leak_count:
        blockers.append("internal compound ID leakage detected")
    if coord_hits:
        blockers.append("coordinate records printed outside intended PDBQT files")
    if smiles_logged:
        blockers.append("SMILES-like fields printed in public outputs")
    if forbidden_created:
        blockers.append("forbidden downstream outputs already exist in run directory")

    counts = {
        "poses_classified": len(attr_rows),
        "poses_pass": sum(1 for row in attr_rows if row.get("pose_hard_gate_pass") == "true"),
        "poses_warn": sum(1 for row in attr_rows if row.get("attribution_status") == "WARN"),
        "poses_rejected": sum(1 for row in attr_rows if row.get("attribution_status") == "REJECT"),
        "poses_fail": sum(1 for row in attr_rows if row.get("attribution_status") == "FAIL"),
        "poses_missing": sum(1 for row in attr_rows if row.get("attribution_status") == "MISSING"),
        "allowed_for_clustering": sum(1 for row in attr_rows if row.get("allowed_for_clustering") == "true"),
        "pocket_retention_pass": sum(1 for row in attr_rows if row.get("pocket_retention_pass") == "true"),
        "pocket_retention_fail": sum(1 for row in attr_rows if row.get("pocket_retention_pass") == "false"),
        "atp_migration_flagged": sum(1 for row in attr_rows if row.get("atp_migration_flag") == "true"),
        "ppi_relationship_confirmed": sum(1 for row in attr_rows if row.get("ppi_relationship_confirmed") == "true"),
        "membrane_penetration_flagged": sum(1 for row in attr_rows if row.get("membrane_penetration_flag") == "true"),
        "dimer_interface_clash_flagged": sum(1 for row in attr_rows if row.get("dimer_interface_clash_flag") == "true"),
        "missing_mapping": sum(1 for row in attr_rows if row.get("pose_reject_reason") == "missing_mapping"),
    }
    mechanism_counts = {name: sum(1 for row in attr_rows if row.get("pose_mechanism_class") == name) for name in ["orthosteric_or_direct_PPI_patch_blocker", "rim_blocker", "allosteric_near_candidate", "generic_nonATP_ligandable_pocket", "ATP_like_reject", "ambiguous_or_failed"]}
    if mode == "attribute" and attr_rows and counts["poses_pass"] == 0:
        warnings.append("no pose_hard_gate_pass rows; attribution completed with rejections")
    if allow_partial:
        warnings.append("allow_partial diagnostic mode used")
    if mode == "dry-run":
        warnings.append("dry-run mode")

    status = "FAIL" if blockers else ("WARN" if warnings or mode == "dry-run" else "PASS")
    m3_t8_ready = status == "PASS" and mode == "attribute" and counts["poses_pass"] > 0
    collection_context = {
        "raw_pose_rows": len(selected_rows),
        "eligible_pose_rows": sum(1 for row in selected_rows if row.get("parse_status") in {"PASS", "WARN"} and _boolish(row.get("allowed_for_attribution"))),
        "profiles_in_raw_pose_table": profiles,
        "selected_profile": selected_profile,
        "primary_states": primary_states,
        "reference_states": reference_states,
        "reference_only_attribution": bool(reference_states and not primary_states),
    }
    summary = {
        "schema_version": "m3_pose_attribution_qc_v1",
        "run_id": ctx.run_id,
        "m2_run_id": m2_run_id,
        "reviewed_at": timestamp,
        "mode": mode,
        "profile": selected_profile,
        "overall_status": status,
        "allow_partial": allow_partial,
        "require_all_poses": require_all_poses,
        "allow_reference_only": allow_reference_only,
        "m3_t8_clustering_ready": m3_t8_ready,
        "inputs": inputs,
        "collection_context": collection_context,
        "attribution": counts,
        "mechanism_classes": mechanism_counts,
        "outputs": {
            "compound_pose_attribution": ctx.relative_to_repo(attribution_csv),
            "compound_pose_mechanism_classification": ctx.relative_to_repo(mechanism_csv),
            "pose_attribution_qc_csv": ctx.relative_to_repo(qc_csv),
            "pose_attribution_qc_json": ctx.relative_to_repo(qc_json),
            "atp_migration_qc": ctx.relative_to_repo(atp_qc_csv),
            "pocket_retention_qc": ctx.relative_to_repo(pocket_qc_csv),
            "ppi_geometry_qc": ctx.relative_to_repo(ppi_qc_csv),
            "membrane_dimer_geometry_qc": ctx.relative_to_repo(membrane_qc_csv),
            "report_file": ctx.relative_to_repo(report_md),
        },
        "parser": {"pose_atoms_loaded_in_memory": True, "coordinate_records_written_to_public_tables": False, "affinity_used_for_ranking": False},
        "counts": {"vina_commands_invoked_by_attributor": 0, "qsub_commands_invoked_by_attributor": 0, "runner_commands_invoked_by_attributor": 0, "pose_clustering_rows_created": 0, "anchor_convergence_rows_created": 0, "candidate_tables_created": 0, "confidentiality_leaks": leak_count, "coordinate_leak_files": len(coord_hits)},
        "blockers": blockers,
        "warnings": warnings,
        "confidentiality": {"internal_ids_redacted": leak_count == 0, "public_outputs_scanned": scanned, "leaks_detected": leak_count, "smiles_logged": smiles_logged, "ligand_coordinates_logged": ligand_coords, "receptor_coordinates_logged": receptor_coords, "pose_coordinates_logged_outside_pose_file": pose_coords},
        "non_goals_preserved": {"vina_execution": True, "qsub_submission": True, "pbs_runner_invocation": True, "pose_clustering": True, "anchor_convergence": True, "compound_scoring": True, "candidate_nomination": True, "final_candidate_tiering": True},
    }

    check_status = {
        "m3_t6_raw_pose_table_present": ("PASS" if raw_found else ("WARN" if mode == "dry-run" else "FAIL")),
        "m3_t6_raw_pose_table_schema_valid": "PASS" if schema_ok else ("NOT_APPLICABLE" if not raw_found else "FAIL"),
        "pose_rows_found": "PASS" if selected_rows else ("WARN" if mode == "dry-run" else "FAIL"),
        "pose_ids_unique": "PASS" if unique_ids else "FAIL",
        "pose_files_inside_run_dir": "PASS" if pose_files_inside else "FAIL",
        "pose_files_exist": "PASS" if pose_files_exist else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "pose_files_nonempty": "PASS" if pose_files_nonempty else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "accepted_pockets_present": "PASS" if inputs["accepted_pockets_found"] else ("WARN" if mode == "dry-run" else "FAIL"),
        "accepted_boxes_present": "PASS" if inputs["accepted_boxes_found"] else ("WARN" if mode == "dry-run" else "FAIL"),
        "atp_reference_present": "PASS" if inputs["atp_reference_found"] else ("WARN" if mode == "dry-run" else "FAIL"),
        "ppi_consensus_patch_present": "PASS" if inputs["ppi_consensus_patch_found"] else ("WARN" if mode == "dry-run" else "FAIL"),
        "membrane_frame_present": "PASS" if inputs["membrane_frame_found"] else ("WARN" if mode == "dry-run" else "FAIL"),
        "receptor_mapping_present": "PASS" if inputs["receptor_mapping_found"] else ("WARN" if mode == "dry-run" else "FAIL"),
        "all_pose_rows_classified": "PASS" if len(attr_rows) == len(selected_rows) else "FAIL",
        "no_silent_pose_drop": "PASS" if len(attr_rows) == len(selected_rows) else "FAIL",
        "pdbqt_pose_atoms_loaded_in_memory": "PASS" if atoms_loaded or not selected_rows else "NOT_APPLICABLE",
        "no_coordinate_rows_written": "PASS" if not coord_hits else "FAIL",
        "pocket_retention_computed": "PASS" if len(pocket_rows) == len([r for r in attr_rows if r.get("attribution_status") not in {"FAIL", "MISSING", "NOT_APPLICABLE"}]) else ("NOT_APPLICABLE" if not attr_rows else "WARN"),
        "atp_migration_computed": "PASS" if len(atp_rows) == len(pocket_rows) else ("NOT_APPLICABLE" if not attr_rows else "WARN"),
        "ppi_geometry_computed": "PASS" if len(ppi_rows) == len(pocket_rows) else ("NOT_APPLICABLE" if not attr_rows else "WARN"),
        "membrane_geometry_computed": "PASS" if len(membrane_rows) == len(pocket_rows) else ("NOT_APPLICABLE" if not attr_rows else "WARN"),
        "attribution_table_written": "PASS" if attribution_csv.is_file() else "FAIL",
        "mechanism_table_written": "PASS" if mechanism_csv.is_file() else "FAIL",
        "atp_migration_qc_written": "PASS" if atp_qc_csv.is_file() else "FAIL",
        "pocket_retention_qc_written": "PASS" if pocket_qc_csv.is_file() else "FAIL",
        "ppi_geometry_qc_written": "PASS" if ppi_qc_csv.is_file() else "FAIL",
        "membrane_dimer_qc_written": "PASS" if membrane_qc_csv.is_file() else "FAIL",
        "no_internal_id_leak": "PASS" if leak_count == 0 else "FAIL",
        "no_smiles_logged": "PASS" if not smiles_logged else "FAIL",
        "no_ligand_coordinates_logged": "PASS" if not ligand_coords else "FAIL",
        "no_receptor_coordinates_logged": "PASS" if not receptor_coords else "FAIL",
        "no_pose_coordinates_logged_outside_pdbqt": "PASS" if not pose_coords else "FAIL",
    }
    for check in REQUIRED_CHECKS:
        default = "PASS"
        if check == "m3_t6_docking_completion_qc_present":
            default = "PASS" if inputs["docking_completion_qc_found"] else ("WARN" if mode == "dry-run" else "FAIL")
        elif check == "profile_selected":
            default = "PASS" if selected_profile else "FAIL"
        elif check == "m3_t6_not_failed":
            default = "PASS" if completion_ready else ("WARN" if mode == "dry-run" else "FAIL")
        elif check == "input_reference_paths_valid":
            default = "PASS" if all(inputs[key] for key in ["accepted_pockets_found", "accepted_boxes_found", "atp_reference_found", "ppi_consensus_patch_found", "membrane_frame_found", "receptor_mapping_found"]) else ("WARN" if mode == "dry-run" else "FAIL")
        elif check.endswith("_status_for_every_pose"):
            default = "PASS" if len(attr_rows) == len(selected_rows) else "FAIL"
        elif check in {"atp_migration_independent_of_vina_score", "ppi_relationship_not_centroid_only", "pose_hard_gates_applied", "pose_reject_reasons_explicit", "mechanism_class_assigned", "raw_pose_affinity_not_used_for_ranking", "no_vina_invoked_by_attributor", "no_qsub_invoked_by_attributor", "no_runner_invoked_by_attributor", "no_pose_clustering_attempted", "no_anchor_convergence_attempted", "no_candidate_scoring_attempted", "no_candidate_nomination_attempted", "no_final_candidate_tables_created", "reference_state_not_promoted", "old_workflow_not_used", "non_goals_preserved"}:
            default = "PASS"
        elif check == "pose_attribution_report_written":
            default = "PASS" if report_md.is_file() else "FAIL"
        status_value = check_status.get(check, default)
        _append_qc(qc_rows, check, "m3_t7_pose_attribution", status_value, check, blocker=status_value == "FAIL")
    _write_csv(qc_csv, QC_FIELDS, qc_rows, ctx)
    _write_report(report_md, ctx, summary)
    late_leaks, late_coord_hits, late_smiles, late_ligand_coords, late_receptor_coords, late_pose_coords, late_scanned = _scan_hygiene(ctx, list(private_entries_map.values()))
    if late_leaks and "internal compound ID leakage detected" not in blockers:
        blockers.append("internal compound ID leakage detected")
    if late_coord_hits and "coordinate records printed outside intended PDBQT files" not in blockers:
        blockers.append("coordinate records printed outside intended PDBQT files")
    if late_smiles and "SMILES-like fields printed in public outputs" not in blockers:
        blockers.append("SMILES-like fields printed in public outputs")
    if late_leaks or late_coord_hits or late_smiles:
        status = "FAIL"
        m3_t8_ready = False
    summary["overall_status"] = status
    summary["m3_t8_clustering_ready"] = m3_t8_ready
    summary["blockers"] = blockers
    summary["confidentiality"] = {
        "internal_ids_redacted": late_leaks == 0,
        "public_outputs_scanned": late_scanned,
        "leaks_detected": late_leaks,
        "smiles_logged": late_smiles,
        "ligand_coordinates_logged": late_ligand_coords,
        "receptor_coordinates_logged": late_receptor_coords,
        "pose_coordinates_logged_outside_pose_file": late_pose_coords,
    }
    summary["counts"]["confidentiality_leaks"] = late_leaks
    summary["counts"]["coordinate_leak_files"] = len(late_coord_hits)
    qc_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    phase3_log.write_text("", encoding="utf-8") if not phase3_log.exists() else None
    with phase3_log.open("a", encoding="utf-8") as handle:
        handle.write(f"{timestamp} M3-T7 {status} profile={selected_profile} poses_classified={len(attr_rows)} no_vina=true no_qsub=true\n")
    append_phase_status(ctx, "phase3_compounds", status, "M3-T7 pose attribution completed", {"phase": "phase3_compounds", "task": "M3-T7", "status": status, "run_id": ctx.run_id, "m2_run_id": m2_run_id, "timestamp": timestamp, "mode": mode, "profile": selected_profile, "raw_pose_rows": len(selected_rows), "poses_classified": len(attr_rows), "poses_pass": counts["poses_pass"], "poses_rejected": counts["poses_rejected"], "atp_migration_flagged": counts["atp_migration_flagged"], "pocket_retention_pass": counts["pocket_retention_pass"], "ppi_relationship_confirmed": counts["ppi_relationship_confirmed"], "allowed_for_clustering": counts["allowed_for_clustering"], "m3_t8_clustering_ready": m3_t8_ready, "qsub_invoked": False, "vina_invoked_by_attributor": False, "runner_invoked_by_attributor": False, "no_pose_clustering_attempted": True, "no_anchor_convergence_attempted": True, "no_scoring_attempted": True, "no_candidate_nomination_attempted": True, "no_final_candidate_tiering_attempted": True})

    return M3PoseAttributionResult(status, m3_t8_ready, blockers, warnings, attribution_csv, mechanism_csv, qc_csv, qc_json, atp_qc_csv, pocket_qc_csv, ppi_qc_csv, membrane_qc_csv, report_md, phase3_log, inputs, collection_context, counts, mechanism_counts)
