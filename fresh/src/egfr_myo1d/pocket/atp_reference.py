"""M2.5 ATP-site reference builder.

This module creates ATP-site reference CSVs for later non-ATP pocket gates. It
does not run fpocket, P2Rank, Vina, PyRosetta, docking, scoring, or candidate
nomination.
"""

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
from egfr_myo1d.structure.pdb_parser import AtomRecord, PDBParseError, parse_pdb


M2_5_PHASE = "m2.5-atp-site-reference"
ATP_REF_ROOT = "phase2_pockets/atp_reference"
ATP_TABLE_ROOT = "phase2_pockets/tables"
REFERENCE_FIELDS = [
    "atp_site_id",
    "reference_source_id",
    "reference_source_path",
    "reference_mode",
    "ligand_resname",
    "ligand_chain_id",
    "ligand_residue_number",
    "uniprot_residue_number",
    "residue_name",
    "atp_site_region",
    "evidence_source",
    "evidence_status",
    "warnings",
]
MAPPING_FIELDS = [
    "state_id",
    "receptor_id",
    "egfr_chain_id",
    "egfr_protomer_id",
    "uniprot_residue_number",
    "reference_residue_number",
    "receptor_residue_number",
    "reference_residue_name",
    "receptor_residue_name",
    "mapping_status",
    "residue_name_match",
    "ca_x",
    "ca_y",
    "ca_z",
    "heavy_atom_count",
    "evidence_status",
    "warnings",
]
CENTROID_FIELDS = [
    "state_id",
    "receptor_id",
    "egfr_chain_id",
    "egfr_protomer_id",
    "centroid_method",
    "atp_centroid_x",
    "atp_centroid_y",
    "atp_centroid_z",
    "residue_count",
    "atom_count",
    "reference_radius_angstrom",
    "evidence_status",
    "warnings",
]
QC_FIELDS = ["check", "status", "details"]
NON_GOALS_CONFIRMED = [
    "no_pyrosetta_import",
    "no_pyrosetta_docking_or_relaxation",
    "no_lightdock_execution",
    "no_vina_execution",
    "no_fpocket_or_p2rank_runtime",
    "no_pocket_discovery",
    "no_compound_docking_or_scoring",
    "no_candidate_nomination",
    "no_scheduler_submission",
    "no_cleanup_deletion",
]


@dataclass
class AtpReferenceReport:
    run_id: str
    mode: str
    profile: str
    status: str
    reference_count: int = 0
    mapping_count: int = 0
    centroid_count: int = 0
    warnings: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    reference_csv: Path | None = None
    residue_mapping_csv: Path | None = None
    centroid_csv: Path | None = None
    mapping_qc_csv: Path | None = None
    tables_reference_csv: Path | None = None
    status_json: Path | None = None
    manifest_path: Path | None = None
    summary_report_path: Path | None = None


def _repo_path(ctx: RunContext, path: Path) -> str:
    return ctx.relative_to_repo(path)


def _resolve_fresh_or_run_path(ctx: RunContext, value: str | Path | None, default_path: Path) -> Path:
    raw = Path(value) if value is not None else default_path
    if raw.is_absolute():
        path = raw.resolve()
    else:
        path = (ctx.repo_root / raw).resolve()
        if not _is_within(path, ctx.fresh_root) and not _is_within(path, ctx.run_dir):
            path = (ctx.run_dir / raw).resolve()
    if not (_is_within(path, ctx.fresh_root) or _is_within(path, ctx.run_dir)):
        raise RunContextError("input path escapes fresh/ and run dir: {0}".format(path))
    return path


def _is_within(child: Path, parent: Path) -> bool:
    try:
        Path(child).resolve().relative_to(Path(parent).resolve())
        return True
    except ValueError:
        return False


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str], ctx: RunContext) -> None:
    safe = ctx.require_within_run_dir(path)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _load_receptor_states(ctx: RunContext) -> dict[str, Any]:
    return load_yaml_config(ctx.repo_root / "fresh" / "configs" / "receptor_states.yaml")


def _state_role_map(ctx: RunContext) -> dict[str, str]:
    data = _load_receptor_states(ctx)
    return {str(k): str(v) for k, v in data.get("state_policy", {}).items()}


def _default_states(ctx: RunContext, include_reference_states: bool) -> list[str]:
    data = _load_receptor_states(ctx)
    states = list(data.get("primary_membrane_validated_states", []))
    if include_reference_states:
        states.extend(data.get("reference_control_states", []))
    return list(dict.fromkeys(states))


def _distance(a: AtomRecord, b: AtomRecord) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def _is_heavy(atom: AtomRecord) -> bool:
    element = (atom.element or atom.atom_name[:1]).upper()
    return element != "H"


def _reference_pdb_from_run_or_config(ctx: RunContext) -> Path | None:
    candidates = [
        ctx.run_dir / "normalized" / "receptors" / "3GT8_raw_full_frame_explicit_AB.pdb",
        ctx.run_dir / "normalized" / "receptors" / "3GT8_raw_dockable_669_1014_explicit_AB.pdb",
        ctx.repo_root / "fresh" / "data" / "raw" / "receptors" / "3GT8_raw.pdb",
    ]
    for candidate in candidates:
        if candidate.is_file() and (_is_within(candidate, ctx.run_dir) or _is_within(candidate, ctx.fresh_root)):
            return candidate.resolve()
    return None


def _reference_mode(mode: str, synthetic_fixture: bool) -> str:
    return "synthetic_fixture" if synthetic_fixture else mode


def _ligand_identity_index(ctx: RunContext) -> dict[tuple[str, int], tuple[int, str]]:
    grouped: dict[tuple[str, int], set[tuple[int, str]]] = {}
    for mapping_path in sorted(ctx.qc_dir.glob("*_receptor_mapping.csv")):
        for row in read_residue_mapping(mapping_path):
            source_resseq = _safe_int(row.get("source_resseq"))
            runtime_resseq = _safe_int(row.get("runtime_resseq"))
            if source_resseq is None:
                continue
            identity = (source_resseq, str(row.get("source_resname", "")))
            source_chain = str(row.get("source_chain", ""))
            runtime_chain = str(row.get("runtime_chain", ""))
            if source_chain:
                grouped.setdefault((source_chain, source_resseq), set()).add(identity)
            if runtime_chain and runtime_resseq is not None:
                grouped.setdefault((runtime_chain, runtime_resseq), set()).add(identity)
    index: dict[tuple[str, int], tuple[int, str]] = {}
    for key, identities in grouped.items():
        source_numbers = {item[0] for item in identities}
        if len(source_numbers) == 1:
            index[key] = sorted(identities)[0]
    return index


def _ligand_reference_rows(
    ctx: RunContext,
    config: dict[str, Any],
    reference_pdb: Path | None,
    synthetic_fixture: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    if reference_pdb is None or not reference_pdb.is_file():
        return [], [], ["no_ligand_bearing_reference_pdb_found"]
    allowed = {str(item).upper() for item in config.get("allowed_ligand_resnames", [])}
    distance_cutoff = float(config.get("ligand_neighbor_distance_angstrom", 5.0))
    try:
        structure = parse_pdb(reference_pdb)
    except PDBParseError as exc:
        return [], [], ["reference_pdb_parse_failed:{0}".format(exc)]
    ligand_atoms = [
        atom
        for atom in structure.atoms
        if atom.record_type == "HETATM" and atom.resname.upper() in allowed and _is_heavy(atom)
    ]
    if not ligand_atoms:
        return [], [], ["no_allowed_atp_like_ligand_found:{0}".format(reference_pdb)]
    ligand_groups: dict[tuple[str, str, int], list[AtomRecord]] = {}
    for atom in ligand_atoms:
        ligand_groups.setdefault((atom.resname, atom.chain_id, atom.residue_number), []).append(atom)
    receptor_atoms = [
        atom for atom in structure.atoms if atom.record_type == "ATOM" and _is_heavy(atom)
    ]
    identity_index = _ligand_identity_index(ctx)
    reference_rows: list[dict[str, Any]] = []
    centroid_rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    unmapped_neighbor_count = 0
    seen: set[tuple[int, str, int]] = set()
    for (ligand_resname, ligand_chain, ligand_residue), atoms in sorted(ligand_groups.items()):
        centroid_rows.append(
            {
                "state_id": "reference_ligand",
                "receptor_id": "3GT8_raw_reference",
                "egfr_chain_id": ligand_chain,
                "egfr_protomer_id": ligand_chain,
                "centroid_method": "ligand_heavy_atom_centroid",
                "atp_centroid_x": _format_float(_mean([atom.x for atom in atoms])),
                "atp_centroid_y": _format_float(_mean([atom.y for atom in atoms])),
                "atp_centroid_z": _format_float(_mean([atom.z for atom in atoms])),
                "residue_count": 0,
                "atom_count": len(atoms),
                "reference_radius_angstrom": config.get("reference_radius_angstrom", 10.0),
                "evidence_status": "REFERENCE_ONLY",
                "warnings": "reference_control_not_primary_state",
            }
        )
        for receptor_atom in receptor_atoms:
            min_distance = min(_distance(receptor_atom, ligand_atom) for ligand_atom in atoms)
            if min_distance > distance_cutoff:
                continue
            identity = identity_index.get((receptor_atom.chain_id, receptor_atom.residue_number))
            if identity is None:
                unmapped_neighbor_count += 1
                continue
            uniprot_residue, residue_name = identity
            key = (uniprot_residue, ligand_chain, ligand_residue)
            if key in seen:
                continue
            seen.add(key)
            reference_rows.append(
                {
                    "atp_site_id": "atp_ref_{0:03d}".format(len(reference_rows) + 1),
                    "reference_source_id": "3GT8_raw_ligand_reference",
                    "reference_source_path": _repo_path(ctx, reference_pdb),
                    "reference_mode": _reference_mode("ligand_bearing_reference", synthetic_fixture),
                    "ligand_resname": ligand_resname,
                    "ligand_chain_id": ligand_chain,
                    "ligand_residue_number": ligand_residue,
                    "uniprot_residue_number": uniprot_residue,
                    "residue_name": residue_name or receptor_atom.resname,
                    "atp_site_region": "unknown_or_configured",
                    "evidence_source": "ligand_neighbor_heavy_atom_cutoff_{0}".format(distance_cutoff),
                    "evidence_status": "PASS",
                    "warnings": "ligand_reference_control_not_primary_state;m1_mapping_restored_ligand_neighbor_identity",
                }
            )
    if unmapped_neighbor_count:
        warnings.append(
            "ligand_neighbor_residues_without_m1_mapping:{0}".format(
                unmapped_neighbor_count
            )
        )
    if not reference_rows:
        warnings.append("no_ligand_neighbor_residues_with_m1_mapping")
    return reference_rows, centroid_rows, warnings


def _configured_minimum(config: dict[str, Any]) -> int:
    section = config.get("configured_residue_set", {}) or {}
    value = section.get("minimum_residue_count", config.get("minimum_configured_residue_count", 1))
    parsed = _safe_int(value)
    return parsed if parsed is not None and parsed > 0 else 1


def _configured_reference_rows(
    ctx: RunContext,
    config: dict[str, Any],
    synthetic_fixture: bool,
) -> tuple[list[dict[str, Any]], list[str]]:
    section = config.get("configured_residue_set", {}) or {}
    residues = section.get("residues", []) or []
    if not residues:
        return [], ["no_configured_atp_residue_set_found"]
    rows: list[dict[str, Any]] = []
    for residue in residues:
        rows.append(
            {
                "atp_site_id": "atp_ref_config_{0:03d}".format(len(rows) + 1),
                "reference_source_id": section.get("source_id", "configured_residue_set"),
                "reference_source_path": _repo_path(ctx, ctx.repo_root / "fresh" / "configs" / "atp_reference.yaml"),
                "reference_mode": _reference_mode("configured_residue_set", synthetic_fixture),
                "ligand_resname": "",
                "ligand_chain_id": "",
                "ligand_residue_number": "",
                "uniprot_residue_number": residue.get("uniprot_residue_number", ""),
                "residue_name": residue.get("residue_name", ""),
                "atp_site_region": residue.get("atp_site_region", "unknown_or_configured"),
                "evidence_source": section.get("evidence_source", "configured_residue_set"),
                "evidence_status": "PASS",
                "warnings": section.get("notes", ""),
            }
        )
    return rows, []


def _merge_reference_rows(ligand_rows: list[dict[str, Any]], configured_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    merged: dict[int, dict[str, Any]] = {}
    for row in ligand_rows + configured_rows:
        try:
            residue = int(str(row.get("uniprot_residue_number", "")).strip())
        except ValueError:
            continue
        if residue not in merged:
            merged[residue] = dict(row)
        else:
            prior = merged[residue]
            prior["evidence_source"] = ";".join(
                sorted({prior.get("evidence_source", ""), row.get("evidence_source", "")} - {""})
            )
            prior["warnings"] = ";".join(
                sorted({prior.get("warnings", ""), row.get("warnings", "")} - {""})
            )
            if not prior.get("residue_name"):
                prior["residue_name"] = row.get("residue_name", "")
            if prior.get("reference_mode") != row.get("reference_mode"):
                modes = {prior.get("reference_mode", ""), row.get("reference_mode", "")}
                if "synthetic_fixture" in modes:
                    prior["reference_mode"] = "synthetic_fixture"
                elif "ligand_bearing_reference" in modes:
                    prior["reference_mode"] = "ligand_bearing_reference"
                else:
                    prior["reference_mode"] = "configured_residue_set"
                prior["warnings"] = ";".join(
                    sorted(
                        {
                            prior.get("warnings", ""),
                            row.get("warnings", ""),
                            "multiple_reference_evidence_modes_available",
                        }
                        - {""}
                    )
                )
    rows = [merged[key] for key in sorted(merged)]
    for idx, row in enumerate(rows, start=1):
        row["atp_site_id"] = "atp_ref_{0:03d}".format(idx)
    return rows


def _format_float(value: float | None) -> str:
    if value is None:
        return ""
    return "{0:.3f}".format(value)


def _mean(values: Iterable[float]) -> float | None:
    items = list(values)
    if not items:
        return None
    return sum(items) / float(len(items))


def _state_receptor_path(ctx: RunContext, state_id: str) -> Path | None:
    candidates = [
        ctx.run_dir / "normalized" / "receptors" / "{0}_dockable_669_1014_explicit_AB.pdb".format(state_id),
        ctx.run_dir / "normalized" / "receptors" / "{0}_full_frame_explicit_AB.pdb".format(state_id),
        ctx.run_dir / "normalized" / "receptors" / "{0}_runtime_offset_receptor_only.pdb".format(state_id),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return ctx.require_within_run_dir(candidate)
    return None


def _atoms_by_identity(structure_atoms: Iterable[AtomRecord]) -> dict[tuple[str, int], list[AtomRecord]]:
    grouped: dict[tuple[str, int], list[AtomRecord]] = {}
    for atom in structure_atoms:
        if atom.record_type == "ATOM" and _is_heavy(atom):
            grouped.setdefault((atom.chain_id, atom.residue_number), []).append(atom)
    return grouped


def _ca_atom(atoms: list[AtomRecord]) -> AtomRecord | None:
    for atom in atoms:
        if atom.atom_name == "CA":
            return atom
    return None


def _map_reference_to_state(
    ctx: RunContext,
    state_id: str,
    state_role: str,
    reference_rows: list[dict[str, Any]],
    config: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    mapping_path = ctx.qc_dir / "{0}_receptor_mapping.csv".format(state_id)
    receptor_path = _state_receptor_path(ctx, state_id)
    warnings: list[str] = []
    if not mapping_path.is_file():
        return [], [], ["missing_receptor_mapping:{0}".format(state_id)]
    mapping_rows = read_residue_mapping(mapping_path)
    receptor_atoms: dict[tuple[str, int], list[AtomRecord]] = {}
    if receptor_path is not None:
        try:
            receptor_atoms = _atoms_by_identity(parse_pdb(receptor_path).atoms)
        except PDBParseError as exc:
            warnings.append("state_receptor_parse_failed:{0}:{1}".format(state_id, exc))
    else:
        warnings.append("missing_state_receptor_pdb:{0}".format(state_id))
    by_source: dict[int, list[dict[str, Any]]] = {}
    for row in mapping_rows:
        try:
            source_residue = int(str(row.get("source_resseq", "")).strip())
        except ValueError:
            continue
        by_source.setdefault(source_residue, []).append(row)
    mapping_output: list[dict[str, Any]] = []
    centroid_atoms: dict[tuple[str, str], list[AtomRecord]] = {}
    residue_counts: dict[tuple[str, str], set[int]] = {}
    missing_count = 0
    mismatch_count = 0
    for ref in reference_rows:
        try:
            residue_number = int(str(ref.get("uniprot_residue_number", "")).strip())
        except ValueError:
            continue
        candidates = by_source.get(residue_number, [])
        if not candidates:
            missing_count += 1
            mapping_output.append(
                _missing_mapping_row(state_id, state_role, ref, "missing_residue_in_receptor_mapping")
            )
            continue
        for item in sorted(candidates, key=lambda x: (x.get("protomer_id", ""), x.get("source_chain", ""))):
            source_chain = str(item.get("source_chain", ""))
            runtime_chain = str(item.get("runtime_chain", ""))
            source_resseq = _safe_int(item.get("source_resseq"))
            runtime_resseq = _safe_int(item.get("runtime_resseq"))
            atoms = receptor_atoms.get((source_chain, source_resseq or -1), [])
            coordinate_source = "source_chain_residue"
            if not atoms and runtime_resseq is not None:
                atoms = receptor_atoms.get((runtime_chain, runtime_resseq), [])
                coordinate_source = "runtime_chain_residue"
            ca = _ca_atom(atoms)
            ref_name = str(ref.get("residue_name", ""))
            receptor_name = str(item.get("source_resname", ""))
            name_match = True if not ref_name else ref_name.upper() == receptor_name.upper()
            row_warnings = []
            if not atoms:
                row_warnings.append("mapped_residue_coordinates_missing")
            if coordinate_source == "runtime_chain_residue":
                row_warnings.append("coordinates_found_by_runtime_numbering")
            if not name_match:
                mismatch_count += 1
                row_warnings.append("residue_name_mismatch")
            status = "PASS" if atoms and name_match else "WARN"
            mapping_output.append(
                {
                    "state_id": state_id,
                    "receptor_id": "EGFR_MYO1D_membrane_compatible_inactive_dimer",
                    "egfr_chain_id": source_chain or runtime_chain,
                    "egfr_protomer_id": str(item.get("protomer_id", "")),
                    "uniprot_residue_number": residue_number,
                    "reference_residue_number": residue_number,
                    "receptor_residue_number": item.get("source_resseq", ""),
                    "reference_residue_name": ref_name,
                    "receptor_residue_name": receptor_name,
                    "mapping_status": status,
                    "residue_name_match": name_match,
                    "ca_x": _format_float(ca.x if ca else None),
                    "ca_y": _format_float(ca.y if ca else None),
                    "ca_z": _format_float(ca.z if ca else None),
                    "heavy_atom_count": len(atoms),
                    "evidence_status": "REFERENCE_CONTROL_ONLY" if "reference_control" in state_role else status,
                    "warnings": ";".join(row_warnings),
                }
            )
            if atoms:
                key = (source_chain or runtime_chain, str(item.get("protomer_id", "")))
                centroid_atoms.setdefault(key, []).extend(atoms)
                residue_counts.setdefault(key, set()).add(residue_number)
    centroid_rows: list[dict[str, Any]] = []
    for (chain_id, protomer_id), atoms in sorted(centroid_atoms.items()):
        centroid_rows.append(
            {
                "state_id": state_id,
                "receptor_id": "EGFR_MYO1D_membrane_compatible_inactive_dimer",
                "egfr_chain_id": chain_id,
                "egfr_protomer_id": protomer_id,
                "centroid_method": "mapped_residue_heavy_atom_centroid",
                "atp_centroid_x": _format_float(_mean([atom.x for atom in atoms])),
                "atp_centroid_y": _format_float(_mean([atom.y for atom in atoms])),
                "atp_centroid_z": _format_float(_mean([atom.z for atom in atoms])),
                "residue_count": len(residue_counts.get((chain_id, protomer_id), set())),
                "atom_count": len(atoms),
                "reference_radius_angstrom": config.get("reference_radius_angstrom", 10.0),
                "evidence_status": "REFERENCE_CONTROL_ONLY" if "reference_control" in state_role else "PASS",
                "warnings": "reference_control_not_primary_state" if "reference_control" in state_role else "",
            }
        )
    if missing_count:
        warnings.append("state_missing_atp_residue_mappings:{0}:{1}".format(state_id, missing_count))
    if mismatch_count:
        warnings.append("state_residue_name_mismatches:{0}:{1}".format(state_id, mismatch_count))
    if not centroid_rows:
        warnings.append("state_atp_centroid_unavailable:{0}".format(state_id))
    return mapping_output, centroid_rows, warnings


def _safe_int(value: Any) -> int | None:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def _missing_mapping_row(state_id: str, state_role: str, ref: dict[str, Any], reason: str) -> dict[str, Any]:
    return {
        "state_id": state_id,
        "receptor_id": "EGFR_MYO1D_membrane_compatible_inactive_dimer",
        "egfr_chain_id": "",
        "egfr_protomer_id": "",
        "uniprot_residue_number": ref.get("uniprot_residue_number", ""),
        "reference_residue_number": ref.get("uniprot_residue_number", ""),
        "receptor_residue_number": "",
        "reference_residue_name": ref.get("residue_name", ""),
        "receptor_residue_name": "",
        "mapping_status": "MISSING",
        "residue_name_match": "",
        "ca_x": "",
        "ca_y": "",
        "ca_z": "",
        "heavy_atom_count": 0,
        "evidence_status": "REFERENCE_CONTROL_ONLY" if "reference_control" in state_role else "WARN",
        "warnings": reason,
    }


def _m2_4_warning(ctx: RunContext) -> list[str]:
    path = ctx.run_dir / "output" / "ppi" / "ppi_consensus_patch.csv"
    if not path.is_file():
        return ["m2_4_consensus_patch_missing_reference_can_still_build"]
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except OSError:
        return ["m2_4_consensus_patch_unreadable_reference_can_still_build"]
    if not rows:
        return ["m2_4_consensus_patch_has_no_runtime_contact_evidence"]
    return []


def _status(blockers: list[str], warnings: list[str], reference_rows: list[dict[str, Any]], centroid_rows: list[dict[str, Any]]) -> str:
    if blockers:
        return "FAIL"
    if not reference_rows:
        return "FAIL"
    if warnings or not centroid_rows:
        return "PASS_WITH_WARNINGS"
    return "PASS"


def _write_summary(ctx: RunContext, report: AtpReferenceReport) -> Path:
    target = ctx.require_within_run_dir(ctx.reports_dir / "m2_5_atp_reference.md")
    lines = [
        "# M2.5 ATP-site Reference",
        "",
        "Status: {0}".format(report.status),
        "Reference rows: {0}".format(report.reference_count),
        "Mapping rows: {0}".format(report.mapping_count),
        "Centroid rows: {0}".format(report.centroid_count),
        "",
        "The ATP reference is a hard-gate reference for later non-ATP pocket work. No pocket discovery or compound scoring was run.",
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


def build_atp_reference(
    ctx: RunContext,
    mode: str = "smoke_input",
    profile: str = "codex_dev",
    states: Iterable[str] | None = None,
    include_reference_states: bool = True,
    config_path: Path | None = None,
    reference_pdb: Path | None = None,
    synthetic_fixture: bool = False,
) -> AtpReferenceReport:
    """Build M2.5 ATP reference outputs under the run directory."""

    ctx.create_directories()
    config_file = _resolve_fresh_or_run_path(
        ctx,
        config_path,
        ctx.repo_root / "fresh" / "configs" / "atp_reference.yaml",
    )
    if not config_file.is_file():
        raise RunContextError("ATP reference config missing: {0}".format(config_file))
    config = load_yaml_config(config_file)
    ref_pdb = (
        _resolve_fresh_or_run_path(ctx, reference_pdb, reference_pdb)
        if reference_pdb is not None
        else _reference_pdb_from_run_or_config(ctx)
    )
    warnings: list[str] = []
    blockers: list[str] = []
    ligand_rows, ligand_centroid_rows, ligand_warnings = _ligand_reference_rows(
        ctx, config, ref_pdb, synthetic_fixture
    )
    warnings.extend(ligand_warnings)
    configured_rows, config_warnings = _configured_reference_rows(ctx, config, synthetic_fixture)
    warnings.extend(config_warnings)
    reference_rows = _merge_reference_rows(ligand_rows, configured_rows)
    if not reference_rows:
        blockers.append("no_atp_reference_evidence_ligand_or_configured_residue_set")

    state_list = list(states) if states else _default_states(ctx, include_reference_states)
    roles = _state_role_map(ctx)
    all_mapping_rows: list[dict[str, Any]] = []
    all_centroid_rows: list[dict[str, Any]] = list(ligand_centroid_rows)
    for state_id in state_list:
        mapping_rows, centroid_rows, state_warnings = _map_reference_to_state(
            ctx, state_id, roles.get(state_id, ""), reference_rows, config
        )
        all_mapping_rows.extend(mapping_rows)
        all_centroid_rows.extend(centroid_rows)
        warnings.extend(state_warnings)

    warnings.extend(_m2_4_warning(ctx))
    primary_states = set(_default_states(ctx, include_reference_states=False))
    primary_centroids = [
        row for row in all_centroid_rows if row.get("state_id") in primary_states
    ]
    if reference_rows and not primary_centroids:
        blockers.append("no_primary_state_atp_centroid_available")

    configured_minimum = _configured_minimum(config)
    if (
        configured_rows
        and not ligand_rows
        and len(configured_rows) < configured_minimum
        and not synthetic_fixture
    ):
        blockers.append(
            "configured_atp_residue_set_insufficient_for_hard_gate:{0}:{1}".format(
                len(configured_rows), configured_minimum
            )
        )
    elif configured_rows and len(configured_rows) < configured_minimum:
        warnings.append(
            "configured_atp_residue_set_below_minimum:{0}:{1}".format(
                len(configured_rows), configured_minimum
            )
        )

    atp_dir = ctx.require_within_run_dir(ctx.run_dir / ATP_REF_ROOT)
    tables_dir = ctx.require_within_run_dir(ctx.run_dir / ATP_TABLE_ROOT)
    atp_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    reference_csv = atp_dir / "atp_site_reference.csv"
    residue_mapping_csv = atp_dir / "atp_site_residue_mapping.csv"
    centroid_csv = atp_dir / "atp_site_centroid_by_state.csv"
    mapping_qc_csv = atp_dir / "atp_site_mapping_qc.csv"
    tables_reference_csv = tables_dir / "atp_site_reference.csv"
    _write_csv(reference_csv, reference_rows, REFERENCE_FIELDS, ctx)
    _write_csv(tables_reference_csv, reference_rows, REFERENCE_FIELDS, ctx)
    _write_csv(residue_mapping_csv, all_mapping_rows, MAPPING_FIELDS, ctx)
    _write_csv(centroid_csv, all_centroid_rows, CENTROID_FIELDS, ctx)
    _write_csv(
        mapping_qc_csv,
        _mapping_qc_rows(all_mapping_rows, state_list),
        QC_FIELDS,
        ctx,
    )

    status = _status(blockers, warnings, reference_rows, primary_centroids)
    report = AtpReferenceReport(
        run_id=ctx.run_id,
        mode=mode,
        profile=profile,
        status=status,
        reference_count=len(reference_rows),
        mapping_count=len(all_mapping_rows),
        centroid_count=len(all_centroid_rows),
        warnings=sorted(set(warnings)),
        blockers=blockers,
        reference_csv=reference_csv,
        residue_mapping_csv=residue_mapping_csv,
        centroid_csv=centroid_csv,
        mapping_qc_csv=mapping_qc_csv,
        tables_reference_csv=tables_reference_csv,
    )
    summary = _write_summary(ctx, report)
    report.summary_report_path = summary
    gate_policy = config.get("gate_policy", {}) or {}
    status_json = atp_dir / "atp_reference_status.json"
    manifest = ctx.manifest_dir / "m2_5_atp_reference_manifest.json"
    payload = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": M2_5_PHASE,
        "mode": mode,
        "profile": profile,
        "status": status,
        "states": state_list,
        "include_reference_states": include_reference_states,
        "config": _repo_path(ctx, config_file),
        "reference_pdb": _repo_path(ctx, ref_pdb) if ref_pdb else None,
        "reference_count": len(reference_rows),
        "mapping_count": len(all_mapping_rows),
        "centroid_count": len(all_centroid_rows),
        "gate_policy": {
            "atp_overlap_fraction_reject_threshold": gate_policy.get("atp_overlap_fraction_reject_threshold", 0.15),
            "atp_centroid_distance_reject_threshold_angstrom": gate_policy.get("atp_centroid_distance_reject_threshold_angstrom", 10.0),
            "future_non_atp_pass_rule": gate_policy.get(
                "rule",
                "future non_atp_pass=false if atp_overlap_fraction >= threshold OR atp_centroid_distance <= threshold",
            ),
        },
        "execution_allowed": False,
        "pyrosetta_imported": False,
        "fpocket_executed": False,
        "p2rank_executed": False,
        "vina_executed": False,
        "pocket_discovery_executed": False,
        "compound_scoring_executed": False,
        "candidate_nomination_executed": False,
        "outputs": _output_records(ctx, reference_csv, residue_mapping_csv, centroid_csv, mapping_qc_csv, tables_reference_csv, summary, status_json),
        "warnings": sorted(set(warnings)),
        "blockers": blockers,
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
    }
    write_json(status_json, payload, ctx)
    write_json(manifest, payload, ctx)
    report.status_json = status_json
    report.manifest_path = manifest

    phase_status = "FAIL" if status == "FAIL" else ("WARN" if status == "PASS_WITH_WARNINGS" else "PASS")
    append_phase_status(
        ctx,
        M2_5_PHASE,
        phase_status,
        "M2.5 ATP reference build {0} with {1} reference residue(s)".format(status, len(reference_rows)),
        {"reference_count": len(reference_rows), "pocket_discovery_executed": False},
    )
    append_job_status(
        ctx,
        "m2_5_atp_reference_local",
        phase_status,
        details={"message": "local ATP reference construction only"},
    )
    return report


def _mapping_qc_rows(rows: list[dict[str, Any]], states: list[str]) -> list[dict[str, Any]]:
    qc_rows: list[dict[str, Any]] = []
    for state in states:
        state_rows = [row for row in rows if row.get("state_id") == state]
        missing = [row for row in state_rows if row.get("mapping_status") == "MISSING"]
        warn = [row for row in state_rows if row.get("mapping_status") == "WARN"]
        status = "WARN" if missing or warn else ("PASS" if state_rows else "WARN")
        qc_rows.append(
            {
                "check": "state_mapping:{0}".format(state),
                "status": status,
                "details": "rows={0};missing={1};warn={2}".format(len(state_rows), len(missing), len(warn)),
            }
        )
    return qc_rows


def _output_records(
    ctx: RunContext,
    reference_csv: Path,
    residue_mapping_csv: Path,
    centroid_csv: Path,
    mapping_qc_csv: Path,
    tables_reference_csv: Path,
    summary: Path,
    status_json: Path,
) -> dict[str, str]:
    return {
        "atp_site_reference_csv": _repo_path(ctx, reference_csv),
        "atp_site_residue_mapping_csv": _repo_path(ctx, residue_mapping_csv),
        "atp_site_centroid_by_state_csv": _repo_path(ctx, centroid_csv),
        "atp_site_mapping_qc_csv": _repo_path(ctx, mapping_qc_csv),
        "tables_atp_site_reference_csv": _repo_path(ctx, tables_reference_csv),
        "summary_report": _repo_path(ctx, summary),
        "status_json": _repo_path(ctx, status_json),
    }
