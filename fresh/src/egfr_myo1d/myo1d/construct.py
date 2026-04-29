"""MYO1D construct preparation and audit helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from egfr_myo1d.myo1d.pdb_writer import select_atoms_by_residue_range, write_pdb_atoms
from egfr_myo1d.structure.myo1d_annotation import expand_multi_range, expand_range
from egfr_myo1d.structure.pdb_parser import AtomRecord, PDBStructure, parse_pdb


CAP_RESNAMES = {"ACE", "NME"}
STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


def residue_annotation(residue_number: int, contract_myo1d: dict[str, Any]) -> str:
    active_face = set(contract_myo1d["active_face"])
    support = set(contract_myo1d["sheet12_support"])
    buffer = set(contract_myo1d["structural_buffer"])
    cap = set(contract_myo1d["contact_monitoring_cap"])
    if residue_number in active_face:
        return "primary_active_face_annotation"
    if residue_number in support:
        return "sheet12_support_annotation"
    if residue_number in buffer:
        return "structural_buffer"
    if residue_number in cap:
        return "contact_monitoring_cap/noise_monitor"
    return "other"


def residue_rows_for_construct(
    structure: PDBStructure,
    construct_id: str,
    source_path: Path,
    output_role: str,
    contract_myo1d: dict[str, Any],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, int, str, str]] = set()
    for residue in structure.residues:
        key = (residue.chain_id, residue.residue_number, residue.insertion_code, residue.resname)
        if key in seen:
            continue
        seen.add(key)
        is_cap = residue.resname in CAP_RESNAMES
        biological = residue.resname in STANDARD_AA
        rows.append(
            {
                "construct_id": construct_id,
                "source_path": str(source_path),
                "chain_id": residue.chain_id,
                "residue_number": residue.residue_number,
                "insertion_code": residue.insertion_code,
                "residue_name": residue.resname,
                "record_type": "/".join(residue.record_types_present),
                "biological": biological,
                "is_cap": is_cap,
                "annotation_class": "cap" if is_cap else residue_annotation(residue.residue_number, contract_myo1d),
                "output_role": output_role,
                "score_bonus_allowed": False,
                "notes": "standard amino acid encoded as HETATM retained" if biological and "HETATM" in residue.record_types_present else "",
            }
        )
    return rows


def detect_terminal_artifact(
    structure: PDBStructure,
    construct_id: str,
    source_path: Path,
    allow_fixture_warning: bool,
) -> dict[str, Any]:
    residue_numbers = sorted({residue.residue_number for residue in structure.residues if residue.resname not in CAP_RESNAMES})
    cap_records = [residue for residue in structure.residues if residue.resname in CAP_RESNAMES]
    first = residue_numbers[0] if residue_numbers else None
    last = residue_numbers[-1] if residue_numbers else None
    has_bad_start = first == 962
    has_tail = any(number >= 1002 for number in residue_numbers)
    terminal_prominent = sum(1 for number in residue_numbers if number >= 998) >= max(3, len(residue_numbers) // 3) if residue_numbers else False
    verdict = "PASS"
    notes = []
    if has_bad_start:
        verdict = "WARN" if allow_fixture_warning else "FAIL"
        notes.append("starts at 962; missing 955-960 structural buffer")
    if has_tail:
        verdict = "WARN" if verdict != "FAIL" else verdict
        notes.append("contains 1002-1006 tail/noise-zone residues")
    if terminal_prominent:
        verdict = "WARN" if verdict != "FAIL" else verdict
        notes.append("terminal 998+ residues are prominent")
    return {
        "construct_id": construct_id,
        "n_terminal_residue": first,
        "c_terminal_residue": last,
        "has_cap_records": bool(cap_records),
        "has_free_artifact_terminal": has_bad_start or (not cap_records and first is not None and first > 955),
        "terminal_contact_region_policy": "998+ contact/noise monitoring; 962-start is fixture-only quarantine",
        "verdict": verdict,
        "notes": "; ".join(notes),
    }


def prepare_myo1d_construct(
    source_path: Path,
    output_path: Path,
    construct_id: str,
    residue_range: tuple[int, int],
    contract_myo1d: dict[str, Any],
    output_role: str,
    include_caps: bool = True,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    structure = parse_pdb(source_path)
    include_resnames = CAP_RESNAMES if include_caps else set()
    selected_atoms = select_atoms_by_residue_range(
        structure.atoms,
        residue_range[0],
        residue_range[1],
        include_resnames=include_resnames,
    )
    write_stats = write_pdb_atoms(output_path, selected_atoms)
    selected_structure = PDBStructure(
        path=output_path,
        atoms=tuple(selected_atoms),
        residues=tuple(
            residue
            for residue in structure.residues
            if residue_range[0] <= residue.residue_number <= residue_range[1] or residue.resname in include_resnames
        ),
    )
    audit_rows = residue_rows_for_construct(
        selected_structure,
        construct_id,
        source_path,
        output_role,
        contract_myo1d,
    )
    terminal_audit = detect_terminal_artifact(
        selected_structure,
        construct_id,
        source_path,
        allow_fixture_warning=output_role != "production_primary",
    )
    output_record = {
        "construct_id": construct_id,
        "source_path": str(source_path),
        "output_path": str(output_path),
        "residue_range": list(residue_range),
        "output_role": output_role,
        "atom_count": write_stats["atom_count"],
        "residue_numbers": sorted({atom.residue_number for atom in selected_atoms}),
        "score_bonus_allowed": False,
    }
    return output_record, audit_rows, terminal_audit


def validate_active_face_presence(prepared_structure: PDBStructure, contract_myo1d: dict[str, Any]) -> tuple[str, list[dict[str, Any]]]:
    residue_numbers = {residue.residue_number for residue in prepared_structure.residues}
    rows: list[dict[str, Any]] = []
    groups = {
        "primary_active_face": contract_myo1d["active_face"],
        "sheet12_support": contract_myo1d["sheet12_support"],
        "structural_buffer": contract_myo1d["structural_buffer"],
        "contact_monitoring_cap": [number for number in contract_myo1d["contact_monitoring_cap"] if number <= 1001],
    }
    status = "PASS"
    for group, residues in groups.items():
        for residue_number in residues:
            present = residue_number in residue_numbers
            if not present:
                status = "FAIL"
            rows.append(
                {
                    "construct_id": contract_myo1d["primary_construct"]["id"],
                    "residue_number": residue_number,
                    "annotation_group": group,
                    "present": present,
                    "score_bonus_allowed": False,
                    "status": "PASS" if present else "FAIL",
                    "notes": "annotation only; no scoring bonus",
                }
            )
    return status, rows

