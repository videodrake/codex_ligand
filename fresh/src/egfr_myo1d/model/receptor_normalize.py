"""M1 Phase 4: EGFR receptor normalization.

Implements the policy from ``milestone1_foundation_codex_handoff_v0_5.md`` §14:

- detects Case A (explicit A/B), Case B (same-chain duplicate dimer),
  Case C (single-chain monomer; not valid for dimer-only main analysis)
- splits Case B into protomer A and protomer B by duplicate-atom-identity
- writes three normalized PDBs per state:
    full_frame_explicit_AB.pdb              (chains A and B preserved)
    dockable_<crop>_explicit_AB.pdb         (cropped to gates.yaml range)
    runtime_offset_receptor_only.pdb        (protomer B residues + 1000)
- writes the canonical residue mapping CSV (handoff §14.3)
- writes the receptor normalization audit CSV
- writes the receptor manifest JSON
- enforces 3GT8_raw role: ``crystallographic_reference_control_not_primary_membrane_state``
- V924R-like residues: WARN only, never mutated
- monomer (Case C): WARN in ``codex_dev``, FAIL in ``hpc_strict``
- ambiguous chain layout: WARN with quarantine note
- HETATM caps preserved; lipid/water/ion HETATM dropped from dockable PDB
- writes are confined to ``ctx.run_dir``
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Iterable

from egfr_myo1d.core.logging_utils import append_phase_status
from egfr_myo1d.core.manifest import load_yaml_config, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError
from egfr_myo1d.io.hashing import describe_file
from egfr_myo1d.io.residue_mapping import MAPPING_CSV_COLUMNS, write_residue_mapping
from egfr_myo1d.model.receptor_qc import (
    CAP_RESNAMES,
    DEFAULT_WARN_MUTATIONS,
    NON_RECEPTOR_HETATM,
    WarnMutation,
    compute_residue_audit_rows,
    detect_normalization_case,
    detect_warn_mutations,
    is_receptor_atom,
    split_duplicate_chain,
)
from egfr_myo1d.myo1d.pdb_writer import write_pdb_atoms
from egfr_myo1d.structure.pdb_parser import AtomRecord, PDBParseError, PDBStructure, parse_pdb


RECEPTOR_AUDIT_CSV_COLUMNS = [
    "state",
    "fixture_role",
    "chain_id",
    "protomer_id",
    "residue_number",
    "insertion_code",
    "residue_name",
    "record_type",
    "atom_count",
    "in_dockable_crop",
    "biological",
    "is_cap",
    "is_warn_mutation",
    "classification",
    "notes",
]


PRIMARY_STATES = ("EGFR_160-185", "EGFR_170-200")
REFERENCE_CONTROL_STATES = ("3GT8_raw",)
KNOWN_STATES = PRIMARY_STATES + REFERENCE_CONTROL_STATES


def resolve_role(state_id: str) -> str:
    if state_id in REFERENCE_CONTROL_STATES:
        return "crystallographic_reference_control_not_primary_membrane_state"
    if state_id in PRIMARY_STATES:
        return "primary_membrane_validated_state"
    return "unknown_state"


@dataclass
class NormalizedReceptor:
    state_id: str
    role: str
    source_file: str
    case: str
    observed_chains: list[str] = field(default_factory=list)
    protomer_count: int = 0
    n_residues_a: int = 0
    n_residues_b: int = 0
    full_frame_pdb: Path | None = None
    dockable_pdb: Path | None = None
    runtime_offset_pdb: Path | None = None
    mapping_csv: Path | None = None
    audit_csv: Path | None = None
    manifest_json: Path | None = None
    v924r_warn: bool = False
    v924r_handled: str = "warn_only_not_mutated"
    warn_mutation_hits: list[dict[str, Any]] = field(default_factory=list)
    dockable_crop: tuple[int, int] = (669, 1014)
    runtime_offset: int = 1000
    warnings: list[str] = field(default_factory=list)
    status: str = "PASS"
    profile: str = "codex_dev"


# ---------------------------------------------------------------------------
# PDB raw-line rewriting helpers
# ---------------------------------------------------------------------------


def _rewrite_atom(
    atom: AtomRecord,
    chain_id: str | None = None,
    residue_number: int | None = None,
) -> AtomRecord:
    """Return a copy of ``atom`` with ``chain_id`` and/or ``residue_number``
    rewritten in both the dataclass fields and the underlying ``raw_line``.
    PDB columns: chain_id at [21], residue_number at [22:26]."""
    line = atom.raw_line
    if len(line) < 26:
        line = line.ljust(26)

    new_chain = atom.chain_id if chain_id is None else chain_id[:1] or "_"
    new_resseq = atom.residue_number if residue_number is None else residue_number

    if chain_id is not None:
        line = line[:21] + new_chain + line[22:]
    if residue_number is not None:
        line = line[:22] + "{0:>4d}".format(new_resseq) + line[26:]

    return replace(
        atom,
        chain_id=new_chain,
        residue_number=new_resseq,
        raw_line=line,
    )


# ---------------------------------------------------------------------------
# Configuration helpers
# ---------------------------------------------------------------------------


def _parse_range_string(text: str) -> tuple[int, int]:
    s, e = text.split("-", 1)
    return int(s), int(e)


def load_receptor_gate(ctx: RunContext) -> dict[str, Any]:
    gates_path = ctx.repo_root / "fresh" / "configs" / "gates.yaml"
    config = load_yaml_config(gates_path)
    return config.get("receptor", {})


# ---------------------------------------------------------------------------
# Mapping row construction
# ---------------------------------------------------------------------------


def _residue_mapping_rows(
    state_id: str,
    role: str,
    source_path: Path,
    a_atoms: list[AtomRecord],
    b_atoms: list[AtomRecord],
    runtime_offset: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for protomer_id, atoms, runtime_chain in (
        ("A", a_atoms, "A"),
        ("B", b_atoms, "B"),
    ):
        grouped: dict[tuple[str, int, str, str], int] = {}
        for atom in atoms:
            key = (atom.chain_id, atom.residue_number, atom.insertion_code, atom.resname)
            grouped[key] = grouped.get(key, 0) + 1
        for (source_chain, source_resseq, source_icode, source_resname), atom_count in grouped.items():
            runtime_resseq = (
                source_resseq + runtime_offset if protomer_id == "B" else source_resseq
            )
            rows.append(
                {
                    "state": state_id,
                    "source_file": str(source_path),
                    "protomer_id": protomer_id,
                    "source_chain": source_chain,
                    "source_resseq": source_resseq,
                    "source_icode": source_icode or "",
                    "source_resname": source_resname,
                    "runtime_chain": runtime_chain,
                    "runtime_resseq": runtime_resseq,
                    "atom_count": atom_count,
                    "role": role,
                }
            )
    rows.sort(
        key=lambda r: (r["protomer_id"], r["source_chain"], int(r["source_resseq"]))
    )
    return rows


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------


def _write_pdb(ctx: RunContext, target: Path, atoms: list[AtomRecord]) -> Path:
    safe = ctx.require_within_run_dir(target)
    safe.parent.mkdir(parents=True, exist_ok=True)
    write_pdb_atoms(safe, atoms)
    return safe


def _write_audit_csv(ctx: RunContext, target: Path, rows: list[dict[str, Any]]) -> Path:
    safe = ctx.require_within_run_dir(target)
    safe.parent.mkdir(parents=True, exist_ok=True)
    with safe.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=RECEPTOR_AUDIT_CSV_COLUMNS, extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            row_to_write = dict(row)
            for col in ("in_dockable_crop", "biological", "is_cap", "is_warn_mutation"):
                if col in row_to_write:
                    row_to_write[col] = str(row_to_write[col]).lower()
            writer.writerow(row_to_write)
    return safe


def _write_manifest(
    ctx: RunContext,
    target: Path,
    report: NormalizedReceptor,
) -> Path:
    safe = ctx.require_within_run_dir(target)

    def desc(path: Path | None) -> dict[str, Any]:
        if path is None:
            return {"present": False, "sha256": None, "size_bytes": None}
        return describe_file(path)

    payload = {
        "state_id": report.state_id,
        "role": report.role,
        "source_file": report.source_file,
        "source_sha256": desc(Path(report.source_file)).get("sha256"),
        "case": report.case,
        "observed_chains": report.observed_chains,
        "protomer_count": report.protomer_count,
        "n_residues_a": report.n_residues_a,
        "n_residues_b": report.n_residues_b,
        "dockable_crop": list(report.dockable_crop),
        "runtime_offset": report.runtime_offset,
        "v924r_warn": report.v924r_warn,
        "v924r_handled": report.v924r_handled,
        "warn_mutation_hits": report.warn_mutation_hits,
        "outputs": {
            "full_frame_pdb": str(report.full_frame_pdb) if report.full_frame_pdb else None,
            "full_frame_sha256": desc(report.full_frame_pdb).get("sha256"),
            "dockable_pdb": str(report.dockable_pdb) if report.dockable_pdb else None,
            "dockable_sha256": desc(report.dockable_pdb).get("sha256"),
            "runtime_offset_pdb": str(report.runtime_offset_pdb) if report.runtime_offset_pdb else None,
            "runtime_offset_sha256": desc(report.runtime_offset_pdb).get("sha256"),
            "mapping_csv": str(report.mapping_csv) if report.mapping_csv else None,
            "audit_csv": str(report.audit_csv) if report.audit_csv else None,
        },
        "warnings": report.warnings,
        "status": report.status,
        "profile": report.profile,
        "score_bonus_allowed": False,
        "timestamp": now_iso(),
    }
    write_json(safe, payload, ctx)
    return safe


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def normalize_receptor(
    ctx: RunContext,
    source_pdb: Path,
    state_id: str,
    profile: str = "codex_dev",
    strict: bool = False,
    dockable_crop: tuple[int, int] | None = None,
    runtime_offset: int | None = None,
    warn_mutations: Iterable[WarnMutation] | None = None,
) -> NormalizedReceptor:
    """Normalize a receptor PDB per handoff §14 and emit canonical M1 outputs.

    Returns a :class:`NormalizedReceptor` report.
    """
    if profile not in {"codex_dev", "hpc_strict"}:
        raise ValueError("profile must be codex_dev or hpc_strict; got: {0}".format(profile))

    role = resolve_role(state_id)
    gate = load_receptor_gate(ctx)
    if dockable_crop is None:
        dockable_crop = _parse_range_string(gate.get("dockable_crop_default", "669-1014"))
    if runtime_offset is None:
        runtime_offset = int(gate.get("runtime_offset_second_protomer", 1000))
    warns = list(warn_mutations) if warn_mutations is not None else list(DEFAULT_WARN_MUTATIONS)

    report = NormalizedReceptor(
        state_id=state_id,
        role=role,
        source_file=str(source_pdb),
        case="unknown",
        dockable_crop=dockable_crop,
        runtime_offset=runtime_offset,
        profile=profile,
    )

    src_path = Path(source_pdb)
    if not src_path.is_file():
        report.warnings.append("source_missing: {0}".format(src_path))
        report.status = "FAIL" if profile == "hpc_strict" or strict else "WARN"
        _finalize(ctx, report)
        return report

    try:
        structure = parse_pdb(src_path)
    except PDBParseError as exc:
        report.warnings.append("parse_error: {0}".format(exc))
        report.status = "FAIL"
        _finalize(ctx, report)
        return report

    report.observed_chains = list(structure.chain_ids)
    case = detect_normalization_case(structure)
    report.case = case

    # V924R-like detection (full structure, before splitting)
    hits = detect_warn_mutations(structure.atoms, warns)
    if hits:
        report.v924r_warn = True
        report.warn_mutation_hits = hits
        report.warnings.append("warn_mutation_detected_not_mutated")
        if report.status == "PASS":
            report.status = "WARN"

    # Receptor-only filter (drop lipid/water/ion HETATM globally for normalization
    # outputs except where caps need to ride along; lipids dropped from dockable below)
    receptor_atoms = [a for a in structure.atoms if is_receptor_atom(a)]

    # Determine protomer split
    if case == "A_explicit_AB":
        a_atoms = [a for a in receptor_atoms if a.chain_id == "A"]
        b_atoms = [a for a in receptor_atoms if a.chain_id == "B"]
        report.protomer_count = 2
    elif case == "B_duplicate_chain":
        a_raw, b_raw = split_duplicate_chain(list(receptor_atoms))
        a_atoms = [_rewrite_atom(a, chain_id="A") for a in a_raw]
        b_atoms = [_rewrite_atom(a, chain_id="B") for a in b_raw]
        report.protomer_count = 2
        report.warnings.append("case_B_duplicate_chain_split_into_A_B")
        if report.status == "PASS":
            report.status = "WARN"
    elif case == "C_monomer":
        a_atoms = list(receptor_atoms)
        b_atoms = []
        report.protomer_count = 1
        report.warnings.append("single_chain_monomer_not_valid_for_dimer_only_main_analysis")
        if profile == "hpc_strict" or strict:
            report.status = "FAIL"
            _finalize(ctx, report)
            return report
        if report.status == "PASS":
            report.status = "WARN"
    else:  # ambiguous
        a_atoms = list(receptor_atoms)
        b_atoms = []
        report.protomer_count = len(report.observed_chains)
        report.warnings.append(
            "ambiguous_chain_layout: chains={0}".format(report.observed_chains)
        )
        if profile == "hpc_strict" or strict:
            report.status = "FAIL"
            _finalize(ctx, report)
            return report
        if report.status == "PASS":
            report.status = "WARN"

    # Compute residue counts per protomer
    a_residue_keys = {(a.chain_id, a.residue_number, a.insertion_code, a.resname) for a in a_atoms}
    b_residue_keys = {(a.chain_id, a.residue_number, a.insertion_code, a.resname) for a in b_atoms}
    report.n_residues_a = len(a_residue_keys)
    report.n_residues_b = len(b_residue_keys)

    out_dir = ctx.run_dir / "normalized" / "receptors"

    # Full frame explicit A/B (chains as labeled, original residue numbers)
    if state_id in REFERENCE_CONTROL_STATES:
        full_frame_path = out_dir / "{0}_explicit_AB.pdb".format(state_id)
    else:
        full_frame_path = out_dir / "{0}_full_frame_explicit_AB.pdb".format(state_id)
    full_atoms = list(a_atoms) + list(b_atoms)
    report.full_frame_pdb = _write_pdb(ctx, full_frame_path, full_atoms)

    # Dockable crop with lipid/water/ion drop
    crop_start, crop_end = dockable_crop
    dockable_atoms = [
        a for a in full_atoms
        if crop_start <= a.residue_number <= crop_end
        and not (a.record_type == "HETATM" and a.resname in NON_RECEPTOR_HETATM)
    ]
    if state_id in REFERENCE_CONTROL_STATES:
        # For 3GT8_raw, omit dockable crop emission per handoff §14.2 (kinase-only crystal,
        # requires_membrane_gate=false). Still leave entry empty in manifest.
        report.dockable_pdb = None
    else:
        dockable_path = out_dir / "{0}_dockable_{1}_{2}_explicit_AB.pdb".format(
            state_id, crop_start, crop_end
        )
        report.dockable_pdb = _write_pdb(ctx, dockable_path, dockable_atoms)

    # Runtime offset receptor only (protomer A unchanged, protomer B residues + offset)
    runtime_atoms = list(a_atoms) + [
        _rewrite_atom(a, residue_number=a.residue_number + runtime_offset) for a in b_atoms
    ]
    runtime_path = out_dir / "{0}_runtime_offset_receptor_only.pdb".format(state_id)
    report.runtime_offset_pdb = _write_pdb(ctx, runtime_path, runtime_atoms)

    # Mapping CSV
    mapping_rows = _residue_mapping_rows(
        state_id=state_id,
        role=role,
        source_path=src_path,
        a_atoms=a_atoms,
        b_atoms=b_atoms,
        runtime_offset=runtime_offset,
    )
    mapping_csv = ctx.qc_dir / "{0}_receptor_mapping.csv".format(state_id)
    write_residue_mapping(mapping_csv, mapping_rows, ctx=ctx)
    report.mapping_csv = mapping_csv

    # Audit CSV
    audit_rows = compute_residue_audit_rows(
        a_atoms,
        state_id=state_id,
        protomer_id="A",
        role=role,
        dockable_crop=dockable_crop,
        warn_mutations=warns,
    ) + compute_residue_audit_rows(
        b_atoms,
        state_id=state_id,
        protomer_id="B",
        role=role,
        dockable_crop=dockable_crop,
        warn_mutations=warns,
    )
    audit_path = ctx.qc_dir / "{0}_receptor_normalization_audit.csv".format(state_id)
    report.audit_csv = _write_audit_csv(ctx, audit_path, audit_rows)

    # Cap presence informational
    has_caps = any(a.resname in CAP_RESNAMES for a in a_atoms + b_atoms)
    if has_caps:
        report.warnings.append("caps_preserved")

    _finalize(ctx, report)
    return report


def _finalize(ctx: RunContext, report: NormalizedReceptor) -> None:
    manifest_path = ctx.manifest_dir / "{0}_receptor_manifest.json".format(report.state_id)
    report.manifest_json = _write_manifest(ctx, manifest_path, report)
    append_phase_status(
        ctx,
        phase="prepare-receptor",
        status=report.status,
        message=(
            "prepare-receptor state={0} case={1} status={2} "
            "protomer_count={3} v924r_warn={4} warnings={5}".format(
                report.state_id,
                report.case,
                report.status,
                report.protomer_count,
                report.v924r_warn,
                len(report.warnings),
            )
        ),
        details={
            "state_id": report.state_id,
            "role": report.role,
            "case": report.case,
            "observed_chains": report.observed_chains,
            "protomer_count": report.protomer_count,
            "n_residues_a": report.n_residues_a,
            "n_residues_b": report.n_residues_b,
            "v924r_warn": report.v924r_warn,
            "v924r_handled": report.v924r_handled,
            "dockable_crop": list(report.dockable_crop),
            "runtime_offset": report.runtime_offset,
            "score_bonus_allowed": False,
            "profile": report.profile,
            "manifest": str(manifest_path),
            "outputs": {
                "full_frame_pdb": str(report.full_frame_pdb) if report.full_frame_pdb else None,
                "dockable_pdb": str(report.dockable_pdb) if report.dockable_pdb else None,
                "runtime_offset_pdb": str(report.runtime_offset_pdb) if report.runtime_offset_pdb else None,
                "mapping_csv": str(report.mapping_csv) if report.mapping_csv else None,
                "audit_csv": str(report.audit_csv) if report.audit_csv else None,
            },
        },
    )
