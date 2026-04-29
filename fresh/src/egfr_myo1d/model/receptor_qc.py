"""Receptor-side QC helpers for M1 Phase 4.

Detects:
- V924R-like crystallization mutation (handoff §13)
- Residues outside the dockable crop range
- Capped HETATM residues vs lipid/water/ion HETATM (drop policy)
- Single-chain vs explicit-AB vs duplicate-chain detection (Cases A/B/C)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from egfr_myo1d.structure.pdb_parser import AtomRecord, PDBStructure


# Lipid/water/ion HETATM resnames that should be dropped from dockable receptor
# but retained in full_frame. CHARMM-GUI POPC + POPS membrane composition.
LIPID_RESNAMES = frozenset(
    {"POPC", "POPS", "POPE", "POPG", "DPPC", "DPPE", "DOPC", "CHL", "CHL1", "CER"}
)
WATER_RESNAMES = frozenset({"HOH", "WAT", "TIP3", "H2O"})
ION_RESNAMES = frozenset({"NA", "CL", "K", "MG", "ZN", "CA", "SOD", "CLA", "POT"})
NON_RECEPTOR_HETATM = LIPID_RESNAMES | WATER_RESNAMES | ION_RESNAMES

# Caps to keep
CAP_RESNAMES = frozenset({"ACE", "NME"})


@dataclass(frozen=True)
class WarnMutation:
    residue_number: int
    expected: str
    warn_if: str
    label: str = ""


# Default V924R-like check from handoff §13 / §3 of phase prompt
DEFAULT_WARN_MUTATIONS = (
    WarnMutation(residue_number=924, expected="VAL", warn_if="ARG", label="3GT8_V924R"),
)


def detect_normalization_case(structure: PDBStructure) -> str:
    """Classify the receptor input shape.

    Returns one of:
    - ``A_explicit_AB``  : chains A and B both present
    - ``B_duplicate_chain``: single chain (or chain pair) with duplicate atom
      identities indicating two protomers stored in one chain ID
    - ``C_monomer``      : single chain, no duplicates → not valid for dimer
    - ``ambiguous``      : multi-chain but neither A nor B present
    """
    chains = set(structure.chain_ids)
    duplicate_count = structure.duplicate_atom_identity_count()

    if {"A", "B"}.issubset(chains):
        return "A_explicit_AB"
    if duplicate_count > 0:
        return "B_duplicate_chain"
    if len(chains) == 1:
        return "C_monomer"
    return "ambiguous"


def split_duplicate_chain(
    atoms: list[AtomRecord],
) -> tuple[list[AtomRecord], list[AtomRecord]]:
    """Split atoms in a same-chain duplicate dimer into two protomers.

    Strategy: track the first occurrence of each
    ``(chain, resseq, icode, resname, atom_name)`` tuple as protomer A; any
    repeat is protomer B. Original order within each protomer is preserved.
    """
    seen: dict[tuple[str, int, str, str, str], int] = {}
    a_atoms: list[AtomRecord] = []
    b_atoms: list[AtomRecord] = []
    for atom in atoms:
        key = (
            atom.chain_id,
            atom.residue_number,
            atom.insertion_code,
            atom.resname,
            atom.atom_name,
        )
        if key in seen:
            b_atoms.append(atom)
        else:
            seen[key] = 1
            a_atoms.append(atom)
    return a_atoms, b_atoms


def is_receptor_atom(atom: AtomRecord) -> bool:
    """True if the atom is part of the receptor (not lipid/water/ion)."""
    if atom.record_type == "ATOM":
        return True
    # HETATM: keep caps and any standard-named amino acid; drop lipids/water/ions
    if atom.resname in CAP_RESNAMES:
        return True
    if atom.resname in NON_RECEPTOR_HETATM:
        return False
    return True  # unknown HETATM resname: keep, recorded as biological


def detect_warn_mutations(
    atoms: Iterable[AtomRecord],
    warn_mutations: Iterable[WarnMutation] = DEFAULT_WARN_MUTATIONS,
) -> list[dict]:
    """Return a list of mutation hits across all chains/protomers.

    Hits do not mutate the receptor. They are reported only.
    """
    hits: list[dict] = []
    seen: set[tuple[str, int, str]] = set()
    for atom in atoms:
        for mut in warn_mutations:
            if atom.residue_number != mut.residue_number:
                continue
            if atom.resname == mut.warn_if:
                key = (atom.chain_id, atom.residue_number, atom.resname)
                if key in seen:
                    continue
                seen.add(key)
                hits.append(
                    {
                        "label": mut.label,
                        "residue_number": atom.residue_number,
                        "chain_id": atom.chain_id,
                        "expected_resname": mut.expected,
                        "observed_resname": atom.resname,
                        "warn_if": mut.warn_if,
                    }
                )
    return hits


def compute_residue_audit_rows(
    atoms: Iterable[AtomRecord],
    state_id: str,
    protomer_id: str,
    role: str,
    dockable_crop: tuple[int, int],
    warn_mutations: Iterable[WarnMutation] = DEFAULT_WARN_MUTATIONS,
) -> list[dict]:
    """Per-residue audit records used by the normalization audit CSV."""
    crop_start, crop_end = dockable_crop
    grouped: dict[tuple[str, int, str, str], list[AtomRecord]] = {}
    for atom in atoms:
        key = (atom.chain_id, atom.residue_number, atom.insertion_code, atom.resname)
        grouped.setdefault(key, []).append(atom)

    mut_targets = {
        (mut.residue_number, mut.warn_if): mut for mut in warn_mutations
    }

    rows: list[dict] = []
    for (chain_id, residue_number, icode, resname), residue_atoms in grouped.items():
        record_types = sorted({a.record_type for a in residue_atoms})
        in_crop = crop_start <= residue_number <= crop_end
        is_cap = resname in CAP_RESNAMES
        is_warn_mut = (residue_number, resname) in mut_targets
        if any(a.record_type == "HETATM" and a.resname in NON_RECEPTOR_HETATM for a in residue_atoms):
            biological = False
            classification = "non_receptor_hetatm"
        elif is_cap:
            biological = False
            classification = "cap"
        else:
            biological = True
            classification = "biological_receptor"
        notes = []
        if is_warn_mut:
            notes.append("warn_mutation:{0}->{1}".format(
                mut_targets[(residue_number, resname)].expected, resname
            ))
        if "HETATM" in record_types and biological and not is_cap:
            notes.append("standard_residue_as_HETATM_kept")
        rows.append(
            {
                "state": state_id,
                "fixture_role": role,
                "chain_id": chain_id,
                "protomer_id": protomer_id,
                "residue_number": residue_number,
                "insertion_code": icode or "_",
                "residue_name": resname,
                "record_type": "/".join(record_types),
                "atom_count": len(residue_atoms),
                "in_dockable_crop": in_crop,
                "biological": biological,
                "is_cap": is_cap,
                "is_warn_mutation": is_warn_mut,
                "classification": classification,
                "notes": ";".join(notes),
            }
        )
    return rows
