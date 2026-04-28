"""Small PDB parser for Task 3 structural input QC.

The parser intentionally preserves PDB identity fields. It does not renumber,
normalize chains, infer biological assemblies, or discard HETATM records.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


class PDBParseError(ValueError):
    """Raised when a PDB file cannot be parsed for QC."""


@dataclass(frozen=True)
class AtomRecord:
    record_type: str
    serial: int | None
    atom_name: str
    altloc: str
    resname: str
    chain_id: str
    residue_number: int
    insertion_code: str
    x: float
    y: float
    z: float
    element: str
    raw_line: str
    source_line_number: int


@dataclass(frozen=True)
class ResidueSummary:
    chain_id: str
    residue_number: int
    insertion_code: str
    resname: str
    record_types_present: tuple[str, ...]
    atom_count: int
    has_backbone_N: bool
    has_backbone_CA: bool
    has_backbone_C: bool
    has_backbone_O: bool


@dataclass(frozen=True)
class PDBStructure:
    path: Path
    atoms: tuple[AtomRecord, ...]
    residues: tuple[ResidueSummary, ...]

    @property
    def chain_ids(self) -> tuple[str, ...]:
        return tuple(sorted({atom.chain_id for atom in self.atoms}))

    def atoms_for_chain(self, chain_id: str) -> tuple[AtomRecord, ...]:
        return tuple(atom for atom in self.atoms if atom.chain_id == chain_id)

    def residues_for_chain(self, chain_id: str) -> tuple[ResidueSummary, ...]:
        return tuple(residue for residue in self.residues if residue.chain_id == chain_id)

    def residue_numbers_for_chain(self, chain_id: str) -> set[int]:
        return {residue.residue_number for residue in self.residues if residue.chain_id == chain_id}

    def residue_resnames(self, residue_number: int) -> set[str]:
        return {
            residue.resname
            for residue in self.residues
            if residue.residue_number == residue_number
        }

    def duplicate_atom_identity_count(self) -> int:
        seen: set[tuple[str, int, str, str, str]] = set()
        duplicates = 0
        for atom in self.atoms:
            key = (
                atom.chain_id,
                atom.residue_number,
                atom.insertion_code,
                atom.resname,
                atom.atom_name,
            )
            if key in seen:
                duplicates += 1
            seen.add(key)
        return duplicates


def _parse_int(value: str) -> int | None:
    text = value.strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def parse_atom_line(line: str, source_line_number: int) -> AtomRecord:
    record_type = line[0:6].strip()
    if record_type not in {"ATOM", "HETATM"}:
        raise PDBParseError(f"line {source_line_number}: unsupported record {record_type}")

    serial = _parse_int(line[6:11])
    atom_name = line[12:16].strip()
    altloc = line[16:17].strip()
    resname = line[17:20].strip()
    chain_id = line[21:22].strip() or "_"
    residue_number = _parse_int(line[22:26])
    insertion_code = line[26:27].strip()

    if residue_number is None:
        raise PDBParseError(f"line {source_line_number}: invalid residue number")

    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError as exc:
        raise PDBParseError(f"line {source_line_number}: invalid non-numeric coordinates") from exc

    element = line[76:78].strip() if len(line) >= 78 else ""
    return AtomRecord(
        record_type=record_type,
        serial=serial,
        atom_name=atom_name,
        altloc=altloc,
        resname=resname,
        chain_id=chain_id,
        residue_number=residue_number,
        insertion_code=insertion_code,
        x=x,
        y=y,
        z=z,
        element=element,
        raw_line=line.rstrip("\n"),
        source_line_number=source_line_number,
    )


def summarize_residues(atoms: Iterable[AtomRecord]) -> tuple[ResidueSummary, ...]:
    grouped: dict[tuple[str, int, str, str], list[AtomRecord]] = {}
    for atom in atoms:
        grouped.setdefault(
            (atom.chain_id, atom.residue_number, atom.insertion_code, atom.resname), []
        ).append(atom)

    residues: list[ResidueSummary] = []
    for (chain_id, residue_number, insertion_code, resname), residue_atoms in grouped.items():
        atom_names = {atom.atom_name.strip() for atom in residue_atoms}
        record_types = tuple(sorted({atom.record_type for atom in residue_atoms}))
        residues.append(
            ResidueSummary(
                chain_id=chain_id,
                residue_number=residue_number,
                insertion_code=insertion_code,
                resname=resname,
                record_types_present=record_types,
                atom_count=len(residue_atoms),
                has_backbone_N="N" in atom_names,
                has_backbone_CA="CA" in atom_names,
                has_backbone_C="C" in atom_names,
                has_backbone_O="O" in atom_names,
            )
        )

    return tuple(
        sorted(
            residues,
            key=lambda residue: (
                residue.chain_id,
                residue.residue_number,
                residue.insertion_code,
                residue.resname,
            ),
        )
    )


def parse_pdb(path: Path) -> PDBStructure:
    pdb_path = Path(path)
    atoms: list[AtomRecord] = []
    with pdb_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(("ATOM  ", "HETATM")):
                atoms.append(parse_atom_line(line, line_number))

    if not atoms:
        raise PDBParseError(f"{pdb_path}: no ATOM or HETATM records found")

    return PDBStructure(path=pdb_path, atoms=tuple(atoms), residues=summarize_residues(atoms))

