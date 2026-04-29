"""PDB writing helpers that preserve parsed identity fields."""

from __future__ import annotations

from pathlib import Path

from egfr_myo1d.structure.pdb_parser import AtomRecord


def write_pdb_atoms(path: Path, atoms: list[AtomRecord] | tuple[AtomRecord, ...]) -> dict[str, int]:
    """Write atom records without renumbering or changing chain/residue identity."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("REMARK   Prepared by fresh EGFR-MYO1D Task 4 QC layer\n")
        for atom in atoms:
            handle.write(atom.raw_line.rstrip("\n") + "\n")
        handle.write("END\n")
    return {"atom_count": len(atoms)}


def select_atoms_by_residue_range(
    atoms: list[AtomRecord] | tuple[AtomRecord, ...],
    start: int,
    end: int,
    include_resnames: set[str] | None = None,
) -> list[AtomRecord]:
    include = include_resnames or set()
    return [
        atom
        for atom in atoms
        if start <= atom.residue_number <= end or atom.resname in include
    ]

