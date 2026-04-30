"""Small Vina PDBQT/log parsers for M3-T6 collection.

The parser returns pose metadata and centers only. It never writes atom
coordinate records to public outputs and does not infer chemistry or mechanism.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class AtomRecord:
    x: float
    y: float
    z: float
    element: str


@dataclass(frozen=True)
class PoseRecord:
    pose_rank: int | None
    pose_model_index: int
    vina_affinity_kcal_mol: float | None
    pose_center_x: float | None
    pose_center_y: float | None
    pose_center_z: float | None
    pose_atom_count: int
    pose_heavy_atom_count: int
    parse_status: str
    parse_notes: str


VINA_RESULT_RE = re.compile(r"REMARK\s+VINA\s+RESULT:\s*([-+]?\d+(?:\.\d+)?)", re.IGNORECASE)
MODEL_RE = re.compile(r"^\s*MODEL\s+(\d+)", re.IGNORECASE)
LOG_ROW_RE = re.compile(r"^\s*(\d+)\s+([-+]?\d+(?:\.\d+)?)\s+[-+]?\d", re.MULTILINE)


def is_coordinate_record(line: str) -> bool:
    return line.startswith("ATOM") or line.startswith("HETATM")


def _parse_atom(line: str) -> AtomRecord | None:
    if not is_coordinate_record(line):
        return None
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        parts = line.split()
        numeric: list[float] = []
        for part in parts:
            try:
                numeric.append(float(part))
            except ValueError:
                continue
        if len(numeric) < 3:
            return None
        x, y, z = numeric[-6:-3] if len(numeric) >= 6 else numeric[-3:]
    element = line[76:78].strip() if len(line) >= 78 else ""
    if not element:
        parts = line.split()
        element = parts[-1] if parts else ""
    return AtomRecord(x=x, y=y, z=z, element=element.upper())


def compute_pose_center(atom_records: list[AtomRecord]) -> tuple[float, float, float]:
    if not atom_records:
        raise ValueError("cannot compute center without atom records")
    count = float(len(atom_records))
    return (
        sum(atom.x for atom in atom_records) / count,
        sum(atom.y for atom in atom_records) / count,
        sum(atom.z for atom in atom_records) / count,
    )


def _pose_from_lines(lines: list[str], model_index: int, fallback_rank: int | None) -> PoseRecord:
    rank = fallback_rank
    atoms: list[AtomRecord] = []
    affinity: float | None = None
    notes: list[str] = []
    for line in lines:
        model_match = MODEL_RE.match(line)
        if model_match:
            try:
                rank = int(model_match.group(1))
            except ValueError:
                notes.append("model_rank_parse_failed")
        result_match = VINA_RESULT_RE.search(line)
        if result_match:
            affinity = float(result_match.group(1))
        atom = _parse_atom(line)
        if atom is not None:
            atoms.append(atom)
    if not atoms:
        return PoseRecord(rank, model_index, affinity, None, None, None, 0, 0, "FAIL", "no_atom_records")
    center = compute_pose_center(atoms)
    heavy = sum(1 for atom in atoms if not atom.element.startswith("H"))
    parse_status = "PASS"
    if affinity is None:
        parse_status = "WARN"
        notes.append("pdbqt_affinity_missing")
    if rank is None:
        rank = model_index
        parse_status = "WARN"
        notes.append("pose_rank_inferred_from_order")
    return PoseRecord(
        pose_rank=rank,
        pose_model_index=model_index,
        vina_affinity_kcal_mol=affinity,
        pose_center_x=center[0],
        pose_center_y=center[1],
        pose_center_z=center[2],
        pose_atom_count=len(atoms),
        pose_heavy_atom_count=heavy,
        parse_status=parse_status,
        parse_notes=";".join(notes) if notes else "parsed",
    )


def parse_vina_pdbqt(path: Path) -> list[PoseRecord]:
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()
    models: list[list[str]] = []
    current: list[str] | None = None
    for line in lines:
        if MODEL_RE.match(line):
            if current:
                models.append(current)
            current = [line]
            continue
        if current is not None:
            current.append(line)
            if line.strip().upper() == "ENDMDL":
                models.append(current)
                current = None
    if current:
        models.append(current)
    if not models:
        models = [lines]
    poses = [_pose_from_lines(model, index + 1, None if len(models) > 1 else 1) for index, model in enumerate(models)]
    return [pose for pose in poses if pose.pose_atom_count > 0]


def _split_models(lines: list[str]) -> list[list[str]]:
    models: list[list[str]] = []
    current: list[str] | None = None
    for line in lines:
        if MODEL_RE.match(line):
            if current:
                models.append(current)
            current = [line]
            continue
        if current is not None:
            current.append(line)
            if line.strip().upper() == "ENDMDL":
                models.append(current)
                current = None
    if current:
        models.append(current)
    if not models:
        models = [lines]
    return models


def extract_pose_atoms(path: Path, pose_model_index: int | None = None, pose_rank: int | None = None) -> list[AtomRecord]:
    """Return atom records for one Vina pose model without logging atom lines."""
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    models = _split_models(text.splitlines())
    selected: list[str] | None = None
    if pose_model_index is not None and 1 <= int(pose_model_index) <= len(models):
        selected = models[int(pose_model_index) - 1]
    if selected is None and pose_rank is not None:
        for model in models:
            rank: int | None = None
            for line in model:
                match = MODEL_RE.match(line)
                if match:
                    rank = int(match.group(1))
                    break
            if rank == int(pose_rank):
                selected = model
                break
    if selected is None and len(models) == 1:
        selected = models[0]
    if selected is None:
        return []
    atoms: list[AtomRecord] = []
    for line in selected:
        atom = _parse_atom(line)
        if atom is not None:
            atoms.append(atom)
    return atoms


def parse_vina_log(path: Path, max_bytes: int = 200000) -> dict[int, float]:
    data = Path(path).read_bytes()[:max_bytes]
    text = data.decode("utf-8", errors="ignore")
    affinities: dict[int, float] = {}
    for match in LOG_ROW_RE.finditer(text):
        rank = int(match.group(1))
        affinity = float(match.group(2))
        affinities.setdefault(rank, affinity)
    return affinities
