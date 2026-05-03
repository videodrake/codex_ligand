"""Scalar geometry helpers for M3-T7 pose attribution.

This module reads coordinate-bearing pose atoms in memory and returns only
distances, fractions, and flags. It does not write atom-level coordinates or
infer compound efficacy.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from egfr_myo1d.compound.pdbqt_parse import AtomRecord, compute_pose_center


DEFAULT_CONFIG: dict[str, Any] = {
    "schema_version": "m3_pose_attribution_config_v1",
    "pocket": {
        "centroid_margin_A": 2.0,
        "atom_fraction_min": 0.50,
        "heavy_atom_fraction_min": 0.50,
        "nearest_pocket_accept_distance_A": 4.0,
    },
    "atp": {
        "centroid_reject_distance_A": 8.0,
        "atom_overlap_reject_fraction": 0.20,
        "heavy_atom_overlap_reject_fraction": 0.20,
        "atom_near_atp_cutoff_A": 4.5,
    },
    "ppi": {
        "contact_cutoff_A": 4.5,
        "near_cutoff_A": 8.0,
        "rim_cutoff_A": 10.0,
        "hotspot_contact_min": 1,
        "relationship_requires_residue_contact_or_rim": True,
    },
    "membrane": {"penetration_margin_A": 1.5, "use_membrane_core_bounds": True},
    "dimer": {"receptor_steric_clash_cutoff_A": 2.0, "central_interface_clash_cutoff_A": 3.0},
    "mechanism": {
        "direct_ppi_contact_cutoff_A": 4.5,
        "rim_blocker_cutoff_A": 8.0,
        "allosteric_near_cutoff_A": 12.0,
    },
}


@dataclass(frozen=True)
class BoxProxy:
    pocket_family_id: str
    box_id: str
    state_id: str
    protomer_id: str
    center: tuple[float, float, float]
    size: tuple[float, float, float]


@dataclass(frozen=True)
class ReferencePoint:
    point_id: str
    xyz: tuple[float, float, float]
    residue_public: str = ""
    state_id: str = ""
    protomer_id: str = ""
    radius: float | None = None


@dataclass(frozen=True)
class MembraneFrame:
    origin: tuple[float, float, float]
    normal: tuple[float, float, float]
    core_z_min: float | None = None
    core_z_max: float | None = None


def _float(value: Any) -> float | None:
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    if math.isfinite(result):
        return result
    return None


def _field(row: dict[str, Any], names: list[str]) -> str:
    lower = {key.lower(): key for key in row}
    for name in names:
        key = lower.get(name.lower())
        if key is not None and str(row.get(key, "")).strip():
            return str(row.get(key, "")).strip()
    return ""


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def heavy_atoms(atoms: list[AtomRecord]) -> list[AtomRecord]:
    return [atom for atom in atoms if not (atom.element or "").upper().startswith("H")]


def centroid(atoms: list[AtomRecord], heavy_only: bool = False) -> tuple[float, float, float]:
    selected = heavy_atoms(atoms) if heavy_only else atoms
    if not selected:
        selected = atoms
    return compute_pose_center(selected)


def inside_box(point: tuple[float, float, float], box: BoxProxy, margin: float = 0.0) -> bool:
    return all(abs(point[i] - box.center[i]) <= (box.size[i] / 2.0 + margin) for i in range(3))


def fraction_inside_box(atoms: list[AtomRecord], box: BoxProxy, margin: float = 0.0, heavy_only: bool = False) -> float:
    selected = heavy_atoms(atoms) if heavy_only else atoms
    if not selected:
        return 0.0
    count = sum(1 for atom in selected if inside_box((atom.x, atom.y, atom.z), box, margin=margin))
    return count / float(len(selected))


def box_from_row(row: dict[str, Any]) -> BoxProxy | None:
    cx = _float(_field(row, ["box_center_x", "center_x", "x"]))
    cy = _float(_field(row, ["box_center_y", "center_y", "y"]))
    cz = _float(_field(row, ["box_center_z", "center_z", "z"]))
    sx = _float(_field(row, ["box_size_x", "size_x"]))
    sy = _float(_field(row, ["box_size_y", "size_y"]))
    sz = _float(_field(row, ["box_size_z", "size_z"]))
    if None in {cx, cy, cz, sx, sy, sz} or sx <= 0 or sy <= 0 or sz <= 0:
        return None
    return BoxProxy(
        pocket_family_id=_field(row, ["pocket_family_id", "pocket_id"]),
        box_id=_field(row, ["box_id"]) or _field(row, ["pocket_family_id", "pocket_id"]),
        state_id=_field(row, ["state_id"]),
        protomer_id=_field(row, ["protomer_id", "protomer", "chain_id"]),
        center=(cx, cy, cz),
        size=(sx, sy, sz),
    )


def reference_point_from_row(row: dict[str, Any], prefix: str, index: int) -> ReferencePoint | None:
    x = _float(_field(row, [f"{prefix}_center_x", f"{prefix}_centroid_x", f"{prefix}_x", "center_x", "centroid_x", "egfr_contact_centroid_x", "x", "box_center_x"]))
    y = _float(_field(row, [f"{prefix}_center_y", f"{prefix}_centroid_y", f"{prefix}_y", "center_y", "centroid_y", "egfr_contact_centroid_y", "y", "box_center_y"]))
    z = _float(_field(row, [f"{prefix}_center_z", f"{prefix}_centroid_z", f"{prefix}_z", "center_z", "centroid_z", "egfr_contact_centroid_z", "z", "box_center_z"]))
    if None in {x, y, z}:
        return None
    point_id = _field(row, [f"{prefix}_id", f"{prefix}_patch_id", "hotspot_id", "residue_public", "egfr_residue_number", "residue_id", "site_id", "id"]) or f"{prefix}_{index}"
    residue_public = _field(row, ["residue_public", "residue_label", "hotspot_label", "egfr_residue_number", "residue_id"])
    return ReferencePoint(
        point_id=point_id,
        xyz=(x, y, z),
        residue_public=residue_public or point_id,
        state_id=_field(row, ["state_id"]),
        protomer_id=_field(row, ["protomer_id", "egfr_protomer_id", "protomer", "chain_id"]),
        radius=_float(_field(row, ["radius", f"{prefix}_radius"])),
    )


def membrane_frame_from_json(data: dict[str, Any]) -> MembraneFrame:
    origin_values = data.get("origin") or data.get("membrane_origin") or [0.0, 0.0, 0.0]
    normal_values = data.get("normal") or data.get("membrane_normal") or data.get("z_axis") or [0.0, 0.0, 1.0]
    origin = tuple(float(origin_values[i]) for i in range(3))
    normal_raw = tuple(float(normal_values[i]) for i in range(3))
    norm = math.sqrt(sum(value * value for value in normal_raw)) or 1.0
    normal = tuple(value / norm for value in normal_raw)
    core_min = _float(data.get("core_z_min") or data.get("membrane_core_z_min") or data.get("z_min"))
    core_max = _float(data.get("core_z_max") or data.get("membrane_core_z_max") or data.get("z_max"))
    return MembraneFrame(origin=origin, normal=normal, core_z_min=core_min, core_z_max=core_max)


def project_membrane_z(point: tuple[float, float, float], frame: MembraneFrame) -> float:
    delta = tuple(point[i] - frame.origin[i] for i in range(3))
    return sum(delta[i] * frame.normal[i] for i in range(3))


def nearest_point(point: tuple[float, float, float], refs: list[ReferencePoint]) -> tuple[ReferencePoint | None, float | None]:
    if not refs:
        return None, None
    best = min(refs, key=lambda ref: distance(point, ref.xyz))
    return best, distance(point, best.xyz)


def nearest_atom_distance(atoms: list[AtomRecord], refs: list[ReferencePoint]) -> float | None:
    if not atoms or not refs:
        return None
    return min(distance((atom.x, atom.y, atom.z), ref.xyz) for atom in atoms for ref in refs)


def atom_overlap_fraction(atoms: list[AtomRecord], refs: list[ReferencePoint], cutoff: float, heavy_only: bool = False) -> float:
    selected = heavy_atoms(atoms) if heavy_only else atoms
    if not selected or not refs:
        return 0.0
    near = 0
    for atom in selected:
        if min(distance((atom.x, atom.y, atom.z), ref.xyz) for ref in refs) <= cutoff:
            near += 1
    return near / float(len(selected))


def ppi_contacts(atoms: list[AtomRecord], refs: list[ReferencePoint], cutoff: float) -> tuple[int, list[str], float | None, str]:
    selected = heavy_atoms(atoms) or atoms
    if not selected or not refs:
        return 0, [], None, ""
    contacted: dict[str, str] = {}
    nearest: tuple[float, ReferencePoint] | None = None
    for ref in refs:
        ref_nearest = min(distance((atom.x, atom.y, atom.z), ref.xyz) for atom in selected)
        if nearest is None or ref_nearest < nearest[0]:
            nearest = (ref_nearest, ref)
        if ref_nearest <= cutoff:
            contacted[ref.point_id] = ref.residue_public or ref.point_id
    residues = sorted(set(contacted.values()))
    if nearest is None:
        return 0, residues, None, ""
    return len(contacted), residues, nearest[0], nearest[1].point_id
