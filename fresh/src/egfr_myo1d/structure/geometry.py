"""Basic geometry helpers for Task 3 QC reports."""

from __future__ import annotations

import math
from typing import Iterable

from egfr_myo1d.structure.pdb_parser import AtomRecord


def centroid(atoms: Iterable[AtomRecord]) -> list[float] | None:
    atom_list = list(atoms)
    if not atom_list:
        return None
    return [
        sum(atom.x for atom in atom_list) / len(atom_list),
        sum(atom.y for atom in atom_list) / len(atom_list),
        sum(atom.z for atom in atom_list) / len(atom_list),
    ]


def vector_norm(vector: list[float]) -> float:
    return math.sqrt(sum(value * value for value in vector))


def normalize_vector(vector: list[float]) -> list[float]:
    norm = vector_norm(vector)
    if norm == 0.0:
        raise ValueError("cannot normalize zero vector")
    return [value / norm for value in vector]


def dot(a: list[float], b: list[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def distance(a: list[float], b: list[float]) -> float:
    return vector_norm([x - y for x, y in zip(a, b)])


def parse_numeric_vector(value: object, field_name: str) -> list[float]:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{field_name} must be a numeric list of length 3")
    vector: list[float] = []
    for item in value:
        if not isinstance(item, (int, float)):
            raise ValueError(f"{field_name} must contain only numeric values")
        vector.append(float(item))
    return vector


def validate_membrane_frame(data: dict[str, object]) -> dict[str, object]:
    messages: list[dict[str, object]] = []
    report: dict[str, object] = {
        "receptor_id": data.get("receptor_id"),
        "coordinate_convention": data.get("coordinate_convention", ""),
        "status": "PASS",
        "messages": messages,
    }

    try:
        membrane_normal = parse_numeric_vector(data.get("membrane_normal"), "membrane_normal")
        dimer_axis = parse_numeric_vector(data.get("dimer_axis"), "dimer_axis")
        normal_norm = vector_norm(membrane_normal)
        axis_norm = vector_norm(dimer_axis)
        if normal_norm == 0.0:
            raise ValueError("membrane_normal has zero norm")
        if axis_norm == 0.0:
            raise ValueError("dimer_axis has zero norm")
        normal_unit = normalize_vector(membrane_normal)
        axis_unit = normalize_vector(dimer_axis)
        report["membrane_normal_unit"] = normal_unit
        report["dimer_axis_unit"] = axis_unit
        report["membrane_normal_norm"] = normal_norm
        report["dimer_axis_norm"] = axis_norm
        if abs(dot(normal_unit, axis_unit)) > 0.90:
            messages.append(
                {
                    "severity": "WARN",
                    "code": "membrane_normal_parallel_to_dimer_axis",
                    "message": "membrane normal and dimer axis are nearly parallel",
                }
            )
    except ValueError as exc:
        report["status"] = "FAIL"
        messages.append({"severity": "FAIL", "code": "invalid_membrane_vector", "message": str(exc)})
        return report

    centroid_a = data.get("protomer_a_centroid")
    centroid_b = data.get("protomer_b_centroid")
    try:
        parsed_a = parse_numeric_vector(centroid_a, "protomer_a_centroid")
        parsed_b = parse_numeric_vector(centroid_b, "protomer_b_centroid")
        report["protomer_centroid_distance"] = distance(parsed_a, parsed_b)
        if distance(parsed_a, parsed_b) == 0.0:
            messages.append(
                {
                    "severity": "WARN",
                    "code": "identical_protomer_centroids",
                    "message": "protomer centroids are identical",
                }
            )
    except ValueError as exc:
        messages.append(
            {
                "severity": "WARN",
                "code": "missing_or_invalid_protomer_centroid",
                "message": str(exc),
            }
        )

    if any(message["severity"] == "WARN" for message in messages):
        report["status"] = "WARN"
    return report

