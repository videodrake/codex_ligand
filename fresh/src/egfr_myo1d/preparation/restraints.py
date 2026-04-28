"""Future docking restraint/mask contract writers for Task 4."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_active_face_contract(contract_myo1d: dict[str, Any]) -> dict[str, Any]:
    return {
        "construct_id": contract_myo1d["primary_construct"]["id"],
        "primary_active_face_residues": contract_myo1d["active_face"],
        "sheet12_support_residues": contract_myo1d["sheet12_support"],
        "structural_buffer_residues": contract_myo1d["structural_buffer"],
        "contact_monitoring_cap_residues": [number for number in contract_myo1d["contact_monitoring_cap"] if number <= 1001],
        "back_face_policy": "detect_and_report_not_score_bonus",
        "active_face_use": "qc_annotation_and_future_sampling_contract_only",
        "score_bonus_allowed": False,
    }


def write_active_face_contract(path: Path, contract_myo1d: dict[str, Any]) -> dict[str, Any]:
    data = build_active_face_contract(contract_myo1d)
    write_json(path, data)
    return data

