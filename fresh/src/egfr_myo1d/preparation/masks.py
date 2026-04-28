"""EGFR future-sampling mask contracts for Task 4."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from egfr_myo1d.structure.geometry import parse_numeric_vector, validate_membrane_frame


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_membrane_frame(path: Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def build_egfr_masks(contract: dict[str, Any], membrane_frame_path: Path) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    receptor = contract["receptor"]
    membrane_data = load_membrane_frame(membrane_frame_path)
    frame_report = validate_membrane_frame(membrane_data)
    atp_site = receptor.get("non_target_masks", {}).get("atp_site_residues", [])
    membrane_proximal: list[int] = []
    dimer_core: list[int] = []
    audit_rows = [
        {
            "receptor_id": receptor["id"],
            "mask_name": "membrane_proximal_residues",
            "source": str(membrane_frame_path),
            "policy": "exclude_or_warn_for_future_ppi_sampling",
            "threshold": "synthetic placeholder: no residues classified without production geometry",
            "residue_count": len(membrane_proximal),
            "status": frame_report["status"],
            "notes": "Task 4 validates frame metadata but does not infer production membrane exclusion",
        },
        {
            "receptor_id": receptor["id"],
            "mask_name": "atp_site_seed_residues",
            "source": "contract",
            "policy": "flag_non_target_overlap",
            "threshold": "contract seed list",
            "residue_count": len(atp_site),
            "status": "PASS",
            "notes": "future docking poses overlapping these residues are non-target for PPI objective",
        },
    ]
    membrane_mask = {
        "receptor_id": receptor["id"],
        "membrane_frame_source": Path(membrane_frame_path).name,
        "membrane_proximal_policy": "exclude_or_warn_for_future_ppi_sampling",
        "atp_site_policy": "flag_non_target_overlap",
        "dimer_interface_policy": "preserve_context_and_report_accessibility",
        "residue_masks": {
            "atp_site_seed_residues": atp_site,
            "membrane_proximal_residues": membrane_proximal,
            "dimer_core_residues": dimer_core,
        },
        "membrane_frame_qc": frame_report,
    }
    non_target_mask = {
        "receptor_id": receptor["id"],
        "atp_site_policy": "flag_non_target_overlap",
        "active_site_is_not_target": True,
        "residue_masks": {
            "atp_site_seed_residues": atp_site,
            "membrane_proximal_residues": membrane_proximal,
            "dimer_core_residues": dimer_core,
        },
    }
    return membrane_mask, non_target_mask, audit_rows

