"""M2.6 raw pocket normalization.

This module converts parsed fpocket evidence into a standardized raw candidate
table. It does not apply ATP, membrane, dimer, or PPI gates.
"""

from __future__ import annotations

from typing import Any


RAW_CANDIDATE_FIELDS = [
    "candidate_pocket_id",
    "state_id",
    "state_role",
    "receptor_id",
    "receptor_path",
    "source_tool",
    "source_pocket_id",
    "egfr_chain_ids",
    "egfr_protomer_ids",
    "multi_protomer",
    "pocket_center_x",
    "pocket_center_y",
    "pocket_center_z",
    "pocket_extent_x",
    "pocket_extent_y",
    "pocket_extent_z",
    "pocket_volume",
    "fpocket_score",
    "fpocket_druggability_score",
    "pocket_residue_count",
    "receptor_residues",
    "mapped_uniprot_residues",
    "mapping_status",
    "raw_candidate_status",
    "evidence_status",
    "warnings",
]


def _split_values(value: Any) -> list[str]:
    return [item for item in str(value or "").split(";") if item]


def _is_multi_protomer(value: Any) -> bool:
    protomers = {item for item in _split_values(value) if item}
    return len(protomers) > 1


def normalize_raw_pockets(raw_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Normalize per-state raw fpocket rows into M2.6 raw candidate rows."""

    rows: list[dict[str, Any]] = []
    for idx, raw in enumerate(
        sorted(
            raw_rows,
            key=lambda row: (
                str(row.get("state_id", "")),
                str(row.get("source_pocket_id", "")),
            ),
        ),
        start=1,
    ):
        receptor_residues = raw.get("raw_residue_ids", "")
        mapped_residues = raw.get("mapped_uniprot_residues", "")
        residue_count = len(set(_split_values(mapped_residues) or _split_values(receptor_residues)))
        warnings = [item for item in _split_values(raw.get("warnings", ""))]
        if raw.get("state_role") == "reference":
            warnings.append("reference_control_raw_evidence_not_primary_state")
        rows.append(
            {
                "candidate_pocket_id": "m2_6_raw_pocket_{0:04d}".format(idx),
                "state_id": raw.get("state_id", ""),
                "state_role": raw.get("state_role", ""),
                "receptor_id": raw.get("receptor_id", ""),
                "receptor_path": raw.get("receptor_path", ""),
                "source_tool": raw.get("source_tool", ""),
                "source_pocket_id": raw.get("source_pocket_id", ""),
                "egfr_chain_ids": raw.get("egfr_chain_ids", ""),
                "egfr_protomer_ids": raw.get("egfr_protomer_ids", ""),
                "multi_protomer": _is_multi_protomer(raw.get("egfr_protomer_ids", "")),
                "pocket_center_x": raw.get("centroid_x", ""),
                "pocket_center_y": raw.get("centroid_y", ""),
                "pocket_center_z": raw.get("centroid_z", ""),
                "pocket_extent_x": raw.get("_extent_x", ""),
                "pocket_extent_y": raw.get("_extent_y", ""),
                "pocket_extent_z": raw.get("_extent_z", ""),
                "pocket_volume": raw.get("fpocket_volume", ""),
                "fpocket_score": raw.get("fpocket_score", ""),
                "fpocket_druggability_score": raw.get("fpocket_druggability_score", ""),
                "pocket_residue_count": residue_count,
                "receptor_residues": receptor_residues,
                "mapped_uniprot_residues": mapped_residues,
                "mapping_status": raw.get("mapping_status", ""),
                "raw_candidate_status": "raw_ungated",
                "evidence_status": raw.get("evidence_status", ""),
                "warnings": ";".join(sorted(set(warnings))),
            }
        )
    return rows
