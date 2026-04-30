"""M2.7 dimer accessibility gate."""

from __future__ import annotations

from typing import Any

from egfr_myo1d.pocket.ppi_relationship import format_float, residue_set, split_values


DIMER_ACCESSIBILITY_FIELDS = [
    "pocket_family_id",
    "representative_candidate_pocket_id",
    "state_ids",
    "state_roles",
    "protomer_ids",
    "interface_overlap_count",
    "interface_residue_count",
    "interface_overlap_fraction",
    "sasa_ratio_dimer_vs_isolated_protomer",
    "sasa_status",
    "dimer_accessibility_pass",
    "dimer_accessibility_reason",
    "warnings",
]


def evaluate_dimer_accessibility(
    family: dict[str, Any],
    interface_residues: set[str],
    threshold: float = 0.25,
) -> dict[str, Any]:
    warnings: list[str] = ["sasa_not_computed_m2_7_optional"]
    pocket_residues = residue_set(family.get("residue_union_uniprot", ""))
    overlap = pocket_residues & interface_residues
    fraction = len(overlap) / float(len(pocket_residues)) if pocket_residues else 0.0
    if not pocket_residues:
        warnings.append("pocket_residues_unavailable")
    if not interface_residues:
        warnings.append("interface_residue_index_unavailable")
    if not pocket_residues:
        passed = False
        reason = "pocket_residues_unavailable"
    elif not interface_residues:
        passed = False
        reason = "interface_residue_index_unavailable"
    else:
        passed = fraction < threshold
        reason = "interface_overlap_below_threshold" if passed else "interface_overlap_fraction_exceeds_threshold"
    return {
        "pocket_family_id": family.get("pocket_family_id", ""),
        "representative_candidate_pocket_id": family.get("representative_candidate_pocket_id", ""),
        "state_ids": family.get("member_state_ids", ""),
        "state_roles": family.get("member_state_roles", ""),
        "protomer_ids": family.get("member_protomer_ids", ""),
        "interface_overlap_count": len(overlap),
        "interface_residue_count": len(interface_residues),
        "interface_overlap_fraction": format_float(fraction),
        "sasa_ratio_dimer_vs_isolated_protomer": "",
        "sasa_status": "not_computed_m2_7_optional",
        "dimer_accessibility_pass": passed,
        "dimer_accessibility_reason": reason,
        "warnings": ";".join(sorted(set(warnings))),
    }
