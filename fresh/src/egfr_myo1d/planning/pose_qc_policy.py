"""Future pose-acceptance policy for EGFR-MYO1D PPI sampling.

Task 6 only serializes deterministic policy. It does not classify real docking
poses and does not run a docking engine.
"""

from __future__ import annotations

from typing import Any


POSE_CLASSES = [
    "accepted_active_8_9_supported_12",
    "accepted_active_8_9_only",
    "review_sheet12_dominant",
    "reject_tail_dominant_artifact",
    "reject_flipped_or_back_face",
    "reject_membrane_proximal",
    "reject_atp_site_overlap",
    "reject_chain_or_residue_reset",
    "reject_quarantined_input",
]


def build_pose_acceptance_policy(contract: dict[str, Any]) -> dict[str, Any]:
    myo = contract.get("myo1d", {})
    receptor = contract.get("receptor", {})
    active_face = list(myo.get("active_face", []))
    sheet12 = list(myo.get("sheet12_support", []))
    monitoring_cap = list(myo.get("contact_monitoring_cap", []))
    atp_site = list(receptor.get("non_target_masks", {}).get("atp_site_residues", []))
    return {
        "policy_version": "task6_spec_only_v0_1",
        "current_task_executes_pose_classification": False,
        "future_pose_classes": POSE_CLASSES,
        "residue_numbering_policy": "validated_original_numbering_required",
        "chain_identity_policy": "validated_receptor_and_partner_chain_ids_required",
        "active_face_residues": active_face,
        "sheet12_support_residues": sheet12,
        "short_c_terminal_cap_residues": [number for number in monitoring_cap if number <= 1001],
        "extended_tail_noise_zone_residues": [number for number in monitoring_cap if number >= 998],
        "score_bonus_allowed": False,
        "key_residue_bonus_weight": 0.0,
        "egfr_masks": {
            "membrane_proximal_policy": "reject_or_warn_for_future_ppi_sampling",
            "atp_site_policy": "reject_or_flag_non_target_overlap",
            "atp_site_seed_residues": atp_site,
            "dimer_context_policy": "preserve_protomer_specific_accessibility",
        },
        "deterministic_thresholds": {
            "minimum_active_face_contact_count": 1,
            "minimum_sheet12_support_contact_count": 1,
            "tail_dominant_fraction_warn": 0.50,
            "tail_dominant_fraction_reject": 0.70,
            "orientation_back_face_reject": True,
            "membrane_proximal_overlap_reject": True,
            "atp_site_overlap_reject": True,
        },
        "class_rules": {
            "accepted_active_8_9_supported_12": "active-face contact evidence and sheet-12 support present, orientation passes, no non-target overlap",
            "accepted_active_8_9_only": "active-face evidence present, orientation passes, sheet-12 absent or weak, no non-target overlap",
            "review_sheet12_dominant": "sheet-12 dominates but active-face evidence is weak or missing",
            "reject_tail_dominant_artifact": "residues 998+ or terminal residues dominate contacts without active-face engagement",
            "reject_flipped_or_back_face": "contact exists but active-face orientation fails",
            "reject_membrane_proximal": "receptor-side contact lies on membrane-facing inaccessible surface",
            "reject_atp_site_overlap": "receptor-side contact overlaps ATP-site mask or ATP-cleft seed residues",
            "reject_chain_or_residue_reset": "pose residue numbering cannot map back to validated receptor/partner numbering",
            "reject_quarantined_input": "pose derives from quarantined receptor or partner input",
        },
        "non_scoring_notice": "Task 6 records future QC thresholds only; no ranking, scoring, or pose classification is performed.",
    }


def audit_pose_policy(policy: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for pose_class in policy["future_pose_classes"]:
        decision = "accept" if pose_class.startswith("accepted_") else ("review" if pose_class.startswith("review_") else "reject")
        rows.append(
            {
                "pose_class": pose_class,
                "decision": decision,
                "threshold_summary": policy["class_rules"].get(pose_class, ""),
                "status": "PASS",
                "notes": "future policy only; no pose classification executed",
            }
        )
    return rows
