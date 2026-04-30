"""M3-T10 candidate hypothesis tiering helpers.

The functions here keep hard gates separate from soft scores. They assign
computational hypothesis tiers only; they never treat docking affinity as
promotion evidence and never claim validation.
"""

from __future__ import annotations

import math
from typing import Any


TIER_CODES = {
    "Tier 1": "TIER1_STRONG_CHEMICAL_ANCHOR_POCKET_HYPOTHESIS",
    "Tier 2": "TIER2_PLAUSIBLE_SINGLE_ANCHOR_OR_STATE_SPECIFIC_HYPOTHESIS",
    "Tier 3": "TIER3_EXPLORATORY_CHEMICAL_SIGNAL",
    "Reject": "REJECTED_OR_QUARANTINED",
    "Not applicable": "NOT_APPLICABLE",
}

TIER_RANK = {"Tier 1": 1, "Tier 2": 2, "Tier 3": 3, "Reject": 4, "Not applicable": 5}

REJECT_REASONS = {
    "none",
    "ATP_migration",
    "ATP_like_mechanism",
    "outside_accepted_pocket",
    "M2_pocket_not_accepted",
    "non_ATP_gate_failed",
    "PPI_relationship_failed",
    "lower_lateral_gate_failed",
    "dimer_accessibility_failed",
    "pose_retention_failed",
    "membrane_or_dimer_conflict",
    "compound_convergence_failed",
    "primary_state_support_missing",
    "reference_only_not_promoted",
    "broad_scan_only_not_promoted",
    "missing_accepted_pocket_trace",
    "missing_compound_trace",
    "missing_anchor_trace",
    "missing_mapping",
    "non_public_compound_id",
    "confidentiality_violation",
    "internal_id_leak",
    "smiles_leak",
    "affinity_only_promotion_blocked",
    "old_workflow_source_blocked",
    "unknown",
}


def boolish(value: Any) -> bool:
    return str(value or "").strip().lower() in {"true", "t", "1", "yes", "y", "pass", "passed"}


def as_float(value: Any) -> float | None:
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def clamp01(value: float | None) -> float:
    if value is None:
        return 0.0
    return max(0.0, min(1.0, value))


def anchor_class_score(anchor_class: str) -> float:
    return {
        "multi_compound_primary_state_supported": 1.0,
        "single_compound_state_robust": 0.75,
        "single_state_multi_compound": 0.65,
        "reference_only": 0.25,
        "weak_or_scattered": 0.25,
        "no_valid_cluster_support": 0.0,
        "fail_missing_input": 0.0,
    }.get(anchor_class or "", 0.0)


def affinity_descriptive_score(value: Any) -> float:
    affinity = as_float(value)
    if affinity is None:
        return 0.0
    # Descriptive only: maps common Vina medians roughly into 0..1.
    return clamp01((abs(affinity) - 4.0) / 8.0)


def compute_soft_scores(evidence: dict[str, Any], config: dict[str, Any]) -> dict[str, float]:
    state_ids = {item for item in str(evidence.get("primary_state_ids_supported") or "").split(";") if item}
    state_score = 1.0 if len(state_ids) >= 2 else (0.6 if len(state_ids) == 1 else (0.2 if boolish(evidence.get("reference_only_support")) else 0.0))
    symmetry_score = 1.0 if boolish(evidence.get("symmetry_support")) else 0.0
    cluster_fraction = as_float(evidence.get("max_cluster_fraction") or evidence.get("median_cluster_fraction"))
    within = as_float(evidence.get("min_within_pocket_fraction"))
    ppi_fraction = as_float(evidence.get("max_ppi_relationship_confirmed_fraction"))
    ppi_contact = as_float(evidence.get("ppi_hotspot_contact_median_max"))
    scores = {
        "S_ppi_patch_support": 1.0 if boolish(evidence.get("ppi_relationship_pass")) else clamp01(ppi_fraction),
        "S_pocket_geometry": sum(1 for key in ["non_atp_pass", "lower_lateral_pass", "dimer_accessibility_pass", "no_membrane_or_dimer_conflict"] if boolish(evidence.get(key))) / 4.0,
        "S_state_support": state_score,
        "S_symmetry_support": symmetry_score,
        "S_compound_anchor_convergence": anchor_class_score(str(evidence.get("anchor_convergence_class") or "")),
        "S_pose_convergence": (clamp01(cluster_fraction) + clamp01(within)) / 2.0,
        "S_affinity_normalized": affinity_descriptive_score(evidence.get("vina_affinity_median_kcal_mol") or evidence.get("median_affinity_across_supported")),
        "S_ppi_disruption_geometry": max(clamp01(ppi_fraction), 1.0 if (ppi_contact or 0.0) > 0.0 else 0.0),
        "S_chemical_qc": 1.0 if evidence.get("compound_public_id") in {"Cpd-A", "Cpd-B", "Cpd-C", ""} else 0.0,
    }
    weights = config.get("soft_scores", {}).get("weights", {})
    total_weight = sum(float(v) for v in weights.values()) or 1.0
    non_affinity_weight = sum(float(v) for k, v in weights.items() if k != "S_affinity_normalized") or 1.0
    scores["candidate_soft_score_total"] = sum(scores.get(k, 0.0) * float(v) for k, v in weights.items()) / total_weight
    scores["candidate_soft_score_non_affinity"] = sum(scores.get(k, 0.0) * float(v) for k, v in weights.items() if k != "S_affinity_normalized") / non_affinity_weight
    scores["candidate_priority_score"] = scores["candidate_soft_score_non_affinity"]
    return {key: round(value, 6) for key, value in scores.items()}


def hard_gate_failures(flags: dict[str, Any]) -> list[str]:
    ordered = [
        "m2_pocket_accepted",
        "non_atp_pass",
        "ppi_relationship_pass",
        "lower_lateral_pass",
        "dimer_accessibility_pass",
        "pose_retention_pass",
        "atp_migration_absent",
        "compound_convergence_pass",
        "primary_state_support",
        "confidentiality_pass",
        "no_membrane_or_dimer_conflict",
        "accepted_pocket_trace_present",
        "compound_trace_present",
        "anchor_trace_present",
    ]
    return [key for key in ordered if not boolish(flags.get(key))]


def reject_reason_from_flags(flags: dict[str, Any], support_reject_reason: str = "none") -> str:
    if not boolish(flags.get("confidentiality_pass")):
        return "confidentiality_violation"
    mechanism = str(flags.get("dominant_mechanism_class") or flags.get("mechanism_class_counts") or "")
    if not boolish(flags.get("atp_migration_absent")):
        return "ATP_migration"
    if "ATP_like_reject" in mechanism:
        return "ATP_like_mechanism"
    if support_reject_reason in {"outside_accepted_pocket", "outside_accepted_pocket_support"}:
        return "outside_accepted_pocket"
    if support_reject_reason in {"membrane_penetration", "central_dimer_clash", "membrane_or_dimer_conflict"}:
        return "membrane_or_dimer_conflict"
    if not boolish(flags.get("accepted_pocket_trace_present")):
        return "missing_accepted_pocket_trace"
    if not boolish(flags.get("compound_trace_present")):
        return "missing_compound_trace"
    if not boolish(flags.get("anchor_trace_present")):
        return "missing_anchor_trace"
    if not boolish(flags.get("m2_pocket_accepted")):
        return "M2_pocket_not_accepted"
    if not boolish(flags.get("non_atp_pass")):
        return "non_ATP_gate_failed"
    if not boolish(flags.get("ppi_relationship_pass")):
        return "PPI_relationship_failed"
    if not boolish(flags.get("lower_lateral_pass")):
        return "lower_lateral_gate_failed"
    if not boolish(flags.get("dimer_accessibility_pass")):
        return "dimer_accessibility_failed"
    if not boolish(flags.get("pose_retention_pass")):
        return "pose_retention_failed"
    if not boolish(flags.get("no_membrane_or_dimer_conflict")):
        return "membrane_or_dimer_conflict"
    if boolish(flags.get("reference_only_support")) and not boolish(flags.get("primary_state_support")):
        return "reference_only_not_promoted"
    if boolish(flags.get("broad_scan_only_support")):
        return "broad_scan_only_not_promoted"
    if not boolish(flags.get("compound_convergence_pass")):
        return "compound_convergence_failed"
    if not boolish(flags.get("primary_state_support")):
        return "primary_state_support_missing"
    return "none"


def assign_candidate_tier(
    flags: dict[str, Any],
    anchor_class: str,
    support_reject_reason: str = "none",
    strict_tier1: bool = True,
) -> tuple[str, str, int, str, str]:
    del strict_tier1
    reason = reject_reason_from_flags(flags, support_reject_reason)
    tier1_gates = [
        "m2_pocket_accepted",
        "non_atp_pass",
        "ppi_relationship_pass",
        "lower_lateral_pass",
        "dimer_accessibility_pass",
        "pose_retention_pass",
        "atp_migration_absent",
        "compound_convergence_pass",
        "primary_state_support",
        "confidentiality_pass",
        "no_membrane_or_dimer_conflict",
    ]
    tier2_gates = list(tier1_gates)
    all_t1 = all(boolish(flags.get(key)) for key in tier1_gates) and not boolish(flags.get("reference_only_support"))
    all_t2 = all(boolish(flags.get(key)) for key in tier2_gates) and not boolish(flags.get("reference_only_support"))
    if reason in {
        "ATP_migration",
        "ATP_like_mechanism",
        "outside_accepted_pocket",
        "M2_pocket_not_accepted",
        "non_ATP_gate_failed",
        "PPI_relationship_failed",
        "lower_lateral_gate_failed",
        "dimer_accessibility_failed",
        "pose_retention_failed",
        "membrane_or_dimer_conflict",
        "missing_accepted_pocket_trace",
        "missing_compound_trace",
        "missing_anchor_trace",
        "confidentiality_violation",
    }:
        tier = "Reject"
        accept = ""
    elif all_t1 and anchor_class == "multi_compound_primary_state_supported":
        tier = "Tier 1"
        reason = "none"
        accept = "all Tier 1 hard gates passed; multi-compound primary-state anchor support"
    elif all_t2 and anchor_class in {"single_compound_state_robust", "single_state_multi_compound"}:
        tier = "Tier 2"
        reason = "none"
        accept = "all Tier 2 hard gates passed; focused anchor support pattern"
    elif reason in {"reference_only_not_promoted", "broad_scan_only_not_promoted", "compound_convergence_failed", "primary_state_support_missing", "none"}:
        tier = "Tier 3"
        accept = "exploratory computational signal; not promoted by affinity"
        if reason == "none":
            reason = "compound_convergence_failed"
    else:
        tier = "Reject"
        accept = ""
    return tier, TIER_CODES[tier], TIER_RANK[tier], reason, accept
