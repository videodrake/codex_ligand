"""Task 8 PPI-guided pocket selection planning.

This module derives deterministic pocket-selection planning records from Task 7
PPI consensus evidence. It does not run pocket discovery, docking, scoring, or
external scientific tools.
"""

import csv


CONSENSUS_REQUIRED_COLUMNS = [
    "ppi_patch_id",
    "receptor_id",
    "receptor_state",
    "protomer_id",
    "support_methods",
    "support_seeds",
    "support_replicates",
    "n_accepted_poses",
    "egfr_consensus_residues",
    "egfr_contact_centroid_x",
    "egfr_contact_centroid_y",
    "egfr_contact_centroid_z",
    "centroid_spread_angstrom",
    "residue_occupancy_summary",
    "active_face_support_fraction",
    "sheet12_support_fraction",
    "tail_noise_fraction",
    "terminal_artifact_fraction",
    "atp_overlap_fraction",
    "membrane_proximal_fraction",
    "membrane_side_classes",
    "evidence_class",
    "warnings",
]

POCKET_FIELDS = [
    "pocket_plan_id",
    "record_semantics",
    "pocket_detector_runtime_executed",
    "detected_pocket_record",
    "compound_docking_or_scoring_executed",
    "candidate_nomination_executed",
    "ppi_patch_id",
    "receptor_id",
    "receptor_state",
    "protomer_id",
    "source_evidence_class",
    "ppi_relation_class",
    "atp_relation_class",
    "dimer_accessibility_class",
    "membrane_accessibility_class",
    "receptor_state_support",
    "egfr_residues",
    "ppi_centroid_x",
    "ppi_centroid_y",
    "ppi_centroid_z",
    "centroid_spread_angstrom",
    "ppi_confidence_class",
    "future_discovery_allowed",
    "future_compound_docking_planned",
    "blocker",
    "warnings",
    "notes",
]

ACCEPTED_PATCH_CLASSES = set(["CONVERGENT_PATCH", "BROAD_PATCH"])
NON_ACCEPTED_PATCH_CLASSES = set([
    "TAIL_DOMINATED",
    "ATP_OVERLAPPING",
    "MEMBRANE_PROXIMAL",
    "DISPERSED",
    "INSUFFICIENT_EVIDENCE",
    "SCHEMA_INVALID",
])


def parse_float(value, default=0.0):
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return default


def load_consensus_patch_table(path):
    with open(str(path), "r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    missing = [column for column in CONSENSUS_REQUIRED_COLUMNS if column not in fieldnames]
    return missing, rows, fieldnames


def validate_consensus_schema(fieldnames):
    found = set(fieldnames or [])
    audit_rows = []
    missing = []
    for column in CONSENSUS_REQUIRED_COLUMNS:
        present = column in found
        if not present:
            missing.append(column)
        audit_rows.append({
            "check": "required_column",
            "field": column,
            "status": "PASS" if present else "FAIL",
            "details": "" if present else "missing required Task 7 consensus column",
        })
    return missing, audit_rows


def is_patch_accepted_for_pocket_planning(row):
    evidence_class = row.get("evidence_class", "")
    tail_fraction = parse_float(row.get("tail_noise_fraction"))
    terminal_fraction = parse_float(row.get("terminal_artifact_fraction"))
    atp_fraction = parse_float(row.get("atp_overlap_fraction"))
    membrane_fraction = parse_float(row.get("membrane_proximal_fraction"))
    if evidence_class not in ACCEPTED_PATCH_CLASSES:
        return False
    if tail_fraction >= 0.5 or terminal_fraction > 0.0:
        return False
    if atp_fraction >= 0.5 or membrane_fraction >= 0.5:
        return False
    return True


def classify_ppi_relation(row):
    evidence_class = row.get("evidence_class", "")
    if evidence_class == "CONVERGENT_PATCH":
        return "ppi_adjacent"
    if evidence_class == "BROAD_PATCH":
        return "rim"
    if evidence_class == "DISPERSED":
        return "distal"
    return "allosteric"


def classify_atp(row):
    atp_fraction = parse_float(row.get("atp_overlap_fraction"))
    if atp_fraction >= 0.5:
        return "ATP-overlap"
    if atp_fraction > 0.0:
        return "ATP-adjacent"
    return "non-ATP"


def classify_membrane(row):
    membrane_fraction = parse_float(row.get("membrane_proximal_fraction"))
    if membrane_fraction >= 0.5:
        return "membrane-proximal-risk"
    if membrane_fraction > 0.0:
        return "ambiguous"
    return "membrane-compatible"


def classify_dimer(row):
    protomer = str(row.get("protomer_id", "")).strip()
    if protomer in ("A", "B"):
        return "dimer-accessible"
    if not protomer or protomer.lower() in ("unknown", "ambiguous", "x"):
        return "ambiguous"
    return "ambiguous"


def ppi_confidence(row):
    evidence_class = row.get("evidence_class", "")
    spread = parse_float(row.get("centroid_spread_angstrom"))
    n_poses = int(parse_float(row.get("n_accepted_poses")))
    if evidence_class == "CONVERGENT_PATCH" and spread <= 8.0 and n_poses >= 2:
        return "moderate_ppi_evidence"
    if evidence_class == "BROAD_PATCH":
        return "broad_ppi_evidence"
    return "weak_or_blocked_ppi_evidence"


def pocket_warning_text(row, accepted):
    warnings = []
    existing = row.get("warnings", "")
    if existing:
        warnings.append(existing)
    if row.get("evidence_class") in NON_ACCEPTED_PATCH_CLASSES:
        warnings.append("ppi_patch_not_accepted_for_pocket_planning")
    if parse_float(row.get("tail_noise_fraction")) > 0.0 or parse_float(row.get("terminal_artifact_fraction")) > 0.0:
        warnings.append("tail_or_terminal_artifact_visible")
    if parse_float(row.get("atp_overlap_fraction")) > 0.0:
        warnings.append("atp_overlap_visible")
    if parse_float(row.get("membrane_proximal_fraction")) > 0.0:
        warnings.append("membrane_proximal_visible")
    if not accepted:
        warnings.append("no_future_compound_docking_from_this_patch")
    return ";".join([item for item in warnings if item])


def build_pocket_records(consensus_rows):
    pocket_rows = []
    selection_audit = []
    atp_audit = []
    membrane_audit = []
    dimer_audit = []
    accepted_count = 0
    for index, row in enumerate(sorted(consensus_rows, key=lambda item: (item.get("receptor_id", ""), item.get("receptor_state", ""), item.get("protomer_id", ""), item.get("ppi_patch_id", ""))), start=1):
        accepted = is_patch_accepted_for_pocket_planning(row)
        if accepted:
            accepted_count += 1
        atp_class = classify_atp(row)
        membrane_class = classify_membrane(row)
        dimer_class = classify_dimer(row)
        blocker = ""
        if not accepted:
            blocker = "ppi_patch_not_accepted"
        if atp_class == "ATP-overlap":
            blocker = "atp_overlap"
        if membrane_class == "membrane-proximal-risk":
            blocker = "membrane_proximal_risk"
        if dimer_class != "dimer-accessible":
            blocker = "dimer_accessibility_ambiguous"
        pocket_id = "pocket_plan_{0:03d}".format(index)
        pocket_row = {
            "pocket_plan_id": pocket_id,
            "record_semantics": "ppi_guided_pocket_plan",
            "pocket_detector_runtime_executed": False,
            "detected_pocket_record": False,
            "compound_docking_or_scoring_executed": False,
            "candidate_nomination_executed": False,
            "ppi_patch_id": row.get("ppi_patch_id", ""),
            "receptor_id": row.get("receptor_id", ""),
            "receptor_state": row.get("receptor_state", ""),
            "protomer_id": row.get("protomer_id", ""),
            "source_evidence_class": row.get("evidence_class", ""),
            "ppi_relation_class": classify_ppi_relation(row),
            "atp_relation_class": atp_class,
            "dimer_accessibility_class": dimer_class,
            "membrane_accessibility_class": membrane_class,
            "receptor_state_support": row.get("receptor_state", ""),
            "egfr_residues": row.get("egfr_consensus_residues", ""),
            "ppi_centroid_x": row.get("egfr_contact_centroid_x", ""),
            "ppi_centroid_y": row.get("egfr_contact_centroid_y", ""),
            "ppi_centroid_z": row.get("egfr_contact_centroid_z", ""),
            "centroid_spread_angstrom": row.get("centroid_spread_angstrom", ""),
            "ppi_confidence_class": ppi_confidence(row),
            "future_discovery_allowed": accepted and blocker == "",
            "future_compound_docking_planned": False,
            "blocker": blocker,
            "warnings": pocket_warning_text(row, accepted),
            "notes": "planning evidence only; no pocket detector or docking executed",
        }
        pocket_rows.append(pocket_row)
        selection_audit.append({
            "pocket_plan_id": pocket_id,
            "ppi_patch_id": row.get("ppi_patch_id", ""),
            "source_evidence_class": row.get("evidence_class", ""),
            "accepted_for_pocket_planning": accepted,
            "future_discovery_allowed": pocket_row["future_discovery_allowed"],
            "future_compound_docking_planned": False,
            "status": "PASS" if pocket_row["future_discovery_allowed"] else "WARN",
            "notes": pocket_row["warnings"],
        })
        atp_audit.append({
            "pocket_plan_id": pocket_id,
            "ppi_patch_id": row.get("ppi_patch_id", ""),
            "atp_overlap_fraction": row.get("atp_overlap_fraction", ""),
            "atp_relation_class": atp_class,
            "status": "WARN" if atp_class != "non-ATP" else "PASS",
            "notes": "ATP-overlap is blocker/strong warning for PPI-disruptive objective" if atp_class != "non-ATP" else "",
        })
        membrane_audit.append({
            "pocket_plan_id": pocket_id,
            "ppi_patch_id": row.get("ppi_patch_id", ""),
            "membrane_proximal_fraction": row.get("membrane_proximal_fraction", ""),
            "membrane_accessibility_class": membrane_class,
            "status": "WARN" if membrane_class != "membrane-compatible" else "PASS",
            "notes": "membrane-proximal evidence must remain visible" if membrane_class != "membrane-compatible" else "",
        })
        dimer_audit.append({
            "pocket_plan_id": pocket_id,
            "ppi_patch_id": row.get("ppi_patch_id", ""),
            "protomer_id": row.get("protomer_id", ""),
            "dimer_accessibility_class": dimer_class,
            "status": "PASS" if dimer_class == "dimer-accessible" else "WARN",
            "notes": "protomer identity preserved from Task 7 consensus patch",
        })
    return pocket_rows, selection_audit, atp_audit, membrane_audit, dimer_audit, accepted_count
