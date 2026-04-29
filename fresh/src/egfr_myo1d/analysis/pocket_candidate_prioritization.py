"""Task 9 pocket candidate intake and prioritization rules.

This module ingests already-provided pocket detector-style records and compares
them with Task 8 PPI-guided pocket planning evidence. It never runs fpocket,
P2Rank, docking, scoring, or compound nomination.
"""

import csv
import json


CANDIDATE_REQUIRED_COLUMNS = [
    "candidate_pocket_id",
    "source_detector",
    "source_detector_version",
    "receptor_id",
    "receptor_state",
    "protomer_id",
    "pocket_centroid_x",
    "pocket_centroid_y",
    "pocket_centroid_z",
    "pocket_residues",
    "pocket_surface_residues",
    "pocket_volume_angstrom3",
    "pocket_area_angstrom2",
    "detector_druggability_score",
    "nearest_ppi_patch_id",
    "nearest_ppi_centroid_distance_angstrom",
    "nearest_ppi_residue_distance_angstrom",
    "atp_overlap_fraction",
    "atp_nearest_residue_distance_angstrom",
    "membrane_proximal_fraction",
    "membrane_side_class",
    "dimer_buried_fraction",
    "notes",
]

PRIORITIZED_FIELDS = [
    "candidate_pocket_id",
    "record_semantics",
    "source_detector",
    "source_detector_version",
    "receptor_id",
    "receptor_state",
    "protomer_id",
    "pocket_residues",
    "pocket_surface_residues",
    "nearest_ppi_patch_id",
    "ppi_relationship_class",
    "atp_relationship_class",
    "membrane_compatibility_class",
    "dimer_accessibility_class",
    "final_priority_class",
    "carry_forward",
    "detector_druggability_score",
    "nearest_ppi_centroid_distance_angstrom",
    "nearest_ppi_residue_distance_angstrom",
    "atp_overlap_fraction",
    "membrane_proximal_fraction",
    "dimer_buried_fraction",
    "ppi_residue_overlap_count",
    "blocker",
    "warnings",
    "external_detector_output_ingested",
    "pocket_detector_runtime_executed",
    "planned_compound_docking_jobs",
    "compound_docking_or_scoring_executed",
    "candidate_nomination_executed",
    "notes",
]

FAMILY_FIELDS = [
    "pocket_family_id",
    "nearest_ppi_patch_id",
    "receptor_id",
    "receptor_state",
    "protomer_id",
    "final_priority_class",
    "candidate_pocket_ids",
    "n_candidates",
    "ppi_relationship_classes",
    "atp_relationship_classes",
    "membrane_compatibility_classes",
    "dimer_accessibility_classes",
    "representative_candidate_id",
    "notes",
]

ATP_CORE_RESIDUES = set([745, 855])
PPI_CARRY_FORWARD_CLASSES = set(["orthosteric_ppi_overlap", "ppi_rim", "allosteric_near"])
HARD_BLOCK_PRIORITY = set(["BLOCKED_ATP", "BLOCKED_MEMBRANE", "BLOCKED_DIMER", "BLOCKED_SCHEMA"])


def parse_float(value, default=None):
    if value is None:
        return default
    text = str(value).strip()
    if text == "":
        return default
    try:
        return float(text)
    except ValueError:
        return default


def parse_residue_list(value):
    residues = []
    invalid = []
    text = str(value or "").strip()
    if not text:
        return residues, invalid
    for raw_token in text.split(";"):
        token = raw_token.strip()
        if not token:
            continue
        chain_id = ""
        residue_text = token
        if ":" in token:
            parts = token.split(":", 1)
            chain_id = parts[0].strip()
            residue_text = parts[1].strip()
        digits = ""
        insertion_code = ""
        for char in residue_text:
            if char.isdigit() or (char == "-" and not digits):
                digits += char
            else:
                insertion_code += char
        if not digits or digits == "-":
            invalid.append(token)
            continue
        try:
            residue_number = int(digits)
        except ValueError:
            invalid.append(token)
            continue
        normalized = "{0}:{1}".format(chain_id, residue_number) if chain_id else str(residue_number)
        residues.append({
            "raw": token,
            "chain_id": chain_id,
            "residue_number": residue_number,
            "insertion_code": insertion_code,
            "normalized": normalized,
        })
    return residues, invalid


def residue_key_set(value):
    residues, invalid = parse_residue_list(value)
    return set([item["normalized"] for item in residues]), invalid


def read_csv_rows(path):
    with open(str(path), "r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    return rows, fieldnames


def validate_candidate_schema(fieldnames):
    found = set(fieldnames or [])
    missing = []
    audit = []
    for column in CANDIDATE_REQUIRED_COLUMNS:
        present = column in found
        if not present:
            missing.append(column)
        audit.append({
            "check": "required_column",
            "field": column,
            "status": "PASS" if present else "FAIL",
            "details": "" if present else "missing required pocket candidate column",
        })
    return missing, audit


def load_pocket_plan(path):
    with open(str(path), "r", encoding="utf-8") as handle:
        return json.load(handle)


def _candidate_csv_from_plan_outputs(plan, plan_path, repo_root):
    outputs = plan.get("outputs", {})
    for key in ["pocket_plan_records_csv", "pockets_csv", "legacy_pockets_csv"]:
        rel = outputs.get(key)
        if not rel:
            continue
        candidate = repo_root / rel
        if candidate.exists():
            return candidate
        candidate = plan_path.parent / rel
        if candidate.exists():
            return candidate
    return None


def ppi_patch_residue_map_from_plan(plan, plan_path, repo_root):
    records = []
    for key in ["pocket_records", "pocket_plan_records", "planned_pockets"]:
        if isinstance(plan.get(key), list):
            records.extend(plan.get(key))
    if not records:
        candidate_csv = _candidate_csv_from_plan_outputs(plan, plan_path, repo_root)
        if candidate_csv is not None:
            records = read_csv_rows(candidate_csv)[0]

    mapping = {}
    for row in records:
        patch_id = row.get("ppi_patch_id", "")
        if not patch_id:
            continue
        residues = row.get("egfr_residues", row.get("egfr_consensus_residues", ""))
        keys, _invalid = residue_key_set(residues)
        mapping.setdefault(patch_id, set()).update(keys)
    return mapping


def classify_ppi_relationship(row, ppi_residue_map):
    patch_id = row.get("nearest_ppi_patch_id", "")
    pocket_keys, invalid_a = residue_key_set(row.get("pocket_residues", ""))
    surface_keys, invalid_b = residue_key_set(row.get("pocket_surface_residues", ""))
    candidate_keys = pocket_keys.union(surface_keys)
    ppi_keys = ppi_residue_map.get(patch_id, set())
    overlap = candidate_keys.intersection(ppi_keys)
    residue_distance = parse_float(row.get("nearest_ppi_residue_distance_angstrom"))
    centroid_distance = parse_float(row.get("nearest_ppi_centroid_distance_angstrom"))

    if overlap:
        relation = "orthosteric_ppi_overlap"
    elif residue_distance is not None and residue_distance <= 6.0:
        relation = "ppi_rim"
    elif residue_distance is not None and residue_distance <= 12.0:
        relation = "allosteric_near"
    elif centroid_distance is not None and centroid_distance <= 25.0:
        relation = "allosteric_far"
    else:
        relation = "not_ppi_connected"
    return relation, len(overlap), invalid_a + invalid_b


def classify_atp_relationship(row):
    overlap = parse_float(row.get("atp_overlap_fraction"))
    distance = parse_float(row.get("atp_nearest_residue_distance_angstrom"))
    pocket_residues, invalid = parse_residue_list(row.get("pocket_residues", ""))
    has_core = any(item["residue_number"] in ATP_CORE_RESIDUES for item in pocket_residues)
    if (overlap is not None and overlap >= 0.25) or has_core:
        return "ATP_overlap", invalid
    if distance is not None and distance <= 6.0:
        return "ATP_adjacent", invalid
    if overlap == 0.0 and distance is not None and distance > 6.0:
        return "non_ATP", invalid
    return "ATP_unknown", invalid


def classify_membrane_compatibility(row):
    fraction = parse_float(row.get("membrane_proximal_fraction"))
    side_class = str(row.get("membrane_side_class", "")).strip().lower()
    if fraction is None or not side_class:
        return "membrane_unknown"
    if fraction >= 0.5:
        return "membrane_blocked"
    if fraction > 0.0:
        return "membrane_proximal_risk"
    if side_class in ("solvent_exposed", "cytosolic", "lateral_accessible", "lower_lobe_accessible"):
        return "membrane_compatible"
    return "membrane_unknown"


def classify_dimer_accessibility(row):
    protomer = str(row.get("protomer_id", "")).strip()
    buried = parse_float(row.get("dimer_buried_fraction"))
    if buried is None or protomer == "" or protomer.lower() in ("x", "unknown", "ambiguous"):
        return "dimer_unknown"
    if buried >= 0.5:
        return "dimer_buried"
    if buried > 0.0:
        return "dimer_rim"
    return "dimer_accessible"


def final_priority(ppi_class, atp_class, membrane_class, dimer_class, detector_score):
    if atp_class == "ATP_overlap":
        return "BLOCKED_ATP"
    if membrane_class == "membrane_blocked":
        return "BLOCKED_MEMBRANE"
    if dimer_class in ("dimer_buried", "dimer_unknown"):
        return "BLOCKED_DIMER"
    if ppi_class in PPI_CARRY_FORWARD_CLASSES and atp_class in ("non_ATP", "ATP_unknown") and dimer_class in ("dimer_accessible", "dimer_rim"):
        if membrane_class in ("membrane_compatible", "membrane_proximal_risk"):
            if detector_score is not None and detector_score >= 0.5 and ppi_class in ("orthosteric_ppi_overlap", "ppi_rim"):
                return "CARRY_FORWARD_STRONG"
            return "CARRY_FORWARD_WEAK"
    return "EXPLORATORY_ONLY"


def row_warnings(row, ppi_class, atp_class, membrane_class, dimer_class, invalid_residues):
    warnings = []
    if invalid_residues:
        warnings.append("malformed_residue_ids:{0}".format(";".join(invalid_residues)))
    if ppi_class in ("allosteric_far", "not_ppi_connected"):
        warnings.append("ppi_residue_or_rim_support_insufficient")
    if atp_class != "non_ATP":
        warnings.append("atp_relationship:{0}".format(atp_class))
    if membrane_class != "membrane_compatible":
        warnings.append("membrane_relationship:{0}".format(membrane_class))
    if dimer_class != "dimer_accessible":
        warnings.append("dimer_relationship:{0}".format(dimer_class))
    notes = row.get("notes", "")
    if notes:
        warnings.append(notes)
    return ";".join(warnings)


def prioritize_candidate_rows(candidate_rows, ppi_residue_map):
    prioritized = []
    input_audit = []
    ppi_audit = []
    atp_audit = []
    membrane_audit = []
    dimer_audit = []
    blockers = []

    sorted_rows = sorted(candidate_rows, key=lambda item: (
        item.get("receptor_id", ""),
        item.get("receptor_state", ""),
        item.get("protomer_id", ""),
        item.get("candidate_pocket_id", ""),
    ))
    for row in sorted_rows:
        candidate_id = row.get("candidate_pocket_id", "")
        ppi_class, overlap_count, invalid_ppi = classify_ppi_relationship(row, ppi_residue_map)
        atp_class, invalid_atp = classify_atp_relationship(row)
        membrane_class = classify_membrane_compatibility(row)
        dimer_class = classify_dimer_accessibility(row)
        detector_score = parse_float(row.get("detector_druggability_score"))
        priority = final_priority(ppi_class, atp_class, membrane_class, dimer_class, detector_score)
        carry_forward = priority in ("CARRY_FORWARD_STRONG", "CARRY_FORWARD_WEAK")
        invalid_residues = sorted(set(invalid_ppi + invalid_atp))
        blocker = ""
        if priority in HARD_BLOCK_PRIORITY:
            blocker = priority
            blockers.append({
                "candidate_pocket_id": candidate_id,
                "blocker": priority,
                "details": "{0};{1};{2}".format(atp_class, membrane_class, dimer_class),
            })
        warnings = row_warnings(row, ppi_class, atp_class, membrane_class, dimer_class, invalid_residues)
        out = {
            "candidate_pocket_id": candidate_id,
            "record_semantics": "provided_detector_pocket_candidate_prioritized",
            "source_detector": row.get("source_detector", ""),
            "source_detector_version": row.get("source_detector_version", ""),
            "receptor_id": row.get("receptor_id", ""),
            "receptor_state": row.get("receptor_state", ""),
            "protomer_id": row.get("protomer_id", ""),
            "pocket_residues": row.get("pocket_residues", ""),
            "pocket_surface_residues": row.get("pocket_surface_residues", ""),
            "nearest_ppi_patch_id": row.get("nearest_ppi_patch_id", ""),
            "ppi_relationship_class": ppi_class,
            "atp_relationship_class": atp_class,
            "membrane_compatibility_class": membrane_class,
            "dimer_accessibility_class": dimer_class,
            "final_priority_class": priority,
            "carry_forward": carry_forward,
            "detector_druggability_score": row.get("detector_druggability_score", ""),
            "nearest_ppi_centroid_distance_angstrom": row.get("nearest_ppi_centroid_distance_angstrom", ""),
            "nearest_ppi_residue_distance_angstrom": row.get("nearest_ppi_residue_distance_angstrom", ""),
            "atp_overlap_fraction": row.get("atp_overlap_fraction", ""),
            "membrane_proximal_fraction": row.get("membrane_proximal_fraction", ""),
            "dimer_buried_fraction": row.get("dimer_buried_fraction", ""),
            "ppi_residue_overlap_count": overlap_count,
            "blocker": blocker,
            "warnings": warnings,
            "external_detector_output_ingested": True,
            "pocket_detector_runtime_executed": False,
            "planned_compound_docking_jobs": 0,
            "compound_docking_or_scoring_executed": False,
            "candidate_nomination_executed": False,
            "notes": "intake/prioritization only; detector output was supplied, not generated",
        }
        prioritized.append(out)
        input_audit.append({
            "candidate_pocket_id": candidate_id,
            "source_detector": row.get("source_detector", ""),
            "status": "PASS",
            "notes": "candidate pocket record ingested without executing detector runtime",
        })
        ppi_audit.append({
            "candidate_pocket_id": candidate_id,
            "nearest_ppi_patch_id": row.get("nearest_ppi_patch_id", ""),
            "ppi_relationship_class": ppi_class,
            "ppi_residue_overlap_count": overlap_count,
            "nearest_ppi_centroid_distance_angstrom": row.get("nearest_ppi_centroid_distance_angstrom", ""),
            "nearest_ppi_residue_distance_angstrom": row.get("nearest_ppi_residue_distance_angstrom", ""),
            "status": "PASS" if ppi_class in PPI_CARRY_FORWARD_CLASSES else "WARN",
            "notes": "centroid-only proximity is not sufficient for carry-forward" if ppi_class == "allosteric_far" else "",
        })
        atp_audit.append({
            "candidate_pocket_id": candidate_id,
            "atp_relationship_class": atp_class,
            "atp_overlap_fraction": row.get("atp_overlap_fraction", ""),
            "status": "FAIL" if atp_class == "ATP_overlap" else "PASS",
            "notes": "ATP-overlap is a blocker for EGFR-MYO1D PPI-disruptive pocket carry-forward" if atp_class == "ATP_overlap" else "",
        })
        membrane_audit.append({
            "candidate_pocket_id": candidate_id,
            "membrane_compatibility_class": membrane_class,
            "membrane_proximal_fraction": row.get("membrane_proximal_fraction", ""),
            "status": "FAIL" if membrane_class == "membrane_blocked" else ("WARN" if membrane_class != "membrane_compatible" else "PASS"),
            "notes": "membrane-proximal evidence remains visible and may block carry-forward",
        })
        dimer_audit.append({
            "candidate_pocket_id": candidate_id,
            "protomer_id": row.get("protomer_id", ""),
            "dimer_accessibility_class": dimer_class,
            "dimer_buried_fraction": row.get("dimer_buried_fraction", ""),
            "status": "FAIL" if dimer_class in ("dimer_buried", "dimer_unknown") else "PASS",
            "notes": "protomer identity and dimer accessibility are preserved",
        })
    return prioritized, input_audit, ppi_audit, atp_audit, membrane_audit, dimer_audit, blockers


def build_candidate_families(prioritized_rows):
    groups = {}
    for row in prioritized_rows:
        key = (
            row.get("nearest_ppi_patch_id", ""),
            row.get("receptor_id", ""),
            row.get("receptor_state", ""),
            row.get("protomer_id", ""),
            row.get("final_priority_class", ""),
        )
        groups.setdefault(key, []).append(row)
    family_rows = []
    for index, key in enumerate(sorted(groups.keys()), start=1):
        rows = sorted(groups[key], key=lambda item: item.get("candidate_pocket_id", ""))
        family_rows.append({
            "pocket_family_id": "pocket_family_{0:03d}".format(index),
            "nearest_ppi_patch_id": key[0],
            "receptor_id": key[1],
            "receptor_state": key[2],
            "protomer_id": key[3],
            "final_priority_class": key[4],
            "candidate_pocket_ids": ";".join([row.get("candidate_pocket_id", "") for row in rows]),
            "n_candidates": len(rows),
            "ppi_relationship_classes": ";".join(sorted(set([row.get("ppi_relationship_class", "") for row in rows]))),
            "atp_relationship_classes": ";".join(sorted(set([row.get("atp_relationship_class", "") for row in rows]))),
            "membrane_compatibility_classes": ";".join(sorted(set([row.get("membrane_compatibility_class", "") for row in rows]))),
            "dimer_accessibility_classes": ";".join(sorted(set([row.get("dimer_accessibility_class", "") for row in rows]))),
            "representative_candidate_id": rows[0].get("candidate_pocket_id", ""),
            "notes": "candidate family is an intake/prioritization grouping, not compound nomination",
        })
    return family_rows

