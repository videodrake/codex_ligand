"""Task 7 PPI consensus patch post-processing.

This module consumes supplied PPI contact records and summarizes EGFR-side
consensus evidence. It does not run docking, scoring, pocket discovery, or any
external scientific engine.
"""

import csv
import math
from collections import Counter, defaultdict


REQUIRED_COLUMNS = [
    "receptor_id",
    "receptor_state",
    "method",
    "seed",
    "replicate",
    "pose_id",
    "protomer_id",
    "egfr_contact_residues",
    "egfr_contact_centroid_x",
    "egfr_contact_centroid_y",
    "egfr_contact_centroid_z",
    "myo1d_contact_residues",
    "myo1d_active_face_contacts",
    "sheet12_support_contacts",
    "tail_noise_contacts",
    "membrane_side_class",
    "terminal_artifact_flag",
    "atp_overlap_flag",
    "membrane_proximal_flag",
    "pose_acceptance_class",
]

CONSENSUS_FIELDS = [
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

ACTIVE_FACE_RESIDUES = set([961, 962, 963, 964, 968, 969, 970, 971, 972])
SHEET12_SUPPORT_RESIDUES = set([993, 994, 995, 996, 997])
TAIL_NOISE_RESIDUES = set([998, 999, 1000, 1001, 1002, 1003, 1004, 1005, 1006])

ACCEPTED_CLASSES = set([
    "accepted_active_8_9_supported_12",
    "accepted_active_8_9_only",
    "review_sheet12_dominant",
])


def parse_bool(value):
    text = str(value or "").strip().lower()
    if text in ("1", "true", "yes", "y"):
        return True
    if text in ("0", "false", "no", "n", ""):
        return False
    raise ValueError("invalid boolean value: {0}".format(value))


def parse_float(value, column_name):
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        raise ValueError("{0} must be numeric: {1}".format(column_name, value))


def split_residue_list(value):
    text = str(value or "").strip()
    if not text:
        return []
    return [item.strip() for item in text.split(";") if item.strip()]


def parse_egfr_residue_token(token):
    if ":" not in token:
        raise ValueError("EGFR residue must preserve chain as chain:number: {0}".format(token))
    chain_id, number_text = token.split(":", 1)
    chain_id = chain_id.strip()
    number_text = number_text.strip()
    if not chain_id:
        raise ValueError("EGFR residue has empty chain ID: {0}".format(token))
    if not number_text.isdigit():
        raise ValueError("EGFR residue number is not numeric: {0}".format(token))
    residue_number = int(number_text)
    return {
        "chain_id": chain_id,
        "residue_number": residue_number,
        "canonical": "{0}:{1}".format(chain_id, residue_number),
    }


def parse_myo1d_residue_token(token):
    raw = token.strip()
    if ":" in raw:
        _chain_id, raw = raw.split(":", 1)
        raw = raw.strip()
    if not raw.isdigit():
        raise ValueError("MYO1D residue number is not numeric: {0}".format(token))
    residue_number = int(raw)
    return {"residue_number": residue_number, "canonical": str(residue_number)}


def parse_residue_list(value, residue_kind):
    residues = []
    errors = []
    parser = parse_egfr_residue_token if residue_kind == "egfr" else parse_myo1d_residue_token
    for token in split_residue_list(value):
        try:
            residues.append(parser(token))
        except ValueError as exc:
            errors.append(str(exc))
    return residues, errors


def validate_schema(fieldnames):
    fieldnames = fieldnames or []
    found = set(fieldnames)
    rows = []
    missing = []
    for column in REQUIRED_COLUMNS:
        present = column in found
        if not present:
            missing.append(column)
        rows.append({
            "column": column,
            "required": True,
            "present": present,
            "status": "PASS" if present else "FAIL",
            "notes": "" if present else "missing required column",
        })
    for column in sorted(found.difference(REQUIRED_COLUMNS)):
        rows.append({
            "column": column,
            "required": False,
            "present": True,
            "status": "WARN",
            "notes": "extra column preserved by input but ignored by Task 7",
        })
    return missing, rows


def load_contact_table(path):
    with open(str(path), "r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing, schema_rows = validate_schema(reader.fieldnames)
        rows = list(reader) if not missing else []
    return missing, schema_rows, rows


def normalize_contact_rows(rows):
    normalized = []
    parsing_rows = []
    active_rows = []
    tail_rows = []
    membrane_rows = []
    errors = []
    for index, row in enumerate(rows, start=1):
        pose_id = row.get("pose_id", "")
        parsed_record = dict(row)
        egfr_residues, egfr_errors = parse_residue_list(row.get("egfr_contact_residues"), "egfr")
        myo_residues, myo_errors = parse_residue_list(row.get("myo1d_contact_residues"), "myo1d")
        active_residues, active_errors = parse_residue_list(row.get("myo1d_active_face_contacts"), "myo1d")
        sheet12_residues, sheet12_errors = parse_residue_list(row.get("sheet12_support_contacts"), "myo1d")
        tail_residues, tail_errors = parse_residue_list(row.get("tail_noise_contacts"), "myo1d")
        parse_errors = egfr_errors + myo_errors + active_errors + sheet12_errors + tail_errors
        for column, parsed, column_errors in [
            ("egfr_contact_residues", egfr_residues, egfr_errors),
            ("myo1d_contact_residues", myo_residues, myo_errors),
            ("myo1d_active_face_contacts", active_residues, active_errors),
            ("sheet12_support_contacts", sheet12_residues, sheet12_errors),
            ("tail_noise_contacts", tail_residues, tail_errors),
        ]:
            parsing_rows.append({
                "row_index": index,
                "pose_id": pose_id,
                "column": column,
                "status": "FAIL" if column_errors else "PASS",
                "parsed_count": len(parsed),
                "errors": "; ".join(column_errors),
            })
        try:
            centroid = [
                parse_float(row.get("egfr_contact_centroid_x"), "egfr_contact_centroid_x"),
                parse_float(row.get("egfr_contact_centroid_y"), "egfr_contact_centroid_y"),
                parse_float(row.get("egfr_contact_centroid_z"), "egfr_contact_centroid_z"),
            ]
        except ValueError as exc:
            centroid = [0.0, 0.0, 0.0]
            parse_errors.append(str(exc))
        try:
            terminal_artifact_flag = parse_bool(row.get("terminal_artifact_flag"))
            atp_overlap_flag = parse_bool(row.get("atp_overlap_flag"))
            membrane_proximal_flag = parse_bool(row.get("membrane_proximal_flag"))
        except ValueError as exc:
            terminal_artifact_flag = False
            atp_overlap_flag = False
            membrane_proximal_flag = False
            parse_errors.append(str(exc))

        myo_numbers = set([item["residue_number"] for item in myo_residues])
        active_numbers = set([item["residue_number"] for item in active_residues])
        sheet12_numbers = set([item["residue_number"] for item in sheet12_residues])
        explicit_tail_numbers = set([item["residue_number"] for item in tail_residues])
        inferred_tail_numbers = myo_numbers.intersection(TAIL_NOISE_RESIDUES)
        tail_numbers = explicit_tail_numbers.union(inferred_tail_numbers)

        active_rows.append({
            "pose_id": pose_id,
            "active_face_contact_count": len(active_numbers.intersection(ACTIVE_FACE_RESIDUES)),
            "sheet12_support_count": len(sheet12_numbers.intersection(SHEET12_SUPPORT_RESIDUES)),
            "score_bonus_allowed": False,
            "key_residue_bonus_weight": 0.0,
            "status": "PASS",
            "notes": "annotation only; no scoring bonus",
        })
        tail_rows.append({
            "pose_id": pose_id,
            "tail_noise_contact_count": len(tail_numbers),
            "terminal_artifact_flag": terminal_artifact_flag,
            "quarantine_recommended": terminal_artifact_flag or len(tail_numbers) > 0,
            "status": "WARN" if terminal_artifact_flag or len(tail_numbers) > 0 else "PASS",
            "notes": "tail/noise residues 998-1006 are not primary binding evidence" if tail_numbers else "",
        })
        membrane_rows.append({
            "pose_id": pose_id,
            "atp_overlap_flag": atp_overlap_flag,
            "membrane_proximal_flag": membrane_proximal_flag,
            "membrane_side_class": row.get("membrane_side_class", ""),
            "status": "WARN" if atp_overlap_flag or membrane_proximal_flag else "PASS",
            "notes": "non-target or membrane-proximal evidence must stay visible" if atp_overlap_flag or membrane_proximal_flag else "",
        })
        if parse_errors:
            errors.extend(["row {0}: {1}".format(index, item) for item in parse_errors])
        parsed_record["_row_index"] = index
        parsed_record["_egfr_residues"] = egfr_residues
        parsed_record["_myo1d_residue_numbers"] = sorted(myo_numbers)
        parsed_record["_active_face_count"] = len(active_numbers.intersection(ACTIVE_FACE_RESIDUES))
        parsed_record["_sheet12_count"] = len(sheet12_numbers.intersection(SHEET12_SUPPORT_RESIDUES))
        parsed_record["_tail_noise_count"] = len(tail_numbers)
        parsed_record["_centroid"] = centroid
        parsed_record["_terminal_artifact_flag"] = terminal_artifact_flag
        parsed_record["_atp_overlap_flag"] = atp_overlap_flag
        parsed_record["_membrane_proximal_flag"] = membrane_proximal_flag
        parsed_record["_parse_errors"] = parse_errors
        normalized.append(parsed_record)
    return normalized, parsing_rows, active_rows, tail_rows, membrane_rows, errors


def mean(values):
    values = list(values)
    return sum(values) / float(len(values)) if values else 0.0


def distance(a, b):
    return math.sqrt(sum((a[i] - b[i]) * (a[i] - b[i]) for i in range(3)))


def fraction(count, total):
    return count / float(total) if total else 0.0


def format_fraction(value):
    return "{0:.3f}".format(value)


def residue_occupancy_summary(counter):
    items = sorted(counter.items(), key=lambda item: (-item[1], item[0]))
    return ";".join(["{0}={1}".format(residue, count) for residue, count in items])


def consensus_residue_list(counter):
    if not counter:
        return ""
    max_count = max(counter.values())
    threshold = 2 if max_count >= 2 else 1
    residues = sorted([residue for residue, count in counter.items() if count >= threshold])
    return ";".join(residues)


def classify_patch(n_poses, residue_counter, centroid_spread, tail_fraction, terminal_fraction, atp_fraction, membrane_fraction):
    if n_poses == 0:
        return "INSUFFICIENT_EVIDENCE"
    if terminal_fraction > 0.0 or tail_fraction >= 0.5:
        return "TAIL_DOMINATED"
    if atp_fraction >= 0.5:
        return "ATP_OVERLAPPING"
    if membrane_fraction >= 0.5:
        return "MEMBRANE_PROXIMAL"
    if centroid_spread > 20.0:
        return "DISPERSED"
    if n_poses < 2:
        return "INSUFFICIENT_EVIDENCE"
    repeated = [residue for residue, count in residue_counter.items() if count >= 2]
    if repeated and centroid_spread <= 8.0:
        return "CONVERGENT_PATCH"
    return "BROAD_PATCH"


def warning_text(evidence_class, tail_fraction, terminal_fraction, atp_fraction, membrane_fraction, centroid_spread):
    warnings = []
    if evidence_class == "TAIL_DOMINATED" or tail_fraction > 0.0 or terminal_fraction > 0.0:
        warnings.append("tail_or_terminal_artifact_visible")
    if atp_fraction > 0.0:
        warnings.append("atp_overlap_visible")
    if membrane_fraction > 0.0:
        warnings.append("membrane_proximal_visible")
    if evidence_class == "DISPERSED" or centroid_spread > 20.0:
        warnings.append("dispersed_centroids")
    return ";".join(warnings)


def build_consensus_patches(records):
    groups = defaultdict(list)
    for record in records:
        key = (
            record.get("receptor_id", ""),
            record.get("receptor_state", ""),
            record.get("protomer_id", ""),
        )
        groups[key].append(record)

    patch_rows = []
    convergence_rows = []
    for index, key in enumerate(sorted(groups.keys()), start=1):
        receptor_id, receptor_state, protomer_id = key
        group = groups[key]
        n_poses = len(group)
        residue_counter = Counter()
        for record in group:
            for residue in record["_egfr_residues"]:
                residue_counter[residue["canonical"]] += 1
        centroid = [
            mean([record["_centroid"][0] for record in group]),
            mean([record["_centroid"][1] for record in group]),
            mean([record["_centroid"][2] for record in group]),
        ]
        centroid_spread = max([distance(record["_centroid"], centroid) for record in group]) if group else 0.0
        active_fraction = fraction(sum(1 for record in group if record["_active_face_count"] > 0), n_poses)
        sheet12_fraction = fraction(sum(1 for record in group if record["_sheet12_count"] > 0), n_poses)
        tail_fraction = fraction(sum(1 for record in group if record["_tail_noise_count"] > 0), n_poses)
        terminal_fraction = fraction(sum(1 for record in group if record["_terminal_artifact_flag"]), n_poses)
        atp_fraction = fraction(sum(1 for record in group if record["_atp_overlap_flag"]), n_poses)
        membrane_fraction = fraction(sum(1 for record in group if record["_membrane_proximal_flag"]), n_poses)
        evidence_class = classify_patch(n_poses, residue_counter, centroid_spread, tail_fraction, terminal_fraction, atp_fraction, membrane_fraction)
        patch_id = "ppi_patch_{0:03d}".format(index)
        methods = sorted(set([record.get("method", "") for record in group]))
        seeds = sorted(set([record.get("seed", "") for record in group]), key=lambda value: str(value))
        replicates = sorted(set([record.get("replicate", "") for record in group]), key=lambda value: str(value))
        membrane_classes = sorted(set([record.get("membrane_side_class", "") for record in group if record.get("membrane_side_class", "")]))
        row = {
            "ppi_patch_id": patch_id,
            "receptor_id": receptor_id,
            "receptor_state": receptor_state,
            "protomer_id": protomer_id,
            "support_methods": ";".join(methods),
            "support_seeds": ";".join([str(item) for item in seeds]),
            "support_replicates": ";".join([str(item) for item in replicates]),
            "n_accepted_poses": n_poses,
            "egfr_consensus_residues": consensus_residue_list(residue_counter),
            "egfr_contact_centroid_x": "{0:.3f}".format(centroid[0]),
            "egfr_contact_centroid_y": "{0:.3f}".format(centroid[1]),
            "egfr_contact_centroid_z": "{0:.3f}".format(centroid[2]),
            "centroid_spread_angstrom": "{0:.3f}".format(centroid_spread),
            "residue_occupancy_summary": residue_occupancy_summary(residue_counter),
            "active_face_support_fraction": format_fraction(active_fraction),
            "sheet12_support_fraction": format_fraction(sheet12_fraction),
            "tail_noise_fraction": format_fraction(tail_fraction),
            "terminal_artifact_fraction": format_fraction(terminal_fraction),
            "atp_overlap_fraction": format_fraction(atp_fraction),
            "membrane_proximal_fraction": format_fraction(membrane_fraction),
            "membrane_side_classes": ";".join(membrane_classes),
            "evidence_class": evidence_class,
            "warnings": warning_text(evidence_class, tail_fraction, terminal_fraction, atp_fraction, membrane_fraction, centroid_spread),
        }
        patch_rows.append(row)
        convergence_rows.append({
            "ppi_patch_id": patch_id,
            "receptor_id": receptor_id,
            "receptor_state": receptor_state,
            "protomer_id": protomer_id,
            "n_accepted_poses": n_poses,
            "unique_egfr_residue_count": len(residue_counter),
            "max_residue_occupancy": max(residue_counter.values()) if residue_counter else 0,
            "centroid_spread_angstrom": row["centroid_spread_angstrom"],
            "evidence_class": evidence_class,
            "status": "WARN" if evidence_class in ("DISPERSED", "TAIL_DOMINATED", "ATP_OVERLAPPING", "MEMBRANE_PROXIMAL") else "PASS",
            "notes": row["warnings"],
        })
    return patch_rows, convergence_rows

