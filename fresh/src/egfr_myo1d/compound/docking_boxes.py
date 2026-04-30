"""M3-T3 docking-box validation helpers.

This module validates M2 accepted pocket boxes only. It does not generate
Vina manifests, run docking, parse poses, score compounds, or nominate
candidates.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from egfr_myo1d.core.run_context import RunContext


BOX_MANIFEST_FIELDS = [
    "pocket_family_id",
    "state_id",
    "state_role",
    "protomer_id",
    "box_id",
    "source_box_file",
    "source_pocket_file",
    "source_pocket_status",
    "original_box_center_x",
    "original_box_center_y",
    "original_box_center_z",
    "original_box_size_x",
    "original_box_size_y",
    "original_box_size_z",
    "box_margin",
    "box_center_x",
    "box_center_y",
    "box_center_z",
    "box_size_x",
    "box_size_y",
    "box_size_z",
    "ppi_relationship_class",
    "geometry_class",
    "non_atp_pass",
    "lower_lateral_pass",
    "dimer_accessibility_pass",
    "atp_overlap_status",
    "receptor_pdbqt_file",
    "receptor_pdbqt_sha256",
    "receptor_mapping_status",
    "membrane_frame_status",
    "traceability_status",
    "box_qc_status",
    "box_qc_notes",
    "allowed_for_vina_smoke",
]

BOX_QC_FIELDS = [
    "pocket_family_id",
    "state_id",
    "protomer_id",
    "box_id",
    "check_id",
    "category",
    "status",
    "severity",
    "details",
    "recommended_fix",
]

REQUIRED_BOX_COLUMNS = [
    "pocket_family_id",
    "state_id",
    "protomer_id",
    "box_center_x",
    "box_center_y",
    "box_center_z",
    "box_size_x",
    "box_size_y",
    "box_size_z",
]

CENTER_COLUMNS = ["box_center_x", "box_center_y", "box_center_z"]
SIZE_COLUMNS = ["box_size_x", "box_size_y", "box_size_z"]
BOX_CHECK_IDS = [
    "accepted_box_input_present",
    "accepted_pocket_family_traceable",
    "box_state_id_present",
    "box_protomer_id_present",
    "box_id_present",
    "box_center_numeric",
    "box_size_numeric",
    "box_size_positive",
    "box_margin_valid",
    "box_state_has_receptor_pdbqt",
    "box_protomer_maps_to_receptor",
    "box_non_atp_pass",
    "box_lower_lateral_pass",
    "box_dimer_accessibility_pass",
    "box_atp_overlap",
    "box_inside_receptor_sanity_bounds",
    "box_allowed_for_vina_smoke",
    "non_goals_preserved",
]


@dataclass(frozen=True)
class PreparedState:
    state_id: str
    state_role: str
    status: str
    receptor_pdbqt_file: str | None
    receptor_pdbqt_sha256: str | None
    receptor_mapping_file: str | None
    membrane_frame_file: str | None
    protomer_ids: set[str]
    coordinate_bounds: tuple[tuple[float, float], tuple[float, float], tuple[float, float]] | None = None


@dataclass
class DockingBoxValidationResult:
    manifest_rows: list[dict[str, Any]]
    qc_rows: list[dict[str, str]]
    blockers: list[str]
    warnings: list[str]
    counts: dict[str, int]
    all_active_boxes_traceable_to_m2: bool
    all_active_boxes_have_receptor: bool
    all_active_boxes_non_atp: bool
    all_active_boxes_valid_dimensions: bool


def read_csv_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def boolish(value: Any) -> bool | None:
    text = str(value or "").strip().lower()
    if text in {"true", "t", "1", "yes", "y", "pass", "passed"}:
        return True
    if text in {"false", "f", "0", "no", "n", "fail", "failed"}:
        return False
    if text == "":
        return None
    return None


def finite_float(value: Any) -> float | None:
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    if math.isnan(number) or math.isinf(number):
        return None
    return number


def format_float(value: float) -> str:
    return ("{0:.6f}".format(value)).rstrip("0").rstrip(".")


def _severity(status: str, blocker: bool = False) -> str:
    if status == "FAIL" and blocker:
        return "BLOCKER"
    if status in {"FAIL", "WARN"}:
        return "MAJOR"
    if status == "NOT_APPLICABLE":
        return "MINOR"
    return "INFO"


def append_box_qc(
    rows: list[dict[str, str]],
    box_row: dict[str, str],
    check_id: str,
    category: str,
    status: str,
    details: str,
    recommended_fix: str = "",
    blocker: bool = False,
) -> None:
    rows.append(
        {
            "pocket_family_id": (box_row.get("pocket_family_id") or "").strip(),
            "state_id": (box_row.get("state_id") or "").strip(),
            "protomer_id": (box_row.get("protomer_id") or "").strip(),
            "box_id": (box_row.get("box_id") or "").strip(),
            "check_id": check_id,
            "category": category,
            "status": status,
            "severity": _severity(status, blocker),
            "details": details,
            "recommended_fix": recommended_fix,
        }
    )


def source_pocket_ids(rows: list[dict[str, str]]) -> set[str]:
    accepted: set[str] = set()
    for row in rows:
        family = (row.get("pocket_family_id") or "").strip()
        if not family:
            continue
        if m3_pocket_authorized(row):
            accepted.add(family)
    return accepted


def _row_text(row: dict[str, str], keys: list[str]) -> str:
    return " ".join(str(row.get(key, "")) for key in keys).strip().lower()


def _is_reference_only(row: dict[str, str]) -> bool:
    if boolish(row.get("reference_only")) is True:
        return True
    return any((row.get(key) or "").strip().lower() == "reference_only" for key in ["export_tier", "acceptance_tier", "gate_class"])


def _has_failure_text(row: dict[str, str]) -> bool:
    status_text = _row_text(
        row,
        [
            "pocket_acceptance_status",
            "box_export_status",
            "box_qc_status",
            "gate_class",
            "acceptance_reason",
            "review_notes",
        ],
    )
    return "fail" in status_text or "reject" in status_text


def m3_pocket_authorized(row: dict[str, str]) -> bool:
    if _is_reference_only(row) or _has_failure_text(row):
        return False
    if "hard_gate_pass" in row and boolish(row.get("hard_gate_pass")) is not True:
        return False
    export_tier = (row.get("export_tier") or row.get("acceptance_tier") or "").strip().lower()
    gate_class = (row.get("gate_class") or "").strip().lower()
    accepted_status = _row_text(row, ["pocket_acceptance_status"])
    primary_allowed = boolish(row.get("m3_primary_allowed")) is True
    review_allowed = boolish(row.get("m3_review_allowed")) is True
    primary_tier = export_tier == "primary" or gate_class == "accepted_primary_pocket"
    review_tier = export_tier == "secondary_review" or gate_class == "accepted_secondary_cryptic"
    legacy_accepted = "accepted" in accepted_status and export_tier not in {"reference_only", "rejected"}
    if primary_tier:
        return primary_allowed if "m3_primary_allowed" in row else legacy_accepted
    if review_tier:
        return review_allowed if "m3_review_allowed" in row else legacy_accepted
    return primary_allowed or review_allowed or legacy_accepted


def m3_box_authorized(row: dict[str, str], pocket_authorized: bool) -> bool:
    if not pocket_authorized or _is_reference_only(row) or _has_failure_text(row):
        return False
    export_tier = (row.get("export_tier") or row.get("acceptance_tier") or "").strip().lower()
    primary_allowed = boolish(row.get("m3_primary_allowed")) is True
    review_allowed = boolish(row.get("m3_review_allowed")) is True
    accepted_status = _row_text(row, ["pocket_acceptance_status"])
    legacy_accepted = "accepted" in accepted_status
    if export_tier == "primary":
        return primary_allowed if "m3_primary_allowed" in row else legacy_accepted
    if export_tier == "secondary_review":
        return review_allowed if "m3_review_allowed" in row else True
    return primary_allowed or review_allowed or legacy_accepted


def validate_docking_boxes(
    ctx: RunContext,
    box_file: Path | None,
    pocket_file: Path | None,
    box_margin: float,
    prepared_states: dict[str, PreparedState],
    mode: str = "prepare",
) -> DockingBoxValidationResult:
    blockers: list[str] = []
    warnings: list[str] = []
    manifest_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, str]] = []
    if box_margin < 0:
        blockers.append("box_margin must be non-negative")

    pocket_rows: list[dict[str, str]] = []
    pocket_ids: set[str] = set()
    if pocket_file and pocket_file.is_file():
        pocket_rows, _ = read_csv_rows(pocket_file)
        pocket_ids = source_pocket_ids(pocket_rows)
    elif mode == "prepare":
        blockers.append("accepted_pockets_for_m3.csv absent in prepare mode")
    else:
        warnings.append("accepted_pockets_for_m3.csv absent in dry-run mode")

    if not box_file or not box_file.is_file():
        pseudo = {"pocket_family_id": "", "state_id": "", "protomer_id": "", "box_id": ""}
        status = "FAIL" if mode == "prepare" else "WARN"
        append_box_qc(
            qc_rows,
            pseudo,
            "accepted_box_input_present",
            "input",
            status,
            "M2 accepted_pocket_boxes.csv not found",
            "Run M2.8 export or provide accepted pocket boxes.",
            mode == "prepare",
        )
        if mode == "prepare":
            blockers.append("M2 accepted_pocket_boxes.csv absent in prepare mode")
        else:
            warnings.append("M2 accepted_pocket_boxes.csv absent in dry-run mode")
        return DockingBoxValidationResult(
            manifest_rows=[],
            qc_rows=qc_rows,
            blockers=blockers,
            warnings=warnings,
            counts={"boxes_input": 0, "boxes_pass": 0, "boxes_warn": 0, "boxes_fail": 0, "boxes_allowed_for_vina_smoke": 0},
            all_active_boxes_traceable_to_m2=False,
            all_active_boxes_have_receptor=False,
            all_active_boxes_non_atp=False,
            all_active_boxes_valid_dimensions=False,
        )

    rows, fields = read_csv_rows(box_file)
    missing_columns = [field for field in REQUIRED_BOX_COLUMNS if field not in fields]
    if missing_columns:
        blockers.append("accepted_pocket_boxes.csv is missing required columns")

    pass_count = 0
    warn_count = 0
    fail_count = 0
    smoke_count = 0
    all_traceable = True
    all_have_receptor = True
    all_non_atp = True
    all_dimensions = True

    for row in rows:
        notes: list[str] = []
        row_fail = bool(missing_columns)
        row_warn = False
        family = (row.get("pocket_family_id") or "").strip()
        state_id = (row.get("state_id") or "").strip()
        protomer_id = (row.get("protomer_id") or "").strip()
        box_id = (row.get("box_id") or "").strip()
        state = prepared_states.get(state_id)

        append_box_qc(qc_rows, row, "accepted_box_input_present", "input", "PASS", "accepted box input row present")
        pocket_authorized = bool(family and family in pocket_ids)
        if family:
            trace_status = "PASS" if pocket_authorized else "FAIL"
        else:
            trace_status = "FAIL"
        append_box_qc(qc_rows, row, "accepted_pocket_family_traceable", "traceability", trace_status, "pocket family traceability checked", "Ensure box pocket_family_id is present in accepted_pockets_for_m3.csv.", trace_status == "FAIL")
        if trace_status == "FAIL":
            row_fail = True
            all_traceable = False
            notes.append("pocket family not traceable to M2 accepted pockets")

        for check_id, value, label in [
            ("box_state_id_present", state_id, "state_id"),
            ("box_protomer_id_present", protomer_id, "protomer_id"),
        ]:
            ok = bool(value)
            append_box_qc(qc_rows, row, check_id, "schema", "PASS" if ok else "FAIL", f"{label} checked", f"Populate {label}.", not ok)
            if not ok:
                row_fail = True
                notes.append(f"missing {label}")
        append_box_qc(qc_rows, row, "box_id_present", "schema", "PASS" if box_id else "WARN", "box_id checked", "Populate box_id for easier traceability.", False)
        if not box_id:
            row_warn = True

        centers = [finite_float(row.get(col)) for col in CENTER_COLUMNS]
        sizes = [finite_float(row.get(col)) for col in SIZE_COLUMNS]
        centers_ok = all(value is not None for value in centers)
        sizes_numeric = all(value is not None for value in sizes)
        sizes_positive = sizes_numeric and all(float(value) > 0 for value in sizes if value is not None)
        append_box_qc(qc_rows, row, "box_center_numeric", "geometry", "PASS" if centers_ok else "FAIL", "center numeric check", "Use finite numeric center coordinates.", not centers_ok)
        append_box_qc(qc_rows, row, "box_size_numeric", "geometry", "PASS" if sizes_numeric else "FAIL", "size numeric check", "Use finite numeric box sizes.", not sizes_numeric)
        append_box_qc(qc_rows, row, "box_size_positive", "geometry", "PASS" if sizes_positive else "FAIL", "positive size check", "Use strictly positive box sizes.", not sizes_positive)
        if not centers_ok or not sizes_numeric or not sizes_positive:
            row_fail = True
            all_dimensions = False
            notes.append("malformed or non-positive box dimensions")
        append_box_qc(qc_rows, row, "box_margin_valid", "geometry", "PASS" if box_margin >= 0 else "FAIL", "box margin checked", "Use a non-negative box margin.", box_margin < 0)

        receptor_ok = bool(state and state.status in {"PASS", "WARN"} and state.receptor_pdbqt_file)
        append_box_qc(qc_rows, row, "box_state_has_receptor_pdbqt", "receptor", "PASS" if receptor_ok else "FAIL", "state receptor PDBQT checked", "Prepare receptor PDBQT for this state.", not receptor_ok)
        if not receptor_ok:
            row_fail = True
            all_have_receptor = False
            notes.append("state has no prepared receptor PDBQT")

        protomer_ok = bool(state and protomer_id and protomer_id in state.protomer_ids)
        append_box_qc(qc_rows, row, "box_protomer_maps_to_receptor", "traceability", "PASS" if protomer_ok else "FAIL", "protomer mapping checked", "Use a protomer_id present in receptor_mapping.csv.", not protomer_ok)
        if not protomer_ok:
            row_fail = True
            notes.append("protomer_id cannot be mapped to receptor")

        non_atp = boolish(row.get("non_atp_pass"))
        lower = boolish(row.get("lower_lateral_pass"))
        dimer = boolish(row.get("dimer_accessibility_pass"))
        exploratory = "explor" in str(row.get("export_tier", "") + row.get("geometry_class", "") + row.get("warnings", "")).lower()
        append_box_qc(qc_rows, row, "box_non_atp_pass", "m2_gate", "PASS" if non_atp is True else "FAIL", "non-ATP gate checked", "Use only M2 non-ATP accepted boxes.", non_atp is not True)
        if non_atp is not True:
            row_fail = True
            all_non_atp = False
            notes.append("non_atp_pass is not true")
        lower_status = "PASS" if lower is True else ("WARN" if exploratory else "FAIL")
        dimer_status = "PASS" if dimer is True else ("WARN" if exploratory else "FAIL")
        append_box_qc(qc_rows, row, "box_lower_lateral_pass", "m2_gate", lower_status, "lower/lateral gate checked", "Use M2 lower/lateral accepted boxes or mark exploratory and exclude.", lower_status == "FAIL")
        append_box_qc(qc_rows, row, "box_dimer_accessibility_pass", "m2_gate", dimer_status, "dimer accessibility gate checked", "Use M2 dimer-accessible boxes or mark exploratory and exclude.", dimer_status == "FAIL")
        if lower_status == "FAIL" or dimer_status == "FAIL":
            row_fail = True
            notes.append("lower/lateral or dimer accessibility gate is not true")
        elif lower_status == "WARN" or dimer_status == "WARN":
            row_warn = True
            notes.append("exploratory box excluded from smoke readiness")

        append_box_qc(qc_rows, row, "box_atp_overlap", "m2_gate", "NOT_APPLICABLE", "ATP overlap re-check not available in M3-T3; M2 non_atp_pass is enforced", "Use M2 ATP reference/gates if re-check is required.")
        bounds_status = "NOT_APPLICABLE"
        if receptor_ok and state and state.coordinate_bounds and centers_ok and sizes_positive:
            bounds_status = "PASS"
            for axis, center, size in zip(state.coordinate_bounds, centers, sizes):
                assert center is not None and size is not None
                if center + size / 2 < axis[0] - 25.0 or center - size / 2 > axis[1] + 25.0:
                    bounds_status = "FAIL"
                    break
        append_box_qc(qc_rows, row, "box_inside_receptor_sanity_bounds", "geometry", bounds_status, "receptor coordinate sanity bounds checked", "Check M2 box coordinates against the normalized receptor.", bounds_status == "FAIL")
        if bounds_status == "FAIL":
            row_fail = True
            notes.append("box outside receptor coordinate sanity bounds")

        box_authorized = m3_box_authorized(row, pocket_authorized)
        if not box_authorized:
            row_fail = True
            notes.append("box is not authorized by M2.8 M3 export flags")

        allowed = not row_fail and lower is True and dimer is True and non_atp is True
        append_box_qc(qc_rows, row, "box_allowed_for_vina_smoke", "readiness", "PASS" if allowed else ("FAIL" if row_fail else "WARN"), "smoke-readiness flag computed", "Fix blocker checks before Vina smoke.", row_fail)
        append_box_qc(qc_rows, row, "non_goals_preserved", "scope", "PASS", "no Vina, docking jobs, pose parsing, scoring, or candidate nomination")

        if row_fail:
            status = "FAIL"
            fail_count += 1
            blockers.append(f"box {box_id or family or '<missing>'}: validation failed")
        elif row_warn:
            status = "WARN"
            warn_count += 1
            warnings.append(f"box {box_id or family}: validation warning")
        else:
            status = "PASS"
            pass_count += 1
        if allowed:
            smoke_count += 1

        final_sizes = sizes
        if not row_fail and sizes_positive and box_margin > 0:
            final_sizes = [float(value) + 2.0 * box_margin for value in sizes if value is not None]
        final_centers = centers
        manifest_rows.append(
            {
                "pocket_family_id": family,
                "state_id": state_id,
                "state_role": row.get("state_role") or (state.state_role if state else "unknown"),
                "protomer_id": protomer_id,
                "box_id": box_id,
                "source_box_file": ctx.relative_to_repo(box_file),
                "source_pocket_file": ctx.relative_to_repo(pocket_file) if pocket_file else "",
                "source_pocket_status": "accepted" if pocket_authorized else "unknown",
                "original_box_center_x": row.get("box_center_x", ""),
                "original_box_center_y": row.get("box_center_y", ""),
                "original_box_center_z": row.get("box_center_z", ""),
                "original_box_size_x": row.get("box_size_x", ""),
                "original_box_size_y": row.get("box_size_y", ""),
                "original_box_size_z": row.get("box_size_z", ""),
                "box_margin": format_float(box_margin),
                "box_center_x": format_float(float(final_centers[0])) if centers_ok else "",
                "box_center_y": format_float(float(final_centers[1])) if centers_ok else "",
                "box_center_z": format_float(float(final_centers[2])) if centers_ok else "",
                "box_size_x": format_float(float(final_sizes[0])) if not row_fail and sizes_positive else "",
                "box_size_y": format_float(float(final_sizes[1])) if not row_fail and sizes_positive else "",
                "box_size_z": format_float(float(final_sizes[2])) if not row_fail and sizes_positive else "",
                "ppi_relationship_class": row.get("ppi_relationship_class", ""),
                "geometry_class": row.get("geometry_class", ""),
                "non_atp_pass": "true" if non_atp is True else "false",
                "lower_lateral_pass": "true" if lower is True else "false",
                "dimer_accessibility_pass": "true" if dimer is True else "false",
                "atp_overlap_status": "NOT_APPLICABLE",
                "receptor_pdbqt_file": state.receptor_pdbqt_file if state and state.receptor_pdbqt_file else "",
                "receptor_pdbqt_sha256": state.receptor_pdbqt_sha256 if state and state.receptor_pdbqt_sha256 else "",
                "receptor_mapping_status": "PASS" if protomer_ok else "FAIL",
                "membrane_frame_status": "PASS" if state and state.membrane_frame_file else "FAIL",
                "traceability_status": "PASS" if trace_status == "PASS" and protomer_ok else "FAIL",
                "box_qc_status": status,
                "box_qc_notes": "; ".join(notes),
                "allowed_for_vina_smoke": "true" if allowed else "false",
            }
        )

    counts = {
        "boxes_input": len(rows),
        "boxes_pass": pass_count,
        "boxes_warn": warn_count,
        "boxes_fail": fail_count,
        "boxes_allowed_for_vina_smoke": smoke_count,
    }
    if rows and smoke_count == 0:
        blockers.append("no accepted boxes after validation")
    return DockingBoxValidationResult(
        manifest_rows=manifest_rows,
        qc_rows=qc_rows,
        blockers=blockers,
        warnings=warnings,
        counts=counts,
        all_active_boxes_traceable_to_m2=all_traceable and bool(rows),
        all_active_boxes_have_receptor=all_have_receptor and bool(rows),
        all_active_boxes_non_atp=all_non_atp and bool(rows),
        all_active_boxes_valid_dimensions=all_dimensions and bool(rows),
    )
