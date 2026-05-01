"""M2.7 membrane lower/lateral geometry gate."""

from __future__ import annotations

import math
from typing import Any

from egfr_myo1d.pocket.ppi_relationship import format_float, point, safe_float, split_values


MEMBRANE_GEOMETRY_FIELDS = [
    "pocket_family_id",
    "representative_candidate_pocket_id",
    "state_ids",
    "state_roles",
    "protomer_ids",
    "membrane_frame_id",
    "membrane_frame_status",
    "family_center_x",
    "family_center_y",
    "family_center_z",
    "pocket_membrane_z",
    "z_membrane_proximity",
    "z_percentile_within_receptor_surface",
    "lower_gate_class",
    "lower_gate_pass",
    "lateral_score",
    "lateral_gate_pass",
    "lower_lateral_pass",
    "membrane_geometry_reason",
    "warnings",
]


def evaluate_membrane_geometry(
    family: dict[str, Any],
    membrane_frame: dict[str, Any] | None,
    receptor_z_ranges: dict[str, tuple[float, float]],
    thresholds: dict[str, float],
    production: bool = False,
) -> dict[str, Any]:
    warnings: list[str] = []
    center = point(family, "family_center_x", "family_center_y", "family_center_z")
    states = split_values(family.get("member_state_ids", ""))
    roles = split_values(family.get("member_state_roles", ""))
    protomers = split_values(family.get("member_protomer_ids", ""))
    state_id, state_frame, selection_warnings = _select_membrane_state_frame(
        membrane_frame,
        states,
        roles,
    )
    warnings.extend(selection_warnings)
    if center is None:
        return _row(family, "missing", "missing", center, None, None, None, "not_evaluable", False, None, False, "pocket_centroid_missing", ["missing_pocket_center"])
    if state_frame is None:
        reason = "missing_membrane_frame"
        warnings.append(reason)
        return _row(family, "missing", "missing", center, None, None, None, "not_evaluable", False, None, False, reason, warnings)
    normal = _vector(state_frame.get("n_membrane"))
    if normal is None:
        reason = "membrane_normal_missing"
        warnings.append(reason)
        return _row(family, state_id, str(state_frame.get("status", "")), center, None, None, None, "not_evaluable", False, None, False, reason, warnings)
    origin = _vector(state_frame.get("p_jm_anchor")) or (0.0, 0.0, 0.0)
    normal = _normalize(normal)
    z_value = _dot(_sub(center, origin), normal)
    z_min, z_max = receptor_z_ranges.get(state_id, (z_value, z_value))
    z_percentile = 0.5 if z_max == z_min else max(0.0, min(1.0, (z_value - z_min) / (z_max - z_min)))
    lower_sign = _lower_direction_sign(membrane_frame, state_frame)
    if lower_sign is None:
        warnings.append("membrane_lower_direction_ambiguous")
        lower_pass = False
        lower_class = "not_evaluable"
    else:
        lower_pass = (z_value * lower_sign) >= thresholds["lower_z_minimum"]
        lower_class = "lower_or_membrane_proximal" if lower_pass else "upper_or_core"
    lateral_score = _best_lateral_score(center, state_frame, protomers, normal, warnings)
    if lateral_score is None:
        warnings.append("lateral_score_unavailable")
    lateral_pass = lateral_score is not None and lateral_score >= thresholds["lateral_score_pass_threshold"]
    lower_lateral_pass = bool(lower_pass and lateral_pass)
    reason = "lower_lateral_pass" if lower_lateral_pass else "lower_or_lateral_gate_failed"
    return _row(
        family,
        state_id,
        str(state_frame.get("status", "")),
        center,
        z_value,
        abs(z_value),
        z_percentile,
        lower_class,
        lower_pass,
        lateral_score,
        lateral_pass,
        reason,
        warnings,
    )


def _state_frame(membrane_frame: dict[str, Any] | None, state_id: str) -> dict[str, Any] | None:
    if not membrane_frame:
        return None
    states = membrane_frame.get("states", {})
    if isinstance(states, dict) and state_id in states:
        return states[state_id]
    return membrane_frame


def _select_membrane_state_frame(
    membrane_frame: dict[str, Any] | None,
    state_ids: list[str],
    state_roles: list[str],
) -> tuple[str, dict[str, Any] | None, list[str]]:
    if not state_ids:
        return "", _state_frame(membrane_frame, ""), []

    candidates: list[tuple[int, str, str, dict[str, Any] | None]] = []
    for index, state_id in enumerate(state_ids):
        role = state_roles[index] if index < len(state_roles) else ""
        candidates.append((index, state_id, role, _state_frame(membrane_frame, state_id)))

    selected: tuple[int, str, str, dict[str, Any] | None] | None = None
    for candidate in candidates:
        _index, state_id, role, frame = candidate
        if _frame_has_membrane_normal(frame) and _is_primary_membrane_frame(state_id, role, frame):
            selected = candidate
            break
    if selected is None:
        for candidate in candidates:
            if _frame_has_membrane_normal(candidate[3]):
                selected = candidate
                break
    if selected is None:
        for candidate in candidates:
            if candidate[3] is not None:
                selected = candidate
                break
    if selected is None:
        return state_ids[0], None, []

    first = candidates[0]
    warnings: list[str] = []
    if selected[0] != first[0]:
        if _is_reference_frame(first[1], first[2], first[3]):
            warnings.append("reference_state_skipped_for_membrane_frame:{0}".format(first[1]))
        else:
            warnings.append("non_membrane_state_skipped_for_membrane_frame:{0}".format(first[1]))
    return selected[1], selected[3], warnings


def _frame_has_membrane_normal(frame: dict[str, Any] | None) -> bool:
    return frame is not None and _vector(frame.get("n_membrane")) is not None


def _is_primary_membrane_frame(state_id: str, role: str, frame: dict[str, Any] | None) -> bool:
    frame_role = str((frame or {}).get("role", "")).lower()
    frame_status = str((frame or {}).get("status", "")).lower()
    role_text = str(role).lower()
    state_text = str(state_id).lower()
    if state_text == "3gt8_raw" or "reference" in frame_role or "reference" in frame_status:
        return False
    return "primary" in role_text or "primary" in frame_role or state_text.startswith("egfr_")


def _is_reference_frame(state_id: str, role: str, frame: dict[str, Any] | None) -> bool:
    frame_role = str((frame or {}).get("role", "")).lower()
    frame_status = str((frame or {}).get("status", "")).lower()
    role_text = str(role).lower()
    state_text = str(state_id).lower()
    return (
        "reference" in frame_role
        or "reference" in frame_status
        or state_text == "3gt8_raw"
        or (not frame_role and "reference" in role_text)
    )


def _vector(value: Any) -> tuple[float, float, float] | None:
    if not isinstance(value, list) or len(value) != 3:
        return None
    try:
        return (float(value[0]), float(value[1]), float(value[2]))
    except (TypeError, ValueError):
        return None


def _normalize(value: tuple[float, float, float]) -> tuple[float, float, float]:
    norm = math.sqrt(sum(item * item for item in value))
    if norm == 0.0:
        raise ValueError("cannot normalize zero membrane vector")
    return (value[0] / norm, value[1] / norm, value[2] / norm)


def _dot(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def _sub(left: tuple[float, float, float], right: tuple[float, float, float]) -> tuple[float, float, float]:
    return (left[0] - right[0], left[1] - right[1], left[2] - right[2])


def _project(value: tuple[float, float, float], normal: tuple[float, float, float]) -> tuple[float, float, float]:
    scale = _dot(value, normal)
    return (value[0] - scale * normal[0], value[1] - scale * normal[1], value[2] - scale * normal[2])


def _lower_direction_sign(membrane_frame: dict[str, Any] | None, state_frame: dict[str, Any]) -> float | None:
    for key in ("lower_direction_sign", "cytosolic_direction_sign"):
        value = safe_float(state_frame.get(key))
        if value is not None:
            return 1.0 if value >= 0 else -1.0
    convention = str((membrane_frame or {}).get("coordinate_convention", "")).lower()
    if "cytosolic" in convention or "intracellular" in convention:
        return 1.0
    return None


def _lateral_score(
    center: tuple[float, float, float],
    state_frame: dict[str, Any],
    protomer_id: str,
    normal: tuple[float, float, float],
) -> float | None:
    a = _vector(state_frame.get("protomer_a_centroid"))
    b = _vector(state_frame.get("protomer_b_centroid"))
    if a is None or b is None:
        return None
    protomer = a if protomer_id == "A" else b if protomer_id == "B" else None
    if protomer is None:
        return None
    dimer_center = ((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, (a[2] + b[2]) / 2.0)
    outward = _project(_sub(protomer, dimer_center), normal)
    pocket = _project(_sub(center, protomer), normal)
    try:
        outward_n = _normalize(outward)
        pocket_n = _normalize(pocket)
    except ValueError:
        return None
    return _dot(outward_n, pocket_n)


def _best_lateral_score(
    center: tuple[float, float, float],
    state_frame: dict[str, Any],
    protomer_ids: list[str],
    normal: tuple[float, float, float],
    warnings: list[str],
) -> float | None:
    candidates = [protomer for protomer in protomer_ids if protomer in {"A", "B"}]
    if not candidates and any(protomer.lower() in {"multi_protomer", "multi"} for protomer in protomer_ids):
        candidates = ["A", "B"]
        warnings.append("multi_protomer_lateral_score_evaluated_against_A_B")
    if not candidates:
        return None
    scores = [
        score
        for score in (_lateral_score(center, state_frame, protomer, normal) for protomer in candidates)
        if score is not None
    ]
    if not scores:
        return None
    if len(candidates) > 1:
        warnings.append("multi_protomer_lateral_score_max_of_A_B")
    return max(scores)


def _row(
    family: dict[str, Any],
    frame_id: str,
    frame_status: str,
    center: tuple[float, float, float] | None,
    z_value: float | None,
    z_proximity: float | None,
    z_percentile: float | None,
    lower_class: str,
    lower_pass: bool,
    lateral_score: float | None,
    lateral_pass: bool,
    reason: str,
    warnings: list[str],
) -> dict[str, Any]:
    return {
        "pocket_family_id": family.get("pocket_family_id", ""),
        "representative_candidate_pocket_id": family.get("representative_candidate_pocket_id", ""),
        "state_ids": family.get("member_state_ids", ""),
        "state_roles": family.get("member_state_roles", ""),
        "protomer_ids": family.get("member_protomer_ids", ""),
        "membrane_frame_id": frame_id,
        "membrane_frame_status": frame_status,
        "family_center_x": format_float(center[0] if center else None),
        "family_center_y": format_float(center[1] if center else None),
        "family_center_z": format_float(center[2] if center else None),
        "pocket_membrane_z": format_float(z_value),
        "z_membrane_proximity": format_float(z_proximity),
        "z_percentile_within_receptor_surface": format_float(z_percentile),
        "lower_gate_class": lower_class,
        "lower_gate_pass": lower_pass,
        "lateral_score": format_float(lateral_score),
        "lateral_gate_pass": lateral_pass,
        "lower_lateral_pass": bool(lower_pass and lateral_pass),
        "membrane_geometry_reason": reason,
        "warnings": ";".join(sorted(set(warnings))),
    }
