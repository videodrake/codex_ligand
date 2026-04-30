"""M2.7 PPI-pocket relationship classifier."""

from __future__ import annotations

import math
from typing import Any


PPI_RELATIONSHIP_FIELDS = [
    "pocket_family_id",
    "representative_candidate_pocket_id",
    "member_candidate_pocket_ids",
    "state_ids",
    "state_roles",
    "protomer_ids",
    "pocket_residue_count",
    "ppi_patch_residue_count",
    "residue_overlap_count",
    "residue_jaccard",
    "minimum_heavy_atom_distance_to_ppi",
    "minimum_ca_distance_to_ppi",
    "pocket_centroid_to_ppi_centroid_distance",
    "pocket_entrance_to_ppi_distance",
    "same_protomer_flag",
    "same_lower_lateral_neighborhood_flag",
    "ppi_relationship_class",
    "ppi_relationship_pass",
    "ppi_relationship_reason",
    "mapping_status",
    "warnings",
]


def split_values(value: Any) -> list[str]:
    return [item for item in str(value or "").split(";") if item]


def residue_set(value: Any) -> set[str]:
    return set(split_values(value))


def safe_float(value: Any) -> float | None:
    try:
        text = str(value).strip()
        if not text:
            return None
        return float(text)
    except (TypeError, ValueError):
        return None


def format_float(value: float | None) -> str:
    return "" if value is None else "{0:.3f}".format(value)


def point(row: dict[str, Any], x_key: str, y_key: str, z_key: str) -> tuple[float, float, float] | None:
    x = safe_float(row.get(x_key))
    y = safe_float(row.get(y_key))
    z = safe_float(row.get(z_key))
    if x is None or y is None or z is None:
        return None
    return (x, y, z)


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def mean_point(points: list[tuple[float, float, float]]) -> tuple[float, float, float] | None:
    if not points:
        return None
    return (
        sum(item[0] for item in points) / len(points),
        sum(item[1] for item in points) / len(points),
        sum(item[2] for item in points) / len(points),
    )


def jaccard(left: set[str], right: set[str]) -> float:
    if not left and not right:
        return 0.0
    return len(left & right) / float(len(left | right))


def classify_ppi_relationship(
    family: dict[str, Any],
    ppi_rows: list[dict[str, Any]],
    residue_geometry: dict[tuple[str, str, str], dict[str, Any]],
    thresholds: dict[str, float],
    same_lower_lateral_neighborhood_flag: bool = False,
    mapping_status: str = "PASS",
) -> dict[str, Any]:
    """Classify one pocket family against PPI patch evidence.

    The classifier intentionally refuses to promote a centroid-only relationship.
    It needs residue overlap or residue/atom distance evidence for `rim` and
    `allosteric_near`.
    """

    warnings: list[str] = []
    pocket_residues = residue_set(family.get("residue_union_uniprot", ""))
    pocket_states = set(split_values(family.get("member_state_ids", "")))
    pocket_protomers = set(split_values(family.get("member_protomer_ids", "")))
    matched_ppi_rows = _state_protomer_matched_ppi_rows(ppi_rows, pocket_states, pocket_protomers, warnings)
    ppi_patch_residues = {
        str(row.get("egfr_residue_number") or row.get("uniprot_residue_number") or "").strip()
        for row in matched_ppi_rows
        if str(row.get("egfr_residue_number") or row.get("uniprot_residue_number") or "").strip()
    }
    ppi_protomers = {
        str(row.get("egfr_protomer_id") or row.get("protomer_id") or "").strip()
        for row in matched_ppi_rows
        if str(row.get("egfr_protomer_id") or row.get("protomer_id") or "").strip()
    }
    overlap = pocket_residues & ppi_patch_residues
    residue_overlap_count = len(overlap)
    residue_jaccard = jaccard(pocket_residues, ppi_patch_residues)
    same_protomer = bool(pocket_protomers & ppi_protomers) if ppi_protomers else False

    pocket_point = point(family, "family_center_x", "family_center_y", "family_center_z")
    ppi_points = [
        point(row, "egfr_contact_centroid_x", "egfr_contact_centroid_y", "egfr_contact_centroid_z")
        for row in matched_ppi_rows
    ]
    ppi_centroid = mean_point([item for item in ppi_points if item is not None])
    centroid_distance = distance(pocket_point, ppi_centroid) if pocket_point and ppi_centroid else None

    min_ca = _minimum_residue_distance(
        pocket_states,
        pocket_protomers,
        pocket_residues,
        matched_ppi_rows,
        residue_geometry,
        atom_key="ca",
    )
    min_heavy = _minimum_residue_distance(
        pocket_states,
        pocket_protomers,
        pocket_residues,
        matched_ppi_rows,
        residue_geometry,
        atom_key="heavy",
    )

    if mapping_status.lower() in {"missing", "no_residues", "fail", "mapping_reject"}:
        relationship_class = "mapping_ambiguous"
        reason = "mapping_failure_prevents_ppi_relationship"
    elif not ppi_rows:
        relationship_class = "not_evaluable"
        reason = "missing_ppi_consensus_patch"
    elif not matched_ppi_rows:
        relationship_class = "not_evaluable"
        reason = "no_state_protomer_matched_ppi_patch_evidence"
    elif not pocket_residues:
        relationship_class = "mapping_ambiguous"
        reason = "pocket_residues_unavailable"
    elif residue_overlap_count >= 1:
        relationship_class = "orthosteric"
        reason = "pocket_residue_overlaps_ppi_patch"
    elif _leq(min_heavy, thresholds["ppi_rim_heavy_atom_cutoff_angstrom"]) or _leq(
        min_ca, thresholds["ppi_rim_ca_cutoff_angstrom"]
    ):
        relationship_class = "rim"
        reason = "pocket_residue_distance_near_ppi_patch_boundary"
    elif same_lower_lateral_neighborhood_flag and (
        _leq(min_heavy, thresholds["ppi_allosteric_heavy_atom_cutoff_angstrom"])
        or _leq(min_ca, thresholds["ppi_allosteric_ca_cutoff_angstrom"])
    ):
        relationship_class = "allosteric_near"
        reason = "nearby_residue_distance_in_same_lower_lateral_neighborhood"
    elif same_lower_lateral_neighborhood_flag and _leq(
        centroid_distance, thresholds["ppi_allosteric_centroid_cutoff_angstrom"]
    ):
        relationship_class = "low_relevance"
        reason = "centroid_only_relationship_not_allowed"
        warnings.append("centroid_distance_not_sufficient_for_ppi_relationship")
    else:
        relationship_class = "low_relevance"
        reason = "no_ppi_residue_overlap_or_near_distance"

    relationship_pass = relationship_class in {"orthosteric", "rim", "allosteric_near"}
    return {
        "pocket_family_id": family.get("pocket_family_id", ""),
        "representative_candidate_pocket_id": family.get("representative_candidate_pocket_id", ""),
        "member_candidate_pocket_ids": family.get("member_candidate_pocket_ids", ""),
        "state_ids": family.get("member_state_ids", ""),
        "state_roles": family.get("member_state_roles", ""),
        "protomer_ids": family.get("member_protomer_ids", ""),
        "pocket_residue_count": len(pocket_residues),
        "ppi_patch_residue_count": len(ppi_patch_residues),
        "residue_overlap_count": residue_overlap_count,
        "residue_jaccard": format_float(residue_jaccard),
        "minimum_heavy_atom_distance_to_ppi": format_float(min_heavy),
        "minimum_ca_distance_to_ppi": format_float(min_ca),
        "pocket_centroid_to_ppi_centroid_distance": format_float(centroid_distance),
        "pocket_entrance_to_ppi_distance": "",
        "same_protomer_flag": same_protomer,
        "same_lower_lateral_neighborhood_flag": same_lower_lateral_neighborhood_flag,
        "ppi_relationship_class": relationship_class,
        "ppi_relationship_pass": relationship_pass,
        "ppi_relationship_reason": reason,
        "mapping_status": mapping_status,
        "warnings": ";".join(sorted(set(warnings))),
    }


def _leq(value: float | None, threshold: float) -> bool:
    return value is not None and value <= threshold


def _state_protomer_matched_ppi_rows(
    ppi_rows: list[dict[str, Any]],
    pocket_states: set[str],
    pocket_protomers: set[str],
    warnings: list[str],
) -> list[dict[str, Any]]:
    matched: list[dict[str, Any]] = []
    missing_state = 0
    missing_protomer = 0
    for row in ppi_rows:
        ppi_state = str(row.get("state_id", "")).strip()
        ppi_protomer = str(row.get("egfr_protomer_id") or row.get("protomer_id") or "").strip()
        if pocket_states:
            if not ppi_state:
                missing_state += 1
                continue
            if ppi_state not in pocket_states:
                continue
        if pocket_protomers:
            if not ppi_protomer:
                missing_protomer += 1
                continue
            if ppi_protomer not in pocket_protomers:
                continue
        matched.append(row)
    if missing_state:
        warnings.append("ppi_rows_missing_state_id_excluded:{0}".format(missing_state))
    if missing_protomer:
        warnings.append("ppi_rows_missing_protomer_id_excluded:{0}".format(missing_protomer))
    return matched


def _minimum_residue_distance(
    pocket_states: set[str],
    pocket_protomers: set[str],
    pocket_residues: set[str],
    ppi_rows: list[dict[str, Any]],
    residue_geometry: dict[tuple[str, str, str], dict[str, Any]],
    atom_key: str,
) -> float | None:
    distances: list[float] = []
    for state_id in pocket_states:
        for protomer_id in pocket_protomers:
            for pocket_residue in pocket_residues:
                pocket_geom = residue_geometry.get((state_id, protomer_id, pocket_residue))
                if not pocket_geom:
                    continue
                pocket_points = pocket_geom.get(atom_key, [])
                for ppi in ppi_rows:
                    ppi_state = str(ppi.get("state_id", "")).strip()
                    ppi_protomer = str(ppi.get("egfr_protomer_id") or ppi.get("protomer_id") or "").strip()
                    ppi_residue = str(ppi.get("egfr_residue_number") or ppi.get("uniprot_residue_number") or "").strip()
                    if ppi_state and ppi_state != state_id:
                        continue
                    if ppi_protomer and ppi_protomer != protomer_id:
                        continue
                    ppi_geom = residue_geometry.get((state_id, protomer_id, ppi_residue))
                    if not ppi_geom:
                        continue
                    for left in pocket_points:
                        for right in ppi_geom.get(atom_key, []):
                            distances.append(distance(left, right))
    return min(distances) if distances else None
