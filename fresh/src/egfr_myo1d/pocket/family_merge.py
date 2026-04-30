"""M2.6 deterministic raw pocket family merge.

Families are raw, ungated evidence groups. M2.7 is responsible for ATP, PPI,
membrane, and dimer gates.
"""

from __future__ import annotations

import math
from statistics import median
from typing import Any


MERGED_FIELDS = [
    "pocket_family_id",
    "representative_candidate_pocket_id",
    "representative_state_id",
    "representative_state_role",
    "representative_protomer_id",
    "member_candidate_pocket_ids",
    "member_state_ids",
    "member_state_roles",
    "member_protomer_ids",
    "primary_state_count",
    "reference_state_count",
    "raw_candidate_count",
    "family_center_x",
    "family_center_y",
    "family_center_z",
    "centroid_spread_angstrom",
    "residue_union_uniprot",
    "residue_intersection_uniprot",
    "max_pairwise_residue_jaccard",
    "median_pairwise_residue_jaccard",
    "merge_method",
    "merge_threshold_centroid_angstrom",
    "merge_threshold_residue_jaccard",
    "merge_status",
    "evidence_status",
    "warnings",
]

MEMBER_FIELDS = [
    "pocket_family_id",
    "candidate_pocket_id",
    "state_id",
    "state_role",
    "egfr_protomer_ids",
    "source_pocket_id",
    "membership_status",
]


def _split(value: Any) -> list[str]:
    return [item for item in str(value or "").split(";") if item]


def _float(value: Any) -> float | None:
    try:
        text = str(value).strip()
        if not text:
            return None
        return float(text)
    except (TypeError, ValueError):
        return None


def _format(value: float | None) -> str:
    return "" if value is None else "{0:.3f}".format(value)


def _centroid(row: dict[str, Any]) -> tuple[float, float, float] | None:
    x = _float(row.get("pocket_center_x"))
    y = _float(row.get("pocket_center_y"))
    z = _float(row.get("pocket_center_z"))
    if x is None or y is None or z is None:
        return None
    return (x, y, z)


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def _residues(row: dict[str, Any]) -> set[str]:
    return set(_split(row.get("mapped_uniprot_residues", "")))


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / float(len(a | b))


def _parent_find(parent: dict[str, str], item: str) -> str:
    while parent[item] != item:
        parent[item] = parent[parent[item]]
        item = parent[item]
    return item


def _parent_union(parent: dict[str, str], left: str, right: str) -> None:
    root_left = _parent_find(parent, left)
    root_right = _parent_find(parent, right)
    if root_left == root_right:
        return
    if root_left < root_right:
        parent[root_right] = root_left
    else:
        parent[root_left] = root_right


def merge_pocket_families(
    raw_candidates: list[dict[str, Any]],
    centroid_threshold_angstrom: float = 6.0,
    residue_jaccard_threshold: float = 0.30,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Merge raw candidates by centroid proximity or mapped residue overlap."""

    rows = sorted(raw_candidates, key=lambda row: str(row.get("candidate_pocket_id", "")))
    parent = {str(row.get("candidate_pocket_id", "")): str(row.get("candidate_pocket_id", "")) for row in rows}
    for i, left in enumerate(rows):
        left_id = str(left.get("candidate_pocket_id", ""))
        for right in rows[i + 1 :]:
            right_id = str(right.get("candidate_pocket_id", ""))
            left_centroid = _centroid(left)
            right_centroid = _centroid(right)
            centroid_match = (
                left_centroid is not None
                and right_centroid is not None
                and _distance(left_centroid, right_centroid) <= centroid_threshold_angstrom
            )
            residue_match = _jaccard(_residues(left), _residues(right)) >= residue_jaccard_threshold
            if centroid_match or residue_match:
                _parent_union(parent, left_id, right_id)

    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        row_id = str(row.get("candidate_pocket_id", ""))
        grouped.setdefault(_parent_find(parent, row_id), []).append(row)

    family_rows: list[dict[str, Any]] = []
    member_rows: list[dict[str, Any]] = []
    for idx, members in enumerate(
        sorted(grouped.values(), key=lambda group: [str(row.get("candidate_pocket_id", "")) for row in group]),
        start=1,
    ):
        family_id = "m2_6_pocket_family_{0:04d}".format(idx)
        representative = sorted(
            members,
            key=lambda row: (str(row.get("state_role", "")) != "primary", str(row.get("candidate_pocket_id", ""))),
        )[0]
        centroids = [centroid for centroid in (_centroid(row) for row in members) if centroid is not None]
        center = None
        spread = None
        if centroids:
            center = (
                sum(item[0] for item in centroids) / len(centroids),
                sum(item[1] for item in centroids) / len(centroids),
                sum(item[2] for item in centroids) / len(centroids),
            )
            spread = max(_distance(item, center) for item in centroids)
        residue_sets = [_residues(row) for row in members]
        union = set().union(*residue_sets) if residue_sets else set()
        intersection = set(residue_sets[0]) if residue_sets else set()
        for residue_set in residue_sets[1:]:
            intersection &= residue_set
        pairwise_jaccards = [
            _jaccard(residue_sets[i], residue_sets[j])
            for i in range(len(residue_sets))
            for j in range(i + 1, len(residue_sets))
        ]
        primary_states = {row.get("state_id", "") for row in members if row.get("state_role") == "primary"}
        reference_states = {row.get("state_id", "") for row in members if row.get("state_role") == "reference"}
        protomers = sorted(set().union(*(set(_split(row.get("egfr_protomer_ids", ""))) for row in members)))
        warnings = sorted(
            set().union(*(set(_split(row.get("warnings", ""))) for row in members))
        )
        if not primary_states and reference_states:
            warnings.append("reference_only_family_not_promotable_in_m2_6")
        if len(protomers) > 1:
            warnings.append("multi_protomer_family")
        family_rows.append(
            {
                "pocket_family_id": family_id,
                "representative_candidate_pocket_id": representative.get("candidate_pocket_id", ""),
                "representative_state_id": representative.get("state_id", ""),
                "representative_state_role": representative.get("state_role", ""),
                "representative_protomer_id": "multi_protomer" if len(protomers) > 1 else (protomers[0] if protomers else ""),
                "member_candidate_pocket_ids": ";".join(str(row.get("candidate_pocket_id", "")) for row in members),
                "member_state_ids": ";".join(sorted({str(row.get("state_id", "")) for row in members})),
                "member_state_roles": ";".join(sorted({str(row.get("state_role", "")) for row in members})),
                "member_protomer_ids": ";".join(protomers),
                "primary_state_count": len(primary_states),
                "reference_state_count": len(reference_states),
                "raw_candidate_count": len(members),
                "family_center_x": _format(center[0] if center else None),
                "family_center_y": _format(center[1] if center else None),
                "family_center_z": _format(center[2] if center else None),
                "centroid_spread_angstrom": _format(spread),
                "residue_union_uniprot": ";".join(sorted(union, key=lambda value: (not value.isdigit(), value))),
                "residue_intersection_uniprot": ";".join(sorted(intersection, key=lambda value: (not value.isdigit(), value))),
                "max_pairwise_residue_jaccard": _format(max(pairwise_jaccards) if pairwise_jaccards else 0.0),
                "median_pairwise_residue_jaccard": _format(median(pairwise_jaccards) if pairwise_jaccards else 0.0),
                "merge_method": "centroid_distance_or_residue_jaccard",
                "merge_threshold_centroid_angstrom": centroid_threshold_angstrom,
                "merge_threshold_residue_jaccard": residue_jaccard_threshold,
                "merge_status": "raw_family_ungated",
                "evidence_status": "REFERENCE_ONLY_RAW_UNGATED" if not primary_states and reference_states else "RAW_UNGATED",
                "warnings": ";".join(sorted(set(warnings))),
            }
        )
        for member in members:
            member_rows.append(
                {
                    "pocket_family_id": family_id,
                    "candidate_pocket_id": member.get("candidate_pocket_id", ""),
                    "state_id": member.get("state_id", ""),
                    "state_role": member.get("state_role", ""),
                    "egfr_protomer_ids": member.get("egfr_protomer_ids", ""),
                    "source_pocket_id": member.get("source_pocket_id", ""),
                    "membership_status": "raw_family_member_ungated",
                }
            )
    return family_rows, member_rows
