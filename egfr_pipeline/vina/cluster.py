#!/usr/bin/env python3
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set

from egfr_pipeline.config import load_config, project_root_from_config
from egfr_pipeline.residue_utils import parse_residue_set


def load_pose_table(path: Path) -> List[dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_pose_table(path: Path, rows: List[dict]) -> Path:
    if not rows:
        return path
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return path


def euclidean_distance(a: List[float], b: List[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def row_centroid(row: dict) -> List[float]:
    return [float(row["centroid_x"]), float(row["centroid_y"]), float(row["centroid_z"])]


def sort_pose_rows(rows: List[dict]) -> List[dict]:
    def sort_key(row: dict):
        affinity = float(row["affinity"]) if row["affinity"] not in ("", None) else float("inf")
        return (
            row["receptor_id"],
            affinity,
            row["ligand_id"],
            int(row["pose_rank"]),
            row["raw_pose_file"],
        )
    return sorted(rows, key=sort_key)


def assign_pockets(rows: List[dict], cutoff: float) -> List[dict]:
    sorted_rows = sort_pose_rows(rows)
    pocket_state: Dict[str, List[dict]] = {}

    for row in sorted_rows:
        receptor_id = row["receptor_id"]
        centroid = row_centroid(row)
        receptor_pockets = pocket_state.setdefault(receptor_id, [])

        selected = None
        min_distance = None
        for pocket in receptor_pockets:
            dist = euclidean_distance(centroid, pocket["center"])
            if dist <= cutoff and (min_distance is None or dist < min_distance):
                selected = pocket
                min_distance = dist

        if selected is None:
            pocket_index = len(receptor_pockets) + 1
            selected = {
                "id": f"P{pocket_index:03d}",
                "center": centroid[:],
                "count": 0,
            }
            receptor_pockets.append(selected)

        selected["count"] += 1
        count = selected["count"]
        selected["center"] = [
            ((selected["center"][axis] * (count - 1)) + centroid[axis]) / count
            for axis in range(3)
        ]
        row["pocket_id"] = selected["id"]

    return sorted_rows


def _union_find_root(parent: Dict[str, str], x: str) -> str:
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def _union_find_merge(parent: Dict[str, str], rank: Dict[str, int], a: str, b: str):
    ra, rb = _union_find_root(parent, a), _union_find_root(parent, b)
    if ra == rb:
        return
    if rank[ra] < rank[rb]:
        ra, rb = rb, ra
    parent[rb] = ra
    if rank[ra] == rank[rb]:
        rank[ra] += 1


def merge_pockets_by_residue(
    rows: List[dict],
    jaccard_threshold: float = 0.3,
    overlap_threshold: float = 0.5,
) -> List[dict]:
    """Post-hoc merge of pockets within each receptor based on residue overlap.

    Two pockets are merged if their contact residue sets satisfy:
      jaccard >= jaccard_threshold  OR  overlap_coeff >= overlap_threshold

    Transitive closure via Union-Find ensures A-B and B-C merges also merge A-C.
    The merged pocket keeps the ID of the pocket with more poses.
    """
    if not rows:
        return rows

    # Group rows by receptor
    by_receptor: Dict[str, Dict[str, List[dict]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        pid = row.get("pocket_id", "")
        if pid:
            by_receptor[row["receptor_id"]][pid].append(row)

    for receptor_id, pockets in by_receptor.items():
        pocket_ids = sorted(pockets.keys())
        if len(pocket_ids) < 2:
            continue

        # Collect residues per pocket
        pocket_residues: Dict[str, Set[str]] = {}
        for pid in pocket_ids:
            residues: Set[str] = set()
            for row in pockets[pid]:
                raw = row.get("contact_residues", "")
                if raw:
                    residues |= parse_residue_set(raw)
            pocket_residues[pid] = residues

        # Union-Find
        parent = {pid: pid for pid in pocket_ids}
        uf_rank = {pid: 0 for pid in pocket_ids}

        for i, pid_a in enumerate(pocket_ids):
            for pid_b in pocket_ids[i + 1:]:
                res_a, res_b = pocket_residues[pid_a], pocket_residues[pid_b]
                if not res_a or not res_b:
                    continue
                intersection = len(res_a & res_b)
                union = len(res_a | res_b)
                j = intersection / union if union else 0.0
                oc = intersection / min(len(res_a), len(res_b))
                if j >= jaccard_threshold or oc >= overlap_threshold:
                    _union_find_merge(parent, uf_rank, pid_a, pid_b)

        # Build merge groups
        groups: Dict[str, List[str]] = defaultdict(list)
        for pid in pocket_ids:
            groups[_union_find_root(parent, pid)].append(pid)

        # For each group, pick the pocket with most poses as canonical ID
        remap: Dict[str, str] = {}
        for root, members in groups.items():
            if len(members) == 1:
                continue
            canonical = max(members, key=lambda p: len(pockets[p]))
            for pid in members:
                if pid != canonical:
                    remap[pid] = canonical

        # Apply remapping
        if remap:
            for row in rows:
                if row["receptor_id"] == receptor_id and row.get("pocket_id", "") in remap:
                    row["pocket_id"] = remap[row["pocket_id"]]

    return rows


def cluster_pose_table(
    config_path: str,
    pose_table_path: Optional[str] = None,
    cutoff: float = 4.0,
    merge_by_residue: bool = False,
    merge_jaccard: float = 0.3,
    merge_overlap: float = 0.5,
) -> Path:
    config = load_config(config_path)
    project_root = project_root_from_config(config)
    target = Path(pose_table_path) if pose_table_path else project_root / "vina_pose_table.csv"
    rows = load_pose_table(target)
    clustered_rows = assign_pockets(rows, cutoff)
    if merge_by_residue:
        clustered_rows = merge_pockets_by_residue(
            clustered_rows,
            jaccard_threshold=merge_jaccard,
            overlap_threshold=merge_overlap,
        )
    return write_pose_table(target, clustered_rows)
