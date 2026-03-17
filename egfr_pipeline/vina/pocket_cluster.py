#!/usr/bin/env python3
"""Pocket clustering for Vina docking poses.

Iterative centroid-based clustering (k-means style) replaces the original
single-pass greedy algorithm to eliminate centroid drift and order dependence.

Algorithm:
  1. Seed phase: scan poses in affinity order; create a new pocket seed when
     no existing seed is within *cutoff*.  Seeds are **fixed** during this scan.
  2. Assign phase: assign every pose to its nearest pocket centroid.
  3. Update phase: recompute centroids from assigned members.
  4. Repeat 2-3 until convergence or *max_iterations*.
  5. Absorb singletons: pockets with fewer than *min_pocket_size* poses are
     merged into their nearest neighbour.

Post-hoc residue-based merging (Union-Find) is available as an additional
safety net against fragmentation.
"""
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from egfr_pipeline.config import load_config
from egfr_pipeline import paths
from egfr_pipeline.residue_utils import parse_residue_set


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def euclidean_distance(a: List[float], b: List[float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def row_centroid(row: dict) -> List[float]:
    return [float(row["centroid_x"]), float(row["centroid_y"]), float(row["centroid_z"])]


def _compute_centroid(points: List[List[float]]) -> List[float]:
    n = len(points)
    if n == 0:
        return [0.0, 0.0, 0.0]
    return [sum(p[axis] for p in points) / n for axis in range(3)]


# ---------------------------------------------------------------------------
# Sort helper
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Core: iterative centroid-based clustering
# ---------------------------------------------------------------------------

def _seed_pockets(
    centroids: List[List[float]],
    cutoff: float,
) -> List[List[float]]:
    """Create pocket seeds by scanning centroids; no seed movement."""
    seeds: List[List[float]] = []
    for centroid in centroids:
        found = False
        for seed in seeds:
            if euclidean_distance(centroid, seed) <= cutoff:
                found = True
                break
        if not found:
            seeds.append(centroid[:])
    return seeds


def _assign_to_nearest(
    centroids: List[List[float]],
    seeds: List[List[float]],
) -> List[int]:
    """Assign each centroid to its nearest seed index."""
    assignments: List[int] = []
    for centroid in centroids:
        best_idx = 0
        best_dist = euclidean_distance(centroid, seeds[0])
        for idx in range(1, len(seeds)):
            dist = euclidean_distance(centroid, seeds[idx])
            if dist < best_dist:
                best_dist = dist
                best_idx = idx
        assignments.append(best_idx)
    return assignments


def _recompute_seeds(
    centroids: List[List[float]],
    assignments: List[int],
    n_seeds: int,
    old_seeds: List[List[float]],
) -> List[List[float]]:
    """Recompute seed positions as mean of assigned centroids."""
    buckets: Dict[int, List[List[float]]] = defaultdict(list)
    for centroid, idx in zip(centroids, assignments):
        buckets[idx].append(centroid)
    new_seeds: List[List[float]] = []
    for i in range(n_seeds):
        members = buckets.get(i, [])
        if members:
            new_seeds.append(_compute_centroid(members))
        else:
            # Empty seed — retain previous position so it doesn't
            # attract poses to the origin
            new_seeds.append(old_seeds[i][:])
    return new_seeds


def assign_pockets(
    rows: List[dict],
    cutoff: float = 8.0,
    max_iterations: int = 10,
    min_pocket_size: int = 2,
) -> List[dict]:
    """Iterative centroid-based pocket assignment (per receptor).

    Args:
        rows: pose table rows (modified in place — pocket_id field set).
        cutoff: distance threshold for seed creation (Angstrom).
        max_iterations: maximum refinement iterations.
        min_pocket_size: pockets smaller than this are absorbed into the
            nearest neighbour.

    Returns:
        The rows list (sorted by receptor/affinity), with pocket_id assigned.
    """
    sorted_rows = sort_pose_rows(rows)

    # Group rows by receptor
    receptor_rows: Dict[str, List[dict]] = defaultdict(list)
    for row in sorted_rows:
        receptor_rows[row["receptor_id"]].append(row)

    for receptor_id, r_rows in receptor_rows.items():
        centroids = [row_centroid(r) for r in r_rows]
        if not centroids:
            continue

        # Step A: seed creation (fixed centroids, no drift)
        seeds = _seed_pockets(centroids, cutoff)
        if not seeds:
            continue

        # Step B: iterative assign → recompute until stable
        assignments = _assign_to_nearest(centroids, seeds)
        for _ in range(max_iterations):
            new_seeds = _recompute_seeds(centroids, assignments, len(seeds), seeds)
            new_assignments = _assign_to_nearest(centroids, new_seeds)
            if new_assignments == assignments:
                seeds = new_seeds
                break
            seeds = new_seeds
            assignments = new_assignments

        # Step C: absorb small pockets
        if min_pocket_size >= 2:
            # Count members per seed
            seed_counts: Dict[int, int] = defaultdict(int)
            for idx in assignments:
                seed_counts[idx] += 1

            small_seeds = {idx for idx, cnt in seed_counts.items() if cnt < min_pocket_size}
            large_seeds = [idx for idx in range(len(seeds)) if idx not in small_seeds]

            if large_seeds and small_seeds:
                for i, idx in enumerate(assignments):
                    if idx in small_seeds:
                        # Reassign to nearest large pocket
                        best_large = large_seeds[0]
                        best_dist = euclidean_distance(centroids[i], seeds[large_seeds[0]])
                        for lidx in large_seeds[1:]:
                            dist = euclidean_distance(centroids[i], seeds[lidx])
                            if dist < best_dist:
                                best_dist = dist
                                best_large = lidx
                        assignments[i] = best_large

        # Build pocket IDs: renumber contiguously sorted by first appearance
        seen_order: List[int] = []
        for idx in assignments:
            if idx not in seen_order:
                seen_order.append(idx)
        seed_to_pocket: Dict[int, str] = {
            seed_idx: f"P{rank + 1:03d}" for rank, seed_idx in enumerate(seen_order)
        }

        for row, idx in zip(r_rows, assignments):
            row["pocket_id"] = seed_to_pocket[idx]

    return sorted_rows


# ---------------------------------------------------------------------------
# Union-Find helpers
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Residue-based pocket merge
# ---------------------------------------------------------------------------

def merge_pockets_by_residue(
    rows: List[dict],
    jaccard_threshold: float = 0.3,
    overlap_threshold: float = 0.5,
    centroid_fallback_cutoff: float = 6.0,
) -> Tuple[List[dict], List[dict]]:
    """Post-hoc merge of pockets within each receptor based on residue overlap.

    Two pockets are merged if their contact residue sets satisfy:
      jaccard >= jaccard_threshold  OR  overlap_coeff >= overlap_threshold

    Additionally, pockets whose centroids are within *centroid_fallback_cutoff*
    are merged even when residue data is empty or too sparse, preventing
    fragments with few contacts from escaping the merge.

    Transitive closure via Union-Find ensures A-B and B-C merges also merge A-C.
    The merged pocket keeps the ID of the pocket with more poses.

    Returns:
        A tuple of (rows, merge_log) where *merge_log* is a list of dicts
        recording every pairwise merge decision (one entry per merged pair).
    """
    merge_log: List[dict] = []

    if not rows:
        return rows, merge_log

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

        # Collect residues and centroids per pocket
        pocket_residues: Dict[str, Set[str]] = {}
        pocket_centroids: Dict[str, List[float]] = {}
        for pid in pocket_ids:
            residues: Set[str] = set()
            points: List[List[float]] = []
            for row in pockets[pid]:
                raw = row.get("contact_residues", "")
                if raw:
                    residues |= parse_residue_set(raw)
                points.append(row_centroid(row))
            pocket_residues[pid] = residues
            pocket_centroids[pid] = _compute_centroid(points) if points else [0.0, 0.0, 0.0]

        # Union-Find
        parent = {pid: pid for pid in pocket_ids}
        uf_rank = {pid: 0 for pid in pocket_ids}

        for i, pid_a in enumerate(pocket_ids):
            for pid_b in pocket_ids[i + 1:]:
                res_a, res_b = pocket_residues[pid_a], pocket_residues[pid_b]

                should_merge = False
                reasons: List[str] = []

                # Compute all metrics for logging
                intersection = 0
                union_size = 0
                j = 0.0
                oc = 0.0
                if res_a and res_b:
                    shared = res_a & res_b
                    intersection = len(shared)
                    union_size = len(res_a | res_b)
                    j = intersection / union_size if union_size else 0.0
                    min_size = min(len(res_a), len(res_b))
                    oc = intersection / min_size if min_size else 0.0
                    if j >= jaccard_threshold:
                        should_merge = True
                        reasons.append("jaccard")
                    if oc >= overlap_threshold:
                        should_merge = True
                        reasons.append("overlap")

                dist = euclidean_distance(
                    pocket_centroids[pid_a], pocket_centroids[pid_b]
                )

                # Centroid fallback: merge if centroids are very close,
                # even when residue data is sparse or absent
                if not should_merge and centroid_fallback_cutoff > 0:
                    if dist <= centroid_fallback_cutoff:
                        should_merge = True
                        reasons.append("centroid")

                if should_merge:
                    _union_find_merge(parent, uf_rank, pid_a, pid_b)
                    merge_log.append({
                        "receptor_id": receptor_id,
                        "pocket_id_a": pid_a,
                        "pocket_id_b": pid_b,
                        "jaccard": round(j, 4),
                        "overlap_coeff": round(oc, 4),
                        "centroid_dist_A": round(dist, 2),
                        "n_residues_a": len(res_a),
                        "n_residues_b": len(res_b),
                        "n_shared_residues": intersection,
                        "merge_reason": "+".join(reasons),
                    })

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

    return rows, merge_log


# ---------------------------------------------------------------------------
# Merge log I/O
# ---------------------------------------------------------------------------

_MERGE_LOG_COLUMNS = [
    "receptor_id",
    "pocket_id_a",
    "pocket_id_b",
    "jaccard",
    "overlap_coeff",
    "centroid_dist_A",
    "n_residues_a",
    "n_residues_b",
    "n_shared_residues",
    "merge_reason",
]


# ---------------------------------------------------------------------------
# Pocket cap: limit poses per pocket with round-robin ligand selection
# ---------------------------------------------------------------------------

_CAP_REPORT_COLUMNS = [
    "receptor_id",
    "pocket_id",
    "total_poses",
    "kept",
    "capped",
    "n_ligands_kept",
    "best_affinity",
    "worst_kept_affinity",
]


def apply_pocket_cap(
    rows: List[dict],
    max_per_pocket: int = 5,
) -> Tuple[List[dict], List[dict]]:
    """Limit poses per pocket using round-robin ligand selection by affinity.

    Within each (receptor_id, pocket_id) group, selects up to *max_per_pocket*
    poses by cycling through ligands in best-affinity order so that each ligand
    is represented at least once (when possible).

    Returns:
        (rows_with_cap_status, cap_report_rows)
    """
    # Group by (receptor_id, pocket_id)
    groups: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    for i, row in enumerate(rows):
        key = (row["receptor_id"], row.get("pocket_id", ""))
        groups[key].append(i)

    cap_report: List[dict] = []

    for (receptor_id, pocket_id), indices in groups.items():
        # Sort each group by affinity (ascending = best first)
        indices.sort(key=lambda i: float(rows[i]["affinity"])
                     if rows[i].get("affinity") not in ("", None)
                     else float("inf"))

        # Bucket by ligand_id, preserving affinity order within each bucket
        ligand_buckets: Dict[str, List[int]] = defaultdict(list)
        for i in indices:
            ligand_buckets[rows[i]["ligand_id"]].append(i)

        # Round-robin selection across ligands (sorted by best affinity per ligand)
        ligand_order = sorted(
            ligand_buckets.keys(),
            key=lambda lid: float(rows[ligand_buckets[lid][0]]["affinity"])
            if rows[ligand_buckets[lid][0]].get("affinity") not in ("", None)
            else float("inf"),
        )
        # Track cursor per ligand
        cursors = {lid: 0 for lid in ligand_order}
        kept_indices: Set[int] = set()

        while len(kept_indices) < max_per_pocket:
            added = False
            for lid in ligand_order:
                if len(kept_indices) >= max_per_pocket:
                    break
                if cursors[lid] < len(ligand_buckets[lid]):
                    kept_indices.add(ligand_buckets[lid][cursors[lid]])
                    cursors[lid] += 1
                    added = True
            if not added:
                break

        # Mark cap_status
        kept_affinities = []
        for i in indices:
            if i in kept_indices:
                rows[i]["cap_status"] = "kept"
                aff = rows[i].get("affinity", "")
                if aff not in ("", None):
                    kept_affinities.append(float(aff))
            else:
                rows[i]["cap_status"] = "capped"

        # Build report row
        ligands_kept = {rows[i]["ligand_id"] for i in kept_indices}
        cap_report.append({
            "receptor_id": receptor_id,
            "pocket_id": pocket_id,
            "total_poses": len(indices),
            "kept": len(kept_indices),
            "capped": len(indices) - len(kept_indices),
            "n_ligands_kept": len(ligands_kept),
            "best_affinity": round(min(kept_affinities), 4) if kept_affinities else "",
            "worst_kept_affinity": round(max(kept_affinities), 4) if kept_affinities else "",
        })

    return rows, cap_report


def write_cap_report(path: Path, rows: List[dict]) -> Path:
    """Write pocket cap report CSV."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_CAP_REPORT_COLUMNS)
        writer.writeheader()
        for entry in rows:
            writer.writerow(entry)
    return path


def write_merge_log(path: Path, merge_log: List[dict]) -> Path:
    """Write merge decision log CSV.  Writes header even when *merge_log* is empty."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_MERGE_LOG_COLUMNS)
        writer.writeheader()
        for entry in merge_log:
            writer.writerow(entry)
    return path


def write_merge_parameters(path: Path, params: dict) -> Path:
    """Write clustering/merge parameters to a JSON file for reproducibility."""
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(params, handle, indent=2)
    return path


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def cluster_pose_table(
    config_path: str,
    pose_table_path: Optional[str] = None,
    cutoff: float = 8.0,
    max_iterations: int = 10,
    min_pocket_size: int = 2,
    merge_by_residue: bool = True,
    merge_jaccard: float = 0.3,
    merge_overlap: float = 0.5,
    merge_centroid_fallback: float = 6.0,
    max_per_pocket: int = 0,
) -> Path:
    config = load_config(config_path)
    postprocess_root = paths.wa_phase4_vina_postprocess(config)
    target = Path(pose_table_path) if pose_table_path else postprocess_root / "vina_pose_table.csv"
    rows = load_pose_table(target)
    clustered_rows = assign_pockets(
        rows, cutoff=cutoff,
        max_iterations=max_iterations,
        min_pocket_size=min_pocket_size,
    )

    merge_log: List[dict] = []
    if merge_by_residue:
        clustered_rows, merge_log = merge_pockets_by_residue(
            clustered_rows,
            jaccard_threshold=merge_jaccard,
            overlap_threshold=merge_overlap,
            centroid_fallback_cutoff=merge_centroid_fallback,
        )

    # Write merge decision log (always, even if empty)
    merge_log_path = target.parent / "vina_clustering_merge_log.csv"
    write_merge_log(merge_log_path, merge_log)

    # Apply pocket cap if configured
    if max_per_pocket > 0:
        clustered_rows, cap_report = apply_pocket_cap(clustered_rows, max_per_pocket)
        cap_report_path = target.parent / "vina_pocket_cap_report.csv"
        write_cap_report(cap_report_path, cap_report)
        n_kept = sum(1 for r in clustered_rows if r.get("cap_status") == "kept")
        print(f"  Pocket cap: {n_kept} kept / {len(clustered_rows) - n_kept} capped")

    # Write clustering parameters for reproducibility
    params = {
        "clustering": {
            "cutoff": cutoff,
            "max_iterations": max_iterations,
            "min_pocket_size": min_pocket_size,
        },
        "post_hoc_merge": {
            "enabled": merge_by_residue,
            "jaccard_threshold": merge_jaccard,
            "overlap_threshold": merge_overlap,
            "centroid_fallback_cutoff": merge_centroid_fallback,
        },
        "pocket_cap": {
            "max_per_pocket": max_per_pocket,
            "selection_method": "round_robin_by_ligand",
        },
    }
    params_path = target.parent / "vina_clustering_parameters.json"
    write_merge_parameters(params_path, params)

    return write_pose_table(target, clustered_rows)
