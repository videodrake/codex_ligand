#!/usr/bin/env python3
"""Phase 2 Task Group 2.2: Candidate Pocket Normalization and Merge Logic.

Merges closely overlapping candidate pockets within each receptor state so that
downstream docking budget is not wasted on duplicate pockets proposed by
different tools.

Tasks:
  2.2.1 - Merge-threshold definition (centroid distance + residue Jaccard)
  2.2.2 - Receptor-local merge operation (never merge across states)
  2.2.3 - Pocket provenance preservation (which raw pockets → merged pocket)
  2.2.4 - Stable candidate pocket IDs (receptor-local, deterministic)

Inputs:
  - output/phase2_pockets/candidate_pockets_raw.csv

Outputs:
  - output/phase2_pockets/candidate_pockets.csv          (merged/final pockets)
  - output/phase2_pockets/candidate_pocket_merge_table.csv (merge decisions)
  - output/phase2_pockets/candidate_pocket_provenance.csv  (raw → merged mapping)

Usage:
    python -m egfr_pipeline.phase2.pocket_merge [--output_dir ...]
"""

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from egfr_pipeline.csv_utils import load_csv
from egfr_pipeline.parsing_utils import safe_float

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE2_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase2_pockets"
RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]

# Merge thresholds
CENTROID_MERGE_DISTANCE_A = 6.0   # Å — pockets with centroids closer merge
RESIDUE_JACCARD_MERGE = 0.30      # Jaccard overlap ≥ 0.30 also triggers merge

# Merged pocket schema
MERGED_POCKET_COLUMNS = [
    "receptor_id",
    "pocket_id",                  # stable ID: {state}_PKT{nn}
    "n_sources",                  # how many raw proposals merged
    "sources",                    # semicolon-delimited source tools
    "raw_pocket_ids",             # semicolon-delimited raw IDs
    "centroid_x",
    "centroid_y",
    "centroid_z",
    "best_proposal_score",
    "best_druggability_score",
    "mean_volume_A3",
    "residue_ids",                # union of all merged pocket residues
    "residue_nums",
    "n_residues",
    "box_center_x",
    "box_center_y",
    "box_center_z",
    "box_size_x",
    "box_size_y",
    "box_size_z",
]

MERGE_TABLE_COLUMNS = [
    "receptor_id",
    "pocket_id_a",
    "pocket_id_b",
    "centroid_dist_A",
    "residue_jaccard",
    "merge_decision",             # merged / kept_separate
    "merge_reason",
]

PROVENANCE_COLUMNS = [
    "receptor_id",
    "merged_pocket_id",
    "raw_pocket_id",
    "proposal_source",
    "pocket_rank",
    "proposal_score",
    "druggability_score",
]


# ---------------------------------------------------------------------------
# Merge logic (Tasks 2.2.1 + 2.2.2)
# ---------------------------------------------------------------------------

def merge_pockets(
    raw_pockets: List[dict],
    centroid_threshold: float = CENTROID_MERGE_DISTANCE_A,
    jaccard_threshold: float = RESIDUE_JACCARD_MERGE,
) -> Tuple[List[dict], List[dict], List[dict]]:
    """Merge overlapping candidate pockets within each receptor state.

    Returns (merged_pockets, merge_table, provenance).
    """
    # Group by receptor state
    by_state: Dict[str, List[dict]] = defaultdict(list)
    for p in raw_pockets:
        by_state[p["receptor_id"]].append(p)

    all_merged = []
    all_merge_table = []
    all_provenance = []

    for state in RECEPTOR_STATES:
        pockets = by_state.get(state, [])
        if not pockets:
            continue

        merged, merge_rows, prov_rows = _merge_state_pockets(
            state, pockets, centroid_threshold, jaccard_threshold
        )
        all_merged.extend(merged)
        all_merge_table.extend(merge_rows)
        all_provenance.extend(prov_rows)

    return all_merged, all_merge_table, all_provenance


def _merge_state_pockets(
    state: str,
    pockets: List[dict],
    centroid_threshold: float,
    jaccard_threshold: float,
) -> Tuple[List[dict], List[dict], List[dict]]:
    """Merge pockets within one receptor state using union-find."""
    n = len(pockets)
    merge_table = []

    # Union-Find
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    # Pairwise comparison
    for i in range(n):
        for j in range(i + 1, n):
            dist = _centroid_distance(pockets[i], pockets[j])
            jaccard = _residue_jaccard(pockets[i], pockets[j])

            should_merge = dist <= centroid_threshold or jaccard >= jaccard_threshold

            reason_parts = []
            if dist <= centroid_threshold:
                reason_parts.append(f"centroid_dist={dist:.1f}A <= {centroid_threshold}A")
            if jaccard >= jaccard_threshold:
                reason_parts.append(f"jaccard={jaccard:.3f} >= {jaccard_threshold}")

            merge_table.append({
                "receptor_id": state,
                "pocket_id_a": pockets[i]["candidate_pocket_id"],
                "pocket_id_b": pockets[j]["candidate_pocket_id"],
                "centroid_dist_A": f"{dist:.2f}",
                "residue_jaccard": f"{jaccard:.3f}",
                "merge_decision": "merged" if should_merge else "kept_separate",
                "merge_reason": "; ".join(reason_parts) if reason_parts else "no_overlap",
            })

            if should_merge:
                union(i, j)

    # Group by merged cluster
    groups: Dict[int, List[int]] = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)

    # Build merged pockets (Task 2.2.4: stable IDs)
    merged_pockets = []
    provenance = []

    # Sort groups by best score for stable ranking
    sorted_groups = sorted(
        groups.values(),
        key=lambda indices: _best_score(pockets, indices),
        reverse=True,
    )

    for rank, indices in enumerate(sorted_groups, start=1):
        pocket_id = f"{state}_PKT{rank:02d}"
        group_pockets = [pockets[i] for i in indices]

        merged = _build_merged_pocket(state, pocket_id, group_pockets)
        merged_pockets.append(merged)

        # Provenance (Task 2.2.3)
        for p in group_pockets:
            provenance.append({
                "receptor_id": state,
                "merged_pocket_id": pocket_id,
                "raw_pocket_id": p["candidate_pocket_id"],
                "proposal_source": p.get("proposal_source", ""),
                "pocket_rank": p.get("pocket_rank", ""),
                "proposal_score": p.get("proposal_score", ""),
                "druggability_score": p.get("druggability_score", ""),
            })

    return merged_pockets, merge_table, provenance


def _build_merged_pocket(
    state: str,
    pocket_id: str,
    group: List[dict],
) -> dict:
    """Build a single merged pocket from a group of raw pockets."""
    # Centroid: weighted average by n_residues (or simple average)
    coords = []
    for p in group:
        cx = safe_float(p.get("centroid_x", ""))
        cy = safe_float(p.get("centroid_y", ""))
        cz = safe_float(p.get("centroid_z", ""))
        if cx is not None and cy is not None and cz is not None:
            coords.append((cx, cy, cz))

    if coords:
        avg_x = sum(c[0] for c in coords) / len(coords)
        avg_y = sum(c[1] for c in coords) / len(coords)
        avg_z = sum(c[2] for c in coords) / len(coords)
    else:
        avg_x = avg_y = avg_z = 0.0

    # Union of residues
    all_res = set()
    for p in group:
        res_str = p.get("residue_ids", "")
        if res_str:
            for r in res_str.split(";"):
                r = r.strip()
                if r:
                    all_res.add(r)

    sorted_res = sorted(all_res, key=lambda r: _resnum_from_id(r))
    res_nums = [str(_resnum_from_id(r)) for r in sorted_res]

    # Best scores
    scores = [safe_float(p.get("proposal_score", "")) for p in group]
    scores = [s for s in scores if s is not None]
    drug_scores = [safe_float(p.get("druggability_score", "")) for p in group]
    drug_scores = [s for s in drug_scores if s is not None]
    volumes = [safe_float(p.get("volume_A3", "")) for p in group]
    volumes = [v for v in volumes if v is not None]

    # Box: bounding box of all pocket atoms (use max box dimensions)
    box_sizes = []
    for dim in ("box_size_x", "box_size_y", "box_size_z"):
        vals = [safe_float(p.get(dim, "")) for p in group]
        vals = [v for v in vals if v is not None]
        box_sizes.append(max(vals) if vals else 0.0)

    # Sources
    sources = sorted(set(p.get("proposal_source", "") for p in group))
    raw_ids = [p["candidate_pocket_id"] for p in group]

    return {
        "receptor_id": state,
        "pocket_id": pocket_id,
        "n_sources": len(set(sources)),
        "sources": ";".join(sources),
        "raw_pocket_ids": ";".join(raw_ids),
        "centroid_x": f"{avg_x:.3f}",
        "centroid_y": f"{avg_y:.3f}",
        "centroid_z": f"{avg_z:.3f}",
        "best_proposal_score": f"{max(scores):.4f}" if scores else "",
        "best_druggability_score": f"{max(drug_scores):.4f}" if drug_scores else "",
        "mean_volume_A3": f"{sum(volumes)/len(volumes):.1f}" if volumes else "",
        "residue_ids": ";".join(sorted_res),
        "residue_nums": ";".join(res_nums),
        "n_residues": len(sorted_res),
        "box_center_x": f"{avg_x:.3f}",
        "box_center_y": f"{avg_y:.3f}",
        "box_center_z": f"{avg_z:.3f}",
        "box_size_x": f"{box_sizes[0]:.1f}",
        "box_size_y": f"{box_sizes[1]:.1f}",
        "box_size_z": f"{box_sizes[2]:.1f}",
    }


# ---------------------------------------------------------------------------
# Distance / overlap metrics (Task 2.2.1)
# ---------------------------------------------------------------------------

def _centroid_distance(a: dict, b: dict) -> float:
    """Euclidean distance between two pocket centroids."""
    ax = safe_float(a.get("centroid_x", "")) or 0.0
    ay = safe_float(a.get("centroid_y", "")) or 0.0
    az = safe_float(a.get("centroid_z", "")) or 0.0
    bx = safe_float(b.get("centroid_x", "")) or 0.0
    by = safe_float(b.get("centroid_y", "")) or 0.0
    bz = safe_float(b.get("centroid_z", "")) or 0.0
    return math.sqrt((ax - bx)**2 + (ay - by)**2 + (az - bz)**2)


def _residue_jaccard(a: dict, b: dict) -> float:
    """Jaccard overlap between two pocket residue sets."""
    set_a = _parse_residue_set(a.get("residue_ids", ""))
    set_b = _parse_residue_set(b.get("residue_ids", ""))
    if not set_a and not set_b:
        return 0.0
    intersection = set_a & set_b
    union = set_a | set_b
    return len(intersection) / len(union) if union else 0.0


def _parse_residue_set(residue_str: str) -> Set[str]:
    """Parse semicolon-delimited residue IDs into a set."""
    if not residue_str:
        return set()
    return {r.strip() for r in residue_str.split(";") if r.strip()}


def _best_score(pockets: List[dict], indices: List[int]) -> float:
    """Best proposal score among indexed pockets."""
    scores = []
    for i in indices:
        s = safe_float(pockets[i].get("proposal_score", ""))
        if s is not None:
            scores.append(s)
    return max(scores) if scores else 0.0


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_pocket_merge(output_dir: Path) -> Tuple[Path, Path, Path]:
    """Run full TG 2.2 pipeline.

    Returns (merged_csv, merge_table_csv, provenance_csv).
    """
    # Load raw pockets
    raw_path = output_dir / "candidate_pockets_raw.csv"
    if not raw_path.exists():
        raise FileNotFoundError(
            f"Raw pockets not found: {raw_path}\n"
            "Run TG 2.1 (pocket_proposal --parse) first."
        )

    raw_pockets = load_csv(raw_path)
    print(f"  Loaded {len(raw_pockets)} raw pockets from {raw_path.name}")

    # Merge
    merged, merge_table, provenance = merge_pockets(raw_pockets)

    n_raw = len(raw_pockets)
    n_merged = len(merged)
    n_pairs_merged = sum(1 for m in merge_table if m["merge_decision"] == "merged")
    print(f"  Merge: {n_raw} raw → {n_merged} merged ({n_pairs_merged} pairs merged)")

    # Write outputs
    merged_path = output_dir / "candidate_pockets.csv"
    _write_csv(merged_path, MERGED_POCKET_COLUMNS, merged)

    merge_table_path = output_dir / "candidate_pocket_merge_table.csv"
    _write_csv(merge_table_path, MERGE_TABLE_COLUMNS, merge_table)

    prov_path = output_dir / "candidate_pocket_provenance.csv"
    _write_csv(prov_path, PROVENANCE_COLUMNS, provenance)

    return merged_path, merge_table_path, prov_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_csv(path: Path, columns: List[str], rows: List[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _resnum_from_id(residue_id: str) -> int:
    m = re.search(r"(\d+)$", residue_id)
    return int(m.group(1)) if m else 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 TG 2.2: Candidate pocket normalization and merge"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE2_OUTPUT_DIR,
        help="Phase 2 output directory",
    )
    args = parser.parse_args()

    print("Phase 2 — Task Group 2.2: Pocket Normalization & Merge")
    print(f"  Output: {args.output_dir}")
    print()

    merged_path, merge_path, prov_path = run_pocket_merge(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 2.2 COMPLETE")
    print(f"{'='*60}")
    print(f"  Merged pockets:  {merged_path}")
    print(f"  Merge table:     {merge_path}")
    print(f"  Provenance:      {prov_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
