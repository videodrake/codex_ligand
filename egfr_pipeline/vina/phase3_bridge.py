#!/usr/bin/env python3
"""Bridge from Tier 1 Vina blind docking pockets to Phase 3 focused docking.

Converts vina_pocket_table.csv (Tier 1 output) into the
phase3_candidate_pocket_reference.csv format expected by Phase 3 TG 3.0.

Priority is derived from blind docking metrics:
  - best_affinity <= primary_cutoff   → primary
  - best_affinity <= secondary_cutoff → secondary
  - best_affinity <= skip_cutoff      → exploratory
  - otherwise                         → skip

Box size is derived from centroid_spread_A with clamping.
"""

import csv
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from egfr_pipeline.config import load_config
from egfr_pipeline import paths


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Phase 2 export schema (what Phase 3 TG 3.0 expects)
REFERENCE_COLUMNS = [
    "receptor_id",
    "pocket_id",
    "centroid_x",
    "centroid_y",
    "centroid_z",
    "box_size_x",
    "box_size_y",
    "box_size_z",
    "relationship_class",
    "hotspot_overlap_count",
    "hotspot_overlap_fraction",
    "centroid_distance_to_patch_A",
    "druggability_confidence",
    "overall_druggability_tier",
    "best_proposal_score",
    "best_druggability_score",
    "normalized_druggability_score",
    "n_sources",
    "sources",
    "multi_source_support",
    "state_class",
    "n_states_matched",
    "n_residues",
    "residue_ids",
    "mean_volume_A3",
    "docking_priority",
    "priority_basis",
]

# Default priority thresholds (kcal/mol, lower = stronger binding)
DEFAULT_PRIMARY_CUTOFF = -7.0
DEFAULT_SECONDARY_CUTOFF = -5.0
DEFAULT_SKIP_CUTOFF = -3.0

# Box size parameters (Angstrom)
DEFAULT_BOX_SIZE = 22.0
BOX_SPREAD_FACTOR = 2.5   # box = spread * factor + padding
BOX_PADDING = 10.0        # Å padding around spread
BOX_MIN = 18.0
BOX_MAX = 30.0


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def _safe_float(val, default=None):
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def classify_priority(
    best_affinity: Optional[float],
    n_pose: int,
    n_ligand: int,
    *,
    primary_cutoff: float = DEFAULT_PRIMARY_CUTOFF,
    secondary_cutoff: float = DEFAULT_SECONDARY_CUTOFF,
    skip_cutoff: float = DEFAULT_SKIP_CUTOFF,
) -> Tuple[str, str]:
    """Classify a pocket's docking priority from blind Vina metrics.

    Returns (priority, basis_reason).
    """
    if best_affinity is None:
        return "skip", "no_affinity_data"

    if best_affinity <= primary_cutoff:
        return "primary", f"best_affinity={best_affinity:.1f}<={primary_cutoff}"
    if best_affinity <= secondary_cutoff:
        return "secondary", f"best_affinity={best_affinity:.1f}<={secondary_cutoff}"
    if best_affinity <= skip_cutoff:
        return "exploratory", f"best_affinity={best_affinity:.1f}<={skip_cutoff}"
    return "skip", f"best_affinity={best_affinity:.1f}>{skip_cutoff}"


def derive_box_size(centroid_spread: Optional[float]) -> float:
    """Derive docking box size from pocket centroid spread."""
    if centroid_spread is None or centroid_spread < 1.0:
        return DEFAULT_BOX_SIZE
    box = centroid_spread * BOX_SPREAD_FACTOR + BOX_PADDING
    return max(BOX_MIN, min(BOX_MAX, round(box, 1)))


def build_phase3_reference(
    pocket_rows: List[dict],
    *,
    primary_cutoff: float = DEFAULT_PRIMARY_CUTOFF,
    secondary_cutoff: float = DEFAULT_SECONDARY_CUTOFF,
    skip_cutoff: float = DEFAULT_SKIP_CUTOFF,
) -> List[dict]:
    """Convert Vina pocket table rows to Phase 3 candidate reference rows.

    Args:
        pocket_rows: Rows from vina_pocket_table.csv.
        primary_cutoff: Best affinity threshold for primary priority.
        secondary_cutoff: Best affinity threshold for secondary priority.
        skip_cutoff: Best affinity threshold below which pockets are skipped.

    Returns:
        List of dicts matching Phase 3 candidate reference schema.
    """
    reference = []
    for row in pocket_rows:
        receptor_id = row["receptor_id"]
        pocket_id_local = row["pocket_id"]
        # Make pocket_id globally unique (Phase 3 requires unique IDs across receptors)
        pocket_id = f"{receptor_id}_{pocket_id_local}"

        best_aff = _safe_float(row.get("best_affinity"))
        n_pose = int(row.get("n_pose", 0))
        n_ligand = int(row.get("n_ligand", 0))
        spread = _safe_float(row.get("centroid_spread_A"))

        priority, basis = classify_priority(
            best_aff, n_pose, n_ligand,
            primary_cutoff=primary_cutoff,
            secondary_cutoff=secondary_cutoff,
            skip_cutoff=skip_cutoff,
        )

        box_size = derive_box_size(spread)

        n_residues = 0
        residue_ids = row.get("union_contact_residues", "")
        if residue_ids:
            n_residues = len([r for r in residue_ids.split(";") if r])

        reference.append({
            "receptor_id": receptor_id,
            "pocket_id": pocket_id,
            "centroid_x": row.get("centroid_x", ""),
            "centroid_y": row.get("centroid_y", ""),
            "centroid_z": row.get("centroid_z", ""),
            "box_size_x": f"{box_size:.1f}",
            "box_size_y": f"{box_size:.1f}",
            "box_size_z": f"{box_size:.1f}",
            "relationship_class": "",
            "hotspot_overlap_count": "",
            "hotspot_overlap_fraction": "",
            "centroid_distance_to_patch_A": "",
            "druggability_confidence": "",
            "overall_druggability_tier": "",
            "best_proposal_score": "",
            "best_druggability_score": row.get("best_affinity", ""),
            "normalized_druggability_score": "",
            "n_sources": "1",
            "sources": "blind_vina",
            "multi_source_support": "False",
            "state_class": "",
            "n_states_matched": "1",
            "n_residues": str(n_residues),
            "residue_ids": residue_ids,
            "mean_volume_A3": "",
            "docking_priority": priority,
            "priority_basis": basis,
        })

    return reference


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def bridge_vina_to_phase3(
    config_path: str,
    *,
    primary_cutoff: float = DEFAULT_PRIMARY_CUTOFF,
    secondary_cutoff: float = DEFAULT_SECONDARY_CUTOFF,
    skip_cutoff: float = DEFAULT_SKIP_CUTOFF,
) -> Path:
    """Read Vina pocket table and write Phase 3 candidate reference.

    Returns path to the generated reference CSV.
    """
    config = load_config(config_path)
    postprocess_root = paths.wa_phase4_vina_postprocess(config)

    pocket_table_path = postprocess_root / "vina_pocket_table.csv"
    if not pocket_table_path.exists():
        raise FileNotFoundError(
            f"Pocket table not found: {pocket_table_path}\n"
            "Run postprocess steps (cluster + summarize) first."
        )

    with open(pocket_table_path, newline="", encoding="utf-8") as f:
        pocket_rows = list(csv.DictReader(f))

    if not pocket_rows:
        raise ValueError("Pocket table is empty — no pockets to bridge.")

    reference = build_phase3_reference(
        pocket_rows,
        primary_cutoff=primary_cutoff,
        secondary_cutoff=secondary_cutoff,
        skip_cutoff=skip_cutoff,
    )

    # Write to Phase 2 export location (where Phase 3 TG 3.0 expects it)
    phase2_dir = paths.wb_phase2_pocket_analysis(config)
    phase2_dir.mkdir(parents=True, exist_ok=True)
    output_path = phase2_dir / "phase3_candidate_pocket_reference.csv"

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=REFERENCE_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(reference)

    # Summary
    by_priority: Dict[str, int] = {}
    for r in reference:
        pri = r["docking_priority"]
        by_priority[pri] = by_priority.get(pri, 0) + 1

    print(f"  Generated Phase 3 reference: {len(reference)} pockets")
    for pri in ["primary", "secondary", "exploratory", "skip"]:
        if pri in by_priority:
            print(f"    {pri}: {by_priority[pri]}")
    print(f"  → {output_path}")

    return output_path
