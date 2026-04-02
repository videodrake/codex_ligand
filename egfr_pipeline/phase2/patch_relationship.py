#!/usr/bin/env python3
"""Phase 2 Task Group 2.3: Patch Relationship Classification.

Classifies each candidate pocket by its structural relationship to the Phase 1
MYO1D receptor-side patch.

Classes:
  - orthosteric_candidate:  Pocket directly overlaps PPI hotspot residues
  - rim_candidate:          Pocket partially overlaps or borders the patch
  - allosteric_candidate:   Pocket is nearby but does not overlap patch residues
  - low_relevance_candidate: Pocket is structurally distant from the patch

Two independent metrics drive classification:
  1. Hotspot residue overlap (pocket lining residues ∩ patch hotspot residues)
  2. Centroid distance (pocket centroid ↔ patch centroid, computed from receptor PDB)

Tasks:
  2.3.1 - Define relationship classes
  2.3.2 - Relationship metrics (centroid distance + residue overlap)
  2.3.3 - Classification transparency (raw metrics preserved)

Inputs:
  - output/phase2_pockets/candidate_pockets.csv
  - output/phase2_pockets/phase2_patch_reference_normalized.csv
  - input/PPI/phase1/receptor_*.pdb (for patch centroid computation)

Outputs:
  - output/phase2_pockets/pocket_patch_relationship.csv
  - output/phase2_pockets/pocket_patch_relationship_metrics.csv
  - output/phase2_pockets/phase2_relationship_classification_note.md

Usage:
    python -m egfr_pipeline.phase2.patch_relationship [--output_dir ...]
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
from egfr_pipeline.parsing_utils import safe_float, safe_int

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE1_INPUT_DIR = PROJECT_ROOT / "input" / "PPI" / "phase1"
PHASE2_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase2_pockets"
RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]

# Classification thresholds
ORTHOSTERIC_HOTSPOT_FRACTION = 0.25   # ≥25% of patch hotspots in pocket → orthosteric
ORTHOSTERIC_MIN_OVERLAP = 2           # minimum 2 shared hotspot residues
RIM_MIN_OVERLAP = 1                   # ≥1 shared residue → at least rim
RIM_MAX_CENTROID_DIST_A = 12.0        # within 12Å centroid + some overlap → rim
ALLOSTERIC_MAX_CENTROID_DIST_A = 20.0 # within 20Å, no overlap → allosteric
# Beyond 20Å and no overlap → low_relevance

# Output schemas
RELATIONSHIP_COLUMNS = [
    "receptor_id",
    "pocket_id",
    "relationship_class",         # orthosteric/rim/allosteric/low_relevance
    "hotspot_overlap_count",      # number of shared residues
    "hotspot_overlap_fraction",   # shared / total patch hotspots
    "centroid_distance_A",        # Å
    "overlapping_residues",       # semicolon-delimited
    "n_pocket_residues",
    "n_patch_hotspots",
    "classification_basis",       # which metric(s) drove the decision
]

METRICS_COLUMNS = [
    "receptor_id",
    "pocket_id",
    "pocket_centroid_x",
    "pocket_centroid_y",
    "pocket_centroid_z",
    "patch_centroid_x",
    "patch_centroid_y",
    "patch_centroid_z",
    "centroid_distance_A",
    "pocket_residue_ids",
    "patch_hotspot_ids",
    "overlap_residue_ids",
    "hotspot_overlap_count",
    "hotspot_overlap_fraction",
    "pocket_residue_jaccard",     # Jaccard(pocket_res, patch_res)
    "relationship_class",
]


# ---------------------------------------------------------------------------
# Patch centroid computation from PDB
# ---------------------------------------------------------------------------

def compute_patch_centroid(
    receptor_pdb: Path,
    patch_residue_nums: Set[int],
) -> Optional[Tuple[float, float, float]]:
    """Compute centroid of patch residues from receptor PDB (CA atoms).

    Returns (x, y, z) or None if no matching residues found.
    """
    coords = []

    if not receptor_pdb.exists():
        return None

    with open(receptor_pdb, encoding="utf-8") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            chain = line[21].strip()
            if chain and chain != "A":
                continue
            resnum_str = line[22:26].strip()
            resnum = safe_int(resnum_str)
            if resnum is not None and resnum in patch_residue_nums:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except (ValueError, IndexError):
                    pass

    if not coords:
        return None

    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)
    return (cx, cy, cz)


# ---------------------------------------------------------------------------
# Classification logic (Tasks 2.3.1 + 2.3.2)
# ---------------------------------------------------------------------------

def classify_pocket_patch_relationship(
    pocket: dict,
    patch_hotspot_ids: Set[str],
    patch_centroid: Optional[Tuple[float, float, float]],
) -> dict:
    """Classify one pocket's relationship to the PPI patch.

    Returns a dict with relationship class and all metrics.
    """
    # Parse pocket residues
    pocket_res_str = pocket.get("residue_ids", "")
    pocket_res = {r.strip() for r in pocket_res_str.split(";") if r.strip()}

    # Compute overlap
    overlap = pocket_res & patch_hotspot_ids
    n_overlap = len(overlap)
    n_patch = len(patch_hotspot_ids)
    overlap_fraction = n_overlap / n_patch if n_patch > 0 else 0.0

    # Jaccard
    union = pocket_res | patch_hotspot_ids
    jaccard = len(overlap) / len(union) if union else 0.0

    # Centroid distance
    pocket_cx = safe_float(pocket.get("centroid_x", ""))
    pocket_cy = safe_float(pocket.get("centroid_y", ""))
    pocket_cz = safe_float(pocket.get("centroid_z", ""))

    if (patch_centroid is not None
            and pocket_cx is not None
            and pocket_cy is not None
            and pocket_cz is not None):
        dist = math.sqrt(
            (pocket_cx - patch_centroid[0])**2 +
            (pocket_cy - patch_centroid[1])**2 +
            (pocket_cz - patch_centroid[2])**2
        )
    else:
        dist = None

    # Classification (Task 2.3.1)
    rel_class, basis = _assign_class(n_overlap, overlap_fraction, dist)

    return {
        "pocket_id": pocket.get("pocket_id", ""),
        "receptor_id": pocket.get("receptor_id", ""),
        "relationship_class": rel_class,
        "hotspot_overlap_count": n_overlap,
        "hotspot_overlap_fraction": f"{overlap_fraction:.3f}",
        "centroid_distance_A": f"{dist:.2f}" if dist is not None else "",
        "overlapping_residues": ";".join(sorted(overlap, key=_resnum_from_id)),
        "n_pocket_residues": len(pocket_res),
        "n_patch_hotspots": n_patch,
        "classification_basis": basis,
        # Extended metrics
        "pocket_centroid_x": pocket.get("centroid_x", ""),
        "pocket_centroid_y": pocket.get("centroid_y", ""),
        "pocket_centroid_z": pocket.get("centroid_z", ""),
        "patch_centroid_x": f"{patch_centroid[0]:.3f}" if patch_centroid else "",
        "patch_centroid_y": f"{patch_centroid[1]:.3f}" if patch_centroid else "",
        "patch_centroid_z": f"{patch_centroid[2]:.3f}" if patch_centroid else "",
        "pocket_residue_ids": pocket.get("residue_ids", ""),
        "patch_hotspot_ids": ";".join(sorted(patch_hotspot_ids, key=_resnum_from_id)),
        "overlap_residue_ids": ";".join(sorted(overlap, key=_resnum_from_id)),
        "pocket_residue_jaccard": f"{jaccard:.3f}",
    }


def _assign_class(
    n_overlap: int,
    overlap_fraction: float,
    centroid_dist: Optional[float],
) -> Tuple[str, str]:
    """Assign relationship class based on metrics.

    Returns (class_label, basis_description).
    """
    # Orthosteric: significant hotspot overlap
    if (n_overlap >= ORTHOSTERIC_MIN_OVERLAP
            and overlap_fraction >= ORTHOSTERIC_HOTSPOT_FRACTION):
        return (
            "orthosteric_candidate",
            f"hotspot_overlap={n_overlap} (frac={overlap_fraction:.2f} >= {ORTHOSTERIC_HOTSPOT_FRACTION})"
        )

    # Rim: partial overlap or close + some overlap
    if n_overlap >= RIM_MIN_OVERLAP:
        return (
            "rim_candidate",
            f"partial_hotspot_overlap={n_overlap} (frac={overlap_fraction:.2f})"
        )

    if (centroid_dist is not None
            and centroid_dist <= RIM_MAX_CENTROID_DIST_A
            and n_overlap > 0):
        return (
            "rim_candidate",
            f"centroid_close={centroid_dist:.1f}A + overlap={n_overlap}"
        )

    # Allosteric: nearby but no residue overlap
    if centroid_dist is not None and centroid_dist <= ALLOSTERIC_MAX_CENTROID_DIST_A:
        return (
            "allosteric_candidate",
            f"centroid_nearby={centroid_dist:.1f}A <= {ALLOSTERIC_MAX_CENTROID_DIST_A}A, no_hotspot_overlap"
        )

    # Low relevance
    if centroid_dist is not None:
        return (
            "low_relevance_candidate",
            f"centroid_distant={centroid_dist:.1f}A > {ALLOSTERIC_MAX_CENTROID_DIST_A}A, no_overlap"
        )

    # Fallback: no centroid available, classify by overlap only
    if n_overlap == 0:
        return (
            "low_relevance_candidate",
            "no_overlap, no_centroid_available"
        )

    return ("rim_candidate", f"overlap={n_overlap}, no_centroid_available")


def _patch_row_applies_to_state(row: dict, state: str) -> bool:
    """Return whether a normalized patch-reference row applies to a receptor state."""
    states_str = (row.get("states", "") or "").strip()
    if not states_str:
        return True

    normalized_tokens = {
        token.strip().lower()
        for token in states_str.replace(",", ";").split(";")
        if token.strip()
    }
    if not normalized_tokens:
        return True

    if normalized_tokens & {"all", "any", "global", "all_states"}:
        return True

    known_states = {known.lower() for known in RECEPTOR_STATES}
    recognized_states = normalized_tokens & known_states
    if not recognized_states:
        return True

    return state.lower() in recognized_states


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def build_classification_note(
    relationships: List[dict],
    patch_hotspot_ids: Set[str],
) -> List[str]:
    """Build the relationship classification markdown note."""
    lines = [
        "# Phase 2 Patch Relationship Classification Note",
        "",
        "## Overview",
        "",
        f"- Patch hotspot residues: {len(patch_hotspot_ids)}",
        f"- Candidate pockets classified: {len(relationships)}",
        "",
    ]

    # Classification thresholds
    lines.extend([
        "## Classification Thresholds",
        "",
        "| Class | Criteria |",
        "|-------|---------|",
        f"| orthosteric | hotspot overlap >= {ORTHOSTERIC_MIN_OVERLAP} residues "
        f"AND fraction >= {ORTHOSTERIC_HOTSPOT_FRACTION} |",
        f"| rim | hotspot overlap >= {RIM_MIN_OVERLAP} residue(s) "
        f"OR centroid <= {RIM_MAX_CENTROID_DIST_A}A + overlap |",
        f"| allosteric | centroid <= {ALLOSTERIC_MAX_CENTROID_DIST_A}A, no overlap |",
        f"| low_relevance | centroid > {ALLOSTERIC_MAX_CENTROID_DIST_A}A, no overlap |",
        "",
    ])

    # Summary by class
    by_class = defaultdict(list)
    for r in relationships:
        by_class[r["relationship_class"]].append(r)

    lines.extend([
        "## Classification Summary",
        "",
        "| Class | Count | Pockets |",
        "|-------|-------|---------|",
    ])
    for cls in ["orthosteric_candidate", "rim_candidate",
                "allosteric_candidate", "low_relevance_candidate"]:
        pockets = by_class.get(cls, [])
        ids = ", ".join(p["pocket_id"] for p in pockets)
        lines.append(f"| {cls} | {len(pockets)} | {ids} |")
    lines.append("")

    # Detail table
    if relationships:
        lines.extend([
            "## Per-Pocket Detail",
            "",
            "| Pocket | Class | Overlap | Fraction | Dist (A) | Basis |",
            "|--------|-------|---------|----------|----------|-------|",
        ])
        for r in relationships:
            lines.append(
                f"| {r['pocket_id']} | {r['relationship_class']} | "
                f"{r['hotspot_overlap_count']} | {r['hotspot_overlap_fraction']} | "
                f"{r['centroid_distance_A']} | {r['classification_basis']} |"
            )
        lines.append("")

    lines.extend([
        "---",
        "",
        "Generated by `egfr_pipeline.phase2.patch_relationship`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_patch_relationship(
    output_dir: Path,
    receptor_dir: Path = PHASE1_INPUT_DIR,
) -> Tuple[Path, Path, Path]:
    """Run full TG 2.3 pipeline.

    Returns (relationship_csv, metrics_csv, note_path).
    """
    # Load inputs
    pockets = load_csv(output_dir / "candidate_pockets.csv")
    patch_ref = load_csv(output_dir / "phase2_patch_reference_normalized.csv")

    if not pockets:
        raise FileNotFoundError("No candidate pockets found. Run TG 2.2 first.")
    if not patch_ref:
        raise FileNotFoundError("No patch reference found. Run TG 2.0 first.")

    print(f"  Loaded {len(pockets)} merged pockets, {len(patch_ref)} patch residues")

    # Build patch hotspot/residue scope per receptor state.
    patch_hotspot_ids_by_state: Dict[str, Set[str]] = {
        state: set() for state in RECEPTOR_STATES
    }
    patch_residue_nums_by_state: Dict[str, Set[int]] = {
        state: set() for state in RECEPTOR_STATES
    }
    patch_all_ids = set()
    patch_hotspot_ids = set()
    for r in patch_ref:
        rid = r.get("patch_residue_id", "")
        rnum = safe_int(r.get("residue_num", ""))
        if rid:
            patch_all_ids.add(rid)
        if _parse_bool(r.get("is_hotspot", "")) and rid:
            patch_hotspot_ids.add(rid)

        for state in RECEPTOR_STATES:
            if not _patch_row_applies_to_state(r, state):
                continue
            if rid:
                if _parse_bool(r.get("is_hotspot", "")):
                    patch_hotspot_ids_by_state[state].add(rid)
            if rnum is not None:
                patch_residue_nums_by_state[state].add(rnum)

    print(f"  Patch: {len(patch_all_ids)} residues, {len(patch_hotspot_ids)} hotspots")

    # Compute patch centroids per state from receptor PDBs
    patch_centroids: Dict[str, Optional[Tuple[float, float, float]]] = {}
    for state in RECEPTOR_STATES:
        pdb_path = receptor_dir / f"receptor_{state}.pdb"
        centroid = compute_patch_centroid(pdb_path, patch_residue_nums_by_state[state])
        patch_centroids[state] = centroid
        if centroid:
            print(f"  Patch centroid ({state}): ({centroid[0]:.1f}, {centroid[1]:.1f}, {centroid[2]:.1f})")
        else:
            print(f"  Patch centroid ({state}): not computed (PDB or residues missing)")

    # Classify each pocket
    relationships = []
    for pocket in pockets:
        state = pocket.get("receptor_id", "")
        centroid = patch_centroids.get(state)
        state_hotspot_ids = patch_hotspot_ids_by_state.get(state, patch_hotspot_ids)
        result = classify_pocket_patch_relationship(
            pocket, state_hotspot_ids, centroid
        )
        relationships.append(result)

    # Summary
    by_class = defaultdict(int)
    for r in relationships:
        by_class[r["relationship_class"]] += 1
    print(f"  Classification: {dict(by_class)}")

    # Write relationship CSV
    rel_path = output_dir / "pocket_patch_relationship.csv"
    _write_csv(rel_path, RELATIONSHIP_COLUMNS, relationships)

    # Write metrics CSV
    metrics_path = output_dir / "pocket_patch_relationship_metrics.csv"
    _write_csv(metrics_path, METRICS_COLUMNS, relationships)

    # Write note
    note_lines = build_classification_note(relationships, patch_hotspot_ids)
    note_path = output_dir / "phase2_relationship_classification_note.md"
    with open(note_path, "w", encoding="utf-8") as f:
        f.write("\n".join(note_lines))

    return rel_path, metrics_path, note_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_csv(path: Path, columns: List[str], rows: List[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _parse_bool(val) -> bool:
    if isinstance(val, bool):
        return val
    return str(val).strip().lower() in ("true", "1", "yes")


def _resnum_from_id(residue_id: str) -> int:
    m = re.search(r"(\d+)$", residue_id)
    return int(m.group(1)) if m else 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 TG 2.3: Patch relationship classification"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE2_OUTPUT_DIR,
        help="Phase 2 output directory",
    )
    parser.add_argument(
        "--receptor_dir", type=Path, default=PHASE1_INPUT_DIR,
        help="Directory containing receptor PDB files",
    )
    args = parser.parse_args()

    print("Phase 2 — Task Group 2.3: Patch Relationship Classification")
    print(f"  Output: {args.output_dir}")
    print()

    rel_path, metrics_path, note_path = run_patch_relationship(
        args.output_dir, args.receptor_dir
    )

    print(f"\n{'='*60}")
    print("TG 2.3 COMPLETE")
    print(f"{'='*60}")
    print(f"  Relationship:  {rel_path}")
    print(f"  Metrics:       {metrics_path}")
    print(f"  Note:          {note_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
