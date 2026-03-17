#!/usr/bin/env python3
"""Phase 4 Task Group 4.5: Final Review-First Output Design.

Creates review-oriented tables with condensed and expanded views,
supporting traceability to upstream evidence.

Tasks:
  4.5.1 - Condensed review table (one row per pocket x ligand)
  4.5.2 - Expanded evidence table (full provenance)
  4.5.3 - Traceability to upstream phases

Inputs:
  - output/phase4_perturbation/phase4_state_interpretation.csv
  - output/phase4_perturbation/perturbation_candidate_table.csv
  - output/phase4_perturbation/phase4_evidence_normalized.csv
  - output/phase4_perturbation/final_candidate_classes.csv

Outputs:
  - output/phase4_perturbation/phase4_final_review_table.csv
  - output/phase4_perturbation/phase4_expanded_evidence_table.csv

Usage:
    python -m egfr_pipeline.phase4.review_output [--output_dir ...]
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
PHASE4_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase4_perturbation"

# ---------------------------------------------------------------------------
# Output schemas
# ---------------------------------------------------------------------------

REVIEW_TABLE_COLUMNS = [
    "rank",
    "pocket_id",
    "receptor_id",
    "ligand_id",
    "mechanistic_class",
    "perturbation_score",
    "adjusted_confidence",
    # Key axes
    "A1_ppi_interface",
    "A2_druggability",
    "A3_perturbation_relevance",
    "A4_state_robustness",
    # Key evidence (condensed)
    "relationship_class",
    "hotspot_overlap_count",
    "overall_druggability_tier",
    "ligand_support_strength",
    "state_class",
    "accessibility_label",
    "state_caveat",
]

EXPANDED_TABLE_COLUMNS = [
    "rank",
    "pocket_id",
    "receptor_id",
    "ligand_id",
    # Mechanistic
    "mechanistic_class",
    "classification_basis",
    "classification_confidence",
    "adjusted_confidence",
    # Scores
    "perturbation_score",
    "A1_ppi_interface",
    "A2_druggability",
    "A3_perturbation_relevance",
    "A4_state_robustness",
    # Phase 1: PPI patch
    "n_hotspot_residues",
    "hotspot_residues",
    "hotspot_overlap_count",
    "hotspot_overlap_fraction",
    "mean_hotspot_confidence",
    "n_robust_hotspots",
    "n_method_agreement_both",
    # Phase 2: Pocket
    "relationship_class",
    "overall_druggability_tier",
    "druggability_confidence",
    "best_proposal_score",
    # Phase 3: Docking
    "pose_support_count",
    "best_affinity_kcal",
    "mean_affinity_kcal",
    "ligand_support_strength",
    "docking_priority",
    "pocket_status",
    "diversity_verdict",
    "coverage_fraction",
    "normalized_entropy",
    # Phase 4: State
    "state_class",
    "n_states_matched",
    "state_interpretation",
    "accessibility_label",
    "confidence_modifier",
    "state_caveat",
    # Provenance flags
    "phase1_loaded",
    "phase2_loaded",
    "phase3_loaded",
]


# ---------------------------------------------------------------------------
# Build review tables
# ---------------------------------------------------------------------------

def build_review_tables(
    interpretation: List[dict],
    candidates: List[dict],
    evidence: List[dict],
    classes: List[dict],
) -> Tuple[List[dict], List[dict]]:
    """Merge all upstream data into condensed and expanded views.

    Returns (review_table, expanded_table).
    """
    # Index by (pocket_id, ligand_id) for merging
    cand_map = _index_by_key(candidates)
    ev_map = _index_by_key(evidence)
    cls_map = _index_by_key(classes)

    review_rows = []
    expanded_rows = []

    for row in interpretation:
        key = (row["pocket_id"], row.get("ligand_id", ""))
        cand = cand_map.get(key, {})
        ev = ev_map.get(key, {})
        cls = cls_map.get(key, {})

        # Merge all sources (interpretation > candidates > classes > evidence)
        merged = {}
        merged.update(ev)
        merged.update(cls)
        merged.update(cand)
        merged.update(row)

        # Condensed review row
        review_row = {col: merged.get(col, "") for col in REVIEW_TABLE_COLUMNS}
        review_rows.append(review_row)

        # Expanded evidence row
        expanded_row = {col: merged.get(col, "") for col in EXPANDED_TABLE_COLUMNS}
        expanded_rows.append(expanded_row)

    return review_rows, expanded_rows


def _index_by_key(rows: List[dict]) -> Dict[tuple, dict]:
    """Index rows by (pocket_id, ligand_id)."""
    result = {}
    for r in rows:
        pid = r.get("pocket_id") or r.get("candidate_pocket_id", "")
        if not pid:
            continue  # Skip rows with no pocket identifier
        key = (pid, r.get("ligand_id", ""))
        result[key] = r
    return result


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_review_output(
    output_dir: Path = PHASE4_OUTPUT_DIR,
) -> Tuple[Path, Path]:
    """Run full TG 4.5 pipeline.

    Returns (review_table_path, expanded_table_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load all upstream tables
    interp = _load_csv(output_dir / "phase4_state_interpretation.csv")
    if not interp:
        raise FileNotFoundError(
            "State interpretation not found. Run TG 4.4 first.")

    candidates = _load_csv(output_dir / "perturbation_candidate_table.csv")
    evidence = _load_csv(output_dir / "phase4_evidence_normalized.csv")
    classes = _load_csv(output_dir / "final_candidate_classes.csv")

    print(f"  Loaded: {len(interp)} interpreted, {len(candidates)} candidates, "
          f"{len(evidence)} evidence, {len(classes)} classes")

    # Build tables
    review_rows, expanded_rows = build_review_tables(
        interp, candidates, evidence, classes)

    print(f"  Review table: {len(review_rows)} rows, "
          f"{len(REVIEW_TABLE_COLUMNS)} columns")
    print(f"  Expanded table: {len(expanded_rows)} rows, "
          f"{len(EXPANDED_TABLE_COLUMNS)} columns")

    # Write outputs
    review_path = output_dir / "phase4_final_review_table.csv"
    _write_csv(review_path, REVIEW_TABLE_COLUMNS, review_rows)

    expanded_path = output_dir / "phase4_expanded_evidence_table.csv"
    _write_csv(expanded_path, EXPANDED_TABLE_COLUMNS, expanded_rows)

    return review_path, expanded_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_csv(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _write_csv(path: Path, columns: List[str], rows: List[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 4 TG 4.5: Final Review-First Output Design"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE4_OUTPUT_DIR,
        help="Phase 4 output directory",
    )
    args = parser.parse_args()

    print("Phase 4 — Task Group 4.5: Final Review-First Output Design")
    print(f"  Output: {args.output_dir}")
    print()

    review_path, expanded_path = run_review_output(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 4.5 COMPLETE")
    print(f"{'='*60}")
    print(f"  Review table:   {review_path}")
    print(f"  Expanded table: {expanded_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
