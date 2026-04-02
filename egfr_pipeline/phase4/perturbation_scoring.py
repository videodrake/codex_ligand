#!/usr/bin/env python3
"""Phase 4 Task Group 4.3: Perturbation-Relevance Scoring Logic.

Combines axis scores into candidate-level ranking, prevents affinity
domination, and preserves evidence provenance.

Tasks:
  4.3.1 - Combine axis scores into candidate-level ranking
  4.3.2 - Prevent affinity domination
  4.3.3 - Preserve evidence provenance

Inputs:
  - output/phase4_perturbation/final_candidate_classes.csv
  - output/phase4_perturbation/phase4_axis_definition_table.csv

Outputs:
  - output/phase4_perturbation/perturbation_axis_scores.csv
  - output/phase4_perturbation/perturbation_candidate_table.csv
  - output/phase4_perturbation/phase4_ranking_method_note.md

Usage:
    python -m egfr_pipeline.phase4.perturbation_scoring [--output_dir ...]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from egfr_pipeline.csv_utils import load_csv
from egfr_pipeline.parsing_utils import safe_float

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
PHASE4_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase4_perturbation"

# ---------------------------------------------------------------------------
# Axis weights (from TG 4.1, reproduced here for scoring)
# ---------------------------------------------------------------------------

DEFAULT_WEIGHTS = {
    "A1_ppi_interface": 0.30,
    "A2_druggability": 0.25,
    "A3_perturbation_relevance": 0.30,
    "A4_state_robustness": 0.15,
}

# Affinity cap: max contribution of affinity-derived signal to final score.
# Prevents a single high-affinity hit at an irrelevant site from dominating.
AFFINITY_CAP_FRACTION = 0.35

# ---------------------------------------------------------------------------
# Output schemas
# ---------------------------------------------------------------------------

AXIS_SCORE_COLUMNS = [
    "pocket_id",
    "receptor_id",
    "ligand_id",
    "mechanistic_class",
    "classification_confidence",
    # Raw axis scores
    "A1_ppi_interface",
    "A2_druggability",
    "A3_perturbation_relevance",
    "A4_state_robustness",
    # Weighted contributions
    "A1_weighted",
    "A2_weighted",
    "A3_weighted",
    "A4_weighted",
    # Combined
    "perturbation_score",
    "rank",
]

CANDIDATE_TABLE_COLUMNS = [
    "rank",
    "pocket_id",
    "receptor_id",
    "ligand_id",
    "mechanistic_class",
    "classification_confidence",
    "perturbation_score",
    # Axis scores
    "A1_ppi_interface",
    "A2_druggability",
    "A3_perturbation_relevance",
    "A4_state_robustness",
    # Key evidence provenance
    "relationship_class",
    "hotspot_overlap_count",
    "hotspot_overlap_fraction",
    "overall_druggability_tier",
    "ligand_support_strength",
    "pose_support_count",
    "best_affinity_kcal",
    "state_class",
    "n_states_matched",
    "classification_basis",
]


# ---------------------------------------------------------------------------
# 4.3.1: Combine axis scores
# ---------------------------------------------------------------------------

def compute_perturbation_score(
    row: dict,
    weights: Dict[str, float] = None,
) -> Tuple[float, Dict[str, float]]:
    """Compute weighted perturbation score from axis values.

    Returns (perturbation_score, {axis: weighted_contribution}).
    """
    if weights is None:
        weights = DEFAULT_WEIGHTS

    contributions = {}
    for axis, weight in weights.items():
        raw = safe_float(row.get(axis, "0"), 0.0)
        contributions[axis] = round(raw * weight, 4)

    total = sum(contributions.values())
    return round(total, 4), contributions


# ---------------------------------------------------------------------------
# 4.3.2: Affinity domination prevention
# ---------------------------------------------------------------------------

def apply_affinity_cap(
    row: dict,
    score: float,
    cap_fraction: float = AFFINITY_CAP_FRACTION,
) -> float:
    """Prevent affinity from dominating the final score.

    If a site is classified as ppi_irrelevant but has high affinity,
    cap its effective score contribution from affinity-derived axes.

    In practice: for irrelevant sites, cap A2+A3 contribution.
    This ensures biological relevance (A1, A4) remains important.
    """
    mech_class = row.get("mechanistic_class", "")

    # Only apply cap to irrelevant or uncertain sites
    if mech_class not in (
        "ligandable_but_ppi_irrelevant_candidate",
        "uncertain_mechanism_candidate",
    ):
        return score

    # For irrelevant sites: cap the affinity-influenced axes (A2, A3)
    a2 = safe_float(row.get("A2_druggability", "0"), 0.0)
    a3 = safe_float(row.get("A3_perturbation_relevance", "0"), 0.0)
    a1 = safe_float(row.get("A1_ppi_interface", "0"), 0.0)
    a4 = safe_float(row.get("A4_state_robustness", "0"), 0.0)

    w = DEFAULT_WEIGHTS
    bio_component = a1 * w["A1_ppi_interface"] + a4 * w["A4_state_robustness"]
    affinity_component = a2 * w["A2_druggability"] + a3 * w["A3_perturbation_relevance"]

    if cap_fraction <= 0:
        return round(bio_component, 4)
    if cap_fraction >= 1:
        return round(score, 4)

    max_affinity = (cap_fraction / (1.0 - cap_fraction)) * bio_component
    capped_affinity = min(affinity_component, max_affinity)

    return round(bio_component + capped_affinity, 4)


# ---------------------------------------------------------------------------
# Score and rank all candidates
# ---------------------------------------------------------------------------

def score_and_rank(
    classified: List[dict],
    weights: Dict[str, float] = None,
) -> Tuple[List[dict], List[dict]]:
    """Score all candidates and produce ranked tables.

    Returns (axis_scores, candidate_table).
    """
    scored = []

    for row in classified:
        raw_score, contributions = compute_perturbation_score(row, weights)
        final_score = apply_affinity_cap(row, raw_score)

        entry = {
            "pocket_id": row["pocket_id"],
            "receptor_id": row["receptor_id"],
            "ligand_id": row.get("ligand_id", ""),
            "mechanistic_class": row["mechanistic_class"],
            "classification_confidence": row["classification_confidence"],
            "A1_ppi_interface": row.get("A1_ppi_interface", ""),
            "A2_druggability": row.get("A2_druggability", ""),
            "A3_perturbation_relevance": row.get("A3_perturbation_relevance", ""),
            "A4_state_robustness": row.get("A4_state_robustness", ""),
            "A1_weighted": contributions.get("A1_ppi_interface", 0),
            "A2_weighted": contributions.get("A2_druggability", 0),
            "A3_weighted": contributions.get("A3_perturbation_relevance", 0),
            "A4_weighted": contributions.get("A4_state_robustness", 0),
            "perturbation_score": final_score,
            "rank": 0,  # filled after sorting
        }
        # Carry forward provenance
        for col in CANDIDATE_TABLE_COLUMNS:
            if col not in entry:
                entry[col] = row.get(col, "")
        scored.append(entry)

    # Sort by perturbation_score descending, then pocket_id, ligand_id
    scored.sort(key=lambda x: (
        -x["perturbation_score"],
        x["pocket_id"],
        x.get("ligand_id", ""),
    ))

    # Assign ranks
    for i, s in enumerate(scored, 1):
        s["rank"] = i

    # Build candidate table (same data, different column subset)
    candidate_table = []
    for s in scored:
        ct = {col: s.get(col, "") for col in CANDIDATE_TABLE_COLUMNS}
        candidate_table.append(ct)

    return scored, candidate_table


# ---------------------------------------------------------------------------
# Ranking method note
# ---------------------------------------------------------------------------

def build_ranking_note(
    scored: List[dict],
    weights: Dict[str, float],
) -> List[str]:
    """Build the ranking method documentation."""
    lines = [
        "# Phase 4 Ranking Method Note",
        "",
        "## 1. Scoring Formula",
        "",
        "The perturbation score is a weighted sum of 4 axis scores:",
        "",
        "```",
        "perturbation_score = sum(axis_score_i * weight_i)",
        "```",
        "",
        "### Axis Weights",
        "",
        "| Axis | Weight | Rationale |",
        "|------|--------|-----------|",
    ]
    rationales = {
        "A1_ppi_interface": "Core biological question: does site overlap MYO1D patch?",
        "A2_druggability": "Practical requirement: can a small molecule bind here?",
        "A3_perturbation_relevance": "Combined biological + chemical evidence",
        "A4_state_robustness": "Supporting evidence: is site consistently accessible?",
    }
    for axis, weight in sorted(weights.items()):
        lines.append(
            f"| {axis} | {weight:.0%} | {rationales.get(axis, '')} |")
    lines.append("")

    lines.extend([
        "## 2. Affinity Domination Prevention",
        "",
        f"For sites classified as `ligandable_but_ppi_irrelevant` or "
        f"`uncertain_mechanism`, the combined contribution of "
        f"affinity-influenced axes (A2 + A3) is capped at "
        f"{AFFINITY_CAP_FRACTION:.0%} of the total score. This ensures "
        f"that a high-affinity hit at a biologically irrelevant site "
        f"cannot outrank a moderate-affinity hit at an orthosteric site.",
        "",
        "## 3. Ranking Results",
        "",
        "| Rank | Pocket | Ligand | Class | Score | A1 | A2 | A3 | A4 |",
        "|------|--------|--------|-------|-------|----|----|----|----|",
    ])
    for s in scored:
        lig = s.get("ligand_id", "") or "(none)"
        short_class = s["mechanistic_class"].replace("_candidate", "")
        lines.append(
            f"| {s['rank']} | {s['pocket_id']} | {lig} | "
            f"{short_class} | {s['perturbation_score']:.4f} | "
            f"{safe_float(s['A1_ppi_interface'], 0.0):.3f} | "
            f"{safe_float(s['A2_druggability'], 0.0):.3f} | "
            f"{safe_float(s['A3_perturbation_relevance'], 0.0):.3f} | "
            f"{safe_float(s['A4_state_robustness'], 0.0):.3f} |"
        )
    lines.append("")

    # Affinity cap impact
    capped = [s for s in scored
              if s["mechanistic_class"] in (
                  "ligandable_but_ppi_irrelevant_candidate",
                  "uncertain_mechanism_candidate")]
    if capped:
        lines.extend([
            "### Affinity Cap Applied To",
            "",
        ])
        for s in capped:
            lines.append(
                f"- {s['pocket_id']}: capped score = {s['perturbation_score']:.4f}")
        lines.append("")

    lines.extend([
        "## 4. Provenance",
        "",
        "Each ranked candidate preserves references to:",
        "- Phase 1: hotspot overlap count/fraction",
        "- Phase 2: relationship class, druggability tier",
        "- Phase 3: ligand support strength, pose count, best affinity",
        "- Phase 4: mechanistic class, classification basis",
        "",
        "These fields are available in `perturbation_candidate_table.csv` "
        "for full traceability.",
        "",
        "---",
        "",
        "Generated by `egfr_pipeline.phase4.perturbation_scoring`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_perturbation_scoring(
    output_dir: Path = PHASE4_OUTPUT_DIR,
    weights: Dict[str, float] = None,
) -> Tuple[Path, Path, Path]:
    """Run full TG 4.3 pipeline.

    Returns (axis_scores_path, candidate_table_path, note_path).
    """
    if weights is None:
        weights = DEFAULT_WEIGHTS

    output_dir.mkdir(parents=True, exist_ok=True)

    # Load classified candidates from TG 4.2
    classes_path = output_dir / "final_candidate_classes.csv"
    if not classes_path.exists():
        raise FileNotFoundError(
            f"Candidate classes not found: {classes_path}\n"
            "Run TG 4.2 (mechanistic_classification) first."
        )

    classified = load_csv(classes_path)
    print(f"  Loaded {len(classified)} classified candidates")

    # Score and rank
    scored, candidate_table = score_and_rank(classified, weights)

    for s in scored:
        lig = s.get("ligand_id", "") or "(none)"
        print(f"    Rank {s['rank']}: {s['pocket_id']} / {lig} — "
              f"score={s['perturbation_score']:.4f} "
              f"({s['mechanistic_class']})")

    # Write outputs
    axis_path = output_dir / "perturbation_axis_scores.csv"
    _write_csv(axis_path, AXIS_SCORE_COLUMNS, scored)

    table_path = output_dir / "perturbation_candidate_table.csv"
    _write_csv(table_path, CANDIDATE_TABLE_COLUMNS, candidate_table)

    note_lines = build_ranking_note(scored, weights)
    note_path = output_dir / "phase4_ranking_method_note.md"
    with open(note_path, "w", encoding="utf-8") as f:
        f.write("\n".join(note_lines))

    return axis_path, table_path, note_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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
        description="Phase 4 TG 4.3: Perturbation-Relevance Scoring"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE4_OUTPUT_DIR,
        help="Phase 4 output directory",
    )
    args = parser.parse_args()

    print("Phase 4 — Task Group 4.3: Perturbation-Relevance Scoring")
    print(f"  Output: {args.output_dir}")
    print()

    axis_path, table_path, note_path = run_perturbation_scoring(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 4.3 COMPLETE")
    print(f"{'='*60}")
    print(f"  Axis scores:       {axis_path}")
    print(f"  Candidate table:   {table_path}")
    print(f"  Ranking method:    {note_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
