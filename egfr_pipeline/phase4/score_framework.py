#!/usr/bin/env python3
"""Phase 4 Task Group 4.1: Multi-Axis Score Framework Definition.

Defines the 4 core scoring axes for perturbation relevance and computes
per-candidate axis scores from the TG 4.0 normalized evidence.

Tasks:
  4.1.1 - Define core evidence axes (PPI, Druggability, Perturbation, State)
  4.1.2 - Define score semantics (clear interpretation per axis)
  4.1.3 - Preserve raw subcomponents (reviewable, not black-box)

Inputs:
  - output/phase4_perturbation/phase4_evidence_normalized.csv

Outputs:
  - output/phase4_perturbation/phase4_axis_definition_table.csv
  - output/phase4_perturbation/phase4_score_framework.md

Usage:
    python -m egfr_pipeline.phase4.score_framework [--output_dir ...]
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Tuple

from egfr_pipeline.csv_utils import load_csv
from egfr_pipeline.parsing_utils import safe_float

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
PHASE4_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase4_perturbation"

# ---------------------------------------------------------------------------
# Axis definitions
# ---------------------------------------------------------------------------

# Each axis: name, weight, description, submetrics
AXIS_DEFINITIONS = [
    {
        "axis_id": "A1_ppi_interface",
        "axis_name": "PPI Interface Confidence",
        "weight": 0.30,
        "description": (
            "How strong and reliable the receptor-side PPI patch evidence is. "
            "Reflects hotspot residue quality, multi-method agreement, "
            "and robustness across receptor states."
        ),
        "submetrics": (
            "hotspot_overlap_fraction, mean_hotspot_confidence, "
            "n_robust_hotspots, n_method_agreement_both"
        ),
        "range": "0.0–1.0",
        "high_means": "Strong PPI patch overlap with validated hotspots",
    },
    {
        "axis_id": "A2_druggability",
        "axis_name": "Druggability Confidence",
        "weight": 0.25,
        "description": (
            "How plausible the pocket is for small-molecule engagement. "
            "Reflects fpocket/P2Rank scoring, pocket volume, "
            "and multi-tool consensus."
        ),
        "submetrics": (
            "overall_druggability_tier, druggability_confidence, "
            "best_proposal_score"
        ),
        "range": "0.0–1.0",
        "high_means": "High-confidence druggable pocket",
    },
    {
        "axis_id": "A3_perturbation_relevance",
        "axis_name": "Perturbation Relevance",
        "weight": 0.30,
        "description": (
            "How directly the site can plausibly alter MYO1D attachment. "
            "Combines PPI relationship class with docking evidence "
            "to assess mechanistic potential."
        ),
        "submetrics": (
            "relationship_class, hotspot_overlap_count, "
            "ligand_support_strength, pose_support_count"
        ),
        "range": "0.0–1.0",
        "high_means": "Orthosteric site with ligand support",
    },
    {
        "axis_id": "A4_state_robustness",
        "axis_name": "State Robustness / Accessibility",
        "weight": 0.15,
        "description": (
            "Whether the site is persistent across receptor states, "
            "conditionally accessible, or strongly state-dependent. "
            "State-specific sites are not penalized but flagged."
        ),
        "submetrics": (
            "state_class, n_states_matched, diversity_verdict, "
            "coverage_fraction"
        ),
        "range": "0.0–1.0",
        "high_means": "Robust pocket present in multiple conformational states",
    },
]

AXIS_TABLE_COLUMNS = [
    "axis_id", "axis_name", "weight", "description",
    "submetrics", "range", "high_means",
]

# Output: per-candidate axis scores
SCORE_COLUMNS = [
    "pocket_id",
    "receptor_id",
    "ligand_id",
    # Raw axis scores
    "A1_ppi_interface",
    "A2_druggability",
    "A3_perturbation_relevance",
    "A4_state_robustness",
    # Submetrics preserved
    "hotspot_overlap_fraction",
    "mean_hotspot_confidence",
    "n_robust_hotspots",
    "n_method_agreement_both",
    "overall_druggability_tier",
    "druggability_confidence",
    "best_proposal_score",
    "relationship_class",
    "hotspot_overlap_count",
    "ligand_support_strength",
    "pose_support_count",
    "state_class",
    "n_states_matched",
    "diversity_verdict",
    "coverage_fraction",
    # Provenance
    "phase1_loaded",
    "phase2_loaded",
    "phase3_loaded",
]


# ---------------------------------------------------------------------------
# 4.1.1 + 4.1.2: Axis score computation
# ---------------------------------------------------------------------------

# -- Axis A1: PPI Interface Confidence --

def compute_a1_ppi_interface(row: dict) -> float:
    """Compute PPI interface confidence score (0-1).

    Submetrics:
    - hotspot_overlap_fraction (0-1): direct overlap weight (40%)
    - mean_hotspot_confidence (1-3 → 0-1): quality of hotspots (25%)
    - n_robust_hotspots (0-7 → 0-1): robustness support (20%)
    - n_method_agreement_both (0-7 → 0-1): multi-method (15%)
    """
    overlap = safe_float(row.get("hotspot_overlap_fraction", "0")) or 0.0

    # mean_hotspot_confidence: 1=low, 2=medium, 3=high → normalize to 0-1
    raw_conf = safe_float(row.get("mean_hotspot_confidence", "0")) or 0.0
    conf_norm = min(max((raw_conf - 1.0) / 2.0, 0.0), 1.0) if raw_conf > 0 else 0.0

    n_total = safe_float(row.get("n_hotspot_residues", "1")) or 1.0
    n_total = max(n_total, 1.0)

    n_robust = safe_float(row.get("n_robust_hotspots", "0")) or 0.0
    robust_frac = min(n_robust / n_total, 1.0)

    n_both = safe_float(row.get("n_method_agreement_both", "0")) or 0.0
    both_frac = min(n_both / n_total, 1.0)

    score = (0.40 * overlap +
             0.25 * conf_norm +
             0.20 * robust_frac +
             0.15 * both_frac)
    return round(min(max(score, 0.0), 1.0), 4)


# -- Axis A2: Druggability Confidence --

TIER_SCORES = {"tier_1": 1.0, "tier_2": 0.6, "tier_3": 0.2}
CONFIDENCE_SCORES = {"high": 1.0, "medium": 0.6, "low": 0.2}


def compute_a2_druggability(row: dict) -> float:
    """Compute druggability confidence score (0-1).

    Submetrics:
    - overall_druggability_tier (tier_1/2/3 → 1.0/0.6/0.2): 40%
    - druggability_confidence (high/medium/low → 1.0/0.6/0.2): 30%
    - best_proposal_score (0-1): raw fpocket/P2Rank score: 30%
    """
    tier = TIER_SCORES.get(row.get("overall_druggability_tier", ""), 0.0)
    conf = CONFIDENCE_SCORES.get(row.get("druggability_confidence", ""), 0.0)
    raw_score = safe_float(row.get("best_proposal_score", "0")) or 0.0

    score = 0.40 * tier + 0.30 * conf + 0.30 * raw_score
    return round(min(max(score, 0.0), 1.0), 4)


# -- Axis A3: Perturbation Relevance --

RELATIONSHIP_SCORES = {
    "orthosteric_candidate": 1.0,
    "rim_candidate": 0.7,
    "allosteric_candidate": 0.5,
    "low_relevance_candidate": 0.1,
}

SUPPORT_SCORES = {
    "strong": 1.0,
    "moderate": 0.6,
    "weak": 0.3,
    "none": 0.0,
    "pending_strong": 0.8,
    "pending_moderate": 0.5,
    "pending_weak": 0.2,
}


def compute_a3_perturbation_relevance(row: dict) -> float:
    """Compute perturbation relevance score (0-1).

    Submetrics:
    - relationship_class (categorical → 0-1): mechanistic proximity: 45%
    - hotspot_overlap_count (0-7 → 0-1): direct residue evidence: 25%
    - ligand_support_strength (categorical → 0-1): docking support: 20%
    - pose_support_count (0+ → capped 0-1): volume of evidence: 10%
    """
    rel = RELATIONSHIP_SCORES.get(row.get("relationship_class", ""), 0.0)

    overlap_count = safe_float(row.get("hotspot_overlap_count", "0")) or 0.0
    n_hotspots = safe_float(row.get("n_hotspot_residues", "1")) or 1.0
    n_hotspots = max(n_hotspots, 1.0)
    overlap_norm = min(overlap_count / n_hotspots, 1.0)

    support = SUPPORT_SCORES.get(
        row.get("ligand_support_strength", "none"), 0.0)

    pose_count = safe_float(row.get("pose_support_count", "0")) or 0.0
    # Cap at 10 poses for normalization
    pose_norm = min(pose_count / 10.0, 1.0)

    score = (0.45 * rel +
             0.25 * overlap_norm +
             0.20 * support +
             0.10 * pose_norm)
    return round(min(max(score, 0.0), 1.0), 4)


# -- Axis A4: State Robustness / Accessibility --

STATE_CLASS_SCORES = {
    "robust_pocket": 1.0,
    "shifted_pocket": 0.7,
    "state_specific_pocket": 0.4,
}


def compute_a4_state_robustness(row: dict) -> float:
    """Compute state robustness score (0-1).

    Submetrics:
    - state_class (categorical → 0-1): cross-state persistence: 50%
    - n_states_matched (1-3 → 0-1): number of states: 25%
    - diversity_verdict (well_distributed=1, else 0.5): search quality: 15%
    - coverage_fraction (0-1): pocket coverage: 10%

    Note: state_specific is scored 0.4 not 0.0 — these sites are flagged
    but not eliminated, per PRD requirement.
    """
    state = STATE_CLASS_SCORES.get(
        row.get("state_class", ""), 0.3)

    n_matched = safe_float(row.get("n_states_matched", "1")) or 1.0
    # Normalize: 1 state = 0.33, 2 = 0.67, 3 = 1.0
    n_matched_norm = min(n_matched / 3.0, 1.0)

    verdict = row.get("diversity_verdict", "")
    if verdict == "well_distributed":
        diversity = 1.0
    elif verdict in ("over_concentrated", "under_sampled", "single_pocket"):
        diversity = 0.4
    elif verdict:
        diversity = 0.5
    else:
        # No docking evidence (skipped pocket)
        diversity = 0.0

    coverage = safe_float(row.get("coverage_fraction", "0")) or 0.0

    score = (0.50 * state +
             0.25 * n_matched_norm +
             0.15 * diversity +
             0.10 * coverage)
    return round(min(max(score, 0.0), 1.0), 4)


# ---------------------------------------------------------------------------
# Score all candidates
# ---------------------------------------------------------------------------

def compute_axis_scores(normalized: List[dict]) -> List[dict]:
    """Compute all 4 axis scores for each normalized evidence row."""
    scored = []
    for row in normalized:
        result = {
            "pocket_id": row["pocket_id"],
            "receptor_id": row["receptor_id"],
            "ligand_id": row.get("ligand_id", ""),
            # Axis scores
            "A1_ppi_interface": compute_a1_ppi_interface(row),
            "A2_druggability": compute_a2_druggability(row),
            "A3_perturbation_relevance": compute_a3_perturbation_relevance(row),
            "A4_state_robustness": compute_a4_state_robustness(row),
        }
        # Preserve submetrics
        for col in SCORE_COLUMNS:
            if col not in result:
                result[col] = row.get(col, "")
        scored.append(result)
    return scored


# ---------------------------------------------------------------------------
# Framework documentation
# ---------------------------------------------------------------------------

def build_framework_doc(
    axis_defs: List[dict],
    scored: List[dict],
) -> List[str]:
    """Build the score framework documentation."""
    lines = [
        "# Phase 4 Score Framework",
        "",
        "## 1. Design Principles",
        "",
        "This framework scores candidate sites on **perturbation relevance** "
        "rather than simple ligand affinity. The goal is to answer: "
        "*Which sites are most plausible for disrupting MYO1D attachment "
        "to EGFR, and by what mechanism?*",
        "",
        "Key design choices:",
        "- **4 independent axes** prevent any single metric from dominating",
        "- **Raw submetrics preserved** for reviewability",
        "- **Affinity is one input, not the score** — high affinity at an "
        "irrelevant site does not outrank moderate affinity at an "
        "orthosteric site",
        "- **State-specific sites are flagged, not eliminated**",
        "",
        "## 2. Axis Definitions",
        "",
    ]

    total_weight = sum(a["weight"] for a in axis_defs)

    for ax in axis_defs:
        lines.extend([
            f"### {ax['axis_id']}: {ax['axis_name']} "
            f"(weight: {ax['weight']:.0%})",
            "",
            f"**Interpretation**: {ax['description']}",
            "",
            f"- Range: {ax['range']}",
            f"- High score means: {ax['high_means']}",
            f"- Submetrics: {ax['submetrics']}",
            "",
        ])

    lines.extend([
        f"**Total weight**: {total_weight:.0%}",
        "",
        "## 3. Scoring Rules",
        "",
        "Each axis is computed as a weighted sum of its submetrics, "
        "each normalized to [0, 1]. Categorical variables are mapped "
        "to ordinal scores (e.g., tier_1 = 1.0, tier_2 = 0.6, "
        "tier_3 = 0.2). The final perturbation score (TG 4.3) will "
        "combine axis scores using the defined weights.",
        "",
        "### Categorical Mappings",
        "",
        "| Variable | Value | Score |",
        "|----------|-------|-------|",
    ])

    for label, score in sorted(TIER_SCORES.items()):
        lines.append(f"| druggability_tier | {label} | {score} |")
    for label, score in sorted(CONFIDENCE_SCORES.items()):
        lines.append(f"| druggability_confidence | {label} | {score} |")
    for label, score in sorted(RELATIONSHIP_SCORES.items()):
        lines.append(f"| relationship_class | {label} | {score} |")
    for label, score in sorted(SUPPORT_SCORES.items()):
        lines.append(f"| ligand_support | {label} | {score} |")
    for label, score in sorted(STATE_CLASS_SCORES.items()):
        lines.append(f"| state_class | {label} | {score} |")

    lines.append("")

    # Current scores preview
    lines.extend([
        "## 4. Current Axis Scores",
        "",
        "| Pocket | Ligand | A1 PPI | A2 Drug | A3 Perturb | A4 State |",
        "|--------|--------|--------|---------|------------|----------|",
    ])
    for s in scored:
        lig = s["ligand_id"] or "(none)"
        lines.append(
            f"| {s['pocket_id']} | {lig} | "
            f"{s['A1_ppi_interface']:.3f} | {s['A2_druggability']:.3f} | "
            f"{s['A3_perturbation_relevance']:.3f} | "
            f"{s['A4_state_robustness']:.3f} |"
        )
    lines.append("")

    # Caveats
    lines.extend([
        "## 5. Caveats",
        "",
        "1. **Pre-execution affinity**: All ligand support levels are "
        "provisional (`pending_*`). Re-score after Vina execution.",
        "2. **Single receptor state**: Only 3GT8_raw has pocket data. "
        "State robustness axis has limited discrimination.",
        "3. **Single pocket tool**: Only fpocket. P2Rank integration will "
        "improve druggability axis.",
        "4. **Weights are initial**: Axis weights can be tuned after "
        "expert review of first-pass results.",
        "",
        "---",
        "",
        "Generated by `egfr_pipeline.phase4.score_framework`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_score_framework(
    output_dir: Path = PHASE4_OUTPUT_DIR,
) -> Tuple[Path, Path]:
    """Run full TG 4.1 pipeline.

    Returns (axis_table_path, framework_doc_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load normalized evidence from TG 4.0
    norm_path = output_dir / "phase4_evidence_normalized.csv"
    if not norm_path.exists():
        raise FileNotFoundError(
            f"Normalized evidence not found: {norm_path}\n"
            "Run TG 4.0 (evidence_ingestion) first."
        )

    normalized = load_csv(norm_path)
    print(f"  Loaded {len(normalized)} evidence rows")

    # Write axis definition table
    axis_path = output_dir / "phase4_axis_definition_table.csv"
    _write_csv(axis_path, AXIS_TABLE_COLUMNS, AXIS_DEFINITIONS)
    print(f"  Wrote axis definitions: {len(AXIS_DEFINITIONS)} axes")

    # Compute axis scores
    scored = compute_axis_scores(normalized)
    print(f"  Computed axis scores for {len(scored)} candidates")

    # Preview scores
    for s in scored:
        lig = s["ligand_id"] or "(none)"
        print(f"    {s['pocket_id']} / {lig}: "
              f"A1={s['A1_ppi_interface']:.3f} "
              f"A2={s['A2_druggability']:.3f} "
              f"A3={s['A3_perturbation_relevance']:.3f} "
              f"A4={s['A4_state_robustness']:.3f}")

    # Write scored table
    score_path = output_dir / "phase4_axis_scores.csv"
    _write_csv(score_path, SCORE_COLUMNS, scored)

    # Build framework documentation
    doc_lines = build_framework_doc(AXIS_DEFINITIONS, scored)
    doc_path = output_dir / "phase4_score_framework.md"
    with open(doc_path, "w", encoding="utf-8") as f:
        f.write("\n".join(doc_lines))

    return axis_path, doc_path


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
        description="Phase 4 TG 4.1: Multi-Axis Score Framework"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE4_OUTPUT_DIR,
        help="Phase 4 output directory",
    )
    args = parser.parse_args()

    print("Phase 4 — Task Group 4.1: Multi-Axis Score Framework")
    print(f"  Output: {args.output_dir}")
    print()

    axis_path, doc_path = run_score_framework(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 4.1 COMPLETE")
    print(f"{'='*60}")
    print(f"  Axis definitions: {axis_path}")
    print(f"  Framework doc:    {doc_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
