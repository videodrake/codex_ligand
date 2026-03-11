#!/usr/bin/env python3
"""Phase 4 Task Group 4.2: Mechanistic Classification Logic.

Assigns mechanistic interpretation labels to candidate sites based on
Phase 2 patch relationships, hotspot overlap, docking evidence, and
state support.

Tasks:
  4.2.1 - Define mechanistic classes
  4.2.2 - Map Phase 2 classes into Phase 4 interpretation
  4.2.3 - Preserve uncertainty (don't force weak evidence into strong labels)

Inputs:
  - output/phase4_perturbation/phase4_axis_scores.csv

Outputs:
  - output/phase4_perturbation/final_candidate_classes.csv
  - output/phase4_perturbation/phase4_mechanistic_classification_note.md

Usage:
    python -m egfr_pipeline.phase4.mechanistic_classification [--output_dir ...]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
PHASE4_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase4_perturbation"

# ---------------------------------------------------------------------------
# Mechanistic class definitions (4.2.1)
# ---------------------------------------------------------------------------

MECHANISTIC_CLASSES = [
    {
        "class_id": "orthosteric_disruptor_candidate",
        "description": (
            "Site directly overlaps the MYO1D attachment patch. "
            "Ligand binding here would sterically block PPI formation."
        ),
        "minimum_evidence": (
            "relationship_class = orthosteric_candidate AND "
            "hotspot_overlap_count >= 2 AND "
            "ligand_support_strength not in (none)"
        ),
    },
    {
        "class_id": "interface_rim_modulator_candidate",
        "description": (
            "Site is at the rim of the PPI interface. Ligand binding "
            "may indirectly weaken MYO1D attachment by perturbing "
            "peripheral contacts."
        ),
        "minimum_evidence": (
            "relationship_class = rim_candidate AND "
            "hotspot_overlap_count >= 1"
        ),
    },
    {
        "class_id": "allosteric_modulator_candidate",
        "description": (
            "Site is spatially distant from the PPI patch but within "
            "the kinase domain. Ligand binding may induce conformational "
            "changes that indirectly affect MYO1D binding."
        ),
        "minimum_evidence": (
            "relationship_class = allosteric_candidate AND "
            "druggability_tier in (tier_1, tier_2)"
        ),
    },
    {
        "class_id": "ligandable_but_ppi_irrelevant_candidate",
        "description": (
            "Site is druggable but has no mechanistic link to MYO1D "
            "attachment. Useful for selectivity controls but not "
            "perturbation candidates."
        ),
        "minimum_evidence": (
            "relationship_class = low_relevance_candidate AND "
            "hotspot_overlap_count = 0"
        ),
    },
    {
        "class_id": "uncertain_mechanism_candidate",
        "description": (
            "Evidence is insufficient or contradictory to assign a "
            "confident mechanistic label. Retained for future evidence."
        ),
        "minimum_evidence": "Default when no other class matches",
    },
]

# Output schema
CLASSIFICATION_COLUMNS = [
    "pocket_id",
    "receptor_id",
    "ligand_id",
    "mechanistic_class",
    "classification_basis",
    "classification_confidence",
    # Axis scores (carried forward)
    "A1_ppi_interface",
    "A2_druggability",
    "A3_perturbation_relevance",
    "A4_state_robustness",
    # Key evidence
    "relationship_class",
    "hotspot_overlap_count",
    "hotspot_overlap_fraction",
    "overall_druggability_tier",
    "ligand_support_strength",
    "state_class",
    "n_states_matched",
    "pose_support_count",
]


# ---------------------------------------------------------------------------
# 4.2.2: Classification logic
# ---------------------------------------------------------------------------

def classify_candidate(row: dict) -> Tuple[str, str, str]:
    """Classify a single candidate into a mechanistic class.

    Returns (mechanistic_class, classification_basis, confidence).
    """
    rel_class = row.get("relationship_class", "")
    overlap = int(_safe_float(row.get("hotspot_overlap_count", "0")))
    overlap_frac = _safe_float(row.get("hotspot_overlap_fraction", "0"))
    tier = row.get("overall_druggability_tier", "")
    support = row.get("ligand_support_strength", "none")
    pose_count = int(_safe_float(row.get("pose_support_count", "0")))

    # Has any docking evidence (even pending)?
    has_docking = support not in ("none", "")

    # --- Orthosteric disruptor ---
    if rel_class == "orthosteric_candidate" and overlap >= 2:
        basis_parts = [f"orthosteric_relationship", f"overlap={overlap}"]
        if has_docking:
            basis_parts.append(f"ligand_support={support}")
            conf = _ortho_confidence(overlap_frac, support, pose_count)
        else:
            conf = "low"
            basis_parts.append("no_docking_evidence")
        return (
            "orthosteric_disruptor_candidate",
            "; ".join(basis_parts),
            conf,
        )

    # --- Interface rim modulator ---
    if rel_class == "rim_candidate" and overlap >= 1:
        basis_parts = [f"rim_relationship", f"overlap={overlap}"]
        if has_docking:
            basis_parts.append(f"ligand_support={support}")
        conf = "medium" if overlap >= 2 else "low"
        return (
            "interface_rim_modulator_candidate",
            "; ".join(basis_parts),
            conf,
        )

    # --- Allosteric modulator ---
    if rel_class == "allosteric_candidate" and tier in ("tier_1", "tier_2"):
        basis_parts = [f"allosteric_relationship", f"tier={tier}"]
        if has_docking:
            basis_parts.append(f"ligand_support={support}")
        conf = "medium" if tier == "tier_1" else "low"
        return (
            "allosteric_modulator_candidate",
            "; ".join(basis_parts),
            conf,
        )

    # --- Ligandable but PPI-irrelevant ---
    if rel_class == "low_relevance_candidate" and overlap == 0:
        basis_parts = [f"low_relevance_relationship", "no_hotspot_overlap"]
        if tier:
            basis_parts.append(f"tier={tier}")
        return (
            "ligandable_but_ppi_irrelevant_candidate",
            "; ".join(basis_parts),
            "high",  # High confidence it's irrelevant
        )

    # --- Edge cases: orthosteric with weak overlap ---
    if rel_class == "orthosteric_candidate" and overlap == 1:
        basis_parts = [
            "orthosteric_relationship",
            f"overlap={overlap} (below threshold=2)",
        ]
        return (
            "uncertain_mechanism_candidate",
            "; ".join(basis_parts),
            "low",
        )

    # --- Default: uncertain ---
    basis_parts = [f"rel={rel_class}", f"overlap={overlap}", f"tier={tier}"]
    return (
        "uncertain_mechanism_candidate",
        "; ".join(basis_parts),
        "low",
    )


def _ortho_confidence(
    overlap_frac: float, support: str, pose_count: int
) -> str:
    """Determine confidence for orthosteric disruptor classification."""
    strong_support = support in ("strong", "pending_strong")
    moderate_support = support in ("moderate", "pending_moderate")

    if overlap_frac >= 0.40 and (strong_support or moderate_support):
        return "high"
    if overlap_frac >= 0.25 or moderate_support:
        return "medium"
    return "low"


# ---------------------------------------------------------------------------
# Classify all candidates
# ---------------------------------------------------------------------------

def classify_all(scored_rows: List[dict]) -> List[dict]:
    """Classify all candidates and produce the classification table."""
    results = []
    for row in scored_rows:
        mech_class, basis, confidence = classify_candidate(row)
        result = {
            "pocket_id": row["pocket_id"],
            "receptor_id": row["receptor_id"],
            "ligand_id": row.get("ligand_id", ""),
            "mechanistic_class": mech_class,
            "classification_basis": basis,
            "classification_confidence": confidence,
        }
        # Carry forward axis scores and key evidence
        for col in CLASSIFICATION_COLUMNS:
            if col not in result:
                result[col] = row.get(col, "")
        results.append(result)
    return results


# ---------------------------------------------------------------------------
# Classification report
# ---------------------------------------------------------------------------

def build_classification_report(
    classified: List[dict],
) -> List[str]:
    """Build the mechanistic classification documentation."""
    lines = [
        "# Phase 4 Mechanistic Classification Note",
        "",
        "## 1. Classification Logic",
        "",
        "Each candidate site is assigned a mechanistic class based on:",
        "- **Patch relationship** (Phase 2): orthosteric / rim / allosteric / low_relevance",
        "- **Hotspot overlap** (Phase 1 x Phase 2): number of PPI hotspot residues in pocket",
        "- **Docking evidence** (Phase 3): ligand support strength",
        "- **Druggability tier** (Phase 2): for allosteric classification",
        "",
        "## 2. Mechanistic Classes",
        "",
        "| Class | Description | Minimum Evidence |",
        "|-------|-------------|-----------------|",
    ]
    for mc in MECHANISTIC_CLASSES:
        lines.append(
            f"| {mc['class_id']} | {mc['description'][:80]}... | "
            f"{mc['minimum_evidence'][:60]}... |"
        )
    lines.append("")

    # Distribution
    by_class = defaultdict(list)
    for c in classified:
        by_class[c["mechanistic_class"]].append(c)

    lines.extend([
        "## 3. Classification Results",
        "",
        "| Mechanistic Class | Count | Pockets | Confidence |",
        "|-------------------|-------|---------|------------|",
    ])
    for mc in MECHANISTIC_CLASSES:
        cid = mc["class_id"]
        rows = by_class.get(cid, [])
        if rows:
            pockets = sorted(set(r["pocket_id"] for r in rows))
            confs = sorted(set(r["classification_confidence"] for r in rows))
            lines.append(
                f"| {cid} | {len(rows)} | "
                f"{', '.join(pockets)} | {', '.join(confs)} |"
            )
    lines.append("")

    # Per-pocket detail
    lines.extend([
        "## 4. Per-Candidate Detail",
        "",
        "| Pocket | Ligand | Class | Confidence | Basis |",
        "|--------|--------|-------|------------|-------|",
    ])
    for c in classified:
        lig = c["ligand_id"] or "(none)"
        lines.append(
            f"| {c['pocket_id']} | {lig} | "
            f"{c['mechanistic_class']} | "
            f"{c['classification_confidence']} | "
            f"{c['classification_basis']} |"
        )
    lines.append("")

    # Uncertainty note (4.2.3)
    uncertain = [c for c in classified
                 if c["mechanistic_class"] == "uncertain_mechanism_candidate"]
    if uncertain:
        lines.extend([
            "## 5. Uncertainty Preservation",
            "",
            f"{len(uncertain)} candidate(s) classified as `uncertain_mechanism_candidate`. "
            "These are retained for future evidence rather than forced into "
            "stronger labels.",
            "",
        ])

    lines.extend([
        "---",
        "",
        "Generated by `egfr_pipeline.phase4.mechanistic_classification`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_mechanistic_classification(
    output_dir: Path = PHASE4_OUTPUT_DIR,
) -> Tuple[Path, Path]:
    """Run full TG 4.2 pipeline.

    Returns (classes_path, note_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load axis scores from TG 4.1
    scores_path = output_dir / "phase4_axis_scores.csv"
    if not scores_path.exists():
        raise FileNotFoundError(
            f"Axis scores not found: {scores_path}\n"
            "Run TG 4.1 (score_framework) first."
        )

    scored = _load_csv(scores_path)
    print(f"  Loaded {len(scored)} scored candidates")

    # Classify
    classified = classify_all(scored)

    # Summary
    by_class = defaultdict(int)
    for c in classified:
        by_class[c["mechanistic_class"]] += 1
    for cls, count in sorted(by_class.items()):
        print(f"    {cls}: {count}")

    # Write outputs
    classes_path = output_dir / "final_candidate_classes.csv"
    _write_csv(classes_path, CLASSIFICATION_COLUMNS, classified)

    report_lines = build_classification_report(classified)
    note_path = output_dir / "phase4_mechanistic_classification_note.md"
    with open(note_path, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    return classes_path, note_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(val: str) -> float:
    try:
        return float(val)
    except (ValueError, TypeError):
        return 0.0


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
        description="Phase 4 TG 4.2: Mechanistic Classification"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE4_OUTPUT_DIR,
        help="Phase 4 output directory",
    )
    args = parser.parse_args()

    print("Phase 4 — Task Group 4.2: Mechanistic Classification")
    print(f"  Output: {args.output_dir}")
    print()

    classes_path, note_path = run_mechanistic_classification(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 4.2 COMPLETE")
    print(f"{'='*60}")
    print(f"  Candidate classes: {classes_path}")
    print(f"  Classification note: {note_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
