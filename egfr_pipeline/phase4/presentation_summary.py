#!/usr/bin/env python3
"""Phase 4 Task Group 4.7: Presentation-Ready Summary Layer.

Generates compact presentation tables and shortlists separated by
mechanism, with confidence levels and rationale notes.

Tasks:
  4.7.1 - Generate compact presentation tables
  4.7.2 - Separate shortlists by mechanism
  4.7.3 - Include confidence levels and rationale notes

Inputs:
  - output/phase4_perturbation/phase4_final_review_table.csv
  - output/phase4_perturbation/phase4_expanded_evidence_table.csv

Outputs:
  - output/phase4_perturbation/phase4_presentation_shortlist.csv
  - output/phase4_perturbation/phase4_top_candidates_brief.md

Usage:
    python -m egfr_pipeline.phase4.presentation_summary [--output_dir ...]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

from egfr_pipeline.csv_utils import load_csv

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
PHASE4_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase4_perturbation"

# ---------------------------------------------------------------------------
# Output schemas
# ---------------------------------------------------------------------------

SHORTLIST_COLUMNS = [
    "rank",
    "pocket_id",
    "receptor_id",
    "mechanistic_class",
    "mechanistic_short",
    "perturbation_score",
    "adjusted_confidence",
    "top_ligand",
    "top_ligand_support",
    "hotspot_overlap",
    "druggability_tier",
    "state_accessibility",
    "rationale",
]


# ---------------------------------------------------------------------------
# Build pocket-level shortlist (deduplicate ligands)
# ---------------------------------------------------------------------------

def build_shortlist(review: List[dict], expanded: List[dict]) -> List[dict]:
    """Build pocket-level shortlist, picking best ligand per pocket."""
    # Index expanded by (pocket, ligand) for extra fields
    exp_map = {}
    for e in expanded:
        exp_map[(e["pocket_id"], e.get("ligand_id", ""))] = e

    # Group by pocket, pick best ligand (highest score, then alphabetical)
    by_pocket: Dict[str, List[dict]] = defaultdict(list)
    for r in review:
        by_pocket[r["pocket_id"]].append(r)

    shortlist = []
    for pid in sorted(by_pocket.keys(),
                      key=lambda p: min(int(r["rank"]) for r in by_pocket[p])):
        rows = by_pocket[pid]
        # Best row = lowest rank
        best = min(rows, key=lambda r: int(r["rank"]))
        # All ligands for this pocket
        ligands_with_support = [
            (r["ligand_id"], r.get("ligand_support_strength", "none"))
            for r in rows if r.get("ligand_id")
        ]

        if ligands_with_support:
            top_lig, top_support = ligands_with_support[0]
            n_ligands = len(ligands_with_support)
        else:
            top_lig, top_support = "", "none"
            n_ligands = 0

        rationale = _build_rationale(best, n_ligands)

        short_class = {
            "orthosteric_disruptor_candidate": "orthosteric",
            "interface_rim_modulator_candidate": "rim",
            "allosteric_modulator_candidate": "allosteric",
            "ligandable_but_ppi_irrelevant_candidate": "irrelevant",
            "uncertain_mechanism_candidate": "uncertain",
        }.get(best.get("mechanistic_class", ""), "unknown")

        shortlist.append({
            "rank": best["rank"],
            "pocket_id": pid,
            "receptor_id": best["receptor_id"],
            "mechanistic_class": best["mechanistic_class"],
            "mechanistic_short": short_class,
            "perturbation_score": best["perturbation_score"],
            "adjusted_confidence": best.get("adjusted_confidence", ""),
            "top_ligand": top_lig,
            "top_ligand_support": top_support,
            "hotspot_overlap": best.get("hotspot_overlap_count", ""),
            "druggability_tier": best.get("overall_druggability_tier", ""),
            "state_accessibility": best.get("accessibility_label", ""),
            "rationale": rationale,
        })

    return shortlist


def _build_rationale(row: dict, n_ligands: int) -> str:
    """Build a concise human-readable rationale string."""
    parts = []

    mech = row.get("mechanistic_class", "")
    overlap = row.get("hotspot_overlap_count", "0")
    tier = row.get("overall_druggability_tier", "")
    support = row.get("ligand_support_strength", "none")

    if "orthosteric" in mech:
        parts.append(f"Overlaps {overlap} PPI hotspot(s)")
        if tier:
            parts.append(f"{tier} druggability")
        if n_ligands > 0:
            parts.append(f"{n_ligands} ligand(s) tested")
    elif "rim" in mech:
        parts.append(f"Rim contact with {overlap} hotspot(s)")
        if n_ligands > 0:
            parts.append(f"{n_ligands} ligand(s)")
    elif "allosteric" in mech:
        parts.append("Distant from PPI patch")
        if tier:
            parts.append(f"{tier} druggability")
    elif "irrelevant" in mech:
        parts.append("No PPI hotspot overlap")
        parts.append("Control site only")
    else:
        parts.append("Insufficient evidence for classification")

    if "pending" in support:
        parts.append("affinity pending")

    return "; ".join(parts)


# ---------------------------------------------------------------------------
# Brief markdown summary
# ---------------------------------------------------------------------------

def build_brief(shortlist: List[dict], review: List[dict]) -> List[str]:
    """Build the top candidates brief."""
    lines = [
        "# Phase 4 Top Candidates Brief",
        "",
        "## Overview",
        "",
    ]

    n_ortho = sum(1 for s in shortlist if "orthosteric" in s["mechanistic_short"])
    n_rim = sum(1 for s in shortlist if s["mechanistic_short"] == "rim")
    n_allo = sum(1 for s in shortlist if s["mechanistic_short"] == "allosteric")
    n_irr = sum(1 for s in shortlist if s["mechanistic_short"] == "irrelevant")
    n_unc = sum(1 for s in shortlist if s["mechanistic_short"] == "uncertain")

    lines.extend([
        f"**{len(shortlist)} candidate pockets** evaluated for MYO1D "
        f"perturbation relevance:",
        "",
    ])
    if n_ortho:
        lines.append(f"- **{n_ortho} orthosteric disruptor(s)** — "
                     "direct PPI interface overlap")
    if n_rim:
        lines.append(f"- **{n_rim} rim modulator(s)** — "
                     "peripheral interface contact")
    if n_allo:
        lines.append(f"- **{n_allo} allosteric modulator(s)** — "
                     "conformational perturbation")
    if n_irr:
        lines.append(f"- **{n_irr} PPI-irrelevant site(s)** — "
                     "druggable but not perturbation candidates")
    if n_unc:
        lines.append(f"- **{n_unc} uncertain** — awaiting more evidence")
    lines.append("")

    # Shortlist by mechanism
    by_mech = defaultdict(list)
    for s in shortlist:
        by_mech[s["mechanistic_short"]].append(s)

    mech_order = ["orthosteric", "rim", "allosteric", "irrelevant", "uncertain"]
    mech_titles = {
        "orthosteric": "Orthosteric Disruptors (Priority Candidates)",
        "rim": "Interface Rim Modulators",
        "allosteric": "Allosteric Modulators",
        "irrelevant": "Ligandable but PPI-Irrelevant (Controls)",
        "uncertain": "Uncertain Mechanism",
    }

    for mech in mech_order:
        rows = by_mech.get(mech, [])
        if not rows:
            continue

        title = mech_titles[mech]
        lines.extend([
            f"## {title}",
            "",
        ])

        for s in rows:
            score = s["perturbation_score"]
            conf = s["adjusted_confidence"]
            lines.extend([
                f"### {s['pocket_id']}",
                "",
                f"- **Score**: {score} | **Confidence**: {conf}",
                f"- **Druggability**: {s['druggability_tier']}",
                f"- **Hotspot overlap**: {s['hotspot_overlap']} residue(s)",
                f"- **Top ligand**: {s['top_ligand'] or 'none'} "
                f"({s['top_ligand_support']})",
                f"- **State**: {s['state_accessibility']}",
                f"- **Rationale**: {s['rationale']}",
                "",
            ])

    # Ligand coverage matrix
    lines.extend([
        "## Ligand Coverage Matrix",
        "",
    ])
    ligands = sorted(set(r["ligand_id"] for r in review if r["ligand_id"]))
    pockets = sorted(set(r["pocket_id"] for r in review))

    if ligands:
        header = "| Pocket | " + " | ".join(ligands) + " |"
        sep = "|" + "|".join(["--------"] + ["---"] * len(ligands)) + "|"
        lines.extend([header, sep])

        review_map = {}
        for r in review:
            review_map[(r["pocket_id"], r["ligand_id"])] = r

        for pid in pockets:
            cells = []
            for lig in ligands:
                r = review_map.get((pid, lig))
                if r:
                    support = r.get("ligand_support_strength", "—")
                    cells.append(support)
                else:
                    cells.append("—")
            lines.append(f"| {pid} | " + " | ".join(cells) + " |")
        lines.append("")

    # Status note
    lines.extend([
        "## Status",
        "",
        "- All ligand support levels are **provisional** (`pending_*`) — "
        "awaiting server-side Vina execution",
        "- State robustness assessment is **preliminary** — only 1 of 3 "
        "receptor states has pocket data",
        "- Scores will be updated after docking execution and multi-state "
        "pocket detection",
        "",
        "---",
        "",
        "Generated by `egfr_pipeline.phase4.presentation_summary`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_presentation_summary(
    output_dir: Path = PHASE4_OUTPUT_DIR,
) -> Tuple[Path, Path]:
    """Run full TG 4.7 pipeline.

    Returns (shortlist_path, brief_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    review = load_csv(output_dir / "phase4_final_review_table.csv")
    expanded = load_csv(output_dir / "phase4_expanded_evidence_table.csv")

    if not review:
        raise FileNotFoundError(
            "Review table not found. Run TG 4.5 first.")

    print(f"  Loaded {len(review)} review rows")

    shortlist = build_shortlist(review, expanded)
    print(f"  Shortlist: {len(shortlist)} pocket(s)")

    for s in shortlist:
        print(f"    {s['pocket_id']}: {s['mechanistic_short']} "
              f"(score={s['perturbation_score']})")

    brief_lines = build_brief(shortlist, review)

    # Write outputs
    shortlist_path = output_dir / "phase4_presentation_shortlist.csv"
    _write_csv(shortlist_path, SHORTLIST_COLUMNS, shortlist)

    brief_path = output_dir / "phase4_top_candidates_brief.md"
    with open(brief_path, "w", encoding="utf-8") as f:
        f.write("\n".join(brief_lines))

    return shortlist_path, brief_path


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
        description="Phase 4 TG 4.7: Presentation-Ready Summary"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE4_OUTPUT_DIR,
        help="Phase 4 output directory",
    )
    args = parser.parse_args()

    print("Phase 4 — Task Group 4.7: Presentation-Ready Summary")
    print(f"  Output: {args.output_dir}")
    print()

    shortlist_path, brief_path = run_presentation_summary(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 4.7 COMPLETE")
    print(f"{'='*60}")
    print(f"  Shortlist: {shortlist_path}")
    print(f"  Brief:     {brief_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
