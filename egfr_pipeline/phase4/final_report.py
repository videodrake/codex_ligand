#!/usr/bin/env python3
"""Phase 4 Task Group 4.6: Final Report and Interpretation Guide.

Generates the integrated Phase 4 report presenting final ranked candidates
by mechanistic class, explaining the scoring methodology, and including
validation/caution sections.

Tasks:
  4.6.1 - Present final ranked candidates by mechanistic class
  4.6.2 - Explain why affinity alone wasn't final criterion
  4.6.3 - Include validation/caution sections

Inputs:
  - output/phase4_perturbation/phase4_final_review_table.csv
  - output/phase4_perturbation/phase4_expanded_evidence_table.csv
  - output/phase4_perturbation/phase4_axis_definition_table.csv
  - output/phase4_perturbation/phase4_state_interpretation.csv
  - output/phase4_perturbation/phase4_evidence_validation.md (verdict)

Outputs:
  - output/phase4_perturbation/integrated_phase4_report.md

Usage:
    python -m egfr_pipeline.phase4.final_report [--output_dir ...]
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
# Report builder
# ---------------------------------------------------------------------------

def build_report(
    review: List[dict],
    expanded: List[dict],
    axes: List[dict],
    interpretation: List[dict],
) -> List[str]:
    """Build the integrated Phase 4 report."""
    lines = []

    # ── Section 1: Executive Summary ──
    pockets = sorted(set(r["pocket_id"] for r in review))
    ligands = sorted(set(r["ligand_id"] for r in review if r["ligand_id"]))
    receptors = sorted(set(r["receptor_id"] for r in review))
    by_class = defaultdict(list)
    for r in review:
        by_class[r["mechanistic_class"]].append(r)

    lines.extend([
        "# Integrated Phase 4 Report: Perturbation Relevance Scoring",
        "",
        "## 1. Executive Summary",
        "",
        f"- **Receptor states**: {len(receptors)} ({', '.join(receptors)})",
        f"- **Candidate pockets**: {len(pockets)}",
        f"- **Ligands evaluated**: {len(ligands)} ({', '.join(ligands)})",
        f"- **Total ranked entries**: {len(review)}",
        "",
        "### Mechanistic Class Distribution",
        "",
        "| Class | Count | Pockets |",
        "|-------|-------|---------|",
    ])
    for cls in sorted(by_class.keys()):
        rows = by_class[cls]
        cls_pockets = sorted(set(r["pocket_id"] for r in rows))
        lines.append(f"| {cls} | {len(rows)} | {', '.join(cls_pockets)} |")
    lines.append("")

    # ── Section 2: Why Not Affinity Alone ──
    lines.extend([
        "## 2. Why Affinity Alone Is Not the Final Criterion",
        "",
        "This pipeline answers a specific biological question: "
        "*Which candidate sites are most plausible for disrupting "
        "MYO1D attachment to EGFR, and by what mechanism?*",
        "",
        "A high-affinity docking score at a site with no connection to "
        "the MYO1D interface is chemically interesting but biologically "
        "irrelevant to the perturbation goal. Conversely, a moderate-affinity "
        "hit at a site that directly overlaps the PPI patch is a stronger "
        "perturbation candidate.",
        "",
        "The 4-axis scoring framework prevents affinity domination by:",
        "",
        "1. **Weighting PPI evidence (A1) and perturbation relevance (A3) "
        "at 60% combined** — biological relevance outweighs chemical metrics",
        "2. **Capping affinity-influenced axes** for irrelevant sites — "
        "high druggability at a low-relevance site cannot outrank an "
        "orthosteric candidate",
        "3. **Preserving mechanistic class labels** — the ranking is "
        "interpretable, not a black-box number",
        "",
    ])

    # ── Section 3: Scoring Framework ──
    lines.extend([
        "## 3. Scoring Framework",
        "",
        "| Axis | Weight | Meaning |",
        "|------|--------|---------|",
    ])
    for ax in axes:
        lines.append(
            f"| {ax['axis_name']} | {float(ax['weight']):.0%} | "
            f"{ax['high_means']} |")
    lines.append("")

    lines.extend([
        "Each axis score is [0, 1]. The perturbation score is their "
        "weighted sum, with an affinity cap applied to irrelevant/uncertain "
        "sites.",
        "",
    ])

    # ── Section 4: Ranked Candidates by Mechanistic Class ──
    lines.extend([
        "## 4. Ranked Candidates by Mechanistic Class",
        "",
    ])

    class_order = [
        "orthosteric_disruptor_candidate",
        "interface_rim_modulator_candidate",
        "allosteric_modulator_candidate",
        "ligandable_but_ppi_irrelevant_candidate",
        "uncertain_mechanism_candidate",
    ]
    class_labels = {
        "orthosteric_disruptor_candidate": "Orthosteric Disruptors",
        "interface_rim_modulator_candidate": "Interface Rim Modulators",
        "allosteric_modulator_candidate": "Allosteric Modulators",
        "ligandable_but_ppi_irrelevant_candidate": "Ligandable but PPI-Irrelevant",
        "uncertain_mechanism_candidate": "Uncertain Mechanism",
    }

    for cls in class_order:
        rows = by_class.get(cls, [])
        if not rows:
            continue

        label = class_labels.get(cls, cls)
        lines.extend([
            f"### 4.{class_order.index(cls)+1}. {label}",
            "",
            "| Rank | Pocket | Ligand | Score | A1 | A2 | A3 | A4 | "
            "Confidence | Overlap |",
            "|------|--------|--------|-------|----|----|----|----|"
            "------------|---------|",
        ])
        for r in sorted(rows, key=lambda x: int(x["rank"])):
            lig = r["ligand_id"] or "(none)"
            lines.append(
                f"| {r['rank']} | {r['pocket_id']} | {lig} | "
                f"{r['perturbation_score']} | "
                f"{_fmt(r.get('A1_ppi_interface'))} | "
                f"{_fmt(r.get('A2_druggability'))} | "
                f"{_fmt(r.get('A3_perturbation_relevance'))} | "
                f"{_fmt(r.get('A4_state_robustness'))} | "
                f"{r.get('adjusted_confidence', '')} | "
                f"{r.get('hotspot_overlap_count', '')} |"
            )
        lines.append("")

    # ── Section 5: State Robustness Assessment ──
    lines.extend([
        "## 5. State Robustness Assessment",
        "",
    ])

    # Deduplicate by pocket
    seen_pockets = set()
    for r in interpretation:
        pid = r["pocket_id"]
        if pid in seen_pockets:
            continue
        seen_pockets.add(pid)
        lines.extend([
            f"**{pid}**",
            f"- State class: {r.get('state_class', 'unknown')}",
            f"- Interpretation: {r.get('state_interpretation', '')}",
            f"- Accessibility: {r.get('accessibility_label', '')}",
            f"- States matched: {r.get('n_states_matched', '')}",
        ])
        caveat = r.get("state_caveat", "")
        if caveat:
            lines.append(f"- Caveat: {caveat}")
        lines.append("")

    # ── Section 6: Evidence Provenance Summary ──
    lines.extend([
        "## 6. Evidence Provenance Summary",
        "",
        "Each candidate's score traces back to concrete upstream evidence:",
        "",
        "| Phase | Data Source | Key Fields |",
        "|-------|-----------|------------|",
        "| Phase 1 | PPI patch reference | hotspot residues, robustness, "
        "method agreement |",
        "| Phase 2 | Pocket proposal | relationship class, druggability tier, "
        "proposal score |",
        "| Phase 3 | Docking evidence | ligand support, pose count, "
        "diversity verdict |",
        "| Phase 4 | Scoring + classification | axis scores, mechanistic class, "
        "state interpretation |",
        "",
        "Full provenance is available in `phase4_expanded_evidence_table.csv` "
        f"({len(expanded)} rows, {len(expanded[0]) if expanded else 0} columns).",
        "",
    ])

    # ── Section 7: Validation and Caveats ──
    lines.extend([
        "## 7. Validation and Caveats",
        "",
        "### Data Completeness",
        "",
    ])

    # Check provenance flags from expanded table
    p1_loaded = all(r.get("phase1_loaded") == "True" for r in expanded)
    p2_loaded = all(r.get("phase2_loaded") == "True" for r in expanded)
    p3_count = sum(1 for r in expanded if r.get("phase3_loaded") == "True")

    lines.extend([
        f"- Phase 1 PPI evidence: {'Complete' if p1_loaded else 'Partial/Missing'}",
        f"- Phase 2 pocket evidence: {'Complete' if p2_loaded else 'Partial/Missing'}",
        f"- Phase 3 docking evidence: {p3_count}/{len(expanded)} candidates "
        f"({'all' if p3_count == len(expanded) else 'partial'})",
        "",
        "### Known Limitations",
        "",
        "1. **Single receptor state**: Only 3GT8_raw has pocket data. "
        "Scores will change when cl38_48 and cl85_100 fpocket/P2Rank results "
        "are integrated. State robustness axis (A4) currently has minimal "
        "discrimination power.",
        "",
        "2. **Pre-execution affinity**: Vina docking jobs are prepared but "
        "not yet executed. All ligand support levels are provisional "
        "(`pending_*`). Re-run Phase 4 after server-side Vina execution.",
        "",
        "3. **Single pocket detection tool**: Only fpocket results available. "
        "P2Rank integration will improve druggability axis (A2) through "
        "multi-tool consensus.",
        "",
        "4. **No FTMap data**: Fragment hotspot mapping not yet available. "
        "When integrated, this will strengthen druggability confidence.",
        "",
        "5. **Axis weights are initial**: The 30/25/30/15 weighting is a "
        "principled starting point but can be tuned after expert review.",
        "",
        "### Recommended Next Steps",
        "",
        "1. Execute Vina docking on HPC server (`run_phase3_docking.sh`)",
        "2. Re-run Phase 3 TG 3.3-3.7 to populate affinity data",
        "3. Re-run Phase 4 to update scores with real affinity evidence",
        "4. Run fpocket/P2Rank on cl38_48 and cl85_100 receptor states",
        "5. Re-run Phase 2-4 with multi-state data for robust state assessment",
        "",
    ])

    # ── Section 8: File Inventory ──
    lines.extend([
        "## 8. File Inventory",
        "",
        "| File | Description |",
        "|------|-------------|",
        "| phase4_evidence_normalized.csv | Merged Phase 1-3 evidence (TG 4.0) |",
        "| phase4_evidence_validation.md | Evidence validation report (TG 4.0) |",
        "| phase4_axis_definition_table.csv | Axis definitions (TG 4.1) |",
        "| phase4_axis_scores.csv | Raw axis scores (TG 4.1) |",
        "| phase4_score_framework.md | Score framework documentation (TG 4.1) |",
        "| final_candidate_classes.csv | Mechanistic classifications (TG 4.2) |",
        "| phase4_mechanistic_classification_note.md | Classification logic (TG 4.2) |",
        "| perturbation_axis_scores.csv | Weighted axis scores + rank (TG 4.3) |",
        "| perturbation_candidate_table.csv | Ranked candidate table (TG 4.3) |",
        "| phase4_ranking_method_note.md | Ranking methodology (TG 4.3) |",
        "| phase4_state_interpretation.csv | State robustness interpretation (TG 4.4) |",
        "| phase4_accessibility_note.md | Accessibility documentation (TG 4.4) |",
        "| phase4_final_review_table.csv | Condensed review table (TG 4.5) |",
        "| phase4_expanded_evidence_table.csv | Full provenance table (TG 4.5) |",
        "| integrated_phase4_report.md | This report (TG 4.6) |",
        "",
        "---",
        "",
        "Generated by `egfr_pipeline.phase4.final_report`",
    ])

    return lines


def _fmt(val) -> str:
    """Format a numeric value for table display."""
    try:
        return f"{float(val):.3f}"
    except (ValueError, TypeError):
        return str(val) if val else ""


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_final_report(
    output_dir: Path = PHASE4_OUTPUT_DIR,
) -> Path:
    """Run full TG 4.6 pipeline. Returns report_path."""
    output_dir.mkdir(parents=True, exist_ok=True)

    review = load_csv(output_dir / "phase4_final_review_table.csv")
    expanded = load_csv(output_dir / "phase4_expanded_evidence_table.csv")
    axes = load_csv(output_dir / "phase4_axis_definition_table.csv")
    interpretation = load_csv(output_dir / "phase4_state_interpretation.csv")

    if not review:
        raise FileNotFoundError(
            "Review table not found. Run TG 4.5 first.")

    print(f"  Loaded: {len(review)} review rows, {len(expanded)} expanded, "
          f"{len(axes)} axes, {len(interpretation)} interpretations")

    report_lines = build_report(review, expanded, axes, interpretation)

    report_path = output_dir / "integrated_phase4_report.md"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    return report_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 4 TG 4.6: Final Report and Interpretation Guide"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE4_OUTPUT_DIR,
        help="Phase 4 output directory",
    )
    args = parser.parse_args()

    print("Phase 4 — Task Group 4.6: Final Report")
    print(f"  Output: {args.output_dir}")
    print()

    report_path = run_final_report(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 4.6 COMPLETE")
    print(f"{'='*60}")
    print(f"  Report: {report_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
