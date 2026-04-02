#!/usr/bin/env python3
"""Phase 2 Task Group 2.7: Phase 2 Review Report.

Generates a readable Phase 2 review package summarizing candidate pockets,
patch relationship classes, druggability confidence, cross-state alignment,
and the resulting Phase 3 docking reference set.

Tasks:
  2.7.1 - Candidate pocket summary report
  2.7.2 - Patch relationship summary
  2.7.3 - Downstream handoff summary

Inputs:
  - output/phase2_pockets/candidate_pockets.csv
  - output/phase2_pockets/candidate_pockets_raw.csv
  - output/phase2_pockets/candidate_pocket_provenance.csv
  - output/phase2_pockets/pocket_patch_relationship.csv
  - output/phase2_pockets/pocket_patch_relationship_metrics.csv
  - output/phase2_pockets/druggability_proposal_summary.csv
  - output/phase2_pockets/candidate_pocket_support_flags.csv
  - output/phase2_pockets/candidate_pocket_state_classes.csv
  - output/phase2_pockets/phase3_candidate_pocket_reference.csv
  - output/phase2_pockets/phase2_patch_reference_normalized.csv

Outputs:
  - output/phase2_pockets/phase2_candidate_pocket_report.md

Usage:
    python -m egfr_pipeline.phase2.review_report [--output_dir ...]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set

from egfr_pipeline.csv_utils import load_csv

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

PHASE2_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase2_pockets"
RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_review_report(output_dir: Path) -> Path:
    """Run full TG 2.7 pipeline. Returns report path."""
    # Load all Phase 2 outputs
    raw_pockets = load_csv(output_dir / "candidate_pockets_raw.csv")
    merged_pockets = load_csv(output_dir / "candidate_pockets.csv")
    provenance = load_csv(output_dir / "candidate_pocket_provenance.csv")
    relationships = load_csv(output_dir / "pocket_patch_relationship.csv")
    rel_metrics = load_csv(output_dir / "pocket_patch_relationship_metrics.csv")
    druggability = load_csv(output_dir / "druggability_proposal_summary.csv")
    support_flags = load_csv(output_dir / "candidate_pocket_support_flags.csv")
    state_classes = load_csv(output_dir / "candidate_pocket_state_classes.csv")
    phase3_ref = load_csv(output_dir / "phase3_candidate_pocket_reference.csv")
    patch_ref = load_csv(output_dir / "phase2_patch_reference_normalized.csv")

    print(f"  Loaded: {len(raw_pockets)} raw, {len(merged_pockets)} merged, "
          f"{len(relationships)} relationships, {len(druggability)} druggability, "
          f"{len(state_classes)} state classes, {len(phase3_ref)} Phase 3 export")

    lines = _build_report(
        raw_pockets, merged_pockets, provenance, relationships,
        rel_metrics, druggability, support_flags, state_classes,
        phase3_ref, patch_ref,
    )

    report_path = output_dir / "phase2_candidate_pocket_report.md"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    return report_path


# ---------------------------------------------------------------------------
# Report builder
# ---------------------------------------------------------------------------

def _build_report(
    raw_pockets, merged_pockets, provenance, relationships,
    rel_metrics, druggability, support_flags, state_classes,
    phase3_ref, patch_ref,
) -> List[str]:
    lines: List[str] = []

    # -----------------------------------------------------------------------
    # Header
    # -----------------------------------------------------------------------
    lines.extend([
        "# Phase 2 Review Report: Pocket Proposal and Druggability Mapping",
        "",
        "## 1. Executive Summary",
        "",
    ])

    states_with_pockets = sorted(set(p["receptor_id"] for p in merged_pockets))
    states_missing = sorted(set(RECEPTOR_STATES) - set(states_with_pockets))
    n_hotspots = len(set(
        r.get("patch_residue_id", "") for r in patch_ref
        if r.get("is_hotspot", "").lower() in ("true", "1", "yes")
    ))

    lines.extend([
        f"- **Receptor states with pocket data:** {len(states_with_pockets)} of {len(RECEPTOR_STATES)}",
        f"  - Present: {', '.join(states_with_pockets)}",
    ])
    if states_missing:
        lines.append(f"  - Missing: {', '.join(states_missing)}")
    lines.extend([
        f"- **Phase 1 patch hotspot residues:** {n_hotspots}",
        f"- **Raw pocket proposals:** {len(raw_pockets)}",
        f"- **Merged candidate pockets:** {len(merged_pockets)}",
    ])

    # Priority summary from Phase 3 export
    by_priority: Dict[str, int] = defaultdict(int)
    for r in phase3_ref:
        by_priority[r.get("docking_priority", "unknown")] += 1
    lines.append(f"- **Phase 3 docking priorities:** "
                 f"primary={by_priority.get('primary',0)}, "
                 f"secondary={by_priority.get('secondary',0)}, "
                 f"exploratory={by_priority.get('exploratory',0)}, "
                 f"skip={by_priority.get('skip',0)}")
    lines.append("")

    if states_missing:
        lines.extend([
            "> **Note:** Cross-state alignment and state robustness are provisional.",
            "> Re-run Phase 2 after fpocket/P2Rank results for all states are available.",
            "",
        ])

    # -----------------------------------------------------------------------
    # 2. Input Summary (Task 2.7.1 — proposal sources)
    # -----------------------------------------------------------------------
    lines.extend([
        "## 2. Pocket Proposal Sources",
        "",
    ])

    # Per-tool counts from raw
    by_tool: Dict[str, int] = defaultdict(int)
    by_state_tool: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for p in raw_pockets:
        tool = p.get("proposal_source", "unknown")
        state = p.get("receptor_id", "unknown")
        by_tool[tool] += 1
        by_state_tool[state][tool] += 1

    lines.extend([
        "| Tool | Total Proposals |",
        "|------|----------------|",
    ])
    for tool, count in sorted(by_tool.items()):
        lines.append(f"| {tool} | {count} |")
    lines.append("")

    if by_state_tool:
        lines.extend([
            "### Per-State Breakdown",
            "",
            "| State | Tool | Proposals |",
            "|-------|------|-----------|",
        ])
        for state in sorted(by_state_tool.keys()):
            for tool, count in sorted(by_state_tool[state].items()):
                lines.append(f"| {state} | {tool} | {count} |")
        lines.append("")

    # Merge summary
    lines.extend([
        "### Merge Summary",
        "",
        f"- Raw proposals: {len(raw_pockets)}",
        f"- After merge: {len(merged_pockets)}",
        f"- Reduction: {len(raw_pockets) - len(merged_pockets)} proposals merged",
        "",
    ])

    # -----------------------------------------------------------------------
    # 3. Merged Pocket Catalog (Task 2.7.1)
    # -----------------------------------------------------------------------
    lines.extend([
        "## 3. Merged Pocket Catalog",
        "",
        "| Pocket | State | Sources | N Res | Volume (A3) | Proposal Score |",
        "|--------|-------|---------|-------|-------------|---------------|",
    ])
    for p in merged_pockets:
        lines.append(
            f"| {p['pocket_id']} | {p['receptor_id']} | {p.get('sources','')} | "
            f"{p.get('n_residues','')} | {p.get('mean_volume_A3','')} | "
            f"{p.get('best_proposal_score','')} |"
        )
    lines.append("")

    # Residue detail
    lines.extend([
        "### Pocket Lining Residues",
        "",
        "| Pocket | Residues |",
        "|--------|----------|",
    ])
    for p in merged_pockets:
        lines.append(f"| {p['pocket_id']} | {p.get('residue_ids','')} |")
    lines.append("")

    # -----------------------------------------------------------------------
    # 4. Patch Relationship Classification (Task 2.7.2)
    # -----------------------------------------------------------------------
    lines.extend([
        "## 4. PPI Patch Relationship Classification",
        "",
    ])

    by_rel: Dict[str, List[str]] = defaultdict(list)
    for r in relationships:
        by_rel[r["relationship_class"]].append(r["pocket_id"])

    lines.extend([
        "| Class | Count | Pockets |",
        "|-------|-------|---------|",
    ])
    for cls in ["orthosteric_candidate", "rim_candidate",
                "allosteric_candidate", "low_relevance_candidate"]:
        pids = by_rel.get(cls, [])
        lines.append(f"| {cls} | {len(pids)} | {', '.join(pids) or '---'} |")
    lines.append("")

    # Detail table
    lines.extend([
        "### Classification Detail",
        "",
        "| Pocket | Class | Hotspot Overlap | Fraction | Centroid Dist (A) | Basis |",
        "|--------|-------|----------------|----------|-------------------|-------|",
    ])
    for r in relationships:
        lines.append(
            f"| {r['pocket_id']} | {r['relationship_class']} | "
            f"{r.get('hotspot_overlap_count','')} | "
            f"{r.get('hotspot_overlap_fraction','')} | "
            f"{r.get('centroid_distance_A','')} | "
            f"{r.get('classification_basis','')} |"
        )
    lines.append("")

    # Uncertainty notes
    low_confidence_rels = [r for r in relationships
                          if r["relationship_class"] == "low_relevance_candidate"]
    if low_confidence_rels:
        lines.extend([
            "### Low-Relevance Pockets",
            "",
            f"{len(low_confidence_rels)} pocket(s) classified as `low_relevance_candidate`.",
            "These are structurally distant from the PPI patch and are not recommended",
            "for PPI-disruption-oriented docking in Phase 3.",
            "",
        ])

    # -----------------------------------------------------------------------
    # 5. Druggability Confidence (Task 2.7.2 extension)
    # -----------------------------------------------------------------------
    lines.extend([
        "## 5. Druggability Confidence",
        "",
    ])

    by_conf: Dict[str, List[str]] = defaultdict(list)
    by_tier: Dict[str, List[str]] = defaultdict(list)
    for d in druggability:
        by_conf[d.get("druggability_confidence", "")].append(d["pocket_id"])
        by_tier[d.get("overall_druggability_tier", "")].append(d["pocket_id"])

    lines.extend([
        "### Druggability Confidence",
        "",
        "| Confidence | Count | Pockets |",
        "|-----------|-------|---------|",
    ])
    for conf in ["high", "medium", "low"]:
        pids = by_conf.get(conf, [])
        lines.append(f"| {conf} | {len(pids)} | {', '.join(pids) or '---'} |")
    lines.append("")

    lines.extend([
        "### Overall Druggability Tier",
        "",
        "| Tier | Count | Pockets |",
        "|------|-------|---------|",
    ])
    for tier in ["tier_1", "tier_2", "tier_3"]:
        pids = by_tier.get(tier, [])
        lines.append(f"| {tier} | {len(pids)} | {', '.join(pids) or '---'} |")
    lines.append("")

    # Multi-source support
    consensus = [s for s in support_flags
                 if s.get("consensus_pocket", "").lower() == "true"]
    single = [s for s in support_flags
              if s.get("consensus_pocket", "").lower() == "false"]
    lines.extend([
        "### Multi-Source Support",
        "",
        f"- Consensus pockets (>=2 tools): {len(consensus)}",
        f"- Single-tool pockets: {len(single)}",
    ])
    if not consensus:
        lines.append("- **Note:** All pockets are currently single-tool (fpocket only).")
        lines.append("  Multi-source consensus will improve after P2Rank results are available.")
    lines.append("")

    # FTMap status
    ftmap_available = any(
        d.get("hotspot_support_available", "").lower() == "true" for d in druggability
    )
    lines.extend([
        "### FTMap Hotspot Support",
        "",
        f"- FTMap data available: {'Yes' if ftmap_available else 'No'}",
    ])
    if not ftmap_available:
        lines.append("- FTMap is the preferred hotspot-support method for PPI-site")
        lines.append("  ligandability evidence. Schema is ready for future integration.")
    lines.append("")

    # -----------------------------------------------------------------------
    # 6. Cross-State Alignment
    # -----------------------------------------------------------------------
    lines.extend([
        "## 6. Cross-State Pocket Alignment",
        "",
    ])

    by_sc: Dict[str, List[str]] = defaultdict(list)
    for sc in state_classes:
        by_sc[sc["state_class"]].append(sc["pocket_id"])

    lines.extend([
        "| State Class | Count | Pockets |",
        "|-------------|-------|---------|",
    ])
    for cls in ["state_robust_pocket", "state_shifted_pocket",
                "state_specific_pocket", "uncertain_alignment"]:
        pids = by_sc.get(cls, [])
        lines.append(f"| {cls} | {len(pids)} | {', '.join(pids) or '---'} |")
    lines.append("")

    if states_missing:
        lines.extend([
            f"**{len(states_missing)} state(s) missing pocket data:** "
            f"{', '.join(states_missing)}",
            "",
            "Cross-state classification is provisional. All pockets are currently",
            "`state_specific_pocket` because only one state has data.",
            "Re-run TG 2.5 and TG 2.6 after all states have pocket proposals.",
            "",
        ])

    # -----------------------------------------------------------------------
    # 7. Phase 3 Docking Recommendation (Task 2.7.3)
    # -----------------------------------------------------------------------
    lines.extend([
        "## 7. Phase 3 Docking Recommendation",
        "",
    ])

    lines.extend([
        "### Docking Priority Summary",
        "",
        "| Priority | Count | Pockets | Action |",
        "|----------|-------|---------|--------|",
    ])
    actions = {
        "primary": "Full docking budget",
        "secondary": "Reduced budget / diversity round",
        "exploratory": "Only if resources allow",
        "skip": "Do not dock",
    }
    for pri in ["primary", "secondary", "exploratory", "skip"]:
        pids = [r["pocket_id"] for r in phase3_ref
                if r.get("docking_priority") == pri]
        lines.append(
            f"| {pri} | {len(pids)} | {', '.join(pids) or '---'} | {actions[pri]} |"
        )
    lines.append("")

    # Comprehensive pocket table
    lines.extend([
        "### Comprehensive Pocket Summary",
        "",
        "| Pocket | Relationship | Tier | Confidence | State Class | Priority |",
        "|--------|-------------|------|-----------|-------------|----------|",
    ])
    drug_map = {d["pocket_id"]: d for d in druggability}
    sc_map = {s["pocket_id"]: s for s in state_classes}
    ref_map = {r["pocket_id"]: r for r in phase3_ref}
    for p in merged_pockets:
        pid = p["pocket_id"]
        d = drug_map.get(pid, {})
        s = sc_map.get(pid, {})
        r = ref_map.get(pid, {})
        lines.append(
            f"| {pid} | {d.get('relationship_class','')} | "
            f"{d.get('overall_druggability_tier','')} | "
            f"{d.get('druggability_confidence','')} | "
            f"{s.get('state_class','')} | "
            f"{r.get('docking_priority','')} |"
        )
    lines.append("")

    # Primary pocket detail
    primary_pockets = [r for r in phase3_ref if r.get("docking_priority") == "primary"]
    if primary_pockets:
        lines.extend([
            "### Primary Docking Targets — Detail",
            "",
        ])
        for r in primary_pockets:
            lines.extend([
                f"**{r['pocket_id']}** ({r['receptor_id']})",
                f"- Relationship: {r.get('relationship_class','')}",
                f"- Druggability: {r.get('druggability_confidence','')} "
                f"(tier {r.get('overall_druggability_tier','')})",
                f"- Hotspot overlap: {r.get('hotspot_overlap_count','')} residues "
                f"({r.get('hotspot_overlap_fraction','')} of patch)",
                f"- Centroid: ({r.get('centroid_x','')}, {r.get('centroid_y','')}, "
                f"{r.get('centroid_z','')})",
                f"- Box: {r.get('box_size_x','')} x {r.get('box_size_y','')} x "
                f"{r.get('box_size_z','')} A",
                f"- Lining residues: {r.get('residue_ids','')}",
                "",
            ])

    # -----------------------------------------------------------------------
    # 8. Limitations and Caveats
    # -----------------------------------------------------------------------
    lines.extend([
        "## 8. Limitations and Caveats",
        "",
    ])

    caveats = []
    if states_missing:
        caveats.append(
            f"Only {len(states_with_pockets)}/{len(RECEPTOR_STATES)} receptor states "
            "have pocket proposals. Cross-state robustness cannot be assessed yet."
        )
    if not consensus:
        caveats.append(
            "All pockets are single-tool (fpocket only). Multi-source consensus "
            "awaits P2Rank results."
        )
    if not ftmap_available:
        caveats.append(
            "FTMap hotspot support is not yet available. Druggability confidence "
            "relies solely on fpocket scores."
        )
    caveats.append(
        "Druggability tier thresholds are calibrated for fpocket's 0-1 score range. "
        "When P2Rank data is added, normalized scores enable cross-tool comparison."
    )

    for i, c in enumerate(caveats, 1):
        lines.append(f"{i}. {c}")
    lines.append("")

    # -----------------------------------------------------------------------
    # 9. Phase 2 Output File Inventory
    # -----------------------------------------------------------------------
    lines.extend([
        "## 9. Phase 2 Output File Inventory",
        "",
        "| File | TG | Description |",
        "|------|-----|------------|",
        "| phase2_patch_reference_normalized.csv | 2.0 | Phase 1 patch reference (normalized) |",
        "| phase2_patch_reference_validation.md | 2.0 | Patch reference validation note |",
        "| candidate_pockets_raw.csv | 2.1 | Raw pocket proposals (all tools) |",
        "| candidate_pocket_source_summary.csv | 2.1 | Per-tool proposal counts |",
        "| candidate_pockets.csv | 2.2 | Merged candidate pockets |",
        "| candidate_pocket_merge_table.csv | 2.2 | Pairwise merge decisions |",
        "| candidate_pocket_provenance.csv | 2.2 | Raw-to-merged mapping |",
        "| pocket_patch_relationship.csv | 2.3 | Patch relationship classes |",
        "| pocket_patch_relationship_metrics.csv | 2.3 | Classification raw metrics |",
        "| druggability_proposal_summary.csv | 2.4 | Druggability confidence + tiers |",
        "| candidate_pocket_support_flags.csv | 2.4 | Per-tool support flags |",
        "| candidate_pocket_cross_state_comparison.csv | 2.5 | Pairwise cross-state comparisons |",
        "| candidate_pocket_state_classes.csv | 2.5 | State robustness classes |",
        "| phase3_candidate_pocket_reference.csv | 2.6 | Phase 3 docking reference |",
        "| phase2_to_phase3_handoff_note.md | 2.6 | Phase 3 handoff note |",
        "| phase2_candidate_pocket_report.md | 2.7 | This report |",
        "",
    ])

    # -----------------------------------------------------------------------
    # Footer
    # -----------------------------------------------------------------------
    lines.extend([
        "---",
        "",
        "Generated by `egfr_pipeline.phase2.review_report`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 TG 2.7: Phase 2 review report"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE2_OUTPUT_DIR,
        help="Phase 2 output directory",
    )
    args = parser.parse_args()

    print("Phase 2 — Task Group 2.7: Review Report")
    print(f"  Output: {args.output_dir}")
    print()

    report_path = run_review_report(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 2.7 COMPLETE")
    print(f"{'='*60}")
    print(f"  Report: {report_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
