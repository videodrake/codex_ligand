#!/usr/bin/env python3
"""Phase 3 Task Group 3.7: Phase 3 Review Report.

Generates a readable review package showing how docking budget was
distributed and whether candidate pockets were explored in a balanced,
biologically useful way.

Tasks:
  3.7.1 - Search-budget summary report
  3.7.2 - Ligand support summary
  3.7.3 - Diversity outcome summary

Inputs:
  - output/phase3_docking/phase3_candidate_reference_normalized.csv
  - output/phase3_docking/phase3_docking_job_table.csv
  - output/phase3_docking/pocket_search_status.csv
  - output/phase3_docking/phase3_pocket_occupancy_summary.csv
  - output/phase3_docking/phase3_diversity_metrics.csv
  - output/phase3_docking/phase4_docking_evidence_reference.csv
  - output/phase3_docking/phase3_run_metadata.json

Outputs:
  - output/phase3_docking/phase3_diverse_docking_report.md

Usage:
    python -m egfr_pipeline.phase3.review_report [--output_dir ...]
"""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

PHASE3_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase3_docking"


# ---------------------------------------------------------------------------
# Report builder
# ---------------------------------------------------------------------------

def build_report(
    normalized: List[dict],
    jobs: List[dict],
    statuses: List[dict],
    occupancy: List[dict],
    diversity: List[dict],
    evidence: List[dict],
    metadata: dict,
) -> List[str]:
    """Build the Phase 3 review report."""
    lines = []

    # ── Section 1: Executive Summary ──
    dockable = [p for p in normalized if p.get("dockable") == "True"]
    skip = [p for p in normalized if p.get("dockable") != "True"]
    receptors = sorted(set(p["receptor_id"] for p in normalized))
    ligands = sorted(set(j["ligand_id"] for j in jobs))

    lines.extend([
        "# Phase 3 Diversity-Aware Docking Report",
        "",
        "## 1. Executive Summary",
        "",
        f"- **Receptor states**: {len(receptors)} ({', '.join(receptors)})",
        f"- **Candidate pockets**: {len(normalized)} total, "
        f"{len(dockable)} dockable, {len(skip)} skipped",
        f"- **Ligands**: {len(ligands)} ({', '.join(ligands)})",
        f"- **Total docking jobs**: {len(jobs)}",
        f"- **Docking mode**: focused (pocket-guided)",
        "",
    ])

    # Execution status
    has_results = any(e.get("best_affinity_kcal") for e in evidence)
    if has_results:
        lines.append("**Execution status**: Post-docking (affinity data available)")
    else:
        lines.append(
            "**Execution status**: Pre-docking (infrastructure validated, "
            "awaiting server-side Vina execution)"
        )
    lines.append("")

    # ── Section 2: Pocket Catalog ──
    lines.extend([
        "## 2. Pocket Catalog",
        "",
        "| Pocket | Receptor | Priority | Tier | Relationship | Box (A) | Seeds |",
        "|--------|----------|----------|------|-------------|---------|-------|",
    ])
    for p in normalized:
        box = f"{p.get('box_size_x','')}x{p.get('box_size_y','')}x{p.get('box_size_z','')}"
        lines.append(
            f"| {p['pocket_id']} | {p['receptor_id']} | "
            f"{p['docking_priority']} | {p.get('overall_druggability_tier','')} | "
            f"{p.get('relationship_class','')} | {box} | "
            f"{p.get('recommended_seeds','')} |"
        )
    lines.append("")

    # ── Section 3: Search Budget Summary (Task 3.7.1) ──
    lines.extend([
        "## 3. Search Budget Summary",
        "",
        "### Budget Policy",
        "",
    ])
    if metadata:
        policy = metadata.get("budget_policy", {})
        lines.extend([
            f"- Max rounds: {policy.get('max_rounds', 'N/A')}",
            f"- Saturation threshold: {policy.get('saturation_pose_threshold', 'N/A')} acceptable poses",
            f"- Affinity window: {policy.get('saturation_affinity_window_kcal', 'N/A')} kcal/mol",
            f"- Max poses per pocket: {policy.get('max_poses_per_pocket', 'N/A')}",
            f"- Reallocation: {'enabled' if policy.get('enable_reallocation') else 'disabled'}",
            "",
        ])

    lines.extend([
        "### Pocket Budget Status",
        "",
        "| Pocket | Priority | Status | Seeds Allocated | Seeds Done | Poses | Saturated? |",
        "|--------|----------|--------|-----------------|------------|-------|------------|",
    ])
    for ps in statuses:
        sat_round = ps.get("round_saturated", "")
        sat_str = f"Round {sat_round}" if sat_round else "No"
        lines.append(
            f"| {ps['pocket_id']} | {ps['docking_priority']} | "
            f"{ps['status']} | {ps.get('seeds_allocated',0)} | "
            f"{ps.get('seeds_completed',0)} | "
            f"{ps.get('total_poses_generated',0)} | {sat_str} |"
        )
    lines.append("")

    # Job distribution by round
    by_round = defaultdict(int)
    for j in jobs:
        try:
            seed_idx = int(j["seed_index"])
        except (ValueError, TypeError):
            continue
        by_round[seed_idx] += 1
    lines.extend([
        "### Job Distribution by Round",
        "",
        "| Round | Jobs |",
        "|-------|------|",
    ])
    for rid in sorted(by_round.keys()):
        lines.append(f"| {rid} | {by_round[rid]} |")
    lines.append("")

    # ── Section 4: Ligand Support Summary (Task 3.7.2) ──
    lines.extend([
        "## 4. Ligand Support Summary",
        "",
        "| Pocket | Ligand | Poses | Best Affinity | Support |",
        "|--------|--------|-------|---------------|---------|",
    ])
    for e in sorted(evidence, key=lambda x: (x["candidate_pocket_id"], x["ligand_id"])):
        best = e.get("best_affinity_kcal", "—")
        if not best:
            best = "pending"
        lines.append(
            f"| {e['candidate_pocket_id']} | {e['ligand_id']} | "
            f"{e['pose_support_count']} | {best} | "
            f"{e['ligand_support_strength']} |"
        )
    lines.append("")

    # Multimodal vs dominant pocket
    by_ligand = defaultdict(set)
    for e in evidence:
        if int(e.get("pose_support_count", 0)) > 0:
            by_ligand[e["ligand_id"]].add(e["candidate_pocket_id"])

    lines.extend(["### Pocket Distribution per Ligand", ""])
    for lig_id in sorted(by_ligand.keys()):
        pockets = by_ligand[lig_id]
        if len(pockets) > 1:
            lines.append(f"- **{lig_id}**: multimodal ({len(pockets)} pockets)")
        elif len(pockets) == 1:
            lines.append(f"- **{lig_id}**: single pocket ({next(iter(pockets))})")
        else:
            lines.append(f"- **{lig_id}**: no poses")
    lines.append("")

    # ── Section 5: Diversity Outcome (Task 3.7.3) ──
    lines.extend([
        "## 5. Diversity Outcome",
        "",
        "| Receptor | Ligand | Pockets Sampled | Coverage | Top-1 Conc. | Entropy | Verdict |",
        "|----------|--------|-----------------|----------|-------------|---------|---------|",
    ])
    for dm in diversity:
        lines.append(
            f"| {dm['receptor_id']} | {dm['ligand_id']} | "
            f"{dm['n_pockets_sampled']}/{dm['n_pockets_total']} | "
            f"{dm['coverage_fraction']} | {dm['top1_concentration']} | "
            f"{dm['normalized_entropy']} | {dm['diversity_verdict']} |"
        )
    lines.append("")

    # Interpretation
    verdicts = [dm["diversity_verdict"] for dm in diversity]
    well = verdicts.count("well_distributed")
    over = verdicts.count("over_concentrated")
    lines.extend([
        "### Interpretation",
        "",
    ])
    if over > 0:
        lines.append(
            f"**Warning**: {over}/{len(verdicts)} ligand(s) show over-concentration "
            "into dominant pockets. Consider increasing budget for under-sampled pockets."
        )
    elif well == len(verdicts):
        lines.append(
            "All ligands show well-distributed search across available pockets. "
            "The diversity-aware budget policy is working as intended."
        )
    else:
        lines.append(
            f"{well}/{len(verdicts)} ligands well-distributed. "
            "Remaining show moderate distribution."
        )
    lines.append("")

    # ── Section 6: Phase 4 Handoff Summary ──
    lines.extend([
        "## 6. Phase 4 Handoff Summary",
        "",
        f"- Evidence file: `phase4_docking_evidence_reference.csv` ({len(evidence)} rows)",
        f"- Pockets with evidence: {len(set(e['candidate_pocket_id'] for e in evidence))}",
        f"- Ligands with evidence: {len(set(e['ligand_id'] for e in evidence))}",
        "",
    ])

    support_dist = defaultdict(int)
    for e in evidence:
        support_dist[e.get("ligand_support_strength", "")] += 1
    lines.extend([
        "### Support Strength Distribution",
        "",
        "| Strength | Count |",
        "|----------|-------|",
    ])
    for strength in ["strong", "moderate", "weak", "none",
                     "pending_strong", "pending_moderate", "pending_weak"]:
        if strength in support_dist:
            lines.append(f"| {strength} | {support_dist[strength]} |")
    lines.append("")

    # ── Section 7: Caveats ──
    lines.extend([
        "## 7. Caveats",
        "",
        "1. **Single receptor state**: Only 3GT8_raw has pocket data. "
        "Re-run Phase 2-3 when cl38_48 and cl85_100 fpocket/P2Rank results "
        "are available.",
        "2. **Single pocket tool**: Only fpocket results available. "
        "P2Rank integration will improve multi-source consensus.",
        "3. **Pre-execution**: Affinity values are pending server-side "
        "Vina execution. All support levels are provisional (`pending_*`).",
        "4. **Diversity metrics**: With only 2 dockable pockets, "
        "maximum entropy = 1.0 bit. More pockets will make diversity "
        "metrics more informative.",
        "",
    ])

    # ── Section 8: File Inventory ──
    lines.extend([
        "## 8. File Inventory",
        "",
        "| File | Description |",
        "|------|-------------|",
        "| phase3_candidate_reference_normalized.csv | Normalized pocket reference |",
        "| phase3_candidate_reference_validation.md | Reference validation report |",
        "| phase3_docking_job_table.csv | Full job table (receptor x pocket x ligand x seed) |",
        "| phase3_job_box_table.csv | Per-pocket box definitions |",
        "| phase3_budget_policy.md | Budget policy documentation |",
        "| pocket_search_status.csv | Per-pocket search status |",
        "| phase3_budget_tracking.csv | Round-by-round budget tracking |",
        "| run_phase3_docking.sh | Server-executable docking script |",
        "| phase3_run_metadata.json | Run metadata and configuration |",
        "| phase3_round_log.csv | Round execution log |",
        "| vina_pose_table.csv | Extended pose table with provenance |",
        "| phase3_pose_provenance_note.md | Pose schema documentation |",
        "| phase3_pocket_occupancy_summary.csv | Per-pocket occupancy |",
        "| phase3_diversity_metrics.csv | Diversity validation metrics |",
        "| phase3_blind_vs_diverse_comparison.csv | Blind vs diverse comparison |",
        "| phase4_docking_evidence_reference.csv | Phase 4 handoff evidence |",
        "| phase3_to_phase4_handoff_note.md | Phase 4 handoff documentation |",
        "| phase3_diverse_docking_report.md | This report |",
        "",
        "---",
        "",
        "Generated by `egfr_pipeline.phase3.review_report`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_review_report(
    output_dir: Path = PHASE3_OUTPUT_DIR,
) -> Path:
    """Run full TG 3.7 pipeline. Returns report_path."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load all inputs
    normalized = _load_csv(output_dir / "phase3_candidate_reference_normalized.csv")
    jobs = _load_csv(output_dir / "phase3_docking_job_table.csv")
    statuses = _load_csv(output_dir / "pocket_search_status.csv")
    occupancy = _load_csv(output_dir / "phase3_pocket_occupancy_summary.csv")
    diversity = _load_csv(output_dir / "phase3_diversity_metrics.csv")
    evidence = _load_csv(output_dir / "phase4_docking_evidence_reference.csv")

    metadata = {}
    meta_path = output_dir / "phase3_run_metadata.json"
    if meta_path.exists():
        with open(meta_path, encoding="utf-8") as f:
            metadata = json.load(f)

    print(f"  Loaded: {len(normalized)} pockets, {len(jobs)} jobs, "
          f"{len(evidence)} evidence rows")

    # Build report
    report_lines = build_report(
        normalized, jobs, statuses, occupancy, diversity, evidence, metadata
    )

    report_path = output_dir / "phase3_diverse_docking_report.md"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    return report_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_csv(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 3 TG 3.7: Phase 3 review report"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE3_OUTPUT_DIR,
        help="Phase 3 output directory",
    )
    args = parser.parse_args()

    print("Phase 3 — Task Group 3.7: Review Report")
    print(f"  Output: {args.output_dir}")
    print()

    report_path = run_review_report(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 3.7 COMPLETE")
    print(f"{'='*60}")
    print(f"  Report: {report_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
