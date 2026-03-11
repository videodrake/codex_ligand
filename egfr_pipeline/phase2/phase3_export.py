#!/usr/bin/env python3
"""Phase 2 Task Group 2.6: Phase 3 Docking Preparation Export.

Exports a clean candidate pocket catalog and metadata package that Phase 3
can use directly for diversity-aware pocket-guided docking.

Tasks:
  2.6.1 - Candidate pocket docking export (unified reference CSV)
  2.6.2 - Phase 3 compatibility checks (schema + ID stability validation)
  2.6.3 - Downstream readiness note (handoff markdown)

Inputs:
  - output/phase2_pockets/candidate_pockets.csv
  - output/phase2_pockets/pocket_patch_relationship.csv
  - output/phase2_pockets/druggability_proposal_summary.csv
  - output/phase2_pockets/candidate_pocket_state_classes.csv

Outputs:
  - output/phase2_pockets/phase3_candidate_pocket_reference.csv
  - output/phase2_pockets/phase2_to_phase3_handoff_note.md

Usage:
    python -m egfr_pipeline.phase2.phase3_export [--output_dir ...]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE2_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase2_pockets"

EXPECTED_RECEPTOR_STATES = {"3GT8_raw", "3GT8_cl38_48", "3GT8_cl85_100"}

# Phase 3 reference schema — the single file Phase 3 needs to consume
EXPORT_COLUMNS = [
    # Identity
    "receptor_id",
    "pocket_id",
    # Docking box geometry
    "centroid_x",
    "centroid_y",
    "centroid_z",
    "box_size_x",
    "box_size_y",
    "box_size_z",
    # PPI patch relationship (from TG 2.3)
    "relationship_class",
    "hotspot_overlap_count",
    "hotspot_overlap_fraction",
    "centroid_distance_to_patch_A",
    # Druggability (from TG 2.4)
    "druggability_confidence",
    "overall_druggability_tier",
    "best_proposal_score",
    "best_druggability_score",
    "normalized_druggability_score",
    # Multi-source support (from TG 2.4)
    "n_sources",
    "sources",
    "multi_source_support",
    # Cross-state (from TG 2.5)
    "state_class",
    "n_states_matched",
    # Pocket geometry
    "n_residues",
    "residue_ids",
    "mean_volume_A3",
    # Phase 3 guidance
    "docking_priority",           # primary / secondary / exploratory / skip
    "priority_basis",
]

# Docking priority rules
# primary:     tier_1 or tier_2 + orthosteric/rim
# secondary:   tier_2 allosteric, or tier_3 + orthosteric/rim
# exploratory: tier_3 allosteric or uncertain
# skip:        tier_3 low_relevance (waste of docking budget)

PRIORITY_CLASSES = ["primary", "secondary", "exploratory", "skip"]


# ---------------------------------------------------------------------------
# Priority assignment (Task 2.6.1)
# ---------------------------------------------------------------------------

def assign_docking_priority(
    tier: str,
    relationship: str,
    state_class: str,
) -> Tuple[str, str]:
    """Assign a docking priority for Phase 3 budget allocation.

    Returns (priority, basis_explanation).
    """
    ppi_relevant = relationship in (
        "orthosteric_candidate", "rim_candidate"
    )
    allosteric = relationship == "allosteric_candidate"
    low_rel = relationship == "low_relevance_candidate"

    if tier == "tier_1":
        return "primary", f"tier_1 + {relationship}"

    if tier == "tier_2":
        if ppi_relevant:
            return "primary", f"tier_2 + {relationship}"
        if allosteric:
            return "secondary", f"tier_2 + allosteric"
        return "exploratory", f"tier_2 + {relationship}"

    # tier_3
    if ppi_relevant:
        return "secondary", f"tier_3 but {relationship}"
    if allosteric:
        return "exploratory", f"tier_3 + allosteric"
    if low_rel:
        return "skip", f"tier_3 + low_relevance (not recommended for docking)"

    return "exploratory", f"tier_3 + {relationship}"


# ---------------------------------------------------------------------------
# Compatibility checks (Task 2.6.2)
# ---------------------------------------------------------------------------

def validate_export(rows: List[dict]) -> List[str]:
    """Run Phase 3 compatibility checks on the export.

    Returns list of warning/error messages (empty = all good).
    """
    issues = []

    if not rows:
        issues.append("ERROR: No pockets to export")
        return issues

    # Check required fields are non-empty
    required = ["receptor_id", "pocket_id", "centroid_x", "centroid_y", "centroid_z"]
    for i, r in enumerate(rows):
        for field in required:
            if not r.get(field, "").strip():
                issues.append(f"WARNING: Row {i+1} ({r.get('pocket_id','?')}): "
                              f"missing required field '{field}'")

    # Check pocket ID uniqueness
    ids = [r["pocket_id"] for r in rows]
    if len(ids) != len(set(ids)):
        dupes = [pid for pid in ids if ids.count(pid) > 1]
        issues.append(f"ERROR: Duplicate pocket IDs: {set(dupes)}")

    # Check all pockets have relationship class
    no_rel = [r["pocket_id"] for r in rows if not r.get("relationship_class", "")]
    if no_rel:
        issues.append(f"WARNING: {len(no_rel)} pocket(s) missing relationship_class")

    # Check docking priority coverage
    no_priority = [r["pocket_id"] for r in rows if not r.get("docking_priority", "")]
    if no_priority:
        issues.append(f"WARNING: {len(no_priority)} pocket(s) missing docking_priority")

    return issues


# ---------------------------------------------------------------------------
# Handoff note (Task 2.6.3)
# ---------------------------------------------------------------------------

def build_handoff_note(
    rows: List[dict],
    validation_issues: List[str],
) -> List[str]:
    """Build the Phase 2 → Phase 3 handoff markdown note."""
    lines = [
        "# Phase 2 → Phase 3 Handoff Note",
        "",
        "## Overview",
        "",
        "This document accompanies `phase3_candidate_pocket_reference.csv`,",
        "the primary input for Phase 3 diversity-aware pocket-guided docking.",
        "",
        f"- Total candidate pockets exported: {len(rows)}",
    ]

    # Per-state counts
    by_state: Dict[str, int] = defaultdict(int)
    for r in rows:
        by_state[r["receptor_id"]] += 1
    for state, count in sorted(by_state.items()):
        lines.append(f"  - {state}: {count} pockets")
    lines.append("")

    # Priority distribution
    by_priority: Dict[str, List[str]] = defaultdict(list)
    for r in rows:
        by_priority[r["docking_priority"]].append(r["pocket_id"])

    lines.extend([
        "## Docking Priority Distribution",
        "",
        "| Priority | Count | Pockets | Phase 3 Action |",
        "|----------|-------|---------|---------------|",
    ])
    actions = {
        "primary": "Full docking budget allocation",
        "secondary": "Reduced budget, include in diversity round",
        "exploratory": "Minimal budget, only if resources allow",
        "skip": "Not recommended for docking",
    }
    for pri in PRIORITY_CLASSES:
        pids = by_priority.get(pri, [])
        lines.append(
            f"| {pri} | {len(pids)} | {', '.join(pids) or '—'} | {actions.get(pri, '')} |"
        )
    lines.append("")

    # Priority definition
    lines.extend([
        "## Docking Priority Definitions",
        "",
        "| Priority | Criteria |",
        "|----------|----------|",
        "| primary | tier_1, or tier_2 + orthosteric/rim |",
        "| secondary | tier_2 + allosteric, or tier_3 + orthosteric/rim |",
        "| exploratory | tier_2/3 + uncertain class, or tier_3 + allosteric |",
        "| skip | tier_3 + low_relevance (docking budget not justified) |",
        "",
    ])

    # Relationship class distribution
    by_rel: Dict[str, int] = defaultdict(int)
    for r in rows:
        by_rel[r.get("relationship_class", "unknown")] += 1

    lines.extend([
        "## PPI Patch Relationship Distribution",
        "",
        "| Class | Count |",
        "|-------|-------|",
    ])
    for cls in ["orthosteric_candidate", "rim_candidate",
                "allosteric_candidate", "low_relevance_candidate"]:
        lines.append(f"| {cls} | {by_rel.get(cls, 0)} |")
    lines.append("")

    # Druggability tier distribution
    by_tier: Dict[str, int] = defaultdict(int)
    for r in rows:
        by_tier[r.get("overall_druggability_tier", "unknown")] += 1

    lines.extend([
        "## Druggability Tier Distribution",
        "",
        "| Tier | Count |",
        "|------|-------|",
    ])
    for tier in ["tier_1", "tier_2", "tier_3"]:
        lines.append(f"| {tier} | {by_tier.get(tier, 0)} |")
    lines.append("")

    # Per-pocket reference table
    lines.extend([
        "## Per-Pocket Reference",
        "",
        "| Pocket | State | Relationship | Tier | Priority | Box (Å) |",
        "|--------|-------|-------------|------|----------|---------|",
    ])
    for r in rows:
        box = f"{r.get('box_size_x','?')}×{r.get('box_size_y','?')}×{r.get('box_size_z','?')}"
        lines.append(
            f"| {r['pocket_id']} | {r['receptor_id']} | "
            f"{r.get('relationship_class','')} | "
            f"{r.get('overall_druggability_tier','')} | "
            f"{r.get('docking_priority','')} | {box} |"
        )
    lines.append("")

    # Validation
    lines.extend([
        "## Export Validation",
        "",
    ])
    if validation_issues:
        lines.append(f"**Issues found: {len(validation_issues)}**")
        lines.append("")
        for issue in validation_issues:
            lines.append(f"- {issue}")
    else:
        lines.append("All compatibility checks passed.")
    lines.append("")

    # Phase 3 consumption instructions
    lines.extend([
        "## Phase 3 Consumption Instructions",
        "",
        "1. Load `phase3_candidate_pocket_reference.csv`",
        "2. Filter by `docking_priority` to allocate docking budget:",
        "   - `primary`: full docking runs with diversity sampling",
        "   - `secondary`: reduced decoy count or single-seed",
        "   - `exploratory`: only if budget allows",
        "   - `skip`: do not dock",
        "3. Use `centroid_x/y/z` and `box_size_x/y/z` for Vina/AutoDock box setup",
        "4. Use `relationship_class` for Phase 4 perturbation relevance context",
        "5. Use `state_class` for cross-state robustness weighting",
        "",
        "## Cross-State Note",
        "",
    ])

    states_present = set(r["receptor_id"] for r in rows)
    missing = EXPECTED_RECEPTOR_STATES - states_present
    if missing:
        lines.extend([
            f"**{len(missing)} receptor state(s) do not yet have pocket data:**",
            f"{', '.join(sorted(missing))}",
            "",
            "Cross-state alignment (`state_class`) is provisional. Re-export",
            "after all states have fpocket/P2Rank results.",
        ])
    else:
        lines.append("All three receptor states have pocket data.")
    lines.append("")

    lines.extend([
        "---",
        "",
        "Generated by `egfr_pipeline.phase2.phase3_export`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_phase3_export(output_dir: Path) -> Tuple[Path, Path]:
    """Run full TG 2.6 pipeline.

    Returns (export_csv, handoff_note_path).
    """
    # Load all upstream data
    pockets = _load_csv(output_dir / "candidate_pockets.csv")
    relationships = _load_csv(output_dir / "pocket_patch_relationship.csv")
    druggability = _load_csv(output_dir / "druggability_proposal_summary.csv")
    state_classes = _load_csv(output_dir / "candidate_pocket_state_classes.csv")

    if not pockets:
        raise FileNotFoundError("No candidate pockets found. Run TG 2.2 first.")

    print(f"  Loaded: {len(pockets)} pockets, {len(relationships)} relationships, "
          f"{len(druggability)} druggability, {len(state_classes)} state classes")

    # Build lookup maps
    rel_map = {r["pocket_id"]: r for r in relationships}
    drug_map = {r["pocket_id"]: r for r in druggability}
    state_map = {r["pocket_id"]: r for r in state_classes}

    # Build export rows (Task 2.6.1)
    export_rows = []
    for p in pockets:
        pid = p["pocket_id"]
        rel = rel_map.get(pid, {})
        drug = drug_map.get(pid, {})
        sc = state_map.get(pid, {})

        tier = drug.get("overall_druggability_tier", "tier_3")
        rel_class = rel.get("relationship_class", "")
        state_cls = sc.get("state_class", "")

        priority, basis = assign_docking_priority(tier, rel_class, state_cls)

        export_rows.append({
            "receptor_id": p["receptor_id"],
            "pocket_id": pid,
            "centroid_x": p.get("centroid_x", ""),
            "centroid_y": p.get("centroid_y", ""),
            "centroid_z": p.get("centroid_z", ""),
            "box_size_x": p.get("box_size_x", ""),
            "box_size_y": p.get("box_size_y", ""),
            "box_size_z": p.get("box_size_z", ""),
            "relationship_class": rel_class,
            "hotspot_overlap_count": rel.get("hotspot_overlap_count", ""),
            "hotspot_overlap_fraction": rel.get("hotspot_overlap_fraction", ""),
            "centroid_distance_to_patch_A": rel.get("centroid_distance_A", ""),
            "druggability_confidence": drug.get("druggability_confidence", ""),
            "overall_druggability_tier": tier,
            "best_proposal_score": drug.get("best_proposal_score", ""),
            "best_druggability_score": drug.get("best_druggability_score", ""),
            "normalized_druggability_score": drug.get("normalized_druggability_score", ""),
            "n_sources": drug.get("n_sources", p.get("n_sources", "")),
            "sources": drug.get("sources", p.get("sources", "")),
            "multi_source_support": drug.get("multi_source_support", ""),
            "state_class": state_cls,
            "n_states_matched": sc.get("n_states_matched", ""),
            "n_residues": p.get("n_residues", ""),
            "residue_ids": p.get("residue_ids", ""),
            "mean_volume_A3": p.get("mean_volume_A3", ""),
            "docking_priority": priority,
            "priority_basis": basis,
        })

    # Task 2.6.2: Validation
    issues = validate_export(export_rows)
    if issues:
        for issue in issues:
            print(f"  {issue}")
    else:
        print("  Validation: all checks passed")

    # Priority distribution
    by_priority: Dict[str, int] = defaultdict(int)
    for r in export_rows:
        by_priority[r["docking_priority"]] += 1
    print(f"  Docking priorities: {dict(by_priority)}")

    # Write export CSV
    export_path = output_dir / "phase3_candidate_pocket_reference.csv"
    _write_csv(export_path, EXPORT_COLUMNS, export_rows)

    # Task 2.6.3: Handoff note
    note_lines = build_handoff_note(export_rows, issues)
    note_path = output_dir / "phase2_to_phase3_handoff_note.md"
    with open(note_path, "w", encoding="utf-8") as f:
        f.write("\n".join(note_lines))

    return export_path, note_path


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
        description="Phase 2 TG 2.6: Phase 3 docking preparation export"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE2_OUTPUT_DIR,
        help="Phase 2 output directory",
    )
    args = parser.parse_args()

    print("Phase 2 — Task Group 2.6: Phase 3 Docking Preparation Export")
    print(f"  Output: {args.output_dir}")
    print()

    export_path, note_path = run_phase3_export(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 2.6 COMPLETE")
    print(f"{'='*60}")
    print(f"  Export:   {export_path}")
    print(f"  Handoff:  {note_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
