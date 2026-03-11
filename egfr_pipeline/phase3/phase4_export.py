#!/usr/bin/env python3
"""Phase 3 Task Group 3.6: Phase 4-Ready Export.

Exports a docking evidence package that preserves ligand support,
pocket provenance, and diversity-aware search history for Phase 4
perturbation scoring.

Tasks:
  3.6.1 - Final docking evidence export
  3.6.2 - Perturbation-relevance handoff preparation
  3.6.3 - Handoff quality checks

Inputs:
  - output/phase3_docking/vina_pose_table.csv
  - output/phase3_docking/phase3_pocket_occupancy_summary.csv
  - output/phase3_docking/phase3_diversity_metrics.csv
  - output/phase3_docking/pocket_search_status.csv
  - output/phase3_docking/phase3_candidate_reference_normalized.csv

Outputs:
  - output/phase3_docking/phase4_docking_evidence_reference.csv
  - output/phase3_docking/phase3_to_phase4_handoff_note.md

Usage:
    python -m egfr_pipeline.phase3.phase4_export [--output_dir ...]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE3_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase3_docking"

# Phase 4 docking evidence schema (Task 3.6.1 + 3.6.2)
EVIDENCE_COLUMNS = [
    # Identity
    "receptor_id",
    "candidate_pocket_id",
    "ligand_id",
    # Docking evidence
    "pose_support_count",
    "best_affinity_kcal",
    "mean_affinity_kcal",
    "n_seeds_completed",
    "docking_mode",
    # Contact residues (from pose table, if available)
    "contact_residues",
    "n_contact_residues",
    # Pocket context (from Phase 2/3 reference)
    "relationship_class",
    "overall_druggability_tier",
    "druggability_confidence",
    "docking_priority",
    # State context
    "state_class",
    "n_states_matched",
    # Search history
    "pocket_status",
    "diversity_verdict",
    "coverage_fraction",
    "normalized_entropy",
    # Ligand support strength (derived)
    "ligand_support_strength",
    # Budget summary
    "total_poses_in_pocket",
    "fraction_of_receptor_poses",
]

# Ligand support strength thresholds
STRONG_SUPPORT_MIN_POSES = 5
STRONG_SUPPORT_MAX_AFFINITY = -6.0  # kcal/mol
MODERATE_SUPPORT_MIN_POSES = 2
MODERATE_SUPPORT_MAX_AFFINITY = -5.0


# ---------------------------------------------------------------------------
# 3.6.1: Final docking evidence export
# ---------------------------------------------------------------------------

def build_evidence_table(
    pose_rows: List[dict],
    occupancy: List[dict],
    diversity: List[dict],
    pocket_statuses: List[dict],
    normalized_pockets: List[dict],
) -> List[dict]:
    """Build the Phase 4 docking evidence reference.

    One row per (receptor, pocket, ligand) combination.
    """
    # Build lookups
    pocket_map = {p["pocket_id"]: p for p in normalized_pockets if p.get("pocket_id")}
    status_map = {ps["pocket_id"]: ps for ps in pocket_statuses if ps.get("pocket_id")}

    # Diversity lookup by (receptor, ligand)
    div_map = {}
    for dm in diversity:
        div_map[(dm["receptor_id"], dm["ligand_id"])] = dm

    # Occupancy lookup by (receptor, ligand, pocket)
    occ_map = {}
    for occ in occupancy:
        occ_map[(occ["receptor_id"], occ["ligand_id"],
                 occ["candidate_pocket_id"])] = occ

    # Aggregate poses by (receptor, ligand, pocket)
    pose_groups = defaultdict(list)
    for row in pose_rows:
        key = (row.get("receptor_id", ""),
               row.get("ligand_id", ""),
               row.get("candidate_pocket_id", ""))
        pose_groups[key].append(row)

    # Build evidence rows for all occupancy entries with dockable pockets
    evidence = []
    for occ in occupancy:
        receptor_id = occ["receptor_id"]
        ligand_id = occ["ligand_id"]
        pocket_id = occ["candidate_pocket_id"]

        if occ.get("docking_priority") == "skip":
            continue

        pocket = pocket_map.get(pocket_id, {})
        ps = status_map.get(pocket_id, {})
        dm = div_map.get((receptor_id, ligand_id), {})

        # Aggregate contact residues from poses
        key = (receptor_id, ligand_id, pocket_id)
        poses = pose_groups.get(key, [])
        all_contacts = set()
        for p in poses:
            cr = p.get("contact_residues", "")
            if cr:
                for r in cr.split(";"):
                    r = r.strip()
                    if r:
                        all_contacts.add(r)

        n_poses = int(occ.get("n_poses", 0))
        best_aff = occ.get("best_affinity_kcal", "")
        mean_aff = occ.get("mean_affinity_kcal", "")

        # Ligand support strength (Task 3.6.2)
        support = _classify_support(n_poses, best_aff)

        evidence.append({
            "receptor_id": receptor_id,
            "candidate_pocket_id": pocket_id,
            "ligand_id": ligand_id,
            "pose_support_count": n_poses,
            "best_affinity_kcal": best_aff,
            "mean_affinity_kcal": mean_aff,
            "n_seeds_completed": occ.get("n_seeds", 0),
            "docking_mode": "focused",
            "contact_residues": ";".join(sorted(all_contacts)) if all_contacts else "",
            "n_contact_residues": len(all_contacts),
            "relationship_class": pocket.get("relationship_class", ""),
            "overall_druggability_tier": pocket.get("overall_druggability_tier", ""),
            "druggability_confidence": pocket.get("druggability_confidence", ""),
            "docking_priority": pocket.get("docking_priority", ""),
            "state_class": pocket.get("state_class", ""),
            "n_states_matched": pocket.get("n_states_matched", ""),
            "pocket_status": ps.get("status", ""),
            "diversity_verdict": dm.get("diversity_verdict", ""),
            "coverage_fraction": dm.get("coverage_fraction", ""),
            "normalized_entropy": dm.get("normalized_entropy", ""),
            "ligand_support_strength": support,
            "total_poses_in_pocket": n_poses,
            "fraction_of_receptor_poses": occ.get("fraction_of_receptor_poses", ""),
        })

    return evidence


def _classify_support(n_poses: int, best_affinity_str: str) -> str:
    """Classify ligand support strength."""
    if n_poses == 0:
        return "none"

    try:
        best = float(best_affinity_str) if best_affinity_str else None
    except (ValueError, TypeError):
        best = None

    # Pre-execution: no affinity data yet
    if best is None:
        if n_poses >= STRONG_SUPPORT_MIN_POSES:
            return "pending_strong"
        elif n_poses >= MODERATE_SUPPORT_MIN_POSES:
            return "pending_moderate"
        return "pending_weak"

    if n_poses >= STRONG_SUPPORT_MIN_POSES and best <= STRONG_SUPPORT_MAX_AFFINITY:
        return "strong"
    elif n_poses >= MODERATE_SUPPORT_MIN_POSES and best <= MODERATE_SUPPORT_MAX_AFFINITY:
        return "moderate"
    elif n_poses > 0:
        return "weak"
    return "none"


# ---------------------------------------------------------------------------
# 3.6.3: Handoff quality checks
# ---------------------------------------------------------------------------

def validate_handoff(evidence: List[dict]) -> List[str]:
    """Validate the Phase 4 export for completeness."""
    messages = []

    if not evidence:
        messages.append("[ERROR] Empty evidence table")
        return messages

    messages.append(f"[OK] {len(evidence)} evidence rows")

    # Required fields
    required = ["receptor_id", "candidate_pocket_id", "ligand_id",
                 "pose_support_count", "relationship_class",
                 "docking_priority", "ligand_support_strength"]
    missing = [f for f in required if f not in evidence[0]]
    if missing:
        messages.append(f"[ERROR] Missing required fields: {missing}")
    else:
        messages.append("[OK] All required fields present")

    # Perturbation-relevance fields (Task 3.6.2)
    relevance_fields = ["relationship_class", "overall_druggability_tier",
                        "state_class", "ligand_support_strength"]
    for field in relevance_fields:
        empty = sum(1 for r in evidence if not r.get(field))
        if empty:
            messages.append(f"[WARN] {field}: {empty}/{len(evidence)} rows empty")
        else:
            messages.append(f"[OK] {field}: all rows populated")

    # Check pocket provenance
    pockets = set(r["candidate_pocket_id"] for r in evidence)
    messages.append(f"[OK] {len(pockets)} distinct pocket(s) in evidence")

    # Check ligand coverage
    ligands = set(r["ligand_id"] for r in evidence)
    messages.append(f"[OK] {len(ligands)} distinct ligand(s) in evidence")

    # Support strength distribution
    support_dist = defaultdict(int)
    for r in evidence:
        support_dist[r.get("ligand_support_strength", "")] += 1
    messages.append(f"[OK] Support distribution: {dict(sorted(support_dist.items()))}")

    return messages


# ---------------------------------------------------------------------------
# Handoff note generation
# ---------------------------------------------------------------------------

def build_handoff_note(
    evidence: List[dict],
    validation: List[str],
) -> List[str]:
    """Generate the Phase 3 → Phase 4 handoff markdown note."""
    lines = [
        "# Phase 3 → Phase 4 Handoff Note",
        "",
        "## Purpose",
        "",
        "This file documents the docking evidence exported from Phase 3 for",
        "consumption by Phase 4 (Perturbation Scoring).",
        "",
        "## Evidence File",
        "",
        "**`phase4_docking_evidence_reference.csv`**",
        "",
        f"- Rows: {len(evidence)}",
        f"- Schema: {len(EVIDENCE_COLUMNS)} columns",
        "",
        "## Schema Description",
        "",
        "| Column | Description | Phase 4 Use |",
        "|--------|-------------|-------------|",
        "| receptor_id | Receptor state | Group by state |",
        "| candidate_pocket_id | Phase 2/3 pocket | Primary key |",
        "| ligand_id | Docked ligand | Ligand support |",
        "| pose_support_count | Poses in this pocket | Support strength |",
        "| best_affinity_kcal | Best docking score | Binding potential |",
        "| mean_affinity_kcal | Mean docking score | Average binding |",
        "| relationship_class | PPI relationship | Orthosteric/rim/allosteric |",
        "| overall_druggability_tier | Druggability | Tier 1/2/3 |",
        "| state_class | State robustness | Robust/shifted/specific |",
        "| ligand_support_strength | Derived support | strong/moderate/weak/none |",
        "| diversity_verdict | Search distribution | Quality flag |",
        "",
        "## Ligand Support Strength Classification",
        "",
        "| Level | Criteria |",
        "|-------|---------|",
        f"| strong | ≥{STRONG_SUPPORT_MIN_POSES} poses AND best affinity ≤ {STRONG_SUPPORT_MAX_AFFINITY} kcal/mol |",
        f"| moderate | ≥{MODERATE_SUPPORT_MIN_POSES} poses AND best affinity ≤ {MODERATE_SUPPORT_MAX_AFFINITY} kcal/mol |",
        "| weak | >0 poses but below moderate threshold |",
        "| none | 0 poses |",
        "| pending_* | Pre-execution (no affinity data yet) |",
        "",
        "## How to Use in Phase 4",
        "",
        "1. Load `phase4_docking_evidence_reference.csv`",
        "2. For each `candidate_pocket_id`:",
        "   - Check `relationship_class` for PPI relevance type",
        "   - Check `ligand_support_strength` for docking support level",
        "   - Check `state_class` for state robustness",
        "3. Combine with Phase 1 PPI hotspot data for perturbation scoring",
        "4. Pockets with `strong` ligand support + `orthosteric_candidate`",
        "   relationship are highest priority for perturbation analysis",
        "",
    ]

    # Evidence summary
    if evidence:
        lines.extend([
            "## Evidence Summary",
            "",
        ])

        by_pocket = defaultdict(list)
        for e in evidence:
            by_pocket[e["candidate_pocket_id"]].append(e)

        lines.extend([
            "| Pocket | Relationship | Tier | Ligands | Support |",
            "|--------|-------------|------|---------|---------|",
        ])
        for pocket_id in sorted(by_pocket.keys()):
            rows = by_pocket[pocket_id]
            rel = rows[0].get("relationship_class", "")
            tier = rows[0].get("overall_druggability_tier", "")
            ligands = sorted(set(r["ligand_id"] for r in rows))
            supports = sorted(set(r["ligand_support_strength"] for r in rows))
            lines.append(
                f"| {pocket_id} | {rel} | {tier} | "
                f"{', '.join(ligands)} | {', '.join(supports)} |"
            )
        lines.append("")

    # Validation
    lines.extend([
        "## Handoff Validation",
        "",
    ])
    for v in validation:
        lines.append(f"- {v}")
    lines.append("")

    # Caveats
    lines.extend([
        "## Caveats",
        "",
        "- Currently only 3GT8_raw state has pocket data (fpocket only)",
        "- Re-run Phase 2-3 when cl38_48 and cl85_100 pocket data available",
        "- Re-run after server-side Vina execution for real affinity values",
        "- `pending_*` support levels indicate pre-execution state",
        "",
        "---",
        "",
        "Generated by `egfr_pipeline.phase3.phase4_export`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_phase4_export(
    output_dir: Path = PHASE3_OUTPUT_DIR,
) -> Tuple[Path, Path]:
    """Run full TG 3.6 pipeline.

    Returns (evidence_path, handoff_note_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load all inputs
    pose_rows = _load_csv(output_dir / "vina_pose_table.csv")
    occupancy = _load_csv(output_dir / "phase3_pocket_occupancy_summary.csv")
    diversity = _load_csv(output_dir / "phase3_diversity_metrics.csv")
    statuses = _load_csv(output_dir / "pocket_search_status.csv")
    normalized = _load_csv(output_dir / "phase3_candidate_reference_normalized.csv")

    print(f"  Loaded: {len(pose_rows)} poses, {len(occupancy)} occupancy, "
          f"{len(diversity)} diversity, {len(normalized)} pockets")

    # 3.6.1: Build evidence table
    evidence = build_evidence_table(
        pose_rows, occupancy, diversity, statuses, normalized
    )
    print(f"  Built {len(evidence)} evidence rows")

    # 3.6.3: Validate handoff
    validation = validate_handoff(evidence)
    for v in validation:
        print(f"  {v}")

    # Write evidence CSV
    evidence_path = output_dir / "phase4_docking_evidence_reference.csv"
    _write_csv(evidence_path, EVIDENCE_COLUMNS, evidence)

    # Write handoff note
    note_lines = build_handoff_note(evidence, validation)
    note_path = output_dir / "phase3_to_phase4_handoff_note.md"
    with open(note_path, "w", encoding="utf-8") as f:
        f.write("\n".join(note_lines))

    return evidence_path, note_path


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
        description="Phase 3 TG 3.6: Phase 4-ready export"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE3_OUTPUT_DIR,
        help="Phase 3 output directory",
    )
    args = parser.parse_args()

    print("Phase 3 — Task Group 3.6: Phase 4-Ready Export")
    print(f"  Output: {args.output_dir}")
    print()

    evidence_path, note_path = run_phase4_export(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 3.6 COMPLETE")
    print(f"{'='*60}")
    print(f"  Evidence: {evidence_path}")
    print(f"  Handoff:  {note_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
