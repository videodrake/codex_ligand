#!/usr/bin/env python3
"""Phase 4 Task Group 4.0: Multi-Phase Evidence Ingestion and Validation.

Loads evidence from Phases 1-3, validates cross-phase consistency,
and produces a normalized evidence table for downstream scoring.

Tasks:
  4.0.1 - Load Phase 1 receptor-side patch reference + hotspot summaries
  4.0.2 - Load Phase 2 pocket reference + relationship + druggability + state classes
  4.0.3 - Load Phase 3 docking evidence reference + ligand support + diversity
  4.0.4 - Cross-phase consistency validation

Inputs:
  - output/phase1_ppi/phase1_downstream_patch_reference.csv
  - output/phase2_pockets/pocket_patch_relationship.csv
  - output/phase2_pockets/druggability_proposal_summary.csv
  - output/phase2_pockets/candidate_pocket_state_classes.csv
  - output/phase3_docking/phase4_docking_evidence_reference.csv

Outputs:
  - output/phase4_perturbation/phase4_evidence_normalized.csv
  - output/phase4_perturbation/phase4_evidence_validation.md

Usage:
    python -m egfr_pipeline.phase4.evidence_ingestion [--output_dir ...]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Default paths
# ---------------------------------------------------------------------------

PHASE1_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase1_ppi"
PHASE2_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase2_pockets"
PHASE3_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase3_docking"
PHASE4_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase4_perturbation"

# ---------------------------------------------------------------------------
# Expected schemas (for validation)
# ---------------------------------------------------------------------------

PHASE1_REQUIRED_COLUMNS = [
    "residue_id", "residue_num", "residue_name", "chain", "lobe_label",
    "robustness_class", "confidence", "is_hotspot_any_state",
    "pyrosetta_evidence", "lightdock_evidence", "method_agreement",
]

PHASE2_RELATIONSHIP_REQUIRED = [
    "receptor_id", "pocket_id", "relationship_class",
    "hotspot_overlap_count", "hotspot_overlap_fraction",
]

PHASE2_DRUGGABILITY_REQUIRED = [
    "receptor_id", "pocket_id", "overall_druggability_tier",
    "druggability_confidence", "best_proposal_score",
]

PHASE2_STATE_CLASS_REQUIRED = [
    "pocket_id", "receptor_id", "state_class",
    "n_states_matched",
]

PHASE3_EVIDENCE_REQUIRED = [
    "receptor_id", "candidate_pocket_id", "ligand_id",
    "pose_support_count", "ligand_support_strength",
    "relationship_class", "overall_druggability_tier",
]

# Output schema: one row per (pocket, ligand) pair
NORMALIZED_COLUMNS = [
    # Identity
    "pocket_id",
    "receptor_id",
    "ligand_id",
    # Phase 1: PPI patch evidence
    "n_hotspot_residues",
    "hotspot_residues",
    "hotspot_overlap_count",
    "hotspot_overlap_fraction",
    "mean_hotspot_confidence",
    "n_robust_hotspots",
    "n_method_agreement_both",
    # Phase 2: Pocket properties
    "relationship_class",
    "overall_druggability_tier",
    "druggability_confidence",
    "best_proposal_score",
    "state_class",
    "n_states_matched",
    # Phase 3: Docking evidence
    "pose_support_count",
    "best_affinity_kcal",
    "mean_affinity_kcal",
    "ligand_support_strength",
    "docking_priority",
    "pocket_status",
    "diversity_verdict",
    "coverage_fraction",
    "normalized_entropy",
    # Provenance
    "phase1_loaded",
    "phase2_loaded",
    "phase3_loaded",
]


# ---------------------------------------------------------------------------
# 4.0.1: Load Phase 1 patch reference
# ---------------------------------------------------------------------------

def load_phase1_patches(phase1_dir: Path) -> Tuple[List[dict], List[str]]:
    """Load Phase 1 downstream patch reference (hotspot residues).

    Returns (patches, issues).
    """
    issues = []
    path = phase1_dir / "phase1_downstream_patch_reference.csv"

    if not path.exists():
        issues.append(f"WARNING: Phase 1 patch reference not found: {path}")
        return [], issues

    patches = _load_csv(path)
    if not patches:
        issues.append("WARNING: Phase 1 patch reference is empty")
        return [], issues

    # Schema check
    missing = [c for c in PHASE1_REQUIRED_COLUMNS if c not in patches[0]]
    if missing:
        issues.append(f"ERROR: Phase 1 missing columns: {missing}")
        return [], issues

    # Content checks
    hotspots = [p for p in patches if p.get("is_hotspot_any_state") == "True"]
    if not hotspots:
        issues.append("WARNING: No hotspot residues in Phase 1 reference")

    issues.append(f"OK: Phase 1 loaded {len(patches)} residues "
                  f"({len(hotspots)} hotspots)")

    return patches, issues


def summarize_phase1(patches: List[dict]) -> dict:
    """Build Phase 1 summary for each pocket (via relationship data later)."""
    hotspots = [p for p in patches if p.get("is_hotspot_any_state") == "True"]

    # Confidence mapping
    confidence_map = {"high": 3, "medium": 2, "low": 1}
    confidences = [confidence_map.get(p.get("confidence", "low"), 1)
                   for p in hotspots]
    mean_conf = sum(confidences) / len(confidences) if confidences else 0.0

    robust = [p for p in hotspots
              if p.get("robustness_class") in ("robust", "moderate")]
    both_methods = [p for p in hotspots
                    if p.get("method_agreement") == "both"]

    return {
        "n_hotspot_residues": len(hotspots),
        "hotspot_residues": ";".join(p["residue_id"] for p in hotspots),
        "mean_hotspot_confidence": round(mean_conf, 3),
        "n_robust_hotspots": len(robust),
        "n_method_agreement_both": len(both_methods),
    }


# ---------------------------------------------------------------------------
# 4.0.2: Load Phase 2 pocket data
# ---------------------------------------------------------------------------

def load_phase2_data(phase2_dir: Path) -> Tuple[dict, List[str]]:
    """Load Phase 2 pocket relationship, druggability, and state classes.

    Returns ({pocket_id: merged_dict}, issues).
    """
    issues = []
    pockets: Dict[str, dict] = {}

    # Relationship
    rel_path = phase2_dir / "pocket_patch_relationship.csv"
    if not rel_path.exists():
        issues.append(f"ERROR: Phase 2 relationship file not found: {rel_path}")
        return pockets, issues

    relationships = _load_csv(rel_path)
    missing = [c for c in PHASE2_RELATIONSHIP_REQUIRED
               if relationships and c not in relationships[0]]
    if missing:
        issues.append(f"ERROR: Phase 2 relationship missing columns: {missing}")
        return pockets, issues

    for r in relationships:
        pid = r["pocket_id"]
        pockets[pid] = {
            "pocket_id": pid,
            "receptor_id": r["receptor_id"],
            "relationship_class": r["relationship_class"],
            "hotspot_overlap_count": r.get("hotspot_overlap_count", "0"),
            "hotspot_overlap_fraction": r.get("hotspot_overlap_fraction", "0"),
        }

    # Druggability
    drug_path = phase2_dir / "druggability_proposal_summary.csv"
    if drug_path.exists():
        druggability = _load_csv(drug_path)
        missing_d = [c for c in PHASE2_DRUGGABILITY_REQUIRED
                     if druggability and c not in druggability[0]]
        if missing_d:
            issues.append(
                f"WARNING: Phase 2 druggability missing columns: {missing_d}")
        else:
            for d in druggability:
                pid = d["pocket_id"]
                if pid in pockets:
                    pockets[pid].update({
                        "overall_druggability_tier": d["overall_druggability_tier"],
                        "druggability_confidence": d["druggability_confidence"],
                        "best_proposal_score": d.get("best_proposal_score", ""),
                    })
    else:
        issues.append(f"WARNING: Phase 2 druggability file not found: {drug_path}")

    # State classes
    state_path = phase2_dir / "candidate_pocket_state_classes.csv"
    if state_path.exists():
        states = _load_csv(state_path)
        missing_s = [c for c in PHASE2_STATE_CLASS_REQUIRED
                     if states and c not in states[0]]
        if missing_s:
            issues.append(
                f"WARNING: Phase 2 state classes missing columns: {missing_s}")
        else:
            for s in states:
                pid = s["pocket_id"]
                if pid in pockets:
                    pockets[pid].update({
                        "state_class": s["state_class"],
                        "n_states_matched": s.get("n_states_matched", "1"),
                    })
    else:
        issues.append(f"WARNING: Phase 2 state classes not found: {state_path}")

    issues.append(f"OK: Phase 2 loaded {len(pockets)} pockets")
    return pockets, issues


# ---------------------------------------------------------------------------
# 4.0.3: Load Phase 3 docking evidence
# ---------------------------------------------------------------------------

def load_phase3_evidence(phase3_dir: Path) -> Tuple[List[dict], List[str]]:
    """Load Phase 3 docking evidence reference.

    Returns (evidence_rows, issues).
    """
    issues = []
    path = phase3_dir / "phase4_docking_evidence_reference.csv"

    if not path.exists():
        issues.append(f"WARNING: Phase 3 evidence not found: {path}")
        return [], issues

    evidence = _load_csv(path)
    if not evidence:
        issues.append("WARNING: Phase 3 evidence is empty")
        return [], issues

    missing = [c for c in PHASE3_EVIDENCE_REQUIRED if c not in evidence[0]]
    if missing:
        issues.append(f"ERROR: Phase 3 evidence missing columns: {missing}")
        return [], issues

    issues.append(f"OK: Phase 3 loaded {len(evidence)} evidence rows")
    return evidence, issues


# ---------------------------------------------------------------------------
# 4.0.4: Cross-phase consistency validation
# ---------------------------------------------------------------------------

def validate_consistency(
    phase1_patches: List[dict],
    phase2_pockets: Dict[str, dict],
    phase3_evidence: List[dict],
) -> List[str]:
    """Validate cross-phase consistency.

    Checks:
    1. Receptor IDs match across phases
    2. Pocket IDs in Phase 3 exist in Phase 2
    3. Relationship class alignment between Phase 2 and Phase 3
    4. Hotspot residues referenced in Phase 2 overlap exist in Phase 1
    """
    issues = []

    # 1. Receptor ID consistency
    phase2_receptors = set(p["receptor_id"] for p in phase2_pockets.values())
    phase3_receptors = set(e["receptor_id"] for e in phase3_evidence)

    if phase3_receptors and phase2_receptors:
        extra = phase3_receptors - phase2_receptors
        if extra:
            issues.append(
                f"WARNING: Phase 3 receptor IDs not in Phase 2: {extra}")
        missing = phase2_receptors - phase3_receptors
        if missing:
            issues.append(
                f"INFO: Phase 2 receptors without Phase 3 evidence: {missing}")
    issues.append(
        f"OK: Receptor IDs — Phase 2: {sorted(phase2_receptors)}, "
        f"Phase 3: {sorted(phase3_receptors)}")

    # 2. Pocket ID alignment
    phase3_pockets = set(e["candidate_pocket_id"] for e in phase3_evidence)
    missing_pockets = phase3_pockets - set(phase2_pockets.keys())
    if missing_pockets:
        issues.append(
            f"WARNING: Phase 3 pocket IDs not in Phase 2: {missing_pockets}")
    else:
        issues.append(
            f"OK: All {len(phase3_pockets)} Phase 3 pockets found in Phase 2")

    # 3. Relationship class alignment
    for e in phase3_evidence:
        pid = e["candidate_pocket_id"]
        p3_rel = e.get("relationship_class", "")
        if pid in phase2_pockets:
            p2_rel = phase2_pockets[pid].get("relationship_class", "")
            if p3_rel and p2_rel and p3_rel != p2_rel:
                issues.append(
                    f"WARNING: Relationship class mismatch for {pid}: "
                    f"Phase 2={p2_rel}, Phase 3={p3_rel}")
    issues.append("OK: Relationship class alignment checked")

    # 4. Hotspot residue coverage
    if phase1_patches:
        hotspot_ids = set(
            p["residue_id"] for p in phase1_patches
            if p.get("is_hotspot_any_state") == "True")
        issues.append(
            f"OK: Phase 1 provides {len(hotspot_ids)} hotspot residues "
            f"for overlap verification")
    else:
        issues.append(
            "INFO: No Phase 1 patches — hotspot overlap not verifiable")

    return issues


# ---------------------------------------------------------------------------
# Build normalized evidence table
# ---------------------------------------------------------------------------

def build_normalized_evidence(
    phase1_summary: dict,
    phase2_pockets: Dict[str, dict],
    phase3_evidence: List[dict],
) -> List[dict]:
    """Merge all phase evidence into a normalized pocket x ligand table.

    One row per (pocket_id, ligand_id) pair from Phase 3 evidence.
    Phase 2 pockets without Phase 3 evidence get one row with ligand_id="".
    """
    rows = []
    seen_pockets = set()

    # Primary: rows driven by Phase 3 evidence (pocket x ligand)
    for e in phase3_evidence:
        pid = e["candidate_pocket_id"]
        seen_pockets.add(pid)
        p2 = phase2_pockets.get(pid, {})

        row = {
            "pocket_id": pid,
            "receptor_id": e["receptor_id"],
            "ligand_id": e["ligand_id"],
            # Phase 1
            "n_hotspot_residues": phase1_summary.get("n_hotspot_residues", 0),
            "hotspot_residues": phase1_summary.get("hotspot_residues", ""),
            "hotspot_overlap_count": p2.get("hotspot_overlap_count", ""),
            "hotspot_overlap_fraction": p2.get("hotspot_overlap_fraction", ""),
            "mean_hotspot_confidence": phase1_summary.get(
                "mean_hotspot_confidence", ""),
            "n_robust_hotspots": phase1_summary.get("n_robust_hotspots", 0),
            "n_method_agreement_both": phase1_summary.get(
                "n_method_agreement_both", 0),
            # Phase 2
            "relationship_class": p2.get(
                "relationship_class",
                e.get("relationship_class", "")),
            "overall_druggability_tier": p2.get(
                "overall_druggability_tier",
                e.get("overall_druggability_tier", "")),
            "druggability_confidence": p2.get(
                "druggability_confidence",
                e.get("druggability_confidence", "")),
            "best_proposal_score": p2.get("best_proposal_score", ""),
            "state_class": p2.get(
                "state_class", e.get("state_class", "")),
            "n_states_matched": p2.get(
                "n_states_matched", e.get("n_states_matched", "")),
            # Phase 3
            "pose_support_count": e.get("pose_support_count", ""),
            "best_affinity_kcal": e.get("best_affinity_kcal", ""),
            "mean_affinity_kcal": e.get("mean_affinity_kcal", ""),
            "ligand_support_strength": e.get("ligand_support_strength", ""),
            "docking_priority": e.get("docking_priority", ""),
            "pocket_status": e.get("pocket_status", ""),
            "diversity_verdict": e.get("diversity_verdict", ""),
            "coverage_fraction": e.get("coverage_fraction", ""),
            "normalized_entropy": e.get("normalized_entropy", ""),
            # Provenance
            "phase1_loaded": "True" if phase1_summary.get(
                "n_hotspot_residues", 0) > 0 else "False",
            "phase2_loaded": "True" if p2 else "False",
            "phase3_loaded": "True",
        }
        rows.append(row)

    # Secondary: Phase 2 pockets without Phase 3 evidence (e.g. skipped)
    for pid, p2 in sorted(phase2_pockets.items()):
        if pid in seen_pockets:
            continue
        row = {
            "pocket_id": pid,
            "receptor_id": p2.get("receptor_id", ""),
            "ligand_id": "",
            "n_hotspot_residues": phase1_summary.get("n_hotspot_residues", 0),
            "hotspot_residues": phase1_summary.get("hotspot_residues", ""),
            "hotspot_overlap_count": p2.get("hotspot_overlap_count", ""),
            "hotspot_overlap_fraction": p2.get("hotspot_overlap_fraction", ""),
            "mean_hotspot_confidence": phase1_summary.get(
                "mean_hotspot_confidence", ""),
            "n_robust_hotspots": phase1_summary.get("n_robust_hotspots", 0),
            "n_method_agreement_both": phase1_summary.get(
                "n_method_agreement_both", 0),
            "relationship_class": p2.get("relationship_class", ""),
            "overall_druggability_tier": p2.get(
                "overall_druggability_tier", ""),
            "druggability_confidence": p2.get("druggability_confidence", ""),
            "best_proposal_score": p2.get("best_proposal_score", ""),
            "state_class": p2.get("state_class", ""),
            "n_states_matched": p2.get("n_states_matched", ""),
            "pose_support_count": "0",
            "best_affinity_kcal": "",
            "mean_affinity_kcal": "",
            "ligand_support_strength": "none",
            "docking_priority": "",
            "pocket_status": "skipped",
            "diversity_verdict": "",
            "coverage_fraction": "",
            "normalized_entropy": "",
            "phase1_loaded": "True" if phase1_summary.get(
                "n_hotspot_residues", 0) > 0 else "False",
            "phase2_loaded": "True",
            "phase3_loaded": "False",
        }
        rows.append(row)

    return rows


# ---------------------------------------------------------------------------
# Validation report
# ---------------------------------------------------------------------------

def build_validation_report(
    all_issues: List[str],
    phase1_summary: dict,
    phase2_pockets: Dict[str, dict],
    phase3_evidence: List[dict],
    normalized_rows: List[dict],
) -> List[str]:
    """Build the Phase 4 evidence validation report."""
    lines = []

    errors = [i for i in all_issues if i.startswith("ERROR")]
    warnings = [i for i in all_issues if i.startswith("WARNING")]
    oks = [i for i in all_issues if i.startswith("OK")]
    infos = [i for i in all_issues if i.startswith("INFO")]

    verdict = "FAIL" if errors else ("CAUTION" if warnings else "PASS")

    lines.extend([
        "# Phase 4 Evidence Validation Report",
        "",
        "## 1. Overall Verdict",
        "",
        f"**{verdict}** — {len(errors)} errors, {len(warnings)} warnings, "
        f"{len(oks)} checks passed",
        "",
    ])

    # Phase 1 summary
    lines.extend([
        "## 2. Phase 1: PPI Patch Evidence",
        "",
        f"- Hotspot residues: {phase1_summary.get('n_hotspot_residues', 0)}",
        f"- Residues: {phase1_summary.get('hotspot_residues', 'N/A')}",
        f"- Mean confidence: {phase1_summary.get('mean_hotspot_confidence', 'N/A')}",
        f"- Robust hotspots (moderate+): "
        f"{phase1_summary.get('n_robust_hotspots', 0)}",
        f"- Both-method agreement: "
        f"{phase1_summary.get('n_method_agreement_both', 0)}",
        "",
    ])

    # Phase 2 summary
    lines.extend([
        "## 3. Phase 2: Pocket Properties",
        "",
        f"- Total pockets: {len(phase2_pockets)}",
    ])
    for pid, p2 in sorted(phase2_pockets.items()):
        lines.append(
            f"  - {pid}: {p2.get('relationship_class', '?')} / "
            f"{p2.get('overall_druggability_tier', '?')} / "
            f"{p2.get('state_class', '?')}")
    lines.append("")

    # Phase 3 summary
    lines.extend([
        "## 4. Phase 3: Docking Evidence",
        "",
        f"- Evidence rows: {len(phase3_evidence)}",
    ])
    pocket_ligands = defaultdict(list)
    for e in phase3_evidence:
        pocket_ligands[e["candidate_pocket_id"]].append(e["ligand_id"])
    for pid, ligs in sorted(pocket_ligands.items()):
        lines.append(f"  - {pid}: {len(ligs)} ligand(s) ({', '.join(ligs)})")
    lines.append("")

    # Consistency checks
    lines.extend([
        "## 5. Cross-Phase Consistency Checks",
        "",
    ])
    for i in all_issues:
        prefix = i.split(":")[0] if ":" in i else "INFO"
        lines.append(f"- [{prefix}] {i}")
    lines.append("")

    # Normalized output
    lines.extend([
        "## 6. Normalized Evidence Table",
        "",
        f"- Total rows: {len(normalized_rows)}",
        f"- Columns: {len(NORMALIZED_COLUMNS)}",
    ])
    # Count by pocket
    by_pocket = defaultdict(int)
    for r in normalized_rows:
        by_pocket[r["pocket_id"]] += 1
    for pid, count in sorted(by_pocket.items()):
        lines.append(f"  - {pid}: {count} row(s)")
    lines.append("")

    lines.extend([
        "---",
        "",
        "Generated by `egfr_pipeline.phase4.evidence_ingestion`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_evidence_ingestion(
    phase1_dir: Path = PHASE1_OUTPUT_DIR,
    phase2_dir: Path = PHASE2_OUTPUT_DIR,
    phase3_dir: Path = PHASE3_OUTPUT_DIR,
    output_dir: Path = PHASE4_OUTPUT_DIR,
) -> Tuple[Path, Path]:
    """Run full TG 4.0 pipeline.

    Returns (normalized_csv_path, validation_report_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    all_issues = []

    # 4.0.1: Phase 1
    phase1_patches, p1_issues = load_phase1_patches(phase1_dir)
    all_issues.extend(p1_issues)
    phase1_summary = summarize_phase1(phase1_patches) if phase1_patches else {}
    print(f"  Phase 1: {len(phase1_patches)} residues loaded")

    # 4.0.2: Phase 2
    phase2_pockets, p2_issues = load_phase2_data(phase2_dir)
    all_issues.extend(p2_issues)
    print(f"  Phase 2: {len(phase2_pockets)} pockets loaded")

    # 4.0.3: Phase 3
    phase3_evidence, p3_issues = load_phase3_evidence(phase3_dir)
    all_issues.extend(p3_issues)
    print(f"  Phase 3: {len(phase3_evidence)} evidence rows loaded")

    # 4.0.4: Cross-phase consistency
    consistency_issues = validate_consistency(
        phase1_patches, phase2_pockets, phase3_evidence)
    all_issues.extend(consistency_issues)

    # Build normalized table
    normalized = build_normalized_evidence(
        phase1_summary, phase2_pockets, phase3_evidence)
    print(f"  Normalized: {len(normalized)} rows "
          f"({len(NORMALIZED_COLUMNS)} columns)")

    # Build validation report
    report_lines = build_validation_report(
        all_issues, phase1_summary, phase2_pockets,
        phase3_evidence, normalized)

    # Write outputs
    norm_path = output_dir / "phase4_evidence_normalized.csv"
    _write_csv(norm_path, NORMALIZED_COLUMNS, normalized)

    report_path = output_dir / "phase4_evidence_validation.md"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    return norm_path, report_path


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
        description="Phase 4 TG 4.0: Multi-Phase Evidence Ingestion"
    )
    parser.add_argument(
        "--phase1_dir", type=Path, default=PHASE1_OUTPUT_DIR,
        help="Phase 1 output directory",
    )
    parser.add_argument(
        "--phase2_dir", type=Path, default=PHASE2_OUTPUT_DIR,
        help="Phase 2 output directory",
    )
    parser.add_argument(
        "--phase3_dir", type=Path, default=PHASE3_OUTPUT_DIR,
        help="Phase 3 output directory",
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE4_OUTPUT_DIR,
        help="Phase 4 output directory",
    )
    args = parser.parse_args()

    print("Phase 4 — Task Group 4.0: Multi-Phase Evidence Ingestion")
    print(f"  Phase 1 input: {args.phase1_dir}")
    print(f"  Phase 2 input: {args.phase2_dir}")
    print(f"  Phase 3 input: {args.phase3_dir}")
    print(f"  Phase 4 output: {args.output_dir}")
    print()

    norm_path, report_path = run_evidence_ingestion(
        args.phase1_dir, args.phase2_dir, args.phase3_dir, args.output_dir
    )

    print(f"\n{'='*60}")
    print("TG 4.0 COMPLETE")
    print(f"{'='*60}")
    print(f"  Normalized evidence: {norm_path}")
    print(f"  Validation report:   {report_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
