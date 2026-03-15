#!/usr/bin/env python3
"""Phase 3 Task Group 3.0: Phase 2 Candidate Pocket Reference Ingestion.

Loads and validates the Phase 2 candidate pocket reference so that all docking
jobs are anchored to explicitly prioritized receptor-local candidate pockets.

Tasks:
  3.0.1 - Candidate pocket reference ingestion
  3.0.2 - Reference validation (completeness + consistency checks)
  3.0.3 - Internal normalization for Phase 3 docking loops

Inputs:
  - output/phase2_pockets/phase3_candidate_pocket_reference.csv

Outputs:
  - output/phase3_docking/phase3_candidate_reference_normalized.csv
  - output/phase3_docking/phase3_candidate_reference_validation.md

Usage:
    python -m egfr_pipeline.phase3.pocket_reference_ingestion [--output_dir ...]
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
PHASE3_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase3_docking"
RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]

# Required fields from Phase 2 export
REQUIRED_FIELDS = [
    "receptor_id",
    "pocket_id",
    "centroid_x",
    "centroid_y",
    "centroid_z",
    "relationship_class",
    "docking_priority",
]

# Fields needed for docking box setup
BOX_FIELDS = [
    "box_size_x",
    "box_size_y",
    "box_size_z",
]

# Valid values
VALID_PRIORITIES = {"primary", "secondary", "exploratory", "skip"}
VALID_RELATIONSHIPS = {
    "orthosteric_candidate", "rim_candidate",
    "allosteric_candidate", "low_relevance_candidate",
}
VALID_TIERS = {"tier_1", "tier_2", "tier_3"}

# Phase 3 normalized schema — adds budget columns
NORMALIZED_COLUMNS = [
    # Identity
    "receptor_id",
    "pocket_id",
    # Docking geometry
    "centroid_x",
    "centroid_y",
    "centroid_z",
    "box_size_x",
    "box_size_y",
    "box_size_z",
    # PPI context
    "relationship_class",
    "hotspot_overlap_count",
    "hotspot_overlap_fraction",
    "centroid_distance_to_patch_A",
    # Druggability
    "druggability_confidence",
    "overall_druggability_tier",
    "best_druggability_score",
    # Support
    "n_sources",
    "multi_source_support",
    # State robustness
    "state_class",
    "n_states_matched",
    # Phase 3 docking guidance
    "docking_priority",
    "priority_basis",
    # Phase 3 budget fields (initialized here, consumed by TG 3.2)
    "recommended_exhaustiveness",
    "recommended_n_poses",
    "recommended_seeds",
    "dockable",                    # True/False — whether to include in docking
    "skip_reason",                 # why excluded, if dockable=False
]

# Default docking parameters by priority
BUDGET_DEFAULTS = {
    "primary":     {"exhaustiveness": 32, "n_poses": 20, "seeds": 3},
    "secondary":   {"exhaustiveness": 16, "n_poses": 10, "seeds": 1},
    "exploratory": {"exhaustiveness": 8,  "n_poses": 5,  "seeds": 1},
    "skip":        {"exhaustiveness": 0,  "n_poses": 0,  "seeds": 0},
}


# ---------------------------------------------------------------------------
# Validation (Task 3.0.2)
# ---------------------------------------------------------------------------

def validate_pocket_reference(rows: List[dict]) -> Tuple[List[str], List[str], List[str]]:
    """Validate Phase 2 pocket reference.

    Returns (passes, warnings, errors).
    """
    passes = []
    warnings = []
    errors = []

    # V1: Non-empty
    if not rows:
        errors.append("V1_EMPTY: No candidate pockets in Phase 2 reference")
        return passes, warnings, errors
    passes.append(f"V1_NON_EMPTY: {len(rows)} candidate pockets loaded")

    # V2: Required fields present
    if rows:
        available = set(rows[0].keys())
        missing = set(REQUIRED_FIELDS) - available
        if missing:
            errors.append(f"V2_MISSING_FIELDS: {sorted(missing)}")
        else:
            passes.append("V2_REQUIRED_FIELDS: All required fields present")

    # V3: Receptor state coverage
    states_present = set(r["receptor_id"] for r in rows)
    states_missing = set(RECEPTOR_STATES) - states_present
    if not states_missing:
        passes.append(f"V3_STATE_COVERAGE: All {len(RECEPTOR_STATES)} states have pockets")
    elif len(states_present) >= 1:
        warnings.append(
            f"V3_PARTIAL_STATES: {len(states_present)}/{len(RECEPTOR_STATES)} states "
            f"have pockets. Missing: {sorted(states_missing)}"
        )
    else:
        errors.append("V3_NO_STATES: No receptor states have pockets")

    # V4: Valid docking priorities
    invalid_priorities = set()
    for r in rows:
        pri = r.get("docking_priority", "")
        if pri not in VALID_PRIORITIES:
            invalid_priorities.add(pri)
    if invalid_priorities:
        errors.append(f"V4_INVALID_PRIORITY: {sorted(invalid_priorities)}")
    else:
        passes.append("V4_VALID_PRIORITIES: All docking priorities valid")

    # V5: Valid relationship classes
    invalid_rels = set()
    for r in rows:
        rel = r.get("relationship_class", "")
        if rel and rel not in VALID_RELATIONSHIPS:
            invalid_rels.add(rel)
    if invalid_rels:
        warnings.append(f"V5_UNKNOWN_RELATIONSHIP: {sorted(invalid_rels)}")
    else:
        passes.append("V5_VALID_RELATIONSHIPS: All relationship classes valid")

    # V6: Centroid coordinates parseable
    bad_centroids = []
    for r in rows:
        for dim in ("centroid_x", "centroid_y", "centroid_z"):
            try:
                float(r.get(dim, ""))
            except (ValueError, TypeError):
                bad_centroids.append(f"{r.get('pocket_id','?')}.{dim}")
    if bad_centroids:
        errors.append(f"V6_BAD_CENTROIDS: {len(bad_centroids)} unparseable coordinates")
    else:
        passes.append("V6_VALID_CENTROIDS: All centroid coordinates parseable")

    # V7: Box sizes present (warning only — can fall back to defaults)
    missing_box = []
    for r in rows:
        if r.get("docking_priority") == "skip":
            continue
        for dim in BOX_FIELDS:
            val = r.get(dim, "").strip()
            if not val:
                missing_box.append(f"{r.get('pocket_id','?')}.{dim}")
                break
    if missing_box:
        warnings.append(
            f"V7_MISSING_BOX: {len(missing_box)} dockable pocket(s) missing box sizes. "
            "Default box will be used."
        )
    else:
        passes.append("V7_BOX_SIZES: All dockable pockets have box sizes")

    # V8: Pocket ID uniqueness
    ids = [r["pocket_id"] for r in rows]
    if len(ids) != len(set(ids)):
        dupes = sorted(set(pid for pid in ids if ids.count(pid) > 1))
        errors.append(f"V8_DUPLICATE_IDS: {dupes}")
    else:
        passes.append(f"V8_UNIQUE_IDS: All {len(ids)} pocket IDs unique")

    # V9: At least one dockable pocket
    dockable = [r for r in rows if r.get("docking_priority") != "skip"]
    if not dockable:
        errors.append("V9_NO_DOCKABLE: All pockets are 'skip' — nothing to dock")
    else:
        passes.append(f"V9_DOCKABLE: {len(dockable)} pocket(s) eligible for docking")

    return passes, warnings, errors


# ---------------------------------------------------------------------------
# Normalization (Task 3.0.3)
# ---------------------------------------------------------------------------

DEFAULT_BOX_SIZE = 22.0  # Å — default per-dimension if missing


def normalize_for_phase3(rows: List[dict]) -> List[dict]:
    """Normalize Phase 2 reference into Phase 3 internal format.

    Adds budget recommendations and dockable flag.
    """
    normalized = []
    for r in rows:
        priority = r.get("docking_priority", "skip")
        budget = BUDGET_DEFAULTS.get(priority, BUDGET_DEFAULTS["skip"])

        dockable = priority != "skip"
        skip_reason = ""
        if not dockable:
            skip_reason = r.get("priority_basis", "skip priority")

        # Box size: use Phase 2 value or default
        box_x = _safe_float(r.get("box_size_x", "")) or DEFAULT_BOX_SIZE
        box_y = _safe_float(r.get("box_size_y", "")) or DEFAULT_BOX_SIZE
        box_z = _safe_float(r.get("box_size_z", "")) or DEFAULT_BOX_SIZE

        normalized.append({
            "receptor_id": r.get("receptor_id", ""),
            "pocket_id": r.get("pocket_id", ""),
            "centroid_x": r.get("centroid_x", ""),
            "centroid_y": r.get("centroid_y", ""),
            "centroid_z": r.get("centroid_z", ""),
            "box_size_x": f"{box_x:.1f}",
            "box_size_y": f"{box_y:.1f}",
            "box_size_z": f"{box_z:.1f}",
            "relationship_class": r.get("relationship_class", ""),
            "hotspot_overlap_count": r.get("hotspot_overlap_count", ""),
            "hotspot_overlap_fraction": r.get("hotspot_overlap_fraction", ""),
            "centroid_distance_to_patch_A": r.get("centroid_distance_to_patch_A", ""),
            "druggability_confidence": r.get("druggability_confidence", ""),
            "overall_druggability_tier": r.get("overall_druggability_tier", ""),
            "best_druggability_score": r.get("best_druggability_score", ""),
            "n_sources": r.get("n_sources", ""),
            "multi_source_support": r.get("multi_source_support", ""),
            "state_class": r.get("state_class", ""),
            "n_states_matched": r.get("n_states_matched", ""),
            "docking_priority": priority,
            "priority_basis": r.get("priority_basis", ""),
            "recommended_exhaustiveness": budget["exhaustiveness"],
            "recommended_n_poses": budget["n_poses"],
            "recommended_seeds": budget["seeds"],
            "dockable": str(dockable),
            "skip_reason": skip_reason,
        })

    return normalized


# ---------------------------------------------------------------------------
# Validation report
# ---------------------------------------------------------------------------

def build_validation_report(
    rows: List[dict],
    normalized: List[dict],
    passes: List[str],
    warnings: List[str],
    errors: List[str],
) -> List[str]:
    """Build the Phase 3 candidate reference validation markdown."""
    lines = [
        "# Phase 3 Candidate Pocket Reference Validation",
        "",
        "## Validation Summary",
        "",
        f"- Input pockets: {len(rows)}",
        f"- Checks passed: {len(passes)}",
        f"- Warnings: {len(warnings)}",
        f"- Errors: {len(errors)}",
        "",
    ]

    # Status
    if errors:
        lines.append("**Status: FAIL** — errors must be resolved before docking")
    elif warnings:
        lines.append("**Status: PASS with warnings** — docking can proceed with caveats")
    else:
        lines.append("**Status: PASS** — ready for Phase 3 docking")
    lines.append("")

    # Checks
    for section, items, icon in [
        ("Passed", passes, "PASS"),
        ("Warnings", warnings, "WARNING"),
        ("Errors", errors, "FAIL"),
    ]:
        if items:
            lines.append(f"### {section}")
            lines.append("")
            for item in items:
                lines.append(f"- [{icon}] {item}")
            lines.append("")

    # Normalized summary
    if normalized:
        by_state: Dict[str, List[dict]] = defaultdict(list)
        for n in normalized:
            by_state[n["receptor_id"]].append(n)

        lines.extend([
            "## Normalized Pocket Summary",
            "",
            "| State | Total | Dockable | Skip | Primary | Secondary | Exploratory |",
            "|-------|-------|----------|------|---------|-----------|-------------|",
        ])
        for state in sorted(by_state.keys()):
            pockets = by_state[state]
            dockable = [p for p in pockets if p["dockable"] == "True"]
            by_pri = defaultdict(int)
            for p in pockets:
                by_pri[p["docking_priority"]] += 1
            lines.append(
                f"| {state} | {len(pockets)} | {len(dockable)} | "
                f"{by_pri.get('skip', 0)} | {by_pri.get('primary', 0)} | "
                f"{by_pri.get('secondary', 0)} | {by_pri.get('exploratory', 0)} |"
            )
        lines.append("")

    # Budget defaults
    lines.extend([
        "## Default Budget Parameters by Priority",
        "",
        "| Priority | Exhaustiveness | N Poses | Seeds |",
        "|----------|---------------|---------|-------|",
    ])
    for pri in ["primary", "secondary", "exploratory", "skip"]:
        b = BUDGET_DEFAULTS[pri]
        lines.append(
            f"| {pri} | {b['exhaustiveness']} | {b['n_poses']} | {b['seeds']} |"
        )
    lines.append("")

    # Per-pocket detail
    if normalized:
        lines.extend([
            "## Per-Pocket Detail",
            "",
            "| Pocket | State | Priority | Dockable | Box (A) | Exhaustiveness |",
            "|--------|-------|----------|----------|---------|---------------|",
        ])
        for n in normalized:
            box = f"{n['box_size_x']}x{n['box_size_y']}x{n['box_size_z']}"
            lines.append(
                f"| {n['pocket_id']} | {n['receptor_id']} | {n['docking_priority']} | "
                f"{n['dockable']} | {box} | {n['recommended_exhaustiveness']} |"
            )
        lines.append("")

    lines.extend([
        "---",
        "",
        "Generated by `egfr_pipeline.phase3.pocket_reference_ingestion`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_pocket_reference_ingestion(
    phase2_dir: Path = PHASE2_OUTPUT_DIR,
    output_dir: Path = PHASE3_OUTPUT_DIR,
) -> Tuple[Path, Path]:
    """Run full TG 3.0 pipeline.

    Returns (normalized_csv, validation_note_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Task 3.0.1: Load Phase 2 export
    ref_path = phase2_dir / "phase3_candidate_pocket_reference.csv"
    if not ref_path.exists():
        raise FileNotFoundError(
            f"Phase 2 export not found: {ref_path}\n"
            "Run Phase 2 TG 2.6 (phase3_export) first."
        )

    rows = _load_csv(ref_path)
    print(f"  Loaded {len(rows)} candidate pockets from {ref_path.name}")

    # Task 3.0.2: Validate
    passes, warnings, errors = validate_pocket_reference(rows)
    print(f"  Validation: {len(passes)} pass, {len(warnings)} warn, {len(errors)} error")
    for w in warnings:
        print(f"    [WARNING] {w}")
    for e in errors:
        print(f"    [ERROR] {e}")

    # Task 3.0.3: Normalize
    normalized = normalize_for_phase3(rows)
    dockable = [n for n in normalized if n["dockable"] == "True"]
    print(f"  Normalized: {len(normalized)} pockets, {len(dockable)} dockable")

    # Write normalized CSV
    norm_path = output_dir / "phase3_candidate_reference_normalized.csv"
    _write_csv(norm_path, NORMALIZED_COLUMNS, normalized)

    # Write validation report
    report_lines = build_validation_report(rows, normalized, passes, warnings, errors)
    report_path = output_dir / "phase3_candidate_reference_validation.md"
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


def _safe_float(val) -> Optional[float]:
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 3 TG 3.0: Candidate pocket reference ingestion"
    )
    parser.add_argument(
        "--phase2_dir", type=Path, default=PHASE2_OUTPUT_DIR,
        help="Phase 2 output directory containing phase3_candidate_pocket_reference.csv",
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE3_OUTPUT_DIR,
        help="Phase 3 output directory",
    )
    args = parser.parse_args()

    print("Phase 3 — Task Group 3.0: Candidate Pocket Reference Ingestion")
    print(f"  Phase 2 input: {args.phase2_dir}")
    print(f"  Phase 3 output: {args.output_dir}")
    print()

    norm_path, report_path = run_pocket_reference_ingestion(
        args.phase2_dir, args.output_dir
    )

    print(f"\n{'='*60}")
    print("TG 3.0 COMPLETE")
    print(f"{'='*60}")
    print(f"  Normalized: {norm_path}")
    print(f"  Validation: {report_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
