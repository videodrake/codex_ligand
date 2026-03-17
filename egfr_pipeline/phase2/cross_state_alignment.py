#!/usr/bin/env python3
"""Phase 2 Task Group 2.5: Cross-State Pocket Proposal Alignment.

Compares candidate pocket proposals across receptor states to identify
pockets that are conserved, shifted, or state-specific.

Tasks:
  2.5.1 - Receptor-state pocket alignment (centroid + residue overlap)
  2.5.2 - Cross-state proposal categories (robust/shifted/specific/uncertain)
  2.5.3 - Preserve raw metrics (all pairwise comparison data kept)

Inputs:
  - output/phase2_pockets/candidate_pockets.csv
  - output/phase2_pockets/druggability_proposal_summary.csv (optional, for tier)

Outputs:
  - output/phase2_pockets/candidate_pocket_cross_state_comparison.csv
  - output/phase2_pockets/candidate_pocket_state_classes.csv
  - output/phase2_pockets/phase2_cross_state_alignment_note.md

Usage:
    python -m egfr_pipeline.phase2.cross_state_alignment [--output_dir ...]
"""

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE2_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase2_pockets"
RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]

# Cross-state alignment thresholds
CENTROID_MATCH_DISTANCE_A = 8.0    # ≤ 8Å → candidate match
CENTROID_SHIFTED_DISTANCE_A = 15.0 # 8–15Å → shifted
RESIDUE_OVERLAP_MATCH = 0.20       # Jaccard ≥ 0.20 → candidate match

# State class definitions
# state_robust_pocket:   matched in all states with data
# state_shifted_pocket:  matched but centroid shifted (8–15Å) in ≥1 state
# state_specific_pocket: found in only one state
# uncertain_alignment:   ambiguous match quality
N_STATES = len(RECEPTOR_STATES)

# Output schemas
COMPARISON_COLUMNS = [
    "pocket_id_a",
    "receptor_id_a",
    "pocket_id_b",
    "receptor_id_b",
    "centroid_distance_A",
    "residue_jaccard",
    "shared_residue_count",
    "shared_residues",
    "match_type",             # match / shifted / no_match
    "match_basis",            # which metric(s) drove the decision
]

STATE_CLASS_COLUMNS = [
    "pocket_id",
    "receptor_id",
    "state_class",            # state_robust / state_shifted / state_specific / uncertain
    "n_states_matched",
    "matched_states",         # semicolon-delimited
    "matched_pocket_ids",     # semicolon-delimited partner pocket IDs
    "best_match_distance_A",
    "best_match_jaccard",
    "state_class_basis",
    # Carry forward from druggability summary if available
    "druggability_confidence",
    "overall_druggability_tier",
    "relationship_class",
]


# ---------------------------------------------------------------------------
# Pairwise cross-state comparison (Task 2.5.1)
# ---------------------------------------------------------------------------

def compare_pockets_across_states(
    pockets: List[dict],
) -> List[dict]:
    """Compare all pocket pairs across different receptor states.

    Returns list of pairwise comparison rows.
    """
    # Group by state
    by_state: Dict[str, List[dict]] = defaultdict(list)
    for p in pockets:
        by_state[p["receptor_id"]].append(p)

    states_with_data = sorted(by_state.keys())
    comparisons = []

    for i, state_a in enumerate(states_with_data):
        for state_b in states_with_data[i + 1:]:
            for pa in by_state[state_a]:
                for pb in by_state[state_b]:
                    comp = _compare_pocket_pair(pa, pb)
                    comparisons.append(comp)

    return comparisons


def _compare_pocket_pair(pa: dict, pb: dict) -> dict:
    """Compare one pair of pockets from different states."""
    # Centroid distance
    dist = _centroid_distance(pa, pb)

    # Residue overlap
    res_a = _parse_residue_nums(pa.get("residue_nums", ""))
    res_b = _parse_residue_nums(pb.get("residue_nums", ""))
    shared = res_a & res_b
    union = res_a | res_b
    jaccard = len(shared) / len(union) if union else 0.0

    # Match type
    match_type, basis = _classify_match(dist, jaccard, len(shared))

    return {
        "pocket_id_a": pa.get("pocket_id", ""),
        "receptor_id_a": pa.get("receptor_id", ""),
        "pocket_id_b": pb.get("pocket_id", ""),
        "receptor_id_b": pb.get("receptor_id", ""),
        "centroid_distance_A": f"{dist:.2f}" if dist is not None else "",
        "residue_jaccard": f"{jaccard:.3f}",
        "shared_residue_count": len(shared),
        "shared_residues": ";".join(str(r) for r in sorted(shared)),
        "match_type": match_type,
        "match_basis": basis,
    }


def _classify_match(
    dist: Optional[float],
    jaccard: float,
    n_shared: int,
) -> Tuple[str, str]:
    """Classify a cross-state pocket pair match."""
    reasons = []

    # Check centroid distance
    centroid_match = False
    centroid_shifted = False
    if dist is not None:
        if dist <= CENTROID_MATCH_DISTANCE_A:
            centroid_match = True
            reasons.append(f"centroid={dist:.1f}A <= {CENTROID_MATCH_DISTANCE_A}A")
        elif dist <= CENTROID_SHIFTED_DISTANCE_A:
            centroid_shifted = True
            reasons.append(
                f"centroid={dist:.1f}A ({CENTROID_MATCH_DISTANCE_A}-{CENTROID_SHIFTED_DISTANCE_A}A)"
            )

    # Check residue overlap
    residue_match = jaccard >= RESIDUE_OVERLAP_MATCH
    if residue_match:
        reasons.append(f"jaccard={jaccard:.3f} >= {RESIDUE_OVERLAP_MATCH}")

    # Decision
    if centroid_match or residue_match:
        return "match", "; ".join(reasons)
    if centroid_shifted and n_shared > 0:
        return "shifted", "; ".join(reasons) + f"; shared={n_shared}"
    if centroid_shifted:
        return "shifted", "; ".join(reasons)

    basis = []
    if dist is not None:
        basis.append(f"centroid={dist:.1f}A > {CENTROID_SHIFTED_DISTANCE_A}A")
    basis.append(f"jaccard={jaccard:.3f} < {RESIDUE_OVERLAP_MATCH}")
    return "no_match", "; ".join(basis)


# ---------------------------------------------------------------------------
# State class assignment (Task 2.5.2)
# ---------------------------------------------------------------------------

def assign_state_classes(
    pockets: List[dict],
    comparisons: List[dict],
    druggability: Dict[str, dict],
) -> List[dict]:
    """Assign cross-state class to each pocket.

    Returns state class rows.
    """
    # Build match graph: pocket_id → list of (matched_pocket_id, match_type)
    matches: Dict[str, List[Tuple[str, str, str, float, float]]] = defaultdict(list)
    # (matched_id, matched_state, match_type, dist, jaccard)

    for c in comparisons:
        mt = c["match_type"]
        if mt in ("match", "shifted"):
            dist = _safe_float(c["centroid_distance_A"])
            jacc = _safe_float(c["residue_jaccard"])
            matches[c["pocket_id_a"]].append((
                c["pocket_id_b"], c["receptor_id_b"], mt,
                dist if dist is not None else 999.0,
                jacc if jacc is not None else 0.0,
            ))
            matches[c["pocket_id_b"]].append((
                c["pocket_id_a"], c["receptor_id_a"], mt,
                dist if dist is not None else 999.0,
                jacc if jacc is not None else 0.0,
            ))

    # Count states with pocket data
    states_present = set(p["receptor_id"] for p in pockets)
    n_states_with_data = len(states_present)

    results = []
    for p in pockets:
        pid = p["pocket_id"]
        state = p["receptor_id"]
        pocket_matches = matches.get(pid, [])

        # Matched states (distinct)
        matched_states = set()
        matched_ids = []
        best_dist = None
        best_jacc = 0.0
        any_shifted = False

        for mid, mstate, mtype, mdist, mjacc in pocket_matches:
            matched_states.add(mstate)
            matched_ids.append(mid)
            if mtype == "shifted":
                any_shifted = True
            if best_dist is None or mdist < best_dist:
                best_dist = mdist
            if mjacc > best_jacc:
                best_jacc = mjacc

        # +1 for the pocket's own state
        n_matched = len(matched_states) + 1

        # State class
        state_class, basis = _assign_state_class(
            n_matched, n_states_with_data, any_shifted, pocket_matches
        )

        # Druggability info if available
        drug_info = druggability.get(pid, {})

        results.append({
            "pocket_id": pid,
            "receptor_id": state,
            "state_class": state_class,
            "n_states_matched": n_matched,
            "matched_states": ";".join(sorted(matched_states)) if matched_states else "",
            "matched_pocket_ids": ";".join(matched_ids) if matched_ids else "",
            "best_match_distance_A": f"{best_dist:.2f}" if best_dist is not None else "",
            "best_match_jaccard": f"{best_jacc:.3f}" if best_jacc > 0 else "",
            "state_class_basis": basis,
            "druggability_confidence": drug_info.get("druggability_confidence", ""),
            "overall_druggability_tier": drug_info.get("overall_druggability_tier", ""),
            "relationship_class": drug_info.get("relationship_class", ""),
        })

    return results


def _assign_state_class(
    n_matched: int,
    n_states_with_data: int,
    any_shifted: bool,
    pocket_matches: list,
) -> Tuple[str, str]:
    """Assign state class label."""
    # Only one state has data — cannot do cross-state comparison
    if n_states_with_data <= 1:
        return (
            "state_specific_pocket",
            f"only {n_states_with_data} state(s) with pocket data; cross-state comparison not yet possible"
        )

    # Matched in all states
    if n_matched >= n_states_with_data:
        if any_shifted:
            return (
                "state_shifted_pocket",
                f"matched in {n_matched}/{n_states_with_data} states but shifted in ≥1"
            )
        return (
            "state_robust_pocket",
            f"matched in all {n_matched}/{n_states_with_data} states"
        )

    # Matched in some states
    if n_matched > 1:
        if any_shifted:
            return (
                "state_shifted_pocket",
                f"matched in {n_matched}/{n_states_with_data} states, shifted"
            )
        # Matched in 2 of 3 — uncertain vs robust
        if n_states_with_data >= 3 and n_matched == 2:
            return (
                "uncertain_alignment",
                f"matched in {n_matched}/{n_states_with_data} states (not all)"
            )
        return (
            "state_robust_pocket",
            f"matched in {n_matched}/{n_states_with_data} states"
        )

    # No matches — state-specific
    return (
        "state_specific_pocket",
        f"no cross-state match found (1/{n_states_with_data} states)"
    )


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def build_alignment_note(
    comparisons: List[dict],
    state_classes: List[dict],
    states_present: Set[str],
) -> List[str]:
    """Build the cross-state alignment markdown note."""
    lines = [
        "# Phase 2 Cross-State Pocket Proposal Alignment Note",
        "",
        "## Overview",
        "",
        f"- Receptor states with pocket data: {len(states_present)}",
        f"  - States: {', '.join(sorted(states_present))}",
        f"  - Missing: {', '.join(sorted(set(RECEPTOR_STATES) - states_present)) or 'none'}",
        f"- Pairwise comparisons: {len(comparisons)}",
        f"- Pockets classified: {len(state_classes)}",
        "",
    ]

    if len(states_present) < 2:
        lines.extend([
            "## Cross-State Comparison Status",
            "",
            "**WARNING:** Only one receptor state has pocket data. Cross-state",
            "comparison requires at least two states. All pockets are currently",
            "classified as `state_specific_pocket` by default.",
            "",
            "Re-run this module after fpocket/P2Rank results are available for",
            "additional receptor states.",
            "",
        ])

    # Alignment thresholds
    lines.extend([
        "## Alignment Thresholds",
        "",
        "| Parameter | Value |",
        "|-----------|-------|",
        f"| Centroid match distance | ≤ {CENTROID_MATCH_DISTANCE_A} Å |",
        f"| Centroid shifted range | {CENTROID_MATCH_DISTANCE_A}–{CENTROID_SHIFTED_DISTANCE_A} Å |",
        f"| Residue Jaccard match | ≥ {RESIDUE_OVERLAP_MATCH} |",
        "",
    ])

    # State class definitions
    lines.extend([
        "## State Class Definitions",
        "",
        "| Class | Criterion |",
        "|-------|-----------|",
        "| state_robust_pocket | Matched in all states with data |",
        "| state_shifted_pocket | Matched but centroid shifted in ≥1 state |",
        "| state_specific_pocket | Found in only one state |",
        "| uncertain_alignment | Matched in some but not all states |",
        "",
    ])

    # State class distribution
    by_class: Dict[str, List[str]] = defaultdict(list)
    for sc in state_classes:
        by_class[sc["state_class"]].append(sc["pocket_id"])

    lines.extend([
        "## State Class Distribution",
        "",
        "| Class | Count | Pockets |",
        "|-------|-------|---------|",
    ])
    for cls in ["state_robust_pocket", "state_shifted_pocket",
                "state_specific_pocket", "uncertain_alignment"]:
        pids = by_class.get(cls, [])
        lines.append(f"| {cls} | {len(pids)} | {', '.join(pids) or '—'} |")
    lines.append("")

    # Per-pocket detail
    if state_classes:
        lines.extend([
            "## Per-Pocket Detail",
            "",
            "| Pocket | State | Class | Matched States | Basis |",
            "|--------|-------|-------|---------------|-------|",
        ])
        for sc in state_classes:
            lines.append(
                f"| {sc['pocket_id']} | {sc['receptor_id']} | {sc['state_class']} | "
                f"{sc['matched_states'] or '—'} | {sc['state_class_basis']} |"
            )
        lines.append("")

    # Pairwise comparisons
    if comparisons:
        lines.extend([
            "## Pairwise Comparisons",
            "",
            "| Pocket A | Pocket B | Dist (Å) | Jaccard | Match |",
            "|----------|----------|----------|---------|-------|",
        ])
        for c in comparisons:
            lines.append(
                f"| {c['pocket_id_a']} | {c['pocket_id_b']} | "
                f"{c['centroid_distance_A']} | {c['residue_jaccard']} | "
                f"{c['match_type']} |"
            )
        lines.append("")

    lines.extend([
        "---",
        "",
        "Generated by `egfr_pipeline.phase2.cross_state_alignment`",
    ])

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_cross_state_alignment(
    output_dir: Path,
) -> Tuple[Path, Path, Path]:
    """Run full TG 2.5 pipeline.

    Returns (comparison_csv, state_class_csv, note_path).
    """
    # Load inputs
    pockets = _load_csv(output_dir / "candidate_pockets.csv")
    if not pockets:
        raise FileNotFoundError("No candidate pockets found. Run TG 2.2 first.")

    # Load druggability summary if available (optional enrichment)
    drug_rows = _load_csv(output_dir / "druggability_proposal_summary.csv")
    drug_map = {r["pocket_id"]: r for r in drug_rows}

    states_present = set(p["receptor_id"] for p in pockets)
    print(f"  Loaded {len(pockets)} pockets from {len(states_present)} state(s)")
    print(f"  States with data: {sorted(states_present)}")
    missing = set(RECEPTOR_STATES) - states_present
    if missing:
        print(f"  States missing pocket data: {sorted(missing)}")

    # Task 2.5.1: Pairwise cross-state comparison
    comparisons = compare_pockets_across_states(pockets)
    n_matches = sum(1 for c in comparisons if c["match_type"] == "match")
    n_shifted = sum(1 for c in comparisons if c["match_type"] == "shifted")
    print(f"  Pairwise comparisons: {len(comparisons)} "
          f"(match={n_matches}, shifted={n_shifted})")

    # Task 2.5.2: State class assignment
    state_classes = assign_state_classes(pockets, comparisons, drug_map)

    by_class: Dict[str, int] = defaultdict(int)
    for sc in state_classes:
        by_class[sc["state_class"]] += 1
    print(f"  State classes: {dict(by_class)}")

    # Write comparison CSV
    comp_path = output_dir / "candidate_pocket_cross_state_comparison.csv"
    _write_csv(comp_path, COMPARISON_COLUMNS, comparisons)

    # Write state class CSV
    class_path = output_dir / "candidate_pocket_state_classes.csv"
    _write_csv(class_path, STATE_CLASS_COLUMNS, state_classes)

    # Write note
    note_lines = build_alignment_note(comparisons, state_classes, states_present)
    note_path = output_dir / "phase2_cross_state_alignment_note.md"
    with open(note_path, "w", encoding="utf-8") as f:
        f.write("\n".join(note_lines))

    return comp_path, class_path, note_path


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


def _centroid_distance(a: dict, b: dict) -> Optional[float]:
    ax = _safe_float(a.get("centroid_x", ""))
    ay = _safe_float(a.get("centroid_y", ""))
    az = _safe_float(a.get("centroid_z", ""))
    bx = _safe_float(b.get("centroid_x", ""))
    by = _safe_float(b.get("centroid_y", ""))
    bz = _safe_float(b.get("centroid_z", ""))
    if any(v is None for v in (ax, ay, az, bx, by, bz)):
        return None
    return math.sqrt((ax - bx)**2 + (ay - by)**2 + (az - bz)**2)


def _parse_residue_nums(nums_str: str) -> Set[int]:
    if not nums_str:
        return set()
    result = set()
    for n in nums_str.split(";"):
        n = n.strip()
        if n:
            try:
                result.add(int(float(n)))
            except ValueError:
                pass
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 TG 2.5: Cross-state pocket proposal alignment"
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE2_OUTPUT_DIR,
        help="Phase 2 output directory",
    )
    args = parser.parse_args()

    print("Phase 2 — Task Group 2.5: Cross-State Pocket Proposal Alignment")
    print(f"  Output: {args.output_dir}")
    print()

    comp_path, class_path, note_path = run_cross_state_alignment(args.output_dir)

    print(f"\n{'='*60}")
    print("TG 2.5 COMPLETE")
    print(f"{'='*60}")
    print(f"  Comparisons:   {comp_path}")
    print(f"  State classes: {class_path}")
    print(f"  Note:          {note_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
