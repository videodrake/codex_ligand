#!/usr/bin/env python3
"""Phase 1 Task 1.5: Multi-state interface patch comparison.

Compares receptor-side interface patches across the three receptor states
to determine which residues/patches are state-robust (consistent across
conformational states) and which are state-specific.

Reads (per state):
  - ppi_cluster_summary.csv
  - ppi_hotspot_residues.csv
  - ppi_interface_patch_table.csv

Produces:
  - ppi_patch_cross_state_comparison.csv  (per-residue across states)
  - ppi_patch_state_robustness.csv        (robustness classification)
  - phase1_interface_comparison_report.md  (readable summary)

Usage:
    python -m egfr_pipeline.phase1.compare_states \
        [--output_dir output/workflow_b/phase1_ppi_analysis]
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from egfr_pipeline import paths

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_PATH_CONFIG = {"output_root": str(PROJECT_ROOT / "output")}
PHASE1_OUTPUT_DIR = paths.wb_phase1_ppi_analysis(DEFAULT_PATH_CONFIG)
RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]
NLOBE_CLOBE_BOUNDARY = 838

# Robustness classification thresholds
ROBUST_MIN_STATES = 3      # Present in all 3 states → robust
MODERATE_MIN_STATES = 2    # Present in 2 states → moderate
# Present in 1 state → state-specific

# Output schemas
CROSS_STATE_COLUMNS = [
    "chain",
    "residue_id",
    "residue_num",
    "residue_name",
    "lobe_label",
    "n_states_present",
    "states_present",               # Semicolon-separated
    "max_occupancy_by_state",       # JSON-like: "3GT8_raw:0.85;EGFR_160-185:0.72"
    "n_clusters_by_state",          # "3GT8_raw:3;EGFR_160-185:2"
    "is_hotspot_in_states",         # "3GT8_raw;EGFR_160-185"
    "robustness_class",             # robust / moderate / state_specific
    "global_max_occupancy",         # Highest occupancy across all states
    "global_mean_occupancy",        # Mean max_occupancy across states where present
]

ROBUSTNESS_COLUMNS = [
    "chain",
    "residue_id",
    "residue_num",
    "residue_name",
    "lobe_label",
    "robustness_class",
    "n_states_present",
    "states_present",
    "global_max_occupancy",
    "is_hotspot_any_state",
    "is_hotspot_all_present_states",
    "construct_type",
    "orientation_validation_status",
]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_patch_table(state_dir: Path) -> List[dict]:
    """Load ppi_interface_patch_table.csv for one state."""
    path = state_dir / "ppi_interface_patch_table.csv"
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_hotspot_residues(state_dir: Path) -> List[dict]:
    """Load ppi_hotspot_residues.csv for one state."""
    path = state_dir / "ppi_hotspot_residues.csv"
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_cluster_summary(state_dir: Path) -> List[dict]:
    """Load ppi_cluster_summary.csv for one state."""
    path = state_dir / "ppi_cluster_summary.csv"
    if not path.exists():
        return []
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------------------
# Cross-state comparison
# ---------------------------------------------------------------------------

def compare_across_states(
    output_base: Path,
    states: List[str] = None,
) -> Tuple[List[dict], List[dict]]:
    """Compare receptor-side patches across receptor states.

    Returns (cross_state_rows, robustness_rows).
    """
    if states is None:
        states = RECEPTOR_STATES

    # --- Collect per-state data ---
    # Key = (chain, residue_id)
    # Value per state = { max_occupancy, n_clusters, is_hotspot, residue_meta }
    residue_state_data = defaultdict(dict)

    warnings = []
    states_with_data = []

    for state in states:
        state_dir = output_base / state
        patches = load_patch_table(state_dir)
        hotspots = load_hotspot_residues(state_dir)

        if not patches:
            warnings.append(f"No patch table for {state}")
            continue

        states_with_data.append(state)

        # Build hotspot set for this state
        hotspot_set = set()
        for h in hotspots:
            if h.get("is_hotspot") in ("True", "true", True):
                hotspot_set.add((h.get("chain", ""), h.get("residue_id", "")))

        for p in patches:
            key = (p.get("chain", ""), p.get("residue_id", ""))
            max_occ = _safe_float(p.get("max_occupancy", "")) or 0.0
            n_clusters = int(p.get("n_clusters_present", 0) or 0)

            residue_state_data[key][state] = {
                "max_occupancy": max_occ,
                "n_clusters": n_clusters,
                "is_hotspot": key in hotspot_set,
                "residue_num": p.get("residue_num", ""),
                "residue_name": p.get("residue_name", ""),
                "lobe_label": p.get("lobe_label", ""),
                "construct_type": p.get("construct_type", ""),
                "orientation_validation_status": p.get(
                    "orientation_validation_status", ""
                ),
            }

    if not states_with_data:
        return [], []

    # --- Build cross-state comparison ---
    cross_state_rows = []
    robustness_rows = []

    for (chain, rid), state_data in sorted(residue_state_data.items()):
        n_states = len(state_data)
        present_states = sorted(state_data.keys())

        # Collect per-state metrics
        occ_parts = []
        cluster_parts = []
        hotspot_states = []
        all_occs = []

        meta = {}
        for st in present_states:
            d = state_data[st]
            occ_parts.append(f"{st}:{d['max_occupancy']:.4f}")
            cluster_parts.append(f"{st}:{d['n_clusters']}")
            if d["is_hotspot"]:
                hotspot_states.append(st)
            all_occs.append(d["max_occupancy"])
            meta = d  # last one for metadata

        global_max = max(all_occs) if all_occs else 0
        global_mean = sum(all_occs) / len(all_occs) if all_occs else 0

        # Robustness classification
        if n_states >= ROBUST_MIN_STATES:
            robustness = "robust"
        elif n_states >= MODERATE_MIN_STATES:
            robustness = "moderate"
        else:
            robustness = "state_specific"

        cross_state_rows.append({
            "chain": chain,
            "residue_id": rid,
            "residue_num": meta.get("residue_num", ""),
            "residue_name": meta.get("residue_name", ""),
            "lobe_label": meta.get("lobe_label", ""),
            "n_states_present": n_states,
            "states_present": ";".join(present_states),
            "max_occupancy_by_state": ";".join(occ_parts),
            "n_clusters_by_state": ";".join(cluster_parts),
            "is_hotspot_in_states": ";".join(hotspot_states),
            "robustness_class": robustness,
            "global_max_occupancy": round(global_max, 4),
            "global_mean_occupancy": round(global_mean, 4),
        })

        is_hotspot_any = len(hotspot_states) > 0
        is_hotspot_all = len(hotspot_states) == n_states
        construct_type = _collapse_metadata_value(
            [state_data[st].get("construct_type", "") for st in present_states],
            mixed_label="mixed_construct_type",
        )
        orientation_status = _collapse_metadata_value(
            [
                state_data[st].get("orientation_validation_status", "")
                for st in present_states
            ],
            mixed_label="mixed_orientation_status",
        )

        robustness_rows.append({
            "chain": chain,
            "residue_id": rid,
            "residue_num": meta.get("residue_num", ""),
            "residue_name": meta.get("residue_name", ""),
            "lobe_label": meta.get("lobe_label", ""),
            "robustness_class": robustness,
            "n_states_present": n_states,
            "states_present": ";".join(present_states),
            "global_max_occupancy": round(global_max, 4),
            "is_hotspot_any_state": is_hotspot_any,
            "is_hotspot_all_present_states": is_hotspot_all,
            "construct_type": construct_type,
            "orientation_validation_status": orientation_status,
        })

    # Sort: receptor first, then by robustness (robust > moderate > specific),
    # then by global_max_occupancy descending
    rank = {"robust": 0, "moderate": 1, "state_specific": 2}
    cross_state_rows.sort(key=lambda r: (
        0 if r["chain"] == "A" else 1,
        rank.get(r["robustness_class"], 3),
        -r.get("global_max_occupancy", 0),
    ))
    robustness_rows.sort(key=lambda r: (
        0 if r["chain"] == "A" else 1,
        rank.get(r["robustness_class"], 3),
        -r.get("global_max_occupancy", 0),
    ))

    return cross_state_rows, robustness_rows


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_comparison_report(
    output_base: Path,
    cross_state_rows: List[dict],
    robustness_rows: List[dict],
    states: List[str],
) -> Path:
    """Generate phase1_interface_comparison_report.md."""
    report_path = output_base / "phase1_interface_comparison_report.md"
    construct_type = _collapse_metadata_value(
        [row.get("construct_type", "") for row in robustness_rows],
        mixed_label="mixed_construct_type",
    )
    orientation_status = _collapse_metadata_value(
        [row.get("orientation_validation_status", "") for row in robustness_rows],
        mixed_label="mixed_orientation_status",
    )
    evidence_line = _describe_orientation_evidence(orientation_status)

    lines = [
        "# Phase 1 Interface Comparison Report",
        "",
        "## Multi-State Receptor-Side Patch Comparison",
        "",
        f"States compared: {', '.join(states)}",
        f"Construct type: {construct_type or 'unknown'}",
        f"Evidence: {evidence_line}",
        "",
    ]

    # Per-state cluster summaries
    lines.append("## Per-State Summary")
    lines.append("")
    lines.append("| State | Clusters | Orient-valid models | Hotspot residues (receptor) |")
    lines.append("|-------|----------|--------------------|-----------------------------|")

    for state in states:
        state_dir = output_base / state
        summaries = load_cluster_summary(state_dir)
        if not summaries:
            lines.append(f"| {state} | — | — | — |")
            continue
        n_clusters = len(summaries)
        n_valid = sum(int(s.get("n_members_orientation_valid", 0) or 0) for s in summaries)
        # Count unique hotspot residues
        all_hs = set()
        for s in summaries:
            hs = s.get("receptor_hotspot_residues", "")
            if hs:
                all_hs.update(hs.split(";"))
        lines.append(f"| {state} | {n_clusters} | {n_valid} | {len(all_hs)} |")

    lines.append("")

    # Robustness summary
    receptor_robust = [r for r in robustness_rows
                       if r["chain"] == "A" and r["robustness_class"] == "robust"]
    receptor_moderate = [r for r in robustness_rows
                         if r["chain"] == "A" and r["robustness_class"] == "moderate"]
    receptor_specific = [r for r in robustness_rows
                         if r["chain"] == "A" and r["robustness_class"] == "state_specific"]

    lines.extend([
        "## Robustness Classification (Receptor-side)",
        "",
        f"- **Robust** (all {len(states)} states): {len(receptor_robust)} residues",
        f"- **Moderate** (2 states): {len(receptor_moderate)} residues",
        f"- **State-specific** (1 state): {len(receptor_specific)} residues",
        "",
    ])

    # Robust residues table
    if receptor_robust:
        lines.extend([
            "### Robust Residues (present in all states)",
            "",
            "| Residue | Lobe | Global max occupancy | Hotspot in all? |",
            "|---------|------|---------------------|-----------------|",
        ])
        for r in receptor_robust:
            lines.append(
                f"| {r['residue_id']} | {r['lobe_label']} | "
                f"{r['global_max_occupancy']} | "
                f"{'Yes' if r['is_hotspot_all_present_states'] else 'No'} |"
            )
        lines.append("")

    # Lobe distribution across robustness classes
    n_robust_nlobe = sum(1 for r in receptor_robust if r["lobe_label"] == "N-lobe")
    n_robust_clobe = sum(1 for r in receptor_robust if r["lobe_label"] == "C-lobe")
    n_moderate_nlobe = sum(1 for r in receptor_moderate if r["lobe_label"] == "N-lobe")
    n_moderate_clobe = sum(1 for r in receptor_moderate if r["lobe_label"] == "C-lobe")

    lines.extend([
        "### Lobe Distribution by Robustness",
        "",
        "| Class | N-lobe | C-lobe | Total |",
        "|-------|--------|--------|-------|",
        f"| Robust | {n_robust_nlobe} | {n_robust_clobe} | {len(receptor_robust)} |",
        f"| Moderate | {n_moderate_nlobe} | {n_moderate_clobe} | {len(receptor_moderate)} |",
        f"| State-specific | — | — | {len(receptor_specific)} |",
        "",
    ])

    # Partner-side summary
    partner_rows = [r for r in robustness_rows if r["chain"] == "B"]
    if partner_rows:
        partner_robust = [r for r in partner_rows if r["robustness_class"] == "robust"]
        lines.extend([
            "## Partner-side (beta-meander) Robustness",
            "",
            f"Total partner residues observed: {len(partner_rows)}",
            f"Robust (all states): {len(partner_robust)}",
            "",
        ])
        if partner_robust:
            lines.append("Robust partner residues: " +
                         ", ".join(r["residue_id"] for r in partner_robust))
            lines.append("")

    # Numbering consistency check
    lines.extend([
        "## Cross-State Numbering Consistency",
        "",
        "Cross-state comparison assumes the receptor inputs were normalized during",
        "TG 1.0 input preparation to a shared 3GT8-derived numbering scheme and",
        "a consistent receptor chain convention.",
        "",
        "Use the Phase 1 input validation outputs as the first authority on",
        "numbering safety before interpreting residue differences across states.",
        "If upstream validation reports a numbering or chain mismatch, treat",
        "cross-state residue differences as potentially ambiguous until resolved.",
        "",
    ])

    # Interpretation guidance
    lines.extend([
        "## Interpretation Notes",
        "",
        "- **Robust residues** are the strongest candidates for the MYO1D",
        "  attachment patch — they appear regardless of EGFR conformational state.",
        "- **State-specific residues** may indicate conformationally gated",
        "  interaction sites that open only in certain MD-sampled states.",
        "- **N-lobe residues** (< 838) are newly detectable in Phase 1 due to",
        "  the full kinase domain receptor. Their presence/absence is informative.",
        f"- Orientation evidence status for this comparison: **{orientation_status or 'unknown'}**.",
        "",
    ])

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    return report_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(val: str) -> Optional[float]:
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def _collapse_metadata_value(values: List[str], mixed_label: str) -> str:
    """Collapse repeated metadata values without hiding disagreements."""
    normalized = [value.strip() for value in values if value and value.strip()]
    if not normalized:
        return ""
    unique_values = sorted(set(normalized))
    if len(unique_values) == 1:
        return unique_values[0]
    return mixed_label


def _describe_orientation_evidence(status: str) -> str:
    """Render orientation evidence text for the comparison report."""
    if status == "orientation_validated":
        return "orientation-validated PyRosetta models only"
    if status == "not_available":
        return "orientation validation not available for the contributing states"
    if status == "mixed_orientation_status":
        return (
            "mixed orientation validation status across contributing states "
            "(see ppi_patch_state_robustness.csv)"
        )
    if status:
        return status
    return "orientation status unknown"


def _write_csv(path: Path, rows: List[dict], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 1 TG 1.5: Multi-state interface patch comparison"
    )
    parser.add_argument("--output_dir", type=Path,
                        default=PHASE1_OUTPUT_DIR,
                        help="Phase 1 output base directory")
    args = parser.parse_args()

    print("Phase 1 — Task 1.5: Multi-State Interface Patch Comparison")
    print(f"Output base: {args.output_dir}")
    print(f"States: {', '.join(RECEPTOR_STATES)}")

    cross_rows, robust_rows = compare_across_states(args.output_dir, RECEPTOR_STATES)

    if not cross_rows:
        print("\nNo cross-state comparison data. Run TG 1.3 for ≥1 state first.")
        return 0

    # Write CSVs
    cross_path = args.output_dir / "ppi_patch_cross_state_comparison.csv"
    robust_path = args.output_dir / "ppi_patch_state_robustness.csv"

    _write_csv(cross_path, cross_rows, CROSS_STATE_COLUMNS)
    _write_csv(robust_path, robust_rows, ROBUSTNESS_COLUMNS)

    # Generate report
    report_path = generate_comparison_report(
        args.output_dir, cross_rows, robust_rows, RECEPTOR_STATES
    )

    # Summary
    receptor_rows = [r for r in robust_rows if r["chain"] == "A"]
    n_robust = sum(1 for r in receptor_rows if r["robustness_class"] == "robust")
    n_moderate = sum(1 for r in receptor_rows if r["robustness_class"] == "moderate")
    n_specific = sum(1 for r in receptor_rows if r["robustness_class"] == "state_specific")

    print(f"\n{'='*60}")
    print("MULTI-STATE COMPARISON COMPLETE")
    print(f"{'='*60}")
    print(f"  Total receptor residues compared: {len(receptor_rows)}")
    print(f"  Robust (all states):    {n_robust}")
    print(f"  Moderate (2 states):    {n_moderate}")
    print(f"  State-specific (1):     {n_specific}")
    print(f"\n  Outputs:")
    print(f"    {cross_path.name}")
    print(f"    {robust_path.name}")
    print(f"    {report_path.name}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
