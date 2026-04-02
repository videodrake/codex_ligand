#!/usr/bin/env python3
"""Phase 1 Task 1.6: Phase 1 Review Report.

Generates a readable Phase 1 review package that summarizes receptor-side
MYO1D attachment evidence before moving into Phase 2 pocket proposal.

Reads all Phase 1 outputs and produces:
  - phase1_interface_report.md           (comprehensive review report)
  - phase1_downstream_patch_reference.csv (Phase 2 handoff file)

The patch reference CSV includes construct_type and orientation_validation_status
fields required by Phase 2 Task Group 2.0 (Patch Reference Ingestion).

Usage:
    python -m egfr_pipeline.phase1.review_report \
        [--output_dir output/workflow_b/phase1_ppi_analysis]
"""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

from egfr_pipeline import paths
from egfr_pipeline.csv_utils import load_csv
from egfr_pipeline.parsing_utils import safe_float

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_PATH_CONFIG = {"output_root": str(PROJECT_ROOT / "output")}
PHASE1_OUTPUT_DIR = paths.wb_phase1_ppi_analysis(DEFAULT_PATH_CONFIG)
RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]
NLOBE_CLOBE_BOUNDARY = 838

# Patch reference schema (Phase 2 handoff)
PATCH_REFERENCE_COLUMNS = [
    "residue_id",
    "residue_num",
    "residue_name",
    "chain",
    "lobe_label",
    "evidence_source",              # primary / secondary / auxiliary
    "construct_type",               # full_kinase_domain
    "orientation_validation_status",# orientation_validated
    "robustness_class",             # robust / moderate / state_specific
    "n_states_present",
    "states_present",
    "global_max_occupancy",
    "is_hotspot_any_state",
    "pyrosetta_evidence",           # True/False
    "lightdock_evidence",           # True/False
    "method_agreement",             # both / pyrosetta_only / lightdock_only / none
    "confidence",                   # high / medium / low
    "notes",
]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_review_report(output_base: Path) -> tuple:
    """Generate Phase 1 review report and downstream patch reference.

    Returns (report_path, patch_ref_path).
    """
    # --- Collect all data ---
    # Per-state cluster summaries
    state_summaries = {}
    for state in RECEPTOR_STATES:
        sd = output_base / state
        state_summaries[state] = {
            "clusters": load_csv(sd / "ppi_cluster_summary.csv"),
            "hotspots": load_csv(sd / "ppi_hotspot_residues.csv"),
            "patches": load_csv(sd / "ppi_interface_patch_table.csv"),
            "orientation_log": load_csv(sd / "orientation_filter_log.csv"),
            "convergence": load_csv(sd / "cross_method_convergence.csv"),
        }

    # Cross-state comparison
    robustness = load_csv(output_base / "ppi_patch_state_robustness.csv")
    cross_state = load_csv(output_base / "ppi_patch_cross_state_comparison.csv")

    # Input metadata
    metadata_dir = _phase1_metadata_dir()
    receptor_meta = load_csv(metadata_dir / "receptor_metadata.csv")
    partner_meta = load_csv(metadata_dir / "partner_metadata.csv")

    # --- Build report ---
    lines = _build_report(state_summaries, robustness, cross_state,
                          receptor_meta, partner_meta)
    report_path = output_base / "phase1_interface_report.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    # --- Build patch reference ---
    patch_rows = _build_patch_reference(robustness, state_summaries)
    patch_path = output_base / "phase1_downstream_patch_reference.csv"
    with open(patch_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=PATCH_REFERENCE_COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(patch_rows)

    return report_path, patch_path


def _phase1_metadata_dir() -> Path:
    runtime_dir = paths.wa_phase2_runtime_inputs(DEFAULT_PATH_CONFIG)
    if runtime_dir.exists():
        return runtime_dir
    legacy_runtime_dir = PROJECT_ROOT / "output" / "phase1_ppi" / "runtime_inputs"
    if legacy_runtime_dir.exists():
        return legacy_runtime_dir
    return PROJECT_ROOT / "input" / "PPI" / "phase1"


def _build_report(
    state_summaries: dict,
    robustness: List[dict],
    cross_state: List[dict],
    receptor_meta: List[dict],
    partner_meta: List[dict],
) -> List[str]:
    """Build the full Phase 1 review report."""
    construct_type = _collapse_metadata_value(
        [row.get("construct_type", "") for row in robustness],
        mixed_label="mixed_construct_type",
    )
    orientation_status = _collapse_metadata_value(
        [row.get("orientation_validation_status", "") for row in robustness],
        mixed_label="mixed_orientation_status",
    )

    lines = [
        "# Phase 1 Interface Report",
        "",
        "## PPI-First Interface Mapping: Full Kinase Domain × Extended Beta-Meander",
        "",
        "---",
        "",
    ]

    # Section 1: Input Summary
    lines.extend(_section_inputs(receptor_meta, partner_meta))

    # Section 2: Per-State Evidence
    lines.extend(_section_per_state(state_summaries))

    # Section 3: Orientation Filtering
    lines.extend(_section_orientation(state_summaries))

    # Section 4: Receptor Hotspot Residues
    lines.extend(_section_hotspots(state_summaries))

    # Section 5: Cross-State Robustness
    lines.extend(_section_robustness(robustness))

    # Section 6: LightDock Convergence
    lines.extend(_section_lightdock(state_summaries))

    # Section 7: Evidence Hierarchy
    lines.extend(_section_evidence_hierarchy(orientation_status))

    # Section 8: Downstream Handoff
    lines.extend(_section_handoff(construct_type, orientation_status))

    lines.extend([
        "---",
        "",
        "Generated by `egfr_pipeline.phase1.review_report`",
    ])

    return lines


def _section_inputs(receptor_meta, partner_meta) -> List[str]:
    lines = [
        "## 1. Input Summary",
        "",
    ]
    if receptor_meta:
        lines.extend([
            "### Receptor States",
            "",
            "| State | Residue range | N residues | Source |",
            "|-------|--------------|------------|--------|",
        ])
        for r in receptor_meta:
            lines.append(
                f"| {r.get('state_name', '')} | "
                f"{r.get('residue_start', '')}-{r.get('residue_end', '')} | "
                f"{r.get('n_residues', '')} | "
                f"{r.get('source_pdb', '')} |"
            )
        lines.append("")

    if partner_meta:
        lines.extend([
            "### Partner",
            "",
        ])
        for p in partner_meta:
            lines.append(
                f"- **{p.get('partner_id', '')}**: residues "
                f"{p.get('residue_start', '')}-{p.get('residue_end', '')}, "
                f"{p.get('n_residues', '')} residues, "
                f"active face: {p.get('active_face_residues', 'N/A')}"
            )
        lines.append("")

    lines.extend([
        "### Key Differences from Legacy",
        "",
        "| Parameter | Legacy (pilot) | Phase 1 |",
        "|-----------|---------------|---------|",
        "| Receptor | C-lobe fragment (~45 res) | Full kinase domain (~309 res) |",
        "| Partner | 962–1006 (truncated) | 955–1006 (extended) |",
        "| Construct | EGFR dimer | EGFR monomer |",
        "| Orientation filter | None | Mandatory |",
        "| N-lobe | Not detectable | Fully covered |",
        "",
    ])

    return lines


def _section_per_state(state_summaries) -> List[str]:
    lines = [
        "## 2. Per-State PyRosetta Evidence",
        "",
        "| State | Clusters | Orient-valid | Receptor hotspots | N-lobe hotspots | C-lobe hotspots |",
        "|-------|----------|-------------|-------------------|-----------------|-----------------|",
    ]

    for state in RECEPTOR_STATES:
        data = state_summaries.get(state, {})
        clusters = data.get("clusters", [])
        if not clusters:
            lines.append(f"| {state} | — | — | — | — | — |")
            continue

        n_clusters = len(clusters)
        n_valid = sum(int(c.get("n_members_orientation_valid", 0) or 0) for c in clusters)
        n_r_hotspots = sum(
            len(c.get("receptor_hotspot_residues", "").split(";"))
            for c in clusters if c.get("receptor_hotspot_residues")
        )
        n_nlobe = sum(int(c.get("n_nlobe_hotspots", 0) or 0) for c in clusters)
        n_clobe = sum(int(c.get("n_clobe_hotspots", 0) or 0) for c in clusters)

        lines.append(f"| {state} | {n_clusters} | {n_valid} | {n_r_hotspots} | {n_nlobe} | {n_clobe} |")

    lines.append("")
    return lines


def _section_orientation(state_summaries) -> List[str]:
    lines = [
        "## 3. Orientation Filtering Summary",
        "",
        "All models were evaluated with the dual-vector PCA orientation test.",
        "Only PASS models (active face toward receptor) contribute to consensus.",
        "",
        "| State | Total models | PASS | FAIL | AMBIGUOUS | Pass rate |",
        "|-------|-------------|------|------|-----------|-----------|",
    ]

    for state in RECEPTOR_STATES:
        data = state_summaries.get(state, {})
        orient = data.get("orientation_log", [])
        if not orient:
            lines.append(f"| {state} | — | — | — | — | — |")
            continue

        n_total = len(orient)
        n_pass = sum(1 for o in orient if o.get("orientation_class") == "pass")
        n_fail = sum(1 for o in orient if o.get("orientation_class") == "fail")
        n_ambig = sum(1 for o in orient if o.get("orientation_class") == "ambiguous")
        rate = f"{100 * n_pass / n_total:.1f}%" if n_total > 0 else "—"

        lines.append(f"| {state} | {n_total} | {n_pass} | {n_fail} | {n_ambig} | {rate} |")

    lines.append("")
    return lines


def _section_hotspots(state_summaries) -> List[str]:
    lines = [
        "## 4. Receptor-Side Hotspot Residues",
        "",
    ]

    for state in RECEPTOR_STATES:
        data = state_summaries.get(state, {})
        hotspots = data.get("hotspots", [])
        receptor_hs = [h for h in hotspots
                       if h.get("chain") == "A"
                       and h.get("is_hotspot") in ("True", "true", True)]
        if not receptor_hs:
            lines.append(f"### {state}: No hotspot data")
            lines.append("")
            continue

        lines.extend([
            f"### {state}",
            "",
            "| Cluster | Residue | Lobe | Occupancy | Mean ΔE |",
            "|---------|---------|------|-----------|---------|",
        ])

        for h in sorted(receptor_hs, key=lambda x: (
            x.get("cluster_id", ""), -(safe_float(x.get("occupancy", "")) or 0)
        )):
            lines.append(
                f"| {h.get('cluster_id', '')} | {h.get('residue_id', '')} | "
                f"{h.get('lobe_label', '')} | {h.get('occupancy', '')} | "
                f"{h.get('mean_delta_e_total', '')} |"
            )
        lines.append("")

    return lines


def _section_robustness(robustness) -> List[str]:
    lines = [
        "## 5. Cross-State Robustness",
        "",
    ]

    receptor_rows = [r for r in robustness if r.get("chain") == "A"]
    if not receptor_rows:
        lines.append("No cross-state comparison data available.")
        lines.append("")
        return lines

    robust = [r for r in receptor_rows if r.get("robustness_class") == "robust"]
    moderate = [r for r in receptor_rows if r.get("robustness_class") == "moderate"]
    specific = [r for r in receptor_rows if r.get("robustness_class") == "state_specific"]

    lines.extend([
        f"- **Robust** (all 3 states): {len(robust)} residues",
        f"- **Moderate** (2 states): {len(moderate)} residues",
        f"- **State-specific** (1 state): {len(specific)} residues",
        "",
    ])

    if robust:
        lines.extend([
            "### Robust Receptor Residues",
            "",
            "| Residue | Lobe | Max occupancy | Hotspot in all states? |",
            "|---------|------|--------------|------------------------|",
        ])
        for r in robust:
            lines.append(
                f"| {r.get('residue_id', '')} | {r.get('lobe_label', '')} | "
                f"{r.get('global_max_occupancy', '')} | "
                f"{'Yes' if r.get('is_hotspot_all_present_states') in ('True', 'true', True) else 'No'} |"
            )
        lines.append("")

    if moderate:
        lines.extend([
            "### Moderate Receptor Residues",
            "",
            "| Residue | Lobe | States | Max occupancy |",
            "|---------|------|--------|--------------|",
        ])
        for r in moderate:
            lines.append(
                f"| {r.get('residue_id', '')} | {r.get('lobe_label', '')} | "
                f"{r.get('states_present', '')} | {r.get('global_max_occupancy', '')} |"
            )
        lines.append("")

    return lines


def _section_lightdock(state_summaries) -> List[str]:
    lines = [
        "## 6. LightDock Secondary Validation",
        "",
        "LightDock remains secondary evidence only.",
        "",
    ]

    # Check if orientation data is available
    has_orientation = False
    for state in RECEPTOR_STATES:
        data = state_summaries.get(state, {})
        conv = data.get("convergence", [])
        for c in conv:
            ovs = c.get("orientation_validation_status", "")
            if ovs and ovs not in ("not_available", ""):
                has_orientation = True
                break

    if has_orientation:
        lines.append(
            "LightDock poses are now filtered with the PDB-based orientation test "
            "(same algorithm as PyRosetta, PyRosetta-free implementation)."
        )
    else:
        lines.append(
            "LightDock raw support tables use `orientation_validation_status = not_available` "
            "until LightDock execution and orientation filtering are performed."
        )
    lines.append("")

    any_convergence = False
    for state in RECEPTOR_STATES:
        data = state_summaries.get(state, {})
        conv = data.get("convergence", [])
        if not conv:
            continue

        any_convergence = True
        receptor_conv = [c for c in conv if c.get("chain") == "A"]
        n_convergent = sum(1 for c in receptor_conv
                          if c.get("convergence_class") == "convergent")
        n_pyro = sum(1 for c in receptor_conv
                     if c.get("convergence_class") == "pyrosetta_only")
        n_light = sum(1 for c in receptor_conv
                      if c.get("convergence_class") == "lightdock_only")
        n_total = len(receptor_conv)
        jaccard = n_convergent / n_total if n_total > 0 else 0

        lines.extend([
            f"### {state}",
            "",
            f"- Convergent: {n_convergent} / {n_total} receptor residues",
            f"- PyRosetta-only: {n_pyro}",
            f"- LightDock-only: {n_light}",
            f"- Jaccard overlap: {jaccard:.3f}",
            "",
        ])

        # List convergent residues
        convergent_res = [c.get("residue_id", "") for c in receptor_conv
                          if c.get("convergence_class") == "convergent"]
        if convergent_res:
            lines.append(f"Convergent residues: {', '.join(convergent_res)}")
            lines.append("")

    if not any_convergence:
        lines.append("No LightDock convergence data available yet.")
        lines.append("Run LightDock after server-side execution, then re-generate this report.")
        lines.append("")

    return lines


def _section_evidence_hierarchy(orientation_status: str) -> List[str]:
    return [
        "## 7. Evidence Hierarchy",
        "",
        "| Level | Source | Role | Status |",
        "|-------|--------|------|--------|",
        "| Primary | PyRosetta PPI docking | Receptor-side patch definition | Core evidence |",
        "| Secondary | LightDock validation | Method-independence check | Independent validation |",
        "| Auxiliary | AlphaFold-Multimer | Structural context (if available) | Optional supplement |",
        "",
        (
            "Primary evidence (PyRosetta; orientation status: "
            f"{_describe_orientation_status(orientation_status)}) is the sole basis for"
        ),
        "receptor-side patch definitions handed to Phase 2. Secondary and auxiliary",
        "evidence increase confidence but do not replace primary evidence.",
        "",
    ]


def _section_handoff(construct_type: str, orientation_status: str) -> List[str]:
    return [
        "## 8. Downstream Handoff to Phase 2",
        "",
        "The patch reference file `phase1_downstream_patch_reference.csv` contains:",
        "",
        "- Receptor-side residues summarized from the current cross-state comparison output",
        "- Robustness classification (robust / moderate / state_specific)",
        f"- `construct_type` preserved from upstream consensus: `{construct_type or 'unknown'}`",
        (
            "- `orientation_validation_status` preserved from upstream consensus: "
            f"`{orientation_status or 'unknown'}`"
        ),
        "- Method agreement status (PyRosetta / LightDock / both)",
        "- Confidence classification (high / medium / low)",
        "",
        "**Phase 2 compatibility note:** The `construct_type` and",
        "`orientation_validation_status` fields are new in v2. Phase 2 Task Group 2.0",
        "(Patch Reference Ingestion) must be updated to consume these fields.",
        "",
    ]


# ---------------------------------------------------------------------------
# Patch reference for Phase 2
# ---------------------------------------------------------------------------

def _build_patch_reference(
    robustness: List[dict],
    state_summaries: dict,
) -> List[dict]:
    """Build the downstream patch reference CSV for Phase 2."""

    # Build convergence index (any state)
    convergence_by_res = {}  # (chain, residue_id) → convergence_class
    for state in RECEPTOR_STATES:
        data = state_summaries.get(state, {})
        for c in data.get("convergence", []):
            key = (c.get("chain", ""), c.get("residue_id", ""))
            cc = c.get("convergence_class", "")
            # convergent > single method
            if cc == "convergent":
                convergence_by_res[key] = "both"
            elif key not in convergence_by_res:
                convergence_by_res[key] = c.get("method_agreement", "single")

    rows = []
    for r in robustness:
        chain = r.get("chain", "")
        rid = r.get("residue_id", "")
        key = (chain, rid)

        # Only receptor-side for Phase 2 handoff
        if chain != "A":
            continue

        rob_class = r.get("robustness_class", "")
        n_states = int(r.get("n_states_present", 0) or 0)
        max_occ = safe_float(r.get("global_max_occupancy", "")) or 0
        is_hotspot = r.get("is_hotspot_any_state") in ("True", "true", True)

        # Method agreement
        agreement = convergence_by_res.get(key, "pyrosetta_only")
        pyro_evidence = True
        ld_evidence = agreement == "both"

        # Confidence classification
        if rob_class == "robust" and agreement == "both":
            confidence = "high"
        elif rob_class == "robust" or (rob_class == "moderate" and agreement == "both"):
            confidence = "medium"
        else:
            confidence = "low"

        # Notes
        notes_parts = []
        if rob_class == "robust" and is_hotspot:
            notes_parts.append("robust hotspot")
        if r.get("lobe_label") == "N-lobe":
            notes_parts.append("N-lobe (newly detectable in Phase 1)")

        rows.append({
            "residue_id": rid,
            "residue_num": r.get("residue_num", ""),
            "residue_name": r.get("residue_name", ""),
            "chain": chain,
            "lobe_label": r.get("lobe_label", ""),
            "evidence_source": "primary",
            "construct_type": r.get("construct_type", ""),
            "orientation_validation_status": r.get("orientation_validation_status", ""),
            "robustness_class": rob_class,
            "n_states_present": n_states,
            "states_present": r.get("states_present", ""),
            "global_max_occupancy": r.get("global_max_occupancy", ""),
            "is_hotspot_any_state": is_hotspot,
            "pyrosetta_evidence": pyro_evidence,
            "lightdock_evidence": ld_evidence,
            "method_agreement": agreement,
            "confidence": confidence,
            "notes": "; ".join(notes_parts) if notes_parts else "",
        })

    # Sort: high confidence first, then by occupancy
    conf_rank = {"high": 0, "medium": 1, "low": 2}
    rows.sort(key=lambda r: (
        conf_rank.get(r.get("confidence", "low"), 3),
        -(safe_float(r.get("global_max_occupancy", "")) or 0),
    ))

    return rows


def _collapse_metadata_value(values: List[str], mixed_label: str) -> str:
    """Collapse repeated metadata values without hiding disagreements."""
    normalized = [value.strip() for value in values if value and value.strip()]
    if not normalized:
        return ""
    unique_values = sorted(set(normalized))
    if len(unique_values) == 1:
        return unique_values[0]
    return mixed_label


def _describe_orientation_status(status: str) -> str:
    """Render a readable orientation-status label for the report."""
    if status == "orientation_validated":
        return "orientation validated"
    if status == "not_available":
        return "not available"
    if status == "mixed_orientation_status":
        return "mixed across contributing states"
    return status or "unknown"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 1 TG 1.6: Phase 1 review report and Phase 2 handoff"
    )
    parser.add_argument("--output_dir", type=Path,
                        default=PHASE1_OUTPUT_DIR,
                        help="Workflow B Phase 1 analysis output directory")
    args = parser.parse_args()

    print("Phase 1 — Task 1.6: Review Report & Phase 2 Handoff")

    report_path, patch_path = generate_review_report(args.output_dir)

    print(f"\n{'='*60}")
    print("PHASE 1 REVIEW REPORT COMPLETE")
    print(f"{'='*60}")
    print(f"  Report:          {report_path}")
    print(f"  Patch reference: {patch_path}")

    # Count patch reference rows
    rows = load_csv(patch_path)
    n_high = sum(1 for r in rows if r.get("confidence") == "high")
    n_medium = sum(1 for r in rows if r.get("confidence") == "medium")
    n_low = sum(1 for r in rows if r.get("confidence") == "low")
    print(f"\n  Patch reference: {len(rows)} receptor residues")
    print(f"    High confidence:   {n_high}")
    print(f"    Medium confidence: {n_medium}")
    print(f"    Low confidence:    {n_low}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
