#!/usr/bin/env python3
"""Phase 2 Task Group 2.0: Phase 1 Patch Reference Ingestion and Validation.

Reads the Phase 1 downstream patch reference CSV, validates schema and content,
standardizes into a Phase 2 internal format, and produces a validation report.

Tasks:
  2.0.1 - Phase 1 patch reference ingestion
  2.0.2 - Patch-reference validation
  2.0.3 - Patch reference standardization

Inputs:
  - output/phase1_ppi/phase1_downstream_patch_reference.csv

Outputs:
  - output/phase2_pockets/phase2_patch_reference_normalized.csv
  - output/phase2_pockets/phase2_patch_reference_validation.md

Usage:
    python -m egfr_pipeline.phase2.patch_ingestion [--phase1_dir ...] [--output_dir ...]
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE1_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase1_ppi"
PHASE2_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase2_pockets"

RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]

# Expected columns from Phase 1 patch reference (review_report.py PATCH_REFERENCE_COLUMNS)
PHASE1_REQUIRED_COLUMNS = [
    "residue_id",
    "residue_num",
    "residue_name",
    "chain",
    "lobe_label",
    "evidence_source",
    "construct_type",
    "orientation_validation_status",
    "robustness_class",
    "n_states_present",
    "states_present",
    "global_max_occupancy",
    "is_hotspot_any_state",
    "pyrosetta_evidence",
    "lightdock_evidence",
    "method_agreement",
    "confidence",
]

# Phase 2 normalized schema — superset with explicit types
PHASE2_NORMALIZED_COLUMNS = [
    "patch_residue_id",        # e.g., ASP855
    "residue_num",             # int, PDB numbering
    "residue_name",            # 3-letter code
    "chain",                   # A (receptor)
    "lobe",                    # N-lobe / C-lobe
    "construct_type",          # full_kinase_domain
    "orientation_validated",   # True/False
    "robustness",              # robust / moderate / state_specific
    "n_states",                # int
    "states",                  # semicolon-delimited
    "max_occupancy",           # float [0,1]
    "is_hotspot",              # True/False
    "pyrosetta_support",       # True/False
    "lightdock_support",       # True/False
    "method_agreement",        # both / pyrosetta_only / lightdock_only / single
    "confidence",              # high / medium / low
    "phase1_evidence_source",  # primary / secondary / auxiliary
    "notes",
]

# Validation thresholds
VALID_CHAINS = {"A"}
VALID_LOBES = {"N-lobe", "C-lobe"}
VALID_ROBUSTNESS = {"robust", "moderate", "state_specific"}
VALID_CONFIDENCE = {"high", "medium", "low"}
VALID_CONSTRUCT_TYPES = {"full_kinase_domain"}
VALID_METHOD_AGREEMENT = {"both", "pyrosetta_only", "lightdock_only", "single", "none"}
RESIDUE_NUM_RANGE = (699, 1007)  # Full kinase domain PDB numbering
MD_GATE_ARTIFACTS = (
    "md_gate_decision.md",
    "phase1_md_validation_report.md",
    "md_stability_metrics.csv",
)


# ---------------------------------------------------------------------------
# Ingestion (Task 2.0.1)
# ---------------------------------------------------------------------------

def load_patch_reference(phase1_dir: Path) -> Tuple[List[dict], Path]:
    """Load Phase 1 downstream patch reference CSV.

    Returns (rows, file_path).
    Raises FileNotFoundError if the file is missing.
    """
    ref_path = phase1_dir / "phase1_downstream_patch_reference.csv"
    if not ref_path.exists():
        raise FileNotFoundError(
            f"Phase 1 patch reference not found: {ref_path}\n"
            "Run Phase 1 TG 1.6 (review_report) first."
        )

    with open(ref_path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    return rows, ref_path


def inspect_md_gate_artifacts(phase1_dir: Path) -> Dict[str, object]:
    """Collect Phase 1 MD gate artifact presence and decision status."""
    artifacts = [_resolve_md_gate_artifact(phase1_dir, name) for name in MD_GATE_ARTIFACTS]
    decision_artifact = next(
        (artifact for artifact in artifacts if artifact["name"] == "md_gate_decision.md"),
        None,
    )
    decision_status = "missing"
    decision_excerpt = ""
    if decision_artifact and decision_artifact["status"] == "found":
        decision_text = decision_artifact["path"].read_text(encoding="utf-8", errors="replace")
        decision_status = _classify_md_gate_decision(decision_text)
        decision_excerpt = _extract_md_gate_excerpt(decision_text)

    return {
        "artifacts": artifacts,
        "decision_status": decision_status,
        "decision_excerpt": decision_excerpt,
    }


def _resolve_md_gate_artifact(phase1_dir: Path, artifact_name: str) -> Dict[str, object]:
    exact_match = phase1_dir / artifact_name
    if exact_match.is_file():
        return {
            "name": artifact_name,
            "status": "found",
            "location": str(exact_match.relative_to(phase1_dir)),
            "path": exact_match,
        }

    recursive_matches = sorted(path for path in phase1_dir.rglob(artifact_name) if path.is_file())
    if len(recursive_matches) == 1:
        match = recursive_matches[0]
        return {
            "name": artifact_name,
            "status": "found",
            "location": str(match.relative_to(phase1_dir)),
            "path": match,
        }

    if recursive_matches:
        return {
            "name": artifact_name,
            "status": "ambiguous",
            "location": f"{len(recursive_matches)} recursive matches",
            "path": None,
        }

    return {
        "name": artifact_name,
        "status": "missing",
        "location": "not found",
        "path": None,
    }


def _classify_md_gate_decision(text: str) -> str:
    normalized = " ".join(text.lower().split())
    if not normalized:
        return "unknown"

    blocked_patterns = (
        r"\bphase 2\b.{0,40}\bblocked\b",
        r"\bblocked\b.{0,40}\bphase 2\b",
        r"\bdo not use this patch\b",
        r"\breturn to phase 1\b",
        r"\bunstable interface\b",
        r"\b(?:overall|final|gate|decision|classification|verdict)\b.{0,40}\bunstable\b",
        r"\bunstable\b.{0,40}\b(?:overall|final|gate|decision|classification|verdict)\b",
    )
    caution_patterns = (
        r"\bmetastable\b",
        r"\bproceed with caution\b",
        r"\breduced confidence\b",
    )
    pass_patterns = (
        r"\bphase 2\b.{0,40}\bcan proceed\b",
        r"\bproceed to phase 2\b",
        r"\bphase 2\b.{0,40}\bproceed\b",
        r"\bstable\b",
    )

    if any(re.search(pattern, normalized) for pattern in blocked_patterns):
        return "blocked"
    if any(re.search(pattern, normalized) for pattern in caution_patterns):
        return "caution"
    if any(re.search(pattern, normalized) for pattern in pass_patterns):
        return "pass"
    return "unknown"


def _extract_md_gate_excerpt(text: str) -> str:
    for raw_line in text.splitlines():
        line = " ".join(raw_line.split())
        if not line:
            continue
        lowered = line.lower()
        if any(token in lowered for token in ("phase 2", "blocked", "stable", "metastable", "unstable")):
            return line[:160]
    return ""


# ---------------------------------------------------------------------------
# Validation (Task 2.0.2)
# ---------------------------------------------------------------------------

def validate_patch_reference(
    rows: List[dict],
    file_path: Path,
    md_gate_trace: Optional[Dict[str, object]] = None,
) -> Dict[str, list]:
    """Validate Phase 1 patch reference content.

    Returns a dict with keys: errors, warnings, info.
    Each value is a list of message strings.
    """
    messages: Dict[str, list] = {"errors": [], "warnings": [], "info": []}

    if not rows:
        messages["errors"].append("Patch reference is empty (0 rows).")
        return messages

    # --- Schema check ---
    actual_cols = set(rows[0].keys())
    missing_cols = set(PHASE1_REQUIRED_COLUMNS) - actual_cols
    if missing_cols:
        messages["errors"].append(
            f"Missing required columns: {sorted(missing_cols)}"
        )

    extra_cols = actual_cols - set(PHASE1_REQUIRED_COLUMNS) - {"notes"}
    if extra_cols:
        messages["info"].append(
            f"Extra columns (preserved but not required): {sorted(extra_cols)}"
        )

    # --- Per-row validation ---
    seen_residues = set()
    states_covered = set()

    for i, row in enumerate(rows, start=1):
        rid = row.get("residue_id", "")
        rnum_str = row.get("residue_num", "")
        chain = row.get("chain", "")
        lobe = row.get("lobe_label", "")
        robustness = row.get("robustness_class", "")
        confidence = row.get("confidence", "")
        construct = row.get("construct_type", "")
        orient_status = row.get("orientation_validation_status", "")
        method_agree = row.get("method_agreement", "")
        states_str = row.get("states_present", "")

        # Residue ID format
        if not rid:
            messages["errors"].append(f"Row {i}: empty residue_id")
        elif rid in seen_residues:
            messages["warnings"].append(f"Row {i}: duplicate residue_id '{rid}'")
        seen_residues.add(rid)

        # Residue number range
        rnum = _safe_int(rnum_str)
        if rnum is None:
            messages["errors"].append(f"Row {i} ({rid}): non-integer residue_num '{rnum_str}'")
        elif not (RESIDUE_NUM_RANGE[0] <= rnum <= RESIDUE_NUM_RANGE[1]):
            messages["warnings"].append(
                f"Row {i} ({rid}): residue_num {rnum} outside expected range "
                f"{RESIDUE_NUM_RANGE[0]}-{RESIDUE_NUM_RANGE[1]}"
            )

        # Chain
        if chain not in VALID_CHAINS:
            messages["errors"].append(
                f"Row {i} ({rid}): unexpected chain '{chain}' (expected {VALID_CHAINS})"
            )

        # Lobe
        if lobe not in VALID_LOBES:
            messages["warnings"].append(
                f"Row {i} ({rid}): unexpected lobe '{lobe}' (expected {VALID_LOBES})"
            )

        # Robustness
        if robustness not in VALID_ROBUSTNESS:
            messages["warnings"].append(
                f"Row {i} ({rid}): unexpected robustness_class '{robustness}'"
            )

        # Confidence
        if confidence not in VALID_CONFIDENCE:
            messages["warnings"].append(
                f"Row {i} ({rid}): unexpected confidence '{confidence}'"
            )

        # Construct type (v2 field)
        if construct and construct not in VALID_CONSTRUCT_TYPES:
            messages["warnings"].append(
                f"Row {i} ({rid}): unexpected construct_type '{construct}'"
            )

        # Orientation validation status (v2 field)
        if orient_status and orient_status != "orientation_validated":
            messages["warnings"].append(
                f"Row {i} ({rid}): orientation_validation_status = '{orient_status}' "
                "(expected 'orientation_validated')"
            )

        # Method agreement
        if method_agree and method_agree not in VALID_METHOD_AGREEMENT:
            messages["warnings"].append(
                f"Row {i} ({rid}): unexpected method_agreement '{method_agree}'"
            )

        # Occupancy range
        occ = _safe_float(row.get("global_max_occupancy", ""))
        if occ is not None and not (0.0 <= occ <= 1.0):
            messages["warnings"].append(
                f"Row {i} ({rid}): global_max_occupancy {occ} outside [0, 1]"
            )

        # Track states
        if states_str:
            for s in states_str.split(";"):
                s = s.strip()
                if s:
                    states_covered.add(s)

    # --- Cross-row checks ---
    n_high = sum(1 for r in rows if r.get("confidence") == "high")
    n_medium = sum(1 for r in rows if r.get("confidence") == "medium")
    n_low = sum(1 for r in rows if r.get("confidence") == "low")

    messages["info"].append(
        f"Total: {len(rows)} residues | high={n_high}, medium={n_medium}, low={n_low}"
    )

    # Check state coverage
    missing_states = set(RECEPTOR_STATES) - states_covered
    if missing_states:
        messages["warnings"].append(
            f"No patch residues found for state(s): {sorted(missing_states)}. "
            "This may be expected if those states had no orientation-valid models."
        )

    # Lobe distribution
    n_nlobe = sum(1 for r in rows if r.get("lobe_label") == "N-lobe")
    n_clobe = sum(1 for r in rows if r.get("lobe_label") == "C-lobe")
    messages["info"].append(f"Lobe distribution: N-lobe={n_nlobe}, C-lobe={n_clobe}")

    if n_nlobe == 0:
        messages["warnings"].append(
            "No N-lobe residues in patch reference. "
            "Full kinase domain should detect N-lobe contacts."
        )

    if md_gate_trace:
        decision_artifact_found = False
        for artifact in md_gate_trace["artifacts"]:
            status = artifact["status"]
            if artifact["name"] == "md_gate_decision.md" and status == "found":
                decision_artifact_found = True
            if status == "found":
                messages["info"].append(
                    f"MD gate artifact found: {artifact['name']} ({artifact['location']})"
                )
            elif status == "ambiguous":
                messages["warnings"].append(
                    f"MD gate artifact ambiguous: {artifact['name']} ({artifact['location']})"
                )
            else:
                messages["warnings"].append(f"MD gate artifact missing: {artifact['name']}")

        decision_status = md_gate_trace["decision_status"]
        decision_excerpt = md_gate_trace.get("decision_excerpt", "")
        if decision_status == "blocked":
            detail = "MD gate decision blocks Phase 2 progression"
            if decision_excerpt:
                detail = f"{detail}: {decision_excerpt}"
            messages["errors"].append(detail)
        elif decision_status == "caution":
            messages["warnings"].append(
                "MD gate decision is METASTABLE/caution; carry reduced confidence into Phase 2."
            )
        elif decision_status == "pass":
            messages["info"].append("MD gate decision indicates Phase 2 may proceed.")
        elif decision_artifact_found:
            messages["warnings"].append(
                "MD gate decision artifact is present but its status could not be classified."
            )

    return messages


# ---------------------------------------------------------------------------
# Standardization (Task 2.0.3)
# ---------------------------------------------------------------------------

def normalize_patch_reference(rows: List[dict]) -> List[dict]:
    """Convert Phase 1 patch reference into Phase 2 normalized format.

    Renames columns, coerces types, adds Phase 2-specific labels.
    """
    normalized = []
    for row in rows:
        orient_status = row.get("orientation_validation_status", "")
        is_orient_validated = orient_status == "orientation_validated"

        normalized.append({
            "patch_residue_id": row.get("residue_id", ""),
            "residue_num": _safe_int(row.get("residue_num", "")) or "",
            "residue_name": row.get("residue_name", ""),
            "chain": row.get("chain", ""),
            "lobe": row.get("lobe_label", ""),
            "construct_type": row.get("construct_type", "full_kinase_domain"),
            "orientation_validated": is_orient_validated,
            "robustness": row.get("robustness_class", ""),
            "n_states": _safe_int(row.get("n_states_present", "")) or 0,
            "states": row.get("states_present", ""),
            "max_occupancy": _safe_float(row.get("global_max_occupancy", "")) or 0.0,
            "is_hotspot": _parse_bool(row.get("is_hotspot_any_state", "")),
            "pyrosetta_support": _parse_bool(row.get("pyrosetta_evidence", "")),
            "lightdock_support": _parse_bool(row.get("lightdock_evidence", "")),
            "method_agreement": row.get("method_agreement", ""),
            "confidence": row.get("confidence", ""),
            "phase1_evidence_source": row.get("evidence_source", "primary"),
            "notes": row.get("notes", ""),
        })

    return normalized


# ---------------------------------------------------------------------------
# Validation report (Deliverable)
# ---------------------------------------------------------------------------

def build_validation_report(
    rows: List[dict],
    normalized: List[dict],
    messages: Dict[str, list],
    file_path: Path,
    md_gate_trace: Optional[Dict[str, object]] = None,
) -> List[str]:
    """Build Phase 2 patch reference validation markdown report."""
    lines = [
        "# Phase 2 Patch Reference Validation Report",
        "",
        "## Source",
        "",
        f"- **File**: `{file_path.name}`",
        f"- **Path**: `{file_path}`",
        f"- **Rows**: {len(rows)}",
        "",
    ]

    # Validation summary
    n_err = len(messages["errors"])
    n_warn = len(messages["warnings"])
    n_info = len(messages["info"])

    if n_err == 0 and n_warn == 0:
        verdict = "PASS — patch reference is valid and ready for Phase 2"
    elif n_err == 0:
        verdict = f"PASS WITH WARNINGS ({n_warn} warnings) — review before proceeding"
    else:
        verdict = f"FAIL ({n_err} errors) — fix Phase 1 outputs and re-run"

    lines.extend([
        "## Validation Result",
        "",
        f"**{verdict}**",
        "",
        f"- Errors: {n_err}",
        f"- Warnings: {n_warn}",
        f"- Info: {n_info}",
        "",
    ])

    # Detail sections
    if messages["errors"]:
        lines.extend(["### Errors", ""])
        for msg in messages["errors"]:
            lines.append(f"- {msg}")
        lines.append("")

    if messages["warnings"]:
        lines.extend(["### Warnings", ""])
        for msg in messages["warnings"]:
            lines.append(f"- {msg}")
        lines.append("")

    if messages["info"]:
        lines.extend(["### Info", ""])
        for msg in messages["info"]:
            lines.append(f"- {msg}")
        lines.append("")

    # Patch summary table
    lines.extend([
        "## Patch Reference Summary",
        "",
        "| Residue | Lobe | Robustness | Confidence | States | Occupancy | Hotspot | Methods |",
        "|---------|------|------------|------------|--------|-----------|---------|---------|",
    ])

    for r in normalized:
        methods = r.get("method_agreement", "")
        lines.append(
            f"| {r['patch_residue_id']} | {r['lobe']} | "
            f"{r['robustness']} | {r['confidence']} | "
            f"{r['n_states']} | {r['max_occupancy']:.2f} | "
            f"{'Yes' if r['is_hotspot'] else 'No'} | {methods} |"
        )
    lines.append("")

    # Lobe distribution
    n_nlobe = sum(1 for r in normalized if r["lobe"] == "N-lobe")
    n_clobe = sum(1 for r in normalized if r["lobe"] == "C-lobe")
    lines.extend([
        "## Lobe Distribution",
        "",
        f"- N-lobe: {n_nlobe} residues",
        f"- C-lobe: {n_clobe} residues",
        "",
    ])

    # State coverage
    all_states = set()
    for r in normalized:
        for s in r.get("states", "").split(";"):
            s = s.strip()
            if s:
                all_states.add(s)

    lines.extend([
        "## State Coverage",
        "",
        "| State | Covered |",
        "|-------|---------|",
    ])
    for st in RECEPTOR_STATES:
        covered = "Yes" if st in all_states else "No"
        lines.append(f"| {st} | {covered} |")
    lines.append("")

    if md_gate_trace:
        lines.extend(_build_md_gate_traceability_section(md_gate_trace))

    # Phase 2 readiness
    lines.extend([
        "## Phase 2 Readiness",
        "",
        f"- Normalized output: `phase2_patch_reference_normalized.csv`",
        f"- Schema: {len(PHASE2_NORMALIZED_COLUMNS)} columns",
        "- All residue numbers in PDB numbering (no offset)",
        "- Orientation validation status preserved",
        "- Construct type preserved (full_kinase_domain)",
        "",
        "---",
        "",
        "Generated by `egfr_pipeline.phase2.patch_ingestion`",
    ])

    return lines


def _build_md_gate_traceability_section(md_gate_trace: Dict[str, object]) -> List[str]:
    decision_status = md_gate_trace.get("decision_status", "missing")
    lines = [
        "## MD Gate Traceability",
        "",
        f"- Decision status: `{decision_status}`",
        "",
        "| Artifact | Status | Location |",
        "|----------|--------|----------|",
    ]
    for artifact in md_gate_trace["artifacts"]:
        lines.append(
            f"| {artifact['name']} | {artifact['status']} | {artifact['location']} |"
        )
    lines.append("")

    decision_excerpt = md_gate_trace.get("decision_excerpt", "")
    if decision_excerpt:
        lines.append(f"- Decision snippet: `{decision_excerpt}`")
        lines.append("")

    return lines


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_patch_ingestion(
    phase1_dir: Path,
    output_dir: Path,
) -> Tuple[Path, Path]:
    """Run full TG 2.0 pipeline.

    Returns (normalized_csv_path, validation_report_path).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # 2.0.1 — Ingestion
    rows, ref_path = load_patch_reference(phase1_dir)
    md_gate_trace = inspect_md_gate_artifacts(phase1_dir)
    print(f"  Loaded {len(rows)} residues from {ref_path.name}")

    # 2.0.2 — Validation
    messages = validate_patch_reference(rows, ref_path, md_gate_trace)

    n_err = len(messages["errors"])
    n_warn = len(messages["warnings"])
    n_md_found = sum(1 for artifact in md_gate_trace["artifacts"] if artifact["status"] == "found")
    print(f"  Validation: {n_err} errors, {n_warn} warnings")
    print(
        "  MD gate: "
        f"{md_gate_trace['decision_status']} "
        f"({n_md_found}/{len(md_gate_trace['artifacts'])} artifacts found)"
    )

    if n_err > 0:
        print("\n  VALIDATION ERRORS:")
        for msg in messages["errors"]:
            print(f"    - {msg}")

    # 2.0.3 — Standardization
    normalized = normalize_patch_reference(rows)
    print(f"  Normalized: {len(normalized)} residues → Phase 2 schema")

    # Write normalized CSV
    norm_path = output_dir / "phase2_patch_reference_normalized.csv"
    with open(norm_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=PHASE2_NORMALIZED_COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(normalized)

    # Write validation report
    report_lines = build_validation_report(
        rows,
        normalized,
        messages,
        ref_path,
        md_gate_trace,
    )
    report_path = output_dir / "phase2_patch_reference_validation.md"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    return norm_path, report_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _safe_float(val: str) -> Optional[float]:
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def _safe_int(val: str) -> Optional[int]:
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return None


def _parse_bool(val) -> bool:
    if isinstance(val, bool):
        return val
    return str(val).strip().lower() in ("true", "1", "yes")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 TG 2.0: Phase 1 patch reference ingestion and validation"
    )
    parser.add_argument(
        "--phase1_dir", type=Path, default=PHASE1_OUTPUT_DIR,
        help="Phase 1 output directory containing patch reference CSV",
    )
    parser.add_argument(
        "--output_dir", type=Path, default=PHASE2_OUTPUT_DIR,
        help="Phase 2 output directory",
    )
    args = parser.parse_args()

    print("Phase 2 — Task Group 2.0: Patch Reference Ingestion")
    print(f"  Phase 1 source: {args.phase1_dir}")
    print(f"  Phase 2 output: {args.output_dir}")
    print()

    norm_path, report_path = run_patch_ingestion(args.phase1_dir, args.output_dir)

    print(f"\n{'='*60}")
    print("TG 2.0 COMPLETE")
    print(f"{'='*60}")
    print(f"  Normalized CSV: {norm_path}")
    print(f"  Validation:     {report_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
