#!/usr/bin/env python3
"""Phase 1 Task 1.1.5: Standardize score extraction from PyRosetta outputs.

Reads the existing pipeline output CSVs (final_ranking.csv, cluster_summary.csv,
all_scored_summary.csv) and produces a standardized pyrosetta_decoy_scores.csv
that includes all required Phase 1 columns:

  - decoy_id, receptor_id, partner_id, construct_type
  - total_score, I_sc (= dG_separated, our best interface score proxy)
  - dG_separated, dSASA, sc_value, packstat
  - nres_int, delta_unsatHbonds, hbonds_int
  - L_RMSD, cluster_id, rank

This utility runs post-docking and consolidates outputs from potentially
multiple seeds into a single standardized table per receptor state.

Usage:
    # Standardize scores for all completed runs
    python -m egfr_pipeline.phase1.standardize_scores

    # Single state
    python -m egfr_pipeline.phase1.standardize_scores --state 3GT8_raw

    # Custom output directory
    python -m egfr_pipeline.phase1.standardize_scores \
        --output_dir output/workflow_b/phase1_ppi_analysis
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

from egfr_pipeline import paths

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_PATH_CONFIG = {"output_root": str(PROJECT_ROOT / "output")}
PHASE1_RUNS_DIR = paths.wa_phase2_ppi_docking(DEFAULT_PATH_CONFIG)
PHASE1_OUTPUT_DIR = paths.wb_phase1_ppi_analysis(DEFAULT_PATH_CONFIG)

RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]

# Standardized output columns (Task 1.1.5)
SCORE_COLUMNS = [
    "decoy_id",
    "receptor_id",
    "partner_id",
    "construct_type",
    "seed_index",
    "run_type",
    "Rank",
    "cluster_id",
    "total_score",
    "I_sc",           # Interface score — mapped from dG_separated
    "dG_separated",
    "dSASA",
    "dG_density",
    "sc_value",
    "packstat",
    "delta_unsatHbonds",
    "nres_int",
    "hbonds_int",
    "L_RMSD",
    "L_RMSD_best",
    "Binding_Residues_A",
    "Binding_Residues_B",
    "key_contact_ratio",
    "source_file",
]


def parse_run_metadata(run_dir: Path) -> Optional[dict]:
    """Read pyrosetta_run_metadata.json from a run directory."""
    meta_path = run_dir / "pyrosetta_run_metadata.json"
    if not meta_path.exists():
        return None
    with open(meta_path, encoding="utf-8") as f:
        return json.load(f)


def find_final_ranking(run_dir: Path) -> Optional[Path]:
    """Find the final_ranking.csv in a pipeline output directory.

    The pipeline creates output in a subdirectory named after the PDB file.
    """
    # Direct path
    direct = run_dir / "final_ranking.csv"
    if direct.exists():
        return direct

    # Search in subdirectories (pipeline creates <pdb_name>/ subdir)
    for subdir in run_dir.iterdir():
        if subdir.is_dir():
            candidate = subdir / "final_ranking.csv"
            if candidate.exists():
                return candidate
            # Also check final_result/
            candidate2 = subdir / "final_result" / "final_ranking.csv"
            if candidate2.exists():
                return candidate2
    return None


def find_all_scored_summary(run_dir: Path) -> Optional[Path]:
    """Find all_scored_summary.csv (energy funnel data)."""
    direct = run_dir / "all_scored_summary.csv"
    if direct.exists():
        return direct
    for subdir in run_dir.iterdir():
        if subdir.is_dir():
            candidate = subdir / "all_scored_summary.csv"
            if candidate.exists():
                return candidate
    return None


def standardize_final_ranking(
    csv_path: Path,
    metadata: dict,
) -> List[dict]:
    """Convert a final_ranking.csv to standardized format.

    Maps existing columns to the Phase 1 standardized schema.
    """
    rows = []
    with open(csv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Parse cluster_id from Parent field (e.g., "C01_M01" → "C01")
            parent = row.get("Parent", "")
            cluster_id = ""
            if parent and "_" in parent:
                cluster_id = parent.split("_")[0]

            standardized = {
                "decoy_id": row.get("File_PDB", ""),
                "receptor_id": metadata.get("receptor_id", ""),
                "partner_id": metadata.get("partner_id", ""),
                "construct_type": metadata.get("construct_type", ""),
                "seed_index": metadata.get("seed_index", 0),
                "run_type": "production" if metadata.get("is_production") else "test",
                "Rank": row.get("Rank", ""),
                "cluster_id": cluster_id,
                "total_score": row.get("Total_Score", ""),
                # I_sc: We use dG_separated as the interface score proxy.
                # True I_sc (InterfaceScoreCalculator) and dG_separated
                # (InterfaceAnalyzerMover) are highly correlated for
                # rigid-body docking. dG_separated is already computed
                # and validated in the existing pipeline.
                "I_sc": row.get("dG_separated", ""),
                "dG_separated": row.get("dG_separated", ""),
                "dSASA": row.get("dSASA", ""),
                "dG_density": row.get("dG_density", ""),
                "sc_value": row.get("sc_value", ""),
                "packstat": row.get("packstat", ""),
                "delta_unsatHbonds": row.get("delta_unsatHbonds", ""),
                "nres_int": row.get("nres_int", ""),
                "hbonds_int": row.get("hbonds_int", ""),
                "L_RMSD": row.get("L_RMSD", ""),
                "L_RMSD_best": row.get("L_RMSD_best", ""),
                "Binding_Residues_A": row.get("Binding_Residues_A", ""),
                "Binding_Residues_B": row.get("Binding_Residues_B", ""),
                "key_contact_ratio": row.get("key_contact_ratio", ""),
                "source_file": str(csv_path.relative_to(PROJECT_ROOT)),
            }
            rows.append(standardized)

    return rows


def consolidate_state(
    state_name: str,
    output_dir: Path,
    runs_dir: Optional[Path] = None,
) -> Optional[Path]:
    """Consolidate all runs for a receptor state into one standardized CSV.

    Merges test + production + multi-seed outputs.
    """
    source_base = runs_dir or output_dir
    source_state_dir = source_base / state_name
    if not source_state_dir.exists():
        return None
    state_dir = output_dir / state_name
    state_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []

    # Scan all run directories
    for run_dir in sorted(source_state_dir.iterdir()):
        if not run_dir.is_dir():
            continue

        metadata = parse_run_metadata(run_dir)
        if metadata is None:
            print(f"    [SKIP] {run_dir.name}: no run metadata")
            continue

        csv_path = find_final_ranking(run_dir)
        if csv_path is None:
            status = metadata.get("status", "unknown")
            if status == "completed":
                print(f"    [WARN] {run_dir.name}: completed but no final_ranking.csv")
            else:
                print(f"    [SKIP] {run_dir.name}: status={status}")
            continue

        rows = standardize_final_ranking(csv_path, metadata)
        all_rows.extend(rows)
        print(f"    [OK] {run_dir.name}: {len(rows)} models")

    if not all_rows:
        print(f"  No completed runs found for {state_name}")
        return None

    # Write consolidated CSV
    out_path = state_dir / "pyrosetta_decoy_scores.csv"
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=SCORE_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"  Consolidated: {out_path} ({len(all_rows)} total models)")
    return out_path


def run_standardize_scores(
    *,
    output_dir: Path = PHASE1_OUTPUT_DIR,
    runs_dir: Path = PHASE1_RUNS_DIR,
    state: Optional[str] = None,
) -> List[Path]:
    states = [state] if state else RECEPTOR_STATES
    results: List[Path] = []
    for state_name in states:
        path = consolidate_state(state_name, output_dir, runs_dir=runs_dir)
        if path:
            results.append(path)
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Phase 1 TG 1.1.5: Standardize PyRosetta score extraction"
    )
    parser.add_argument("--state", choices=RECEPTOR_STATES,
                        help="Process a single receptor state")
    parser.add_argument("--runs_dir", type=Path,
                        default=PHASE1_RUNS_DIR,
                        help="Workflow A Phase 2 PyRosetta run directory")
    parser.add_argument("--output_dir", type=Path,
                        default=PHASE1_OUTPUT_DIR,
                        help="Workflow B Phase 1 analysis output directory")
    args = parser.parse_args()

    print("Phase 1 — Task 1.1.5: Score Standardization")
    print(f"Runs base:   {args.runs_dir}")
    print(f"Output base: {args.output_dir}")

    states = [args.state] if args.state else RECEPTOR_STATES
    results = []

    for state in states:
        print(f"\nProcessing {state}:")
        path = consolidate_state(state, args.output_dir, runs_dir=args.runs_dir)
        if path:
            results.append((state, path))

    if not results:
        print("\nNo completed runs found to standardize.")
        print("Run docking first:")
        print("  python -m egfr_pipeline.phase1.launch_docking --test --dry-run")
        return 0

    print(f"\n{'='*60}")
    print(f"STANDARDIZATION COMPLETE")
    print(f"{'='*60}")
    for state, path in results:
        print(f"  {state}: {path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
