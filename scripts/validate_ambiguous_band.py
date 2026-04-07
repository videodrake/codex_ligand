#!/usr/bin/env python3
"""Retroactive orientation dot-product validation on Workflow A PPI models.

Purpose: Compute orientation_score (dot product) for all Workflow A PPI
docking models to validate the AMBIGUOUS_BAND = 0.15 threshold before
starting Workflow B.

Usage (HPC only — requires PyRosetta):
    cd /work4/hwang/onepack/my_second_project/codex_ligand2
    python scripts/validate_ambiguous_band.py

Output:
    output/workflow_a/orientation_band_validation.csv
    + console summary with distribution statistics
"""

import csv
import glob
import os
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]
N_SEEDS = 10
PPI_DOCKING_ROOT = "output/workflow_a/phase2_ppi_docking"
OUTPUT_CSV = "output/workflow_a/orientation_band_validation.csv"
CURRENT_BAND = 0.15

# ---------------------------------------------------------------------------
# Late import of PyRosetta (HPC only)
# ---------------------------------------------------------------------------


def _init_pyrosetta():
    """Initialize PyRosetta with minimal verbosity."""
    import pyrosetta  # noqa: F811
    pyrosetta.init("-mute all", silent=True)
    return pyrosetta


def main():
    pyrosetta = _init_pyrosetta()

    # Import orientation functions from existing code
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from egfr_pipeline.phase1.orientation_filter import compute_orientation_score

    # Discover all PPI model PDBs
    rows = []
    n_total = 0
    n_error = 0

    for state in STATES:
        for seed in range(N_SEEDS):
            # Actual layout: {root}/{state}/prod_seed{n}/docking_{state}_ext_beta_meander/final_result/
            seed_dir = os.path.join(PPI_DOCKING_ROOT, state, f"prod_seed{seed}")
            docking_subdir = f"docking_{state}_ext_beta_meander"
            restored_dir = os.path.join(seed_dir, docking_subdir, "restored", "final_result")
            raw_dir = os.path.join(seed_dir, docking_subdir, "final_result")

            # Fallback: glob for any docking_* subdirectory
            pdb_dir = None
            if os.path.isdir(restored_dir):
                pdb_dir = restored_dir
            elif os.path.isdir(raw_dir):
                pdb_dir = raw_dir
            else:
                # Try to discover docking subdirectory dynamically
                pattern = os.path.join(seed_dir, "docking_*", "final_result")
                candidates = sorted(glob.glob(pattern))
                if candidates:
                    pdb_dir = candidates[0]

            if pdb_dir is None or not os.path.isdir(pdb_dir):
                print(f"  SKIP {state}/seed{seed}: directory not found")
                continue

            pdbs = sorted(glob.glob(os.path.join(pdb_dir, "Rank*.pdb")))
            if not pdbs:
                print(f"  SKIP {state}/seed{seed}: no Rank*.pdb files")
                continue

            print(f"  Processing {state}/seed{seed}: {len(pdbs)} models ...", end="", flush=True)
            for pdb_path in pdbs:
                n_total += 1
                try:
                    pose = pyrosetta.pose_from_pdb(pdb_path)
                    result = compute_orientation_score(
                        pose, partner_chain="B", receptor_chain="A"
                    )
                    rows.append({
                        "state": state,
                        "seed": seed,
                        "model": os.path.basename(pdb_path),
                        "orientation_score": result["orientation_score"],
                        "orientation_class": result["orientation_class"],
                        "n_active_face_ca": result.get("n_active_face_ca", ""),
                        "n_receptor_contact_ca": result.get("n_receptor_contact_ca", ""),
                        "error": result.get("error", ""),
                    })
                except Exception as e:
                    n_error += 1
                    rows.append({
                        "state": state,
                        "seed": seed,
                        "model": os.path.basename(pdb_path),
                        "orientation_score": "",
                        "orientation_class": "load_error",
                        "n_active_face_ca": "",
                        "n_receptor_contact_ca": "",
                        "error": str(e),
                    })
            print(" done")

    # Write CSV
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    fieldnames = [
        "state", "seed", "model", "orientation_score", "orientation_class",
        "n_active_face_ca", "n_receptor_contact_ca", "error",
    ]
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {OUTPUT_CSV}  ({len(rows)} rows, {n_error} errors)")

    # Summary statistics
    scores = [
        float(r["orientation_score"])
        for r in rows
        if r["orientation_score"] != "" and r["orientation_class"] not in ("load_error",)
    ]
    if not scores:
        print("No valid scores to summarize.")
        return

    import numpy as np

    scores = np.array(scores)
    abs_scores = np.abs(scores)

    print(f"\n{'='*60}")
    print(f"ORIENTATION DOT-PRODUCT DISTRIBUTION  (n={len(scores)})")
    print(f"{'='*60}")
    print(f"  mean     = {scores.mean():.4f}")
    print(f"  median   = {np.median(scores):.4f}")
    print(f"  std      = {scores.std():.4f}")
    print(f"  min/max  = {scores.min():.4f} / {scores.max():.4f}")

    print(f"\n  AMBIGUOUS_BAND threshold analysis:")
    for band in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]:
        n_amb = int(np.sum(abs_scores < band))
        pct = 100.0 * n_amb / len(scores)
        marker = " ← current" if abs(band - CURRENT_BAND) < 0.001 else ""
        print(f"    |dot| < {band:.2f}:  {n_amb:4d} / {len(scores)}  ({pct:5.1f}%){marker}")

    # Class distribution
    classes = {}
    for r in rows:
        c = r["orientation_class"]
        classes[c] = classes.get(c, 0) + 1
    print(f"\n  Class distribution:")
    for c in ["pass", "fail", "ambiguous", "insufficient_data", "load_error"]:
        if c in classes:
            pct = 100.0 * classes[c] / len(rows)
            print(f"    {c:20s}: {classes[c]:4d}  ({pct:5.1f}%)")

    # Per-state breakdown
    print(f"\n  Per-state pass rate:")
    for state in STATES:
        st_rows = [r for r in rows if r["state"] == state and r["orientation_class"] in ("pass", "fail", "ambiguous")]
        if st_rows:
            n_pass = sum(1 for r in st_rows if r["orientation_class"] == "pass")
            print(f"    {state:20s}: {n_pass}/{len(st_rows)} pass ({100*n_pass/len(st_rows):.1f}%)")

    print(f"\n{'='*60}")
    print("Recommendation: If >10% models fall in ambiguous band,")
    print("consider narrowing AMBIGUOUS_BAND. If <2%, consider widening.")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
