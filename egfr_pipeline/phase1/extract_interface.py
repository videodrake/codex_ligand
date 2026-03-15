#!/usr/bin/env python3
"""Phase 1 Task 1.2: Receptor-side interface residue extraction.

Reads PyRosetta pipeline outputs (final_ranking.csv, InterfaceEnergies CSVs)
and produces structured, per-model long-form tables that preserve:
  - Receptor-side and partner-side residues (separate columns)
  - N-lobe / C-lobe annotation per receptor residue
  - Per-residue interface energy (delta E) when available
  - Model/decoy ID for traceability back to score tables

Interface extraction rule (Task 1.2.1):
  - All-atom NeighborhoodResidueSelector, 6 Å cutoff (same as pipeline)
  - Already computed during pipeline execution and stored in
    Binding_Residues_A / Binding_Residues_B columns of final_ranking.csv
  - This module parses those columns — no re-extraction from PDB needed

Deliverables:
  - pyrosetta_interface_residue_table.csv  (long-form: one row per model × residue)
  - pyrosetta_interface_models.csv         (per-model summary)

Usage:
    python -m egfr_pipeline.phase1.extract_interface [--output_dir output/phase1_ppi]
"""

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE1_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase1_ppi"
RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]

# N-lobe / C-lobe boundary (PDB numbering)
# The hinge region connecting N-lobe and C-lobe is at ~residue 838
NLOBE_CLOBE_BOUNDARY = 838

# Interface extraction rule documentation
INTERFACE_RULE = {
    "method": "NeighborhoodResidueSelector (all-atom)",
    "cutoff_angstrom": 6.0,
    "source": "PyRosetta InterfaceAnalyzer pipeline (analysis.py)",
    "description": (
        "All-atom distance-based interface detection. A residue is "
        "classified as an interface residue if any of its atoms is within "
        "6.0 Å of any atom in the partner chain. This is computed during "
        "pipeline execution and stored in Binding_Residues_A/B columns."
    ),
}

# Residue name normalization
_RESNAME_NORMALIZE = {
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "CYX": "CYS", "HIE": "HIS", "HID": "HIS", "HIP": "HIS",
}

# Output schemas
RESIDUE_TABLE_COLUMNS = [
    "model_id",
    "receptor_id",
    "partner_id",
    "construct_type",
    "seed_index",
    "rank",
    "cluster_id",
    "chain",           # A (receptor) or B (partner)
    "residue_id",      # e.g., "LEU819"
    "residue_num",     # e.g., 819
    "residue_name",    # e.g., "LEU"
    "lobe_label",      # N-lobe / C-lobe / partner / unknown
    "delta_e_total",   # Per-residue interface energy (from InterfaceEnergies CSV)
    "delta_e_fa_atr",
    "delta_e_fa_rep",
    "delta_e_fa_sol",
    "delta_e_fa_elec",
    "source_file",
]

MODEL_TABLE_COLUMNS = [
    "model_id",
    "receptor_id",
    "partner_id",
    "construct_type",
    "seed_index",
    "rank",
    "cluster_id",
    "dG_separated",
    "dSASA",
    "sc_value",
    "packstat",
    "nres_int",
    "n_receptor_interface_residues",
    "n_partner_interface_residues",
    "n_nlobe_interface_residues",
    "n_clobe_interface_residues",
    "receptor_interface_residues",    # Semicolon-separated
    "partner_interface_residues",     # Semicolon-separated
    "source_file",
]


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def normalize_resname(name: str) -> str:
    """Normalize residue name (HSD→HIS, etc.)."""
    return _RESNAME_NORMALIZE.get(name, name)


def parse_residue_ref(ref: str) -> Tuple[str, str, Optional[int]]:
    """Parse 'A:LEU819' → (chain='A', residue_id='LEU819', resnum=819).

    Also handles 'LEU819' (no chain prefix).
    """
    chain = ""
    if ":" in ref:
        chain, ref = ref.split(":", 1)

    # Split into name + number
    match = re.match(r"([A-Z]+)(-?\d+)", ref)
    if match:
        resname = normalize_resname(match.group(1))
        resnum = int(match.group(2))
        return chain, f"{resname}{resnum}", resnum
    return chain, ref, None


def parse_binding_residues_column(raw: str) -> List[Tuple[str, str, Optional[int]]]:
    """Parse 'A:LEU819,A:ALA822' → [(chain, residue_id, resnum), ...]."""
    if not raw or raw in ("None", "No_Chain_2", "Analysis_Failed", ""):
        return []
    results = []
    for token in raw.split(","):
        token = token.strip()
        if token:
            results.append(parse_residue_ref(token))
    return results


def get_lobe_label(chain: str, resnum: Optional[int]) -> str:
    """Assign lobe label based on chain and residue number."""
    if chain == "B" or chain == "partner":
        return "partner"
    if resnum is None:
        return "unknown"
    if resnum < NLOBE_CLOBE_BOUNDARY:
        return "N-lobe"
    else:
        return "C-lobe"


# ---------------------------------------------------------------------------
# InterfaceEnergies CSV parsing
# ---------------------------------------------------------------------------

def load_interface_energies(csv_path: Path) -> Dict[str, dict]:
    """Load per-residue interface energies from InterfaceEnergies CSV.

    Returns {residue_key: {delta_e_total, delta_e_fa_atr, ...}}.
    The residue_key is "A:LEU819" or "B:VAL962".
    """
    energies = {}
    if not csv_path.exists():
        return energies

    with open(csv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            resid = row.get("Residue_ID", "")
            resname = row.get("Residue_Name", "")
            chain = row.get("Chain", "")
            key = f"{chain}:{normalize_resname(resname)}{resid}"
            energies[key] = {
                "delta_e_total": row.get("DeltaE_total", ""),
                "delta_e_fa_atr": row.get("DeltaE_fa_atr", ""),
                "delta_e_fa_rep": row.get("DeltaE_fa_rep", ""),
                "delta_e_fa_sol": row.get("DeltaE_fa_sol", ""),
                "delta_e_fa_elec": row.get("DeltaE_fa_elec", ""),
            }

    return energies


# ---------------------------------------------------------------------------
# Per-run extraction
# ---------------------------------------------------------------------------

def extract_run(
    run_dir: Path,
    metadata: dict,
) -> Tuple[List[dict], List[dict]]:
    """Extract interface residues from a single pipeline run.

    Returns (residue_rows, model_rows).
    """
    receptor_id = metadata.get("receptor_id", "")
    partner_id = metadata.get("partner_id", "")
    construct_type = metadata.get("construct_type", "")
    seed_index = metadata.get("seed_index", 0)

    residue_rows = []
    model_rows = []

    # Find final_ranking.csv (may be in a pipeline-created subdirectory)
    final_csv = _find_csv(run_dir, "final_ranking.csv")
    if final_csv is None:
        return residue_rows, model_rows

    # Find the directory containing InterfaceEnergies CSVs
    result_dir = final_csv.parent

    with open(final_csv, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rank = row.get("Rank", "")
            parent = row.get("Parent", "")
            cluster_id = parent.split("_")[0] if parent and "_" in parent else parent
            file_pdb = row.get("File_PDB", "")
            model_id = file_pdb if file_pdb else f"Rank{rank}_{parent}"

            # Parse receptor-side residues
            binding_a_raw = row.get("Binding_Residues_A", "") or row.get("Binding_Residues", "")
            receptor_residues = parse_binding_residues_column(binding_a_raw)

            # Parse partner-side residues
            binding_b_raw = row.get("Binding_Residues_B", "")
            partner_residues = parse_binding_residues_column(binding_b_raw)

            # Load per-residue interface energies if available
            iface_energies = {}
            file_csv = row.get("File_CSV", "")
            if file_csv:
                iface_csv_name = Path(file_csv).name.replace("_Energies.csv", "_InterfaceEnergies.csv")
                iface_csv_path = result_dir / iface_csv_name
                iface_energies = load_interface_energies(iface_csv_path)

            # --- Build long-form residue rows ---
            # Receptor-side
            for chain, res_id, resnum in receptor_residues:
                if not chain:
                    chain = "A"
                lobe = get_lobe_label(chain, resnum)
                resname_match = re.match(r"([A-Z]+)", res_id)
                resname = resname_match.group(1) if resname_match else ""

                energy_key = f"{chain}:{res_id}"
                energy = iface_energies.get(energy_key, {})

                residue_rows.append({
                    "model_id": model_id,
                    "receptor_id": receptor_id,
                    "partner_id": partner_id,
                    "construct_type": construct_type,
                    "seed_index": seed_index,
                    "rank": rank,
                    "cluster_id": cluster_id,
                    "chain": chain,
                    "residue_id": res_id,
                    "residue_num": resnum if resnum is not None else "",
                    "residue_name": resname,
                    "lobe_label": lobe,
                    "delta_e_total": energy.get("delta_e_total", ""),
                    "delta_e_fa_atr": energy.get("delta_e_fa_atr", ""),
                    "delta_e_fa_rep": energy.get("delta_e_fa_rep", ""),
                    "delta_e_fa_sol": energy.get("delta_e_fa_sol", ""),
                    "delta_e_fa_elec": energy.get("delta_e_fa_elec", ""),
                    "source_file": str(final_csv.relative_to(PROJECT_ROOT)),
                })

            # Partner-side
            for chain, res_id, resnum in partner_residues:
                if not chain:
                    chain = "B"
                resname_match = re.match(r"([A-Z]+)", res_id)
                resname = resname_match.group(1) if resname_match else ""

                energy_key = f"{chain}:{res_id}"
                energy = iface_energies.get(energy_key, {})

                residue_rows.append({
                    "model_id": model_id,
                    "receptor_id": receptor_id,
                    "partner_id": partner_id,
                    "construct_type": construct_type,
                    "seed_index": seed_index,
                    "rank": rank,
                    "cluster_id": cluster_id,
                    "chain": chain,
                    "residue_id": res_id,
                    "residue_num": resnum if resnum is not None else "",
                    "residue_name": resname,
                    "lobe_label": "partner",
                    "delta_e_total": energy.get("delta_e_total", ""),
                    "delta_e_fa_atr": energy.get("delta_e_fa_atr", ""),
                    "delta_e_fa_rep": energy.get("delta_e_fa_rep", ""),
                    "delta_e_fa_sol": energy.get("delta_e_fa_sol", ""),
                    "delta_e_fa_elec": energy.get("delta_e_fa_elec", ""),
                    "source_file": str(final_csv.relative_to(PROJECT_ROOT)),
                })

            # --- Build per-model summary row ---
            n_receptor = len(receptor_residues)
            n_partner = len(partner_residues)
            n_nlobe = sum(1 for _, _, rn in receptor_residues
                          if rn is not None and rn < NLOBE_CLOBE_BOUNDARY)
            n_clobe = sum(1 for _, _, rn in receptor_residues
                          if rn is not None and rn >= NLOBE_CLOBE_BOUNDARY)

            receptor_res_str = ";".join(r[1] for r in receptor_residues)
            partner_res_str = ";".join(r[1] for r in partner_residues)

            model_rows.append({
                "model_id": model_id,
                "receptor_id": receptor_id,
                "partner_id": partner_id,
                "construct_type": construct_type,
                "seed_index": seed_index,
                "rank": rank,
                "cluster_id": cluster_id,
                "dG_separated": row.get("dG_separated", ""),
                "dSASA": row.get("dSASA", ""),
                "sc_value": row.get("sc_value", ""),
                "packstat": row.get("packstat", ""),
                "nres_int": row.get("nres_int", ""),
                "n_receptor_interface_residues": n_receptor,
                "n_partner_interface_residues": n_partner,
                "n_nlobe_interface_residues": n_nlobe,
                "n_clobe_interface_residues": n_clobe,
                "receptor_interface_residues": receptor_res_str,
                "partner_interface_residues": partner_res_str,
                "source_file": str(final_csv.relative_to(PROJECT_ROOT)),
            })

    return residue_rows, model_rows


def _find_csv(run_dir: Path, filename: str) -> Optional[Path]:
    """Find a CSV file in a run directory or its subdirectories."""
    # Direct
    direct = run_dir / filename
    if direct.exists():
        return direct
    # In pipeline output subdirectory
    for subdir in run_dir.iterdir():
        if subdir.is_dir():
            candidate = subdir / filename
            if candidate.exists():
                return candidate
            # final_result subdirectory
            candidate2 = subdir / "final_result" / filename
            if candidate2.exists():
                return candidate2
    return None


# ---------------------------------------------------------------------------
# State-level consolidation
# ---------------------------------------------------------------------------

def process_state(
    state_name: str,
    output_base: Path,
) -> Tuple[Optional[Path], Optional[Path]]:
    """Process all runs for a receptor state.

    Returns (residue_table_path, model_table_path) or (None, None).
    """
    state_dir = output_base / state_name
    if not state_dir.exists():
        print(f"  State directory not found: {state_dir}")
        return None, None

    all_residue_rows = []
    all_model_rows = []

    for run_dir in sorted(state_dir.iterdir()):
        if not run_dir.is_dir():
            continue

        # Load run metadata
        meta_path = run_dir / "pyrosetta_run_metadata.json"
        if not meta_path.exists():
            continue

        with open(meta_path, encoding="utf-8") as f:
            metadata = json.load(f)

        if metadata.get("dry_run"):
            print(f"    [SKIP] {run_dir.name}: dry run")
            continue

        res_rows, mod_rows = extract_run(run_dir, metadata)
        if res_rows:
            all_residue_rows.extend(res_rows)
            all_model_rows.extend(mod_rows)
            print(f"    [OK] {run_dir.name}: {len(mod_rows)} models, {len(res_rows)} residue entries")
        else:
            print(f"    [SKIP] {run_dir.name}: no final_ranking.csv")

    if not all_residue_rows:
        print(f"  No data for {state_name}")
        return None, None

    # Write residue table
    res_path = state_dir / "pyrosetta_interface_residue_table.csv"
    _write_csv(res_path, all_residue_rows, RESIDUE_TABLE_COLUMNS)

    # Write model table
    mod_path = state_dir / "pyrosetta_interface_models.csv"
    _write_csv(mod_path, all_model_rows, MODEL_TABLE_COLUMNS)

    # Print summary
    n_models = len(all_model_rows)
    n_receptor_entries = sum(1 for r in all_residue_rows if r["chain"] == "A")
    n_partner_entries = sum(1 for r in all_residue_rows if r["chain"] == "B")
    unique_receptor_res = len(set(r["residue_id"] for r in all_residue_rows if r["chain"] == "A"))
    unique_partner_res = len(set(r["residue_id"] for r in all_residue_rows if r["chain"] == "B"))

    # Lobe distribution
    n_nlobe = sum(1 for r in all_residue_rows if r["lobe_label"] == "N-lobe")
    n_clobe = sum(1 for r in all_residue_rows if r["lobe_label"] == "C-lobe")

    print(f"  Summary for {state_name}:")
    print(f"    Models: {n_models}")
    print(f"    Receptor residue entries: {n_receptor_entries} ({unique_receptor_res} unique)")
    print(f"    Partner residue entries:  {n_partner_entries} ({unique_partner_res} unique)")
    print(f"    N-lobe entries: {n_nlobe}, C-lobe entries: {n_clobe}")
    print(f"    Residue table: {res_path}")
    print(f"    Model table:   {mod_path}")

    return res_path, mod_path


def _write_csv(path: Path, rows: List[dict], columns: List[str]) -> None:
    """Write rows to CSV with specified columns."""
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
        description="Phase 1 TG 1.2: Extract receptor-side interface residues"
    )
    parser.add_argument("--state", choices=RECEPTOR_STATES,
                        help="Process a single receptor state")
    parser.add_argument("--output_dir", type=Path,
                        default=PHASE1_OUTPUT_DIR,
                        help="Phase 1 output base directory")
    args = parser.parse_args()

    print("Phase 1 — Task 1.2: Interface Residue Extraction")
    print(f"Output base: {args.output_dir}")
    print(f"\nInterface extraction rule:")
    print(f"  Method: {INTERFACE_RULE['method']}")
    print(f"  Cutoff: {INTERFACE_RULE['cutoff_angstrom']} Å")
    print(f"  N-lobe/C-lobe boundary: residue {NLOBE_CLOBE_BOUNDARY} (PDB numbering)")

    states = [args.state] if args.state else RECEPTOR_STATES
    results = []

    for state in states:
        print(f"\nProcessing {state}:")
        res_path, mod_path = process_state(state, args.output_dir)
        if res_path:
            results.append((state, res_path, mod_path))

    if not results:
        print("\nNo completed docking runs found to extract from.")
        print("Run docking first, then re-run this extraction.")
        return 0

    print(f"\n{'='*60}")
    print("EXTRACTION COMPLETE")
    print(f"{'='*60}")
    for state, res_path, mod_path in results:
        print(f"  {state}:")
        print(f"    Residue table: {res_path.name}")
        print(f"    Model table:   {mod_path.name}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
