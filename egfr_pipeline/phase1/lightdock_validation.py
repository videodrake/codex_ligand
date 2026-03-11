#!/usr/bin/env python3
"""Phase 1 Task 1.4: LightDock secondary validation.

Independent secondary validation using LightDock to provide method-independence
evidence for receptor-side patch identification.

Sub-modules:
  1.4.1 - LightDock setup and execution (server-side)
  1.4.2 - LightDock interface extraction (from top-ranked swarm poses)
  1.4.3 - Cross-method convergence analysis (PyRosetta vs LightDock)

LightDock overview:
  - Swarm-based docking: distributes starting poses (swarms) across the
    receptor surface, each swarm runs an independent optimization
  - Scoring: DFIRE2 or fastdfire (knowledge-based, fast)
  - Output: top poses per swarm, ranked by scoring function
  - Interface residues extracted by contact distance (CA-CA < 10 Å)

Dependencies:
  - lightdock3 (pip install lightdock3) — server-side only
  - PDB files from Phase 1 TG 1.0 inputs
  - PyRosetta consensus from TG 1.3 (for convergence comparison)

Usage:
    # Setup (generates run script)
    python -m egfr_pipeline.phase1.lightdock_validation --setup --state 3GT8_raw

    # Extract interfaces (after LightDock completes)
    python -m egfr_pipeline.phase1.lightdock_validation --extract --state 3GT8_raw

    # Cross-method convergence
    python -m egfr_pipeline.phase1.lightdock_validation --convergence --state 3GT8_raw
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
PHASE1_INPUT_DIR = PROJECT_ROOT / "input" / "PPI" / "phase1"
RECEPTOR_STATES = ["3GT8_raw", "3GT8_cl38_48", "3GT8_cl85_100"]

# LightDock parameters
LIGHTDOCK_DEFAULTS = {
    "scoring_function": "fastdfire",
    "n_swarms": 400,        # Number of starting swarm positions
    "n_glowworms": 200,     # Glowworms per swarm (population size)
    "n_steps": 100,         # Optimization steps per swarm
    "n_top_poses": 10,      # Top poses to extract per swarm
    "contact_cutoff_A": 10.0,   # CA-CA cutoff for interface residues
    "receptor_chain": "A",
    "partner_chain": "B",
}

# Output schemas
LIGHTDOCK_INTERFACE_COLUMNS = [
    "model_id",
    "receptor_id",
    "swarm_id",
    "pose_rank",
    "scoring_value",
    "chain",
    "residue_id",
    "residue_num",
    "residue_name",
    "lobe_label",
    "source",           # "lightdock"
]

LIGHTDOCK_MODEL_SUMMARY_COLUMNS = [
    "model_id",
    "receptor_id",
    "swarm_id",
    "pose_rank",
    "scoring_value",
    "n_receptor_interface_residues",
    "n_partner_interface_residues",
    "n_nlobe_interface_residues",
    "n_clobe_interface_residues",
    "receptor_interface_residues",   # Semicolon-separated
    "partner_interface_residues",    # Semicolon-separated
    "source",
]

CONVERGENCE_COLUMNS = [
    "receptor_id",
    "chain",
    "residue_id",
    "residue_num",
    "residue_name",
    "lobe_label",
    "in_pyrosetta",         # True/False
    "in_lightdock",         # True/False
    "pyrosetta_max_occupancy",
    "lightdock_frequency",  # Fraction of top LightDock models with this residue
    "convergence_class",    # convergent / pyrosetta_only / lightdock_only
    "method_agreement",     # both / single
]

NLOBE_CLOBE_BOUNDARY = 838


# ---------------------------------------------------------------------------
# 1.4.1 LightDock setup
# ---------------------------------------------------------------------------

def generate_lightdock_setup(
    state_name: str,
    output_base: Path,
    params: dict = None,
) -> Path:
    """Generate LightDock run configuration and shell script.

    Returns path to run script.
    """
    if params is None:
        params = dict(LIGHTDOCK_DEFAULTS)

    state_dir = output_base / state_name / "lightdock"
    state_dir.mkdir(parents=True, exist_ok=True)

    # Input PDB
    input_pdb = PHASE1_INPUT_DIR / f"docking_{state_name}_ext_beta_meander.pdb"

    # Metadata
    metadata = {
        "receptor_id": state_name,
        "partner_id": "extended_beta_meander_955_1006",
        "construct_type": "full_kinase_domain",
        "method": "LightDock",
        "scoring_function": params["scoring_function"],
        "n_swarms": params["n_swarms"],
        "n_glowworms": params["n_glowworms"],
        "n_steps": params["n_steps"],
        "n_top_poses": params["n_top_poses"],
        "contact_cutoff_A": params["contact_cutoff_A"],
        "input_pdb": str(input_pdb),
        "phase": "Phase 1: PPI-First Interface Mapping",
        "task_group": "TG 1.4: LightDock Secondary Validation",
        "output_dir": str(state_dir),
        "role": "secondary_validation (independent method, not primary evidence)",
    }

    meta_path = state_dir / "lightdock_run_metadata.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    # Generate run script
    # LightDock expects separate receptor and partner PDB files.
    # The docking pair PDB has chain A (receptor) and chain B (partner).
    # We need to split them before running LightDock.
    script_path = state_dir / f"run_lightdock_{state_name}.sh"
    script_content = f"""#!/bin/bash
# LightDock run script for {state_name}
# Generated by egfr_pipeline.phase1.lightdock_validation
#
# Prerequisites:
#   pip install lightdock3
#   or: conda install -c bioconda lightdock
#
# Usage:
#   cd {state_dir}
#   bash run_lightdock_{state_name}.sh

set -e

INPUT_PDB="{input_pdb}"
WORKDIR="{state_dir}"
cd "$WORKDIR"

echo "=== LightDock Secondary Validation: {state_name} ==="
echo "Input: $INPUT_PDB"
echo "Scoring: {params['scoring_function']}"
echo "Swarms: {params['n_swarms']}, Glowworms: {params['n_glowworms']}, Steps: {params['n_steps']}"

# Step 1: Split PDB into receptor (chain A) and partner (chain B)
echo "Splitting input PDB..."
python3 -c "
with open('$INPUT_PDB') as f:
    lines = f.readlines()
with open('receptor.pdb', 'w') as r, open('partner.pdb', 'w') as p:
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            chain = line[21]
            if chain == 'A':
                r.write(line)
            elif chain == 'B':
                p.write(line)
    r.write('END\\n')
    p.write('END\\n')
"

# Step 2: LightDock setup
echo "Running lightdock3_setup..."
lightdock3_setup.py receptor.pdb partner.pdb \\
    -s {params['n_swarms']} \\
    -g {params['n_glowworms']} \\
    --noxt --noh

# Step 3: LightDock run
echo "Running lightdock3..."
lightdock3.py setup.json {params['n_steps']} \\
    -s {params['scoring_function']} \\
    -c {os.cpu_count() or 16}

# Step 4: Generate top poses
echo "Generating top poses..."
lgd_generate_conformations.py receptor.pdb partner.pdb \\
    --top {params['n_top_poses']}

# Step 5: Rank and cluster
echo "Ranking results..."
lgd_cluster_bsas.py setup.json

echo "=== LightDock complete for {state_name} ==="
echo "Run extraction next:"
echo "  python -m egfr_pipeline.phase1.lightdock_validation --extract --state {state_name}"
"""

    with open(script_path, "w") as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)

    print(f"  LightDock setup for {state_name}:")
    print(f"    Metadata: {meta_path}")
    print(f"    Run script: {script_path}")
    print(f"    Swarms: {params['n_swarms']}, Steps: {params['n_steps']}")
    print(f"    Scoring: {params['scoring_function']}")

    return script_path


# ---------------------------------------------------------------------------
# 1.4.2 Interface extraction from LightDock output
# ---------------------------------------------------------------------------

def extract_lightdock_interfaces(
    state_name: str,
    output_base: Path,
    contact_cutoff: float = 10.0,
) -> Tuple[Optional[Path], Optional[Path]]:
    """Extract interface residues from LightDock top-ranked poses.

    LightDock output structure:
      swarm_0/
        lightdock_0.pdb, lightdock_1.pdb, ...
      swarm_1/
        ...
      rank_by_scoring.list  (global ranking)

    Returns (interface_table_path, model_summary_path) or (None, None).
    """
    state_dir = output_base / state_name / "lightdock"
    if not state_dir.exists():
        print(f"  LightDock directory not found: {state_dir}")
        return None, None

    # Find ranked poses
    rank_file = state_dir / "rank_by_scoring.list"
    if not rank_file.exists():
        # Try alternative: look for generated conformations
        rank_file = state_dir / "rank_by_cluster.list"
    if not rank_file.exists():
        print(f"  No ranking file found in {state_dir}")
        print(f"  Run LightDock first, then re-run extraction.")
        return None, None

    # Parse ranking file
    # Format: "swarm_N/lightdock_M.pdb   score"
    ranked_poses = _parse_lightdock_ranking(rank_file)
    if not ranked_poses:
        print(f"  No ranked poses found in {rank_file}")
        return None, None

    print(f"  Found {len(ranked_poses)} ranked poses")

    interface_rows = []
    model_rows = []

    for pose_rank, (pdb_path, score, swarm_id) in enumerate(ranked_poses, 1):
        full_path = state_dir / pdb_path
        if not full_path.exists():
            continue

        model_id = f"swarm{swarm_id}_rank{pose_rank}"

        # Extract interface residues by CA-CA distance
        receptor_res, partner_res = _extract_contacts_from_pdb(
            full_path, "A", "B", contact_cutoff
        )

        # Build per-residue rows
        n_nlobe = 0
        n_clobe = 0
        for rid, rnum, rname in receptor_res:
            lobe = "N-lobe" if rnum < NLOBE_CLOBE_BOUNDARY else "C-lobe"
            if lobe == "N-lobe":
                n_nlobe += 1
            else:
                n_clobe += 1
            interface_rows.append({
                "model_id": model_id,
                "receptor_id": state_name,
                "swarm_id": swarm_id,
                "pose_rank": pose_rank,
                "scoring_value": score,
                "chain": "A",
                "residue_id": rid,
                "residue_num": rnum,
                "residue_name": rname,
                "lobe_label": lobe,
                "source": "lightdock",
            })

        for rid, rnum, rname in partner_res:
            interface_rows.append({
                "model_id": model_id,
                "receptor_id": state_name,
                "swarm_id": swarm_id,
                "pose_rank": pose_rank,
                "scoring_value": score,
                "chain": "B",
                "residue_id": rid,
                "residue_num": rnum,
                "residue_name": rname,
                "lobe_label": "partner",
                "source": "lightdock",
            })

        model_rows.append({
            "model_id": model_id,
            "receptor_id": state_name,
            "swarm_id": swarm_id,
            "pose_rank": pose_rank,
            "scoring_value": score,
            "n_receptor_interface_residues": len(receptor_res),
            "n_partner_interface_residues": len(partner_res),
            "n_nlobe_interface_residues": n_nlobe,
            "n_clobe_interface_residues": n_clobe,
            "receptor_interface_residues": ";".join(r[0] for r in receptor_res),
            "partner_interface_residues": ";".join(r[0] for r in partner_res),
            "source": "lightdock",
        })

    if not interface_rows:
        print(f"  No interface data extracted")
        return None, None

    # Write outputs
    iface_path = state_dir / "lightdock_interface_support_table.csv"
    model_path = state_dir / "lightdock_model_summary.csv"

    _write_csv(iface_path, interface_rows, LIGHTDOCK_INTERFACE_COLUMNS)
    _write_csv(model_path, model_rows, LIGHTDOCK_MODEL_SUMMARY_COLUMNS)

    print(f"  Extracted {len(model_rows)} models, {len(interface_rows)} residue entries")
    print(f"    {iface_path.name}")
    print(f"    {model_path.name}")

    return iface_path, model_path


def _parse_lightdock_ranking(rank_file: Path) -> List[Tuple[str, float, int]]:
    """Parse LightDock ranking file.

    Returns [(pdb_path, score, swarm_id), ...] sorted by score (descending).
    """
    results = []
    with open(rank_file, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            pdb_path = parts[0]
            try:
                score = float(parts[1])
            except ValueError:
                continue

            # Extract swarm ID from path (e.g., "swarm_42/lightdock_3.pdb")
            swarm_match = re.search(r"swarm_(\d+)", pdb_path)
            swarm_id = int(swarm_match.group(1)) if swarm_match else 0

            results.append((pdb_path, score, swarm_id))

    # Sort by score descending (higher = better in LightDock)
    results.sort(key=lambda x: -x[1])
    return results


def _extract_contacts_from_pdb(
    pdb_path: Path,
    receptor_chain: str,
    partner_chain: str,
    cutoff: float,
) -> Tuple[List[Tuple[str, int, str]], List[Tuple[str, int, str]]]:
    """Extract interface residues by CA-CA distance from a PDB file.

    Returns (receptor_residues, partner_residues) as lists of
    (residue_id, residue_num, residue_name).
    """
    # Parse CA atoms
    receptor_cas = []  # [(resnum, resname, x, y, z)]
    partner_cas = []

    with open(pdb_path, encoding="utf-8") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            chain = line[21]
            resname = line[17:20].strip()
            try:
                resnum = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            if chain == receptor_chain:
                receptor_cas.append((resnum, resname, x, y, z))
            elif chain == partner_chain:
                partner_cas.append((resnum, resname, x, y, z))

    # Find contacts
    receptor_contact_set = set()
    partner_contact_set = set()
    cutoff_sq = cutoff * cutoff

    for r_num, r_name, rx, ry, rz in receptor_cas:
        for p_num, p_name, px, py, pz in partner_cas:
            dx = rx - px
            dy = ry - py
            dz = rz - pz
            dist_sq = dx*dx + dy*dy + dz*dz
            if dist_sq < cutoff_sq:
                r_name_norm = _normalize_resname(r_name)
                p_name_norm = _normalize_resname(p_name)
                receptor_contact_set.add((f"{r_name_norm}{r_num}", r_num, r_name_norm))
                partner_contact_set.add((f"{p_name_norm}{p_num}", p_num, p_name_norm))

    receptor_res = sorted(receptor_contact_set, key=lambda x: x[1])
    partner_res = sorted(partner_contact_set, key=lambda x: x[1])

    return receptor_res, partner_res


_RESNAME_MAP = {"HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "CYX": "CYS",
                "HIE": "HIS", "HID": "HIS", "HIP": "HIS"}

def _normalize_resname(name: str) -> str:
    return _RESNAME_MAP.get(name, name)


# ---------------------------------------------------------------------------
# 1.4.3 Cross-method convergence analysis
# ---------------------------------------------------------------------------

def compute_cross_method_convergence(
    state_name: str,
    output_base: Path,
) -> Optional[Path]:
    """Compare receptor-side interface residues between PyRosetta and LightDock.

    Returns path to cross_method_convergence.csv, or None.
    """
    state_dir = output_base / state_name

    # Load PyRosetta patch table (from TG 1.3)
    pyrosetta_patches = _load_pyrosetta_patches(state_dir)

    # Load LightDock interface data (from TG 1.4.2)
    lightdock_dir = state_dir / "lightdock"
    lightdock_residues = _load_lightdock_residues(lightdock_dir)

    if not pyrosetta_patches and not lightdock_residues:
        print(f"  No data for convergence analysis ({state_name})")
        return None

    # Build unified residue set
    all_residues = set()  # (chain, residue_id)
    all_residues.update(pyrosetta_patches.keys())
    all_residues.update(lightdock_residues.keys())

    convergence_rows = []
    for (chain, rid) in sorted(all_residues):
        pyro = pyrosetta_patches.get((chain, rid))
        light = lightdock_residues.get((chain, rid))

        in_pyrosetta = pyro is not None
        in_lightdock = light is not None

        if in_pyrosetta and in_lightdock:
            convergence_class = "convergent"
            method_agreement = "both"
        elif in_pyrosetta:
            convergence_class = "pyrosetta_only"
            method_agreement = "single"
        else:
            convergence_class = "lightdock_only"
            method_agreement = "single"

        # Get metadata from whichever source has it
        meta = pyro or light or {}

        convergence_rows.append({
            "receptor_id": state_name,
            "chain": chain,
            "residue_id": rid,
            "residue_num": meta.get("residue_num", ""),
            "residue_name": meta.get("residue_name", ""),
            "lobe_label": meta.get("lobe_label", ""),
            "in_pyrosetta": in_pyrosetta,
            "in_lightdock": in_lightdock,
            "pyrosetta_max_occupancy": pyro.get("max_occupancy", "") if pyro else "",
            "lightdock_frequency": light.get("frequency", "") if light else "",
            "convergence_class": convergence_class,
            "method_agreement": method_agreement,
        })

    if not convergence_rows:
        return None

    # Sort: convergent first, then by chain
    rank = {"convergent": 0, "pyrosetta_only": 1, "lightdock_only": 2}
    convergence_rows.sort(key=lambda r: (
        0 if r["chain"] == "A" else 1,
        rank.get(r["convergence_class"], 3),
    ))

    out_path = state_dir / "cross_method_convergence.csv"
    _write_csv(out_path, convergence_rows, CONVERGENCE_COLUMNS)

    # Summary
    receptor_rows = [r for r in convergence_rows if r["chain"] == "A"]
    n_convergent = sum(1 for r in receptor_rows if r["convergence_class"] == "convergent")
    n_pyro_only = sum(1 for r in receptor_rows if r["convergence_class"] == "pyrosetta_only")
    n_light_only = sum(1 for r in receptor_rows if r["convergence_class"] == "lightdock_only")
    n_total = len(receptor_rows)

    jaccard = (n_convergent / n_total) if n_total > 0 else 0

    print(f"\n  Cross-method convergence for {state_name} (receptor-side):")
    print(f"    Convergent:       {n_convergent}")
    print(f"    PyRosetta-only:   {n_pyro_only}")
    print(f"    LightDock-only:   {n_light_only}")
    print(f"    Jaccard overlap:  {jaccard:.3f}")
    print(f"    Output: {out_path}")

    return out_path


def _load_pyrosetta_patches(state_dir: Path) -> Dict[Tuple[str, str], dict]:
    """Load PyRosetta patch table keyed by (chain, residue_id)."""
    path = state_dir / "ppi_interface_patch_table.csv"
    if not path.exists():
        return {}
    result = {}
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            key = (row.get("chain", ""), row.get("residue_id", ""))
            result[key] = {
                "residue_num": row.get("residue_num", ""),
                "residue_name": row.get("residue_name", ""),
                "lobe_label": row.get("lobe_label", ""),
                "max_occupancy": row.get("max_occupancy", ""),
            }
    return result


def _load_lightdock_residues(lightdock_dir: Path) -> Dict[Tuple[str, str], dict]:
    """Load LightDock interface residues and compute per-residue frequency."""
    path = lightdock_dir / "lightdock_interface_support_table.csv"
    if not path.exists():
        return {}

    # Count models and per-residue occurrence
    model_ids = set()
    residue_counts = Counter()
    residue_meta = {}

    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            mid = row.get("model_id", "")
            key = (row.get("chain", ""), row.get("residue_id", ""))
            model_ids.add(mid)
            residue_counts[key] += 1
            if key not in residue_meta:
                residue_meta[key] = {
                    "residue_num": row.get("residue_num", ""),
                    "residue_name": row.get("residue_name", ""),
                    "lobe_label": row.get("lobe_label", ""),
                }

    n_models = len(model_ids)
    result = {}
    for key, count in residue_counts.items():
        meta = residue_meta.get(key, {})
        freq = count / n_models if n_models > 0 else 0
        result[key] = {
            "residue_num": meta.get("residue_num", ""),
            "residue_name": meta.get("residue_name", ""),
            "lobe_label": meta.get("lobe_label", ""),
            "frequency": round(freq, 4),
        }

    return result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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
        description="Phase 1 TG 1.4: LightDock secondary validation"
    )
    parser.add_argument("--state", choices=RECEPTOR_STATES,
                        help="Process a single receptor state")
    parser.add_argument("--output_dir", type=Path,
                        default=PHASE1_OUTPUT_DIR,
                        help="Phase 1 output base directory")
    parser.add_argument("--setup", action="store_true",
                        help="Generate LightDock setup and run script")
    parser.add_argument("--extract", action="store_true",
                        help="Extract interface residues from LightDock results")
    parser.add_argument("--convergence", action="store_true",
                        help="Run cross-method convergence analysis")
    parser.add_argument("--all", action="store_true",
                        help="Run all steps (setup + extract + convergence)")
    parser.add_argument("--contact_cutoff", type=float, default=10.0,
                        help="CA-CA contact cutoff for interface extraction (default: 10.0 Å)")
    args = parser.parse_args()

    states = [args.state] if args.state else RECEPTOR_STATES

    if args.setup or args.all:
        print("Phase 1 — Task 1.4.1: LightDock Setup")
        for state in states:
            print(f"\n{state}:")
            generate_lightdock_setup(state, args.output_dir)

    if args.extract or args.all:
        print("\nPhase 1 — Task 1.4.2: LightDock Interface Extraction")
        for state in states:
            print(f"\n{state}:")
            extract_lightdock_interfaces(state, args.output_dir, args.contact_cutoff)

    if args.convergence or args.all:
        print("\nPhase 1 — Task 1.4.3: Cross-Method Convergence Analysis")
        for state in states:
            print(f"\n{state}:")
            compute_cross_method_convergence(state, args.output_dir)

    if not (args.setup or args.extract or args.convergence or args.all):
        print("Specify --setup, --extract, --convergence, or --all")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
