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
    # Setup (generates run script for all states)
    python -m egfr_pipeline.phase1.lightdock_validation --setup

    # Execute LightDock (requires lightdock3 on PATH)
    python -m egfr_pipeline.phase1.lightdock_validation --run --state 3GT8_raw

    # Extract interfaces + orientation filter (after execution)
    python -m egfr_pipeline.phase1.lightdock_validation --extract --state 3GT8_raw

    # Cross-method convergence
    python -m egfr_pipeline.phase1.lightdock_validation --convergence --state 3GT8_raw

    # All non-execution steps (setup + extract + convergence + note)
    python -m egfr_pipeline.phase1.lightdock_validation --all

    # Check LightDock availability
    python -m egfr_pipeline.phase1.lightdock_validation --check
"""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False

# Import orientation constants (no PyRosetta dependency at import time).
# These are plain Python lists/ints/floats — no heavy dependencies.
try:
    from egfr_pipeline.phase1.orientation_filter import (
        ACTIVE_FACE_RESIDUES,
        BACK_FACE_RESIDUES,
        PROBE_RESIDUES,
        PROBE_RESIDUE_FALLBACK,
        AMBIGUOUS_BAND,
    )
except ImportError:
    # Inline fallback so module remains importable even without orientation_filter
    ACTIVE_FACE_RESIDUES = [961, 962, 963, 964, 968, 969, 970, 971, 972]
    BACK_FACE_RESIDUES = [977, 978, 979, 980, 985, 986, 987, 988, 993, 994, 995, 996, 997]
    PROBE_RESIDUES = [962, 964, 971]
    PROBE_RESIDUE_FALLBACK = 962
    AMBIGUOUS_BAND = 0.15

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHASE1_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase1_ppi"
PHASE1_INPUT_DIR = PROJECT_ROOT / "input" / "PPI" / "phase1"
RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]

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
    "construct_type",
    "orientation_validation_status",
    "orientation_score",
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
    "construct_type",
    "orientation_validation_status",
    "orientation_score",
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
    "construct_type",
    "orientation_validation_status",
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
    "method_agreement",     # both / pyrosetta_only / lightdock_only
]

NLOBE_CLOBE_BOUNDARY = 838

N_CPUS_DEFAULT = 16

LIGHTDOCK_EXECUTABLES = [
    "lightdock3_setup.py",
    "lightdock3.py",
    "lgd_generate_conformations.py",
    "lgd_rank.py",
    "lgd_cluster_bsas.py",
]


def _bash_path(path: Path) -> str:
    resolved = Path(path).resolve()
    posix = resolved.as_posix()
    if len(posix) >= 2 and posix[1] == ":":
        return f"/{posix[0].lower()}{posix[2:]}"
    return posix


# ---------------------------------------------------------------------------
# LightDock availability check
# ---------------------------------------------------------------------------

def check_lightdock_available() -> Tuple[bool, str]:
    """Check if LightDock command-line tools are available on PATH."""
    missing = [exe for exe in LIGHTDOCK_EXECUTABLES
               if shutil.which(exe) is None]
    if not missing:
        return True, "All LightDock executables found on PATH"
    return False, (
        f"Missing LightDock executables: {', '.join(missing)}. "
        f"Install via: pip install lightdock3"
    )


# ---------------------------------------------------------------------------
# PDB atom parsing (PyRosetta-free)
# ---------------------------------------------------------------------------

def _parse_pdb_atoms(
    pdb_path: Path,
    atom_names: Tuple[str, ...] = ("CA", "CB"),
) -> Dict[Tuple[str, int, str], Tuple[float, float, float, str]]:
    """Parse atom coordinates from a PDB file.

    Returns dict keyed by (chain, resnum, atom_name) -> (x, y, z, resname).
    """
    atoms = {}
    with open(pdb_path, encoding="utf-8") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name not in atom_names:
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
            atoms[(chain, resnum, atom_name)] = (x, y, z, resname)
    return atoms


def compute_orientation_score_from_pdb(
    pdb_path: Path,
    partner_chain: str = "B",
    receptor_chain: str = "A",
    contact_dist: float = 10.0,
) -> Dict:
    """Compute orientation score from a PDB file (no PyRosetta dependency).

    Mirrors orientation_filter.compute_orientation_score() algorithm using
    pure PDB parsing + numpy.
    """
    if not _HAS_NUMPY:
        return {
            "orientation_score": 0.0,
            "orientation_class": "not_available",
            "n_active_face_ca": 0,
            "error": "numpy not available",
        }

    try:
        atoms = _parse_pdb_atoms(pdb_path, ("CA", "CB"))
    except Exception as e:
        return {
            "orientation_score": 0.0,
            "orientation_class": "parse_error",
            "n_active_face_ca": 0,
            "error": str(e),
        }

    # 1. Collect active-face CA coordinates
    active_cas = []
    for rn in ACTIVE_FACE_RESIDUES:
        key = (partner_chain, rn, "CA")
        if key in atoms:
            active_cas.append(list(atoms[key][:3]))

    n_active = len(active_cas)
    if n_active < 3:
        return {
            "orientation_score": 0.0,
            "orientation_class": "insufficient_data",
            "n_active_face_ca": n_active,
            "n_back_face_ca": 0,
            "n_receptor_contact_ca": 0,
            "error": f"Only {n_active} active-face CA found (need >=3)",
        }

    active_cas = np.array(active_cas)

    # 2. Sheet centroid
    sheet_centroid = active_cas.mean(axis=0)

    # 3. Sheet-plane normal via PCA
    centered = active_cas - sheet_centroid
    _, _s, vt = np.linalg.svd(centered, full_matrices=False)
    normal = vt[-1]
    normal = normal / np.linalg.norm(normal)

    # 4. Orient normal toward active face using multi-probe CA->CB consensus
    votes = 0
    n_valid = 0
    for rn in PROBE_RESIDUES:
        ca_key = (partner_chain, rn, "CA")
        cb_key = (partner_chain, rn, "CB")
        if cb_key not in atoms:
            cb_key = ca_key  # GLY fallback
        if ca_key not in atoms or cb_key not in atoms:
            continue
        ca = np.array(atoms[ca_key][:3])
        cb = np.array(atoms[cb_key][:3])
        ca_cb = cb - ca
        norm_cacb = np.linalg.norm(ca_cb)
        if norm_cacb < 1e-6:
            continue
        ca_cb = ca_cb / norm_cacb
        dot = np.dot(normal, ca_cb)
        votes += 1 if dot > 0 else -1
        n_valid += 1

    if n_valid == 0:
        # Single-probe fallback
        ca_key = (partner_chain, PROBE_RESIDUE_FALLBACK, "CA")
        cb_key = (partner_chain, PROBE_RESIDUE_FALLBACK, "CB")
        if cb_key not in atoms:
            cb_key = ca_key
        if ca_key not in atoms:
            return {
                "orientation_score": 0.0,
                "orientation_class": "probe_error",
                "n_active_face_ca": n_active,
                "error": "No valid probe residues found",
            }
        ca = np.array(atoms[ca_key][:3])
        cb = np.array(atoms[cb_key][:3])
        ca_cb = cb - ca
        norm_cacb = np.linalg.norm(ca_cb)
        if norm_cacb < 1e-6:
            return {
                "orientation_score": 0.0,
                "orientation_class": "probe_error",
                "n_active_face_ca": n_active,
                "error": "Degenerate probe CA-CB vector",
            }
        ca_cb = ca_cb / norm_cacb
        if np.dot(normal, ca_cb) < 0:
            normal = -normal
    elif votes < 0:
        normal = -normal

    # 5. Receptor-direction vector
    receptor_contact_cas = []
    cutoff_sq = contact_dist * contact_dist
    for (ch, rn, aname), (x, y, z, _rn) in atoms.items():
        if ch != receptor_chain or aname != "CA":
            continue
        r_ca = np.array([x, y, z])
        dists_sq = np.sum((active_cas - r_ca) ** 2, axis=1)
        if dists_sq.min() < cutoff_sq:
            receptor_contact_cas.append(r_ca)

    if not receptor_contact_cas:
        # Fallback: all receptor CA atoms
        for (ch, rn, aname), (x, y, z, _rn) in atoms.items():
            if ch == receptor_chain and aname == "CA":
                receptor_contact_cas.append(np.array([x, y, z]))

    if not receptor_contact_cas:
        return {
            "orientation_score": 0.0,
            "orientation_class": "receptor_error",
            "n_active_face_ca": n_active,
            "error": "No receptor CA atoms found",
        }

    receptor_centroid = np.mean(receptor_contact_cas, axis=0)
    receptor_dir = receptor_centroid - sheet_centroid
    receptor_dir_norm = np.linalg.norm(receptor_dir)
    if receptor_dir_norm < 1e-6:
        return {
            "orientation_score": 0.0,
            "orientation_class": "degenerate_geometry",
            "n_active_face_ca": n_active,
            "error": "Sheet centroid coincides with receptor centroid",
        }
    receptor_dir = receptor_dir / receptor_dir_norm

    # 6. Dot product → classify
    dot = float(np.dot(normal, receptor_dir))
    if abs(dot) < AMBIGUOUS_BAND:
        orientation_class = "ambiguous"
    elif dot > 0:
        orientation_class = "pass"
    else:
        orientation_class = "fail"

    n_back = sum(1 for rn in BACK_FACE_RESIDUES
                 if (partner_chain, rn, "CA") in atoms)

    return {
        "orientation_score": round(dot, 4),
        "orientation_class": orientation_class,
        "n_active_face_ca": n_active,
        "n_back_face_ca": n_back,
        "n_receptor_contact_ca": len(receptor_contact_cas),
    }


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
    merged_params = dict(LIGHTDOCK_DEFAULTS)
    if params is not None:
        merged_params.update(params)
    params = merged_params

    state_dir = (output_base / state_name / "lightdock").resolve()
    state_dir.mkdir(parents=True, exist_ok=True)

    # Input PDB (resolve to absolute)
    input_pdb = (PHASE1_INPUT_DIR / f"docking_{state_name}_ext_beta_meander.pdb").resolve()

    n_cpus = params.get("n_cpus", N_CPUS_DEFAULT)
    n_swarms = params["n_swarms"]
    n_steps = params["n_steps"]
    n_top = params["n_top_poses"]
    scoring = params["scoring_function"]

    # Metadata
    metadata = {
        "receptor_id": state_name,
        "partner_id": "extended_beta_meander_955_1006",
        "construct_type": "full_kinase_domain",
        "method": "LightDock",
        "scoring_function": scoring,
        "n_swarms": n_swarms,
        "n_glowworms": params["n_glowworms"],
        "n_steps": n_steps,
        "n_top_poses": n_top,
        "n_cpus": n_cpus,
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

    # Generate run script with corrected LightDock workflow
    script_path = state_dir / f"run_lightdock_{state_name}.sh"
    input_pdb_script = _bash_path(input_pdb)
    workdir_script = _bash_path(state_dir)
    script_path_script = _bash_path(script_path)
    script_content = f"""#!/bin/bash
# LightDock run script for {state_name}
# Generated by egfr_pipeline.phase1.lightdock_validation
#
# Prerequisites:
#   pip install lightdock3
#   or: conda install -c bioconda lightdock
#
# Usage:
#   bash {script_path_script}

set -e

INPUT_PDB="{input_pdb_script}"
WORKDIR="{workdir_script}"
N_SWARMS={n_swarms}
N_STEPS={n_steps}
N_TOP={n_top}
N_CPUS={n_cpus}

cd "$WORKDIR"

# Pre-flight: verify LightDock is installed
for cmd in lightdock3_setup.py lightdock3.py lgd_generate_conformations.py lgd_rank.py lgd_cluster_bsas.py; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: $cmd not found on PATH. Install lightdock3 first."
        exit 1
    fi
done

echo "=== LightDock Secondary Validation: {state_name} ==="
echo "Input: $INPUT_PDB"
echo "Scoring: {scoring}"
echo "Swarms: $N_SWARMS, Glowworms: {params['n_glowworms']}, Steps: $N_STEPS"
echo "CPUs: $N_CPUS"
echo "Started: $(date)"

# Step 1: Split PDB into receptor (chain A) and partner (chain B)
echo ""
echo "[Step 1/5] Splitting input PDB..."
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
        elif line.startswith('TER'):
            chain = line[21] if len(line) > 21 else ''
            if chain == 'A':
                r.write(line)
            elif chain == 'B':
                p.write(line)
    r.write('END\\n')
    p.write('END\\n')
"

# Step 2: LightDock setup
echo "[Step 2/5] Running lightdock3_setup..."
lightdock3_setup.py receptor.pdb partner.pdb \\
    -s $N_SWARMS \\
    -g {params['n_glowworms']} \\
    --noxt --noh

# Step 3: LightDock run
echo "[Step 3/5] Running lightdock3 ($N_CPUS cores)..."
lightdock3.py setup.json $N_STEPS \\
    -s {scoring} \\
    -c $N_CPUS

# Step 4: Generate top poses per swarm
echo "[Step 4/5] Generating top poses per swarm..."
for swarm_dir in swarm_*/; do
    cd "$WORKDIR/$swarm_dir"
    lgd_generate_conformations.py ../receptor.pdb ../partner.pdb \\
        gso_${{N_STEPS}}.out $N_TOP
    cd "$WORKDIR"
done

# Step 5: Rank and cluster
echo "[Step 5/5] Ranking and clustering..."
lgd_rank.py $N_SWARMS $N_STEPS
lgd_cluster_bsas.py setup.json

# Completion marker
touch "$WORKDIR/.lightdock_complete"

echo ""
echo "=== LightDock complete for {state_name} ==="
echo "Finished: $(date)"
echo ""
echo "Next step - extract interfaces:"
echo "  python -m egfr_pipeline.phase1.lightdock_validation --extract --state {state_name}"
"""

    with open(script_path, "w", encoding="utf-8") as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)

    print(f"  LightDock setup for {state_name}:")
    print(f"    Metadata: {meta_path}")
    print(f"    Run script: {script_path}")
    print(f"    Swarms: {n_swarms}, Steps: {n_steps}, CPUs: {n_cpus}")
    print(f"    Scoring: {scoring}")

    return script_path


def run_lightdock_scripts(
    state_name: str,
    output_base: Path,
) -> bool:
    """Execute generated LightDock run script via subprocess.

    Returns True if execution succeeded.
    """
    state_dir = (output_base / state_name / "lightdock").resolve()
    script_path = state_dir / f"run_lightdock_{state_name}.sh"

    if not script_path.exists():
        print(f"  Run script not found: {script_path}")
        print(f"  Run --setup first to generate it.")
        return False

    available, msg = check_lightdock_available()
    if not available:
        print(f"  {msg}")
        print(f"  Transfer the script to a server with LightDock installed:")
        print(f"    {script_path}")
        return False

    print(f"  Executing LightDock for {state_name}...")
    print(f"  Script: {script_path}")

    result = subprocess.run(
        ["bash", str(script_path)],
        cwd=str(state_dir),
    )

    if result.returncode != 0:
        print(f"  LightDock execution failed (exit code {result.returncode})")
        return False

    print(f"  LightDock execution completed for {state_name}")
    return True


# ---------------------------------------------------------------------------
# 1.4.2 Interface extraction from LightDock output
# ---------------------------------------------------------------------------

def extract_lightdock_interfaces(
    state_name: str,
    output_base: Path,
    contact_cutoff: float = 10.0,
    max_poses: int = 200,
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

    # Check for completion marker
    if not (state_dir / ".lightdock_complete").exists():
        rank_check = (state_dir / "rank_by_scoring.list").exists()
        if not rank_check:
            print(f"  LightDock has not been run yet for {state_name}.")
            print(f"  Execute the generated script first, or use --run.")
            return None, None

    # Find ranked poses
    rank_file = state_dir / "rank_by_scoring.list"
    if not rank_file.exists():
        rank_file = state_dir / "rank_by_cluster.list"
    if not rank_file.exists():
        print(f"  No ranking file found in {state_dir}")
        return None, None

    ranked_poses = _parse_lightdock_ranking(rank_file, max_poses=max_poses)
    if not ranked_poses:
        print(f"  No ranked poses found in {rank_file}")
        return None, None

    lightdock_metadata = _load_lightdock_metadata(state_dir)
    construct_type = lightdock_metadata.get("construct_type", "")

    print(f"  Found {len(ranked_poses)} ranked poses (max_poses={max_poses})")

    interface_rows = []
    model_rows = []
    n_orient_pass = n_orient_fail = n_orient_other = 0

    for pose_rank, (pdb_path, score, swarm_id) in enumerate(ranked_poses, 1):
        full_path = state_dir / pdb_path
        if not full_path.exists():
            continue

        model_id = f"swarm{swarm_id}_rank{pose_rank}"

        # Orientation scoring (PDB-based, no PyRosetta)
        orient = compute_orientation_score_from_pdb(
            full_path, partner_chain="B", receptor_chain="A",
            contact_dist=contact_cutoff,
        )
        orientation_class = orient.get("orientation_class", "error")
        orientation_score = orient.get("orientation_score", 0.0)

        if orientation_class == "pass":
            n_orient_pass += 1
        elif orientation_class == "fail":
            n_orient_fail += 1
        else:
            n_orient_other += 1

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
                "construct_type": construct_type,
                "orientation_validation_status": orientation_class,
                "orientation_score": orientation_score,
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
                "construct_type": construct_type,
                "orientation_validation_status": orientation_class,
                "orientation_score": orientation_score,
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
            "construct_type": construct_type,
            "orientation_validation_status": orientation_class,
            "orientation_score": orientation_score,
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
    print(f"  Orientation: P:{n_orient_pass} F:{n_orient_fail} Other:{n_orient_other}")
    print(f"    {iface_path.name}")
    print(f"    {model_path.name}")

    return iface_path, model_path


def _parse_lightdock_ranking(
    rank_file: Path,
    max_poses: int = 200,
) -> List[Tuple[str, float, int]]:
    """Parse LightDock ranking file.

    Returns [(pdb_path, score, swarm_id), ...] sorted by score (descending),
    limited to top max_poses entries.
    """
    results = []
    with open(rank_file, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("PDB"):
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
    if max_poses > 0:
        results = results[:max_poses]
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

    if not receptor_contact_set and receptor_cas and partner_cas:
        nearest = None
        nearest_dist_sq = None
        for r_num, r_name, rx, ry, rz in receptor_cas:
            for p_num, p_name, px, py, pz in partner_cas:
                dx = rx - px
                dy = ry - py
                dz = rz - pz
                dist_sq = dx * dx + dy * dy + dz * dz
                if nearest_dist_sq is None or dist_sq < nearest_dist_sq:
                    nearest_dist_sq = dist_sq
                    nearest = ((r_num, r_name), (p_num, p_name))
        if nearest is not None:
            (r_num, r_name), (p_num, p_name) = nearest
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
    lightdock_metadata = _load_lightdock_metadata(lightdock_dir)

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
            method_agreement = "pyrosetta_only"
        else:
            convergence_class = "lightdock_only"
            method_agreement = "lightdock_only"

        # Get metadata from whichever source has it
        meta = pyro or light or {}
        construct_type = (
            meta.get("construct_type", "")
            or lightdock_metadata.get("construct_type", "")
        )
        orientation_validation_status = meta.get(
            "orientation_validation_status", ""
        ) or "not_available"

        convergence_rows.append({
            "receptor_id": state_name,
            "construct_type": construct_type,
            "orientation_validation_status": orientation_validation_status,
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

    # Per-lobe Jaccard
    nlobe_rows = [r for r in receptor_rows if r.get("lobe_label") == "N-lobe"]
    clobe_rows = [r for r in receptor_rows if r.get("lobe_label") == "C-lobe"]
    n_nlobe_conv = sum(1 for r in nlobe_rows if r["convergence_class"] == "convergent")
    n_clobe_conv = sum(1 for r in clobe_rows if r["convergence_class"] == "convergent")
    jaccard_nlobe = (n_nlobe_conv / len(nlobe_rows)) if nlobe_rows else 0
    jaccard_clobe = (n_clobe_conv / len(clobe_rows)) if clobe_rows else 0

    # Write summary JSON
    summary = {
        "receptor_id": state_name,
        "n_convergent": n_convergent,
        "n_pyrosetta_only": n_pyro_only,
        "n_lightdock_only": n_light_only,
        "n_total": n_total,
        "jaccard_overlap": round(jaccard, 4),
        "jaccard_nlobe": round(jaccard_nlobe, 4),
        "jaccard_clobe": round(jaccard_clobe, 4),
    }
    summary_path = state_dir / "cross_method_convergence_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"\n  Cross-method convergence for {state_name} (receptor-side):")
    print(f"    Convergent:       {n_convergent}")
    print(f"    PyRosetta-only:   {n_pyro_only}")
    print(f"    LightDock-only:   {n_light_only}")
    print(f"    Jaccard overlap:  {jaccard:.3f} (N-lobe: {jaccard_nlobe:.3f}, C-lobe: {jaccard_clobe:.3f})")
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
                "construct_type": row.get("construct_type", ""),
                "orientation_validation_status": row.get(
                    "orientation_validation_status", ""
                ),
            }
    return result


def _load_lightdock_metadata(lightdock_dir: Path) -> Dict[str, str]:
    """Load LightDock run metadata when available."""
    path = lightdock_dir / "lightdock_run_metadata.json"
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return {
        "construct_type": str(data.get("construct_type", "") or ""),
    }


def _load_lightdock_residues(
    lightdock_dir: Path,
    orientation_filter: bool = True,
) -> Dict[Tuple[str, str], dict]:
    """Load LightDock interface residues and compute per-residue frequency.

    When orientation_filter is True, only residues from orientation-validated
    models (orientation_validation_status == 'pass') are included, mirroring
    how PyRosetta cluster_consensus only uses orientation-validated models.
    """
    path = lightdock_dir / "lightdock_interface_support_table.csv"
    if not path.exists():
        return {}

    # Count models and per-residue occurrence
    model_ids = set()
    residue_counts = Counter()
    residue_meta = {}

    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            # Orientation filter: skip non-pass models
            if orientation_filter:
                ovs = row.get("orientation_validation_status", "")
                if ovs and ovs != "pass" and ovs != "not_available":
                    continue

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
# Validation note generation
# ---------------------------------------------------------------------------

def generate_lightdock_note(output_base: Path) -> Optional[Path]:
    """Generate a dynamic LightDock validation note from actual results.

    Reads per-state metadata and convergence summaries, writes a markdown
    note to output_base/phase1_lightdock_validation_results.md.
    """
    lines = [
        "# LightDock Validation Results",
        "",
        "Auto-generated summary of LightDock secondary validation execution.",
        "",
    ]

    any_data = False
    for state in RECEPTOR_STATES:
        state_dir = output_base / state
        ld_dir = state_dir / "lightdock"
        meta = _load_lightdock_metadata(ld_dir)
        if not meta:
            continue

        # Load metadata JSON for run params
        meta_path = ld_dir / "lightdock_run_metadata.json"
        run_params = {}
        if meta_path.exists():
            with open(meta_path, encoding="utf-8") as f:
                run_params = json.load(f)

        lines.extend([f"## {state}", ""])

        if run_params:
            lines.extend([
                f"- Scoring: {run_params.get('scoring_function', 'N/A')}",
                f"- Swarms: {run_params.get('n_swarms', 'N/A')}",
                f"- Steps: {run_params.get('n_steps', 'N/A')}",
                f"- Contact cutoff: {run_params.get('contact_cutoff_A', 'N/A')} A",
                "",
            ])

        # Model summary
        model_path = ld_dir / "lightdock_model_summary.csv"
        if model_path.exists():
            any_data = True
            models = []
            with open(model_path, encoding="utf-8") as f:
                models = list(csv.DictReader(f))
            n_models = len(models)
            n_pass = sum(1 for m in models
                         if m.get("orientation_validation_status") == "pass")
            n_fail = sum(1 for m in models
                         if m.get("orientation_validation_status") == "fail")
            lines.extend([
                f"- Models extracted: {n_models}",
                f"- Orientation: PASS={n_pass}, FAIL={n_fail}, "
                f"Other={n_models - n_pass - n_fail}",
                "",
            ])

        # Convergence summary
        summary_path = state_dir / "cross_method_convergence_summary.json"
        if summary_path.exists():
            with open(summary_path, encoding="utf-8") as f:
                summary = json.load(f)
            lines.extend([
                "### Cross-method convergence",
                "",
                f"- Convergent: {summary.get('n_convergent', 0)}",
                f"- PyRosetta-only: {summary.get('n_pyrosetta_only', 0)}",
                f"- LightDock-only: {summary.get('n_lightdock_only', 0)}",
                f"- Jaccard: {summary.get('jaccard_overlap', 0):.3f}",
                f"  - N-lobe: {summary.get('jaccard_nlobe', 0):.3f}",
                f"  - C-lobe: {summary.get('jaccard_clobe', 0):.3f}",
                "",
            ])

    if not any_data:
        lines.append("No LightDock results available yet.")
        lines.append("")

    lines.extend([
        "## Score comparison note",
        "",
        "PyRosetta uses REU (dG_separated), LightDock uses DFIRE2 scoring.",
        "These scoring functions are not directly comparable.",
        "Cross-method validation relies on residue-level Jaccard overlap,",
        "which is unit-agnostic.",
        "",
        "---",
        "",
        "Generated by `egfr_pipeline.phase1.lightdock_validation`",
    ])

    note_path = output_base / "phase1_lightdock_validation_results.md"
    note_path.parent.mkdir(parents=True, exist_ok=True)
    with open(note_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"  Validation note: {note_path}")
    return note_path


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
    parser.add_argument("--run", action="store_true",
                        help="Execute generated LightDock run scripts (requires LightDock on PATH)")
    parser.add_argument("--extract", action="store_true",
                        help="Extract interface residues from LightDock results")
    parser.add_argument("--convergence", action="store_true",
                        help="Run cross-method convergence analysis")
    parser.add_argument("--note", action="store_true",
                        help="Generate validation results note")
    parser.add_argument("--all", action="store_true",
                        help="Run all steps (setup + extract + convergence + note)")
    parser.add_argument("--contact_cutoff", type=float, default=10.0,
                        help="CA-CA contact cutoff for interface extraction (default: 10.0 A)")
    parser.add_argument("--max_poses", type=int, default=200,
                        help="Maximum poses to extract from ranking (default: 200)")
    parser.add_argument("--n-cpus", type=int, default=N_CPUS_DEFAULT,
                        help=f"CPU count for generated run scripts (default: {N_CPUS_DEFAULT})")
    parser.add_argument("--check", action="store_true",
                        help="Check LightDock availability on PATH")
    args = parser.parse_args()

    if args.check:
        available, msg = check_lightdock_available()
        print(msg)
        return 0 if available else 1

    states = [args.state] if args.state else RECEPTOR_STATES

    if args.setup or args.all:
        print("Phase 1 — Task 1.4.1: LightDock Setup")
        for state in states:
            print(f"\n{state}:")
            generate_lightdock_setup(
                state,
                args.output_dir,
                params={"n_cpus": args.n_cpus},
            )

    if args.run:
        print("\nPhase 1 — Task 1.4.1: LightDock Execution")
        for state in states:
            print(f"\n{state}:")
            run_lightdock_scripts(state, args.output_dir)

    if args.extract or args.all:
        print("\nPhase 1 — Task 1.4.2: LightDock Interface Extraction")
        for state in states:
            print(f"\n{state}:")
            extract_lightdock_interfaces(
                state, args.output_dir, args.contact_cutoff,
                max_poses=args.max_poses,
            )

    if args.convergence or args.all:
        print("\nPhase 1 — Task 1.4.3: Cross-Method Convergence Analysis")
        for state in states:
            print(f"\n{state}:")
            compute_cross_method_convergence(state, args.output_dir)

    if args.note or args.all:
        print("\nPhase 1 — LightDock Validation Note")
        generate_lightdock_note(args.output_dir)

    has_action = any([args.setup, args.run, args.extract,
                      args.convergence, args.note, args.all])
    if not has_action:
        print("Specify --setup, --run, --extract, --convergence, --note, --all, or --check")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
