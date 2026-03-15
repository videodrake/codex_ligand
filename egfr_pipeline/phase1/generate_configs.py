#!/usr/bin/env python3
"""Phase 1 Task 1.1: Generate PyRosetta docking config files for each receptor state.

Creates properly configured INI files adapted for the Phase 1 full-kinase-domain
monomer setup (not the legacy dimer setup).

Key differences from legacy configs:
  - Input PDB: Phase 1 docking pairs (monomer receptor + extended beta-meander)
  - Excluded residues: Chain A monomer only (no +1000 offset dimer residues)
  - Output path: Explicit receptor state labeling
  - Metadata: construct_type = full_kinase_domain

Usage:
    python -m egfr_pipeline.phase1.generate_configs [--output_dir config/phase1]
"""

import argparse
import configparser
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]

# Membrane-proximal residues — monomer chain A only (no dimer offset)
EXCLUDED_RESIDUES_A = "709-720,724-731,736-739,747,783-785,799-805,871-873,917-921"

# Jura 2009: EGFR asymmetric dimer interface (C-lobe surface)
# αC helix(741-756), αE helix(847-859), activation loop(831-852)
KNOWN_BINDING_REGION_A = "741-756,831-859"

# Multi-seed strategy: 5 seeds × 20K each = 100K total per state
# This provides better sampling diversity than a single 100K run
PRODUCTION_MODELS_PER_SEED = 20000
PRODUCTION_N_SEEDS = 5
TEST_MODELS = 1000


def generate_config(
    state_name: str,
    is_production: bool,
    seed_index: int = 0,
    output_dir: Path = None,
) -> Path:
    """Generate a single INI config file for one receptor state.

    Args:
        state_name: e.g., "3GT8_raw"
        is_production: True for production (20K models), False for test (1K)
        seed_index: For multi-seed production runs (0-4)
        output_dir: Directory for config output
    """
    config = configparser.ConfigParser()

    # ---- Path ----
    config["Path"] = {
        "input_pdb_name": f"input/PPI/phase1/docking_{state_name}_ext_beta_meander.pdb",
    }

    # ---- System ----
    if is_production:
        seed_value = str(42 + seed_index * 1000)  # deterministic but distinct seeds
    else:
        seed_value = "auto"

    config["System"] = {
        "n_cpus": "16",  # Safe default for shared HPC (32 cores / 2)
        "random_seed": seed_value,
    }

    # ---- Docking ----
    n_models = PRODUCTION_MODELS_PER_SEED if is_production else TEST_MODELS
    config["Docking"] = {
        "total_global_models": str(n_models),
    }

    # ---- Filter (v1.0 fallback, kept for backward compat) ----
    config["Filter"] = {
        "filter_percentile": "15",
        "min_dSASA": "100.0",
        "contact_distance": "6.0",
        "min_sc_value": "0.3",
        "min_survivors": "50",
    }

    # ---- Cluster ----
    config["Cluster"] = {
        "threshold": "auto",
        "cluster_top_n": "auto",
        "members_per_cluster": "5",
        "member_diversity_dist": "auto",
        "site_group_jaccard": "0.3",
    }

    # ---- Refinement ----
    config["Refinement"] = {
        "refine_per_struct": "10",
        "perturb_trans": "0.1",
        "perturb_rot": "1.0",
    }

    # ---- Constraints ----
    # Key difference: monomer-only excluded residues (no dimer +1000 offset)
    config["Constraints"] = {
        "excluded_residues_A": EXCLUDED_RESIDUES_A,
        "key_residues_B": "",
        "exclusion_contact_dist": "10.0",
        "enable_early_rejection": "true",
        "max_excluded_contacts": "0",
        "key_residue_bonus_weight": "0.0",
        "key_residue_weights": "",
    }

    # ---- MiniRefinement (v2.0) ----
    config["MiniRefinement"] = {
        "enabled": "true",
        "mode": "pack_only",
        "n_rounds": "3",
    }

    # ---- FilterStage1 (v2.0) ----
    config["FilterStage1"] = {
        "total_score_percentile": "10",
        "max_dG_separated": "0.0",
        "enabled": "true",
    }

    # ---- FilterStage2 (v2.0) ----
    config["FilterStage2"] = {
        "min_dSASA": "500.0",
        "max_dG_density": "-1.0",
        "min_sc_value": "0.50",
        "min_packstat": "0.55",
        "max_delta_unsatHbonds": "8",
        "min_nres_int": "10",
        "min_hbonds_int": "0",
        "enable_expensive_metrics": "true",
        "min_survivors": "50",
    }

    # ---- Output ----
    config["Output"] = {
        "save_filter_max": "200",
        "save_top_n": "20",
    }

    # ---- ExperimentalData ----
    config["ExperimentalData"] = {
        "critical_residues_B": "",
        "non_binding_residues_B": "",
        "known_binding_region_A": KNOWN_BINDING_REGION_A,
        "confidence": "low",
    }

    # ---- Write ----
    if output_dir is None:
        output_dir = PROJECT_ROOT / "config" / "phase1"
    output_dir.mkdir(parents=True, exist_ok=True)

    if is_production:
        filename = f"phase1_prod_{state_name}_seed{seed_index}.ini"
    else:
        filename = f"phase1_test_{state_name}.ini"

    output_path = output_dir / filename

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f"# Phase 1 PyRosetta PPI Docking Config\n")
        f.write(f"# Receptor: {state_name} (full kinase domain monomer)\n")
        f.write(f"# Partner: extended beta-meander (955-1006)\n")
        f.write(f"# Construct: full_kinase_domain (NOT legacy dimer)\n")
        if is_production:
            f.write(f"# Run type: production (seed {seed_index}, {n_models} models)\n")
            f.write(f"# Multi-seed strategy: {PRODUCTION_N_SEEDS} seeds × {PRODUCTION_MODELS_PER_SEED} = "
                    f"{PRODUCTION_N_SEEDS * PRODUCTION_MODELS_PER_SEED} total\n")
        else:
            f.write(f"# Run type: test ({n_models} models)\n")
        f.write(f"#\n")
        f.write(f"# Key differences from legacy (dimer) configs:\n")
        f.write(f"#   - Receptor is monomer (chain A only, ~309 residues)\n")
        f.write(f"#   - No +1000 offset dimer residues in excluded_residues_A\n")
        f.write(f"#   - Partner is extended beta-meander (955-1006, not 960-1006)\n")
        f.write(f"#   - n_cpus = 16 (safe for shared HPC with 32 cores)\n")
        f.write(f"\n")
        config.write(f)

    return output_path


def main():
    parser = argparse.ArgumentParser(
        description="Phase 1 TG 1.1: Generate PyRosetta config files for all receptor states"
    )
    parser.add_argument("--output_dir", type=Path,
                        default=PROJECT_ROOT / "config" / "phase1",
                        help="Output directory for config files")
    parser.add_argument("--production_only", action="store_true",
                        help="Only generate production configs (skip test)")
    parser.add_argument("--test_only", action="store_true",
                        help="Only generate test configs (skip production)")
    args = parser.parse_args()

    print("Phase 1 — Task 1.1: Generating PyRosetta config files")
    configs = []

    for state in RECEPTOR_STATES:
        # Test config
        if not args.production_only:
            path = generate_config(state, is_production=False, output_dir=args.output_dir)
            configs.append(path)
            print(f"  Test:  {path.name}")

        # Production configs (multi-seed)
        if not args.test_only:
            for seed_idx in range(PRODUCTION_N_SEEDS):
                path = generate_config(state, is_production=True,
                                       seed_index=seed_idx, output_dir=args.output_dir)
                configs.append(path)
                if seed_idx == 0:
                    print(f"  Prod:  {path.name} ... (×{PRODUCTION_N_SEEDS} seeds)")

    print(f"\nGenerated {len(configs)} config files in {args.output_dir}/")
    print(f"\nMulti-seed production strategy:")
    print(f"  {PRODUCTION_N_SEEDS} seeds × {PRODUCTION_MODELS_PER_SEED} models = "
          f"{PRODUCTION_N_SEEDS * PRODUCTION_MODELS_PER_SEED} total per state")
    print(f"  {len(RECEPTOR_STATES)} states × {PRODUCTION_N_SEEDS * PRODUCTION_MODELS_PER_SEED} = "
          f"{len(RECEPTOR_STATES) * PRODUCTION_N_SEEDS * PRODUCTION_MODELS_PER_SEED} total models")

    return 0


if __name__ == "__main__":
    sys.exit(main())
