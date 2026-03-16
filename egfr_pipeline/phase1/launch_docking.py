#!/usr/bin/env python3
"""Phase 1 Task 1.1: Docking launcher with run metadata generation.

Wraps the existing PipelineManager to:
  1. Validate Phase 1 inputs before launching
  2. Generate pyrosetta_run_metadata.json with construct_type, receptor_id, etc.
  3. Standardize output directory placement by receptor state
  4. Support multi-seed production runs

This script is designed for server-side execution only.
In the Codex workspace, use --dry-run to validate configs without running.

Usage:
    # Dry run (validate configs, generate metadata, no docking)
    python -m egfr_pipeline.phase1.launch_docking --dry-run

    # Test run (single state, 1K models)
    python -m egfr_pipeline.phase1.launch_docking --test --state 3GT8_raw

    # Production run (all states, all seeds)
    python -m egfr_pipeline.phase1.launch_docking --production

    # Single production seed
    python -m egfr_pipeline.phase1.launch_docking --config config/phase1/phase1_prod_3GT8_raw_seed0.ini
"""

import argparse
import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
from egfr_pipeline.phase1.prepare_inputs import (
    PHASE1_RUNTIME_INPUT_DIR,
    prepare_phase1_inputs,
)

# Phase 1 output base directory
PHASE1_OUTPUT_DIR = PROJECT_ROOT / "output" / "phase1_ppi"

# Input validation
PHASE1_INPUT_DIR = PHASE1_RUNTIME_INPUT_DIR
PHASE1_CONFIG_DIR = PROJECT_ROOT / "config" / "phase1"

RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]


def validate_phase1_inputs() -> bool:
    """Verify all Phase 1 inputs from TG 1.0 are in place."""
    prepare_phase1_inputs(PHASE1_INPUT_DIR)
    ok = True
    checks = []

    # Check docking pair PDBs
    for state in RECEPTOR_STATES:
        pdb = PHASE1_INPUT_DIR / f"docking_{state}_ext_beta_meander.pdb"
        exists = pdb.exists()
        checks.append(("docking_pair", state, exists))
        if not exists:
            ok = False

    # Check metadata
    for fname in ["receptor_metadata.csv", "partner_metadata.csv",
                   "docking_pair_metadata.csv"]:
        fpath = PHASE1_INPUT_DIR / fname
        exists = fpath.exists()
        checks.append(("metadata", fname, exists))
        if not exists:
            ok = False

    # Check validation report
    report = PHASE1_INPUT_DIR / "phase1_input_validation_report.md"
    checks.append(("validation", "report", report.exists()))

    for category, name, status in checks:
        mark = "OK" if status else "MISSING"
        print(f"  [{mark}] {category}: {name}")

    return ok


def generate_run_metadata(
    config_path: Path,
    state_name: str,
    output_dir: Path,
    seed_index: int = 0,
    is_production: bool = False,
    dry_run: bool = False,
) -> dict:
    """Generate pyrosetta_run_metadata.json for a single run.

    This is the key Task 1.1.2 deliverable: receptor_id, partner_id,
    construct_type, and run parameters are all recorded.
    """
    import configparser
    config = configparser.ConfigParser()
    config.read(config_path)

    total_models = config.getint("Docking", "total_global_models")
    n_cpus = config.getint("System", "n_cpus", fallback=16)
    seed_raw = config.get("System", "random_seed", fallback="auto")

    metadata = {
        # Task 1.1.2: Normalized execution inputs
        "receptor_id": state_name,
        "partner_id": "extended_beta_meander_955_1006",
        "construct_type": "full_kinase_domain",
        "construct_note": "Phase 1 v2: monomer receptor (not dimer), extended partner (955 not 960/962)",

        # Run parameters
        "config_file": str(config_path.relative_to(PROJECT_ROOT)),
        "input_pdb": config.get("Path", "input_pdb_name"),
        "total_global_models": total_models,
        "n_cpus": n_cpus,
        "random_seed": seed_raw,
        "seed_index": seed_index,
        "is_production": is_production,

        # Filter version
        "filter_version": "v2.0" if config.has_section("FilterStage1") else "v1.0",
        "mini_refinement_enabled": config.getboolean("MiniRefinement", "enabled", fallback=False),

        # Constraints
        "excluded_residues_A": config.get("Constraints", "excluded_residues_A", fallback=""),
        "enable_early_rejection": config.getboolean("Constraints", "enable_early_rejection", fallback=False),

        # Phase 1 context
        "phase": "Phase 1: PPI-First Interface Mapping",
        "task_group": "TG 1.1: PyRosetta Global Docking Standardization",
        "output_dir": str(output_dir.relative_to(PROJECT_ROOT)),

        # Timestamps
        "generated_at": datetime.now().isoformat(),
        "dry_run": dry_run,
    }

    # Write metadata JSON
    output_dir.mkdir(parents=True, exist_ok=True)
    meta_path = output_dir / "pyrosetta_run_metadata.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)

    return metadata


def get_output_dir(state_name: str, seed_index: int = 0, is_production: bool = False) -> Path:
    """Compute standardized output directory for a run.

    Task 1.1.4: Separate outputs by receptor state and partner construct.
    """
    run_type = "prod" if is_production else "test"
    return PHASE1_OUTPUT_DIR / state_name / f"{run_type}_seed{seed_index}"


def run_single(
    config_path: Path,
    state_name: str,
    seed_index: int = 0,
    is_production: bool = False,
    dry_run: bool = False,
) -> dict:
    """Execute a single docking run (or dry-run validation).

    Returns run metadata dict.
    """
    output_dir = get_output_dir(state_name, seed_index, is_production)

    print(f"\n{'='*60}")
    print(f"Run: {state_name} | seed={seed_index} | "
          f"{'production' if is_production else 'test'}")
    print(f"  Config: {config_path}")
    print(f"  Output: {output_dir}")

    # Generate metadata
    meta = generate_run_metadata(
        config_path, state_name, output_dir,
        seed_index, is_production, dry_run
    )
    print(f"  Metadata: {output_dir / 'pyrosetta_run_metadata.json'}")
    print(f"  Models: {meta['total_global_models']}")
    print(f"  Filter: {meta['filter_version']}")

    if dry_run:
        print(f"  [DRY RUN] Skipping actual docking execution")
        return meta

    # Actual execution — import pipeline manager
    print(f"  Starting docking pipeline...")
    try:
        # Change to output directory for pipeline manager
        original_cwd = os.getcwd()
        os.chdir(str(output_dir))

        # The pipeline manager uses the config's input_pdb_name relative to CWD
        # We need to provide an absolute path
        abs_config = str(config_path.resolve())
        abs_pdb = str((PROJECT_ROOT / meta["input_pdb"]).resolve())

        from egfr_pipeline.pyrosetta_docking.pipeline_manager import PipelineManager
        pm = PipelineManager(abs_config, override_input_pdb=abs_pdb)
        pm.execute()

        os.chdir(original_cwd)

        meta["status"] = "completed"
        meta["completed_at"] = datetime.now().isoformat()

    except Exception as e:
        os.chdir(original_cwd)
        meta["status"] = "failed"
        meta["error"] = str(e)
        print(f"  ERROR: {e}")

    # Update metadata with final status
    meta_path = output_dir / "pyrosetta_run_metadata.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)

    return meta


def main():
    parser = argparse.ArgumentParser(
        description="Phase 1 TG 1.1: Launch PyRosetta docking with standardized metadata"
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--test", action="store_true",
                       help="Run test configs (1K models per state)")
    group.add_argument("--production", action="store_true",
                       help="Run all production configs (multi-seed)")
    group.add_argument("--config", type=Path,
                       help="Run a single specific config file")

    parser.add_argument("--state", choices=RECEPTOR_STATES,
                        help="Limit to a single receptor state")
    parser.add_argument("--seed", type=int, default=None,
                        help="Limit to a single seed index (production only)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate configs and generate metadata only (no docking)")
    args = parser.parse_args()

    print("Phase 1 — Task 1.1: PyRosetta Docking Launch")
    print(f"Output base: {PHASE1_OUTPUT_DIR}")

    # Validate inputs
    print("\nValidating Phase 1 inputs (TG 1.0):")
    if not validate_phase1_inputs():
        print("\nERROR: Phase 1 inputs incomplete. Run TG 1.0 first:")
        print("  python -m egfr_pipeline.phase1.prepare_inputs")
        return 1

    # Single config mode
    if args.config:
        if not args.config.exists():
            print(f"ERROR: Config not found: {args.config}")
            return 1
        # Parse state name from filename
        stem = args.config.stem
        state = None
        for s in RECEPTOR_STATES:
            if s in stem:
                state = s
                break
        if state is None:
            print(f"ERROR: Cannot determine receptor state from filename: {stem}")
            return 1
        seed_idx = 0
        for part in stem.split("_"):
            if part.startswith("seed"):
                seed_idx = int(part[4:])
        is_prod = "prod" in stem
        run_single(args.config, state, seed_idx, is_prod, args.dry_run)
        return 0

    # Determine which states to run
    states = [args.state] if args.state else RECEPTOR_STATES

    runs = []

    if args.test or (not args.test and not args.production):
        # Test configs
        for state in states:
            cfg = PHASE1_CONFIG_DIR / f"phase1_test_{state}.ini"
            if not cfg.exists():
                print(f"  WARNING: Config not found: {cfg}")
                print(f"  Run: python -m egfr_pipeline.phase1.generate_configs")
                continue
            meta = run_single(cfg, state, 0, False, args.dry_run)
            runs.append(meta)

    if args.production:
        # Production configs (multi-seed)
        seeds = [args.seed] if args.seed is not None else range(5)
        for state in states:
            for seed_idx in seeds:
                cfg = PHASE1_CONFIG_DIR / f"phase1_prod_{state}_seed{seed_idx}.ini"
                if not cfg.exists():
                    print(f"  WARNING: Config not found: {cfg}")
                    continue
                meta = run_single(cfg, state, seed_idx, True, args.dry_run)
                runs.append(meta)

    # Summary
    print(f"\n{'='*60}")
    print(f"LAUNCH SUMMARY")
    print(f"{'='*60}")
    print(f"  Total runs: {len(runs)}")
    for meta in runs:
        status = "DRY_RUN" if meta.get("dry_run") else meta.get("status", "unknown")
        print(f"  {meta['receptor_id']} seed{meta['seed_index']}: {status} "
              f"({meta['total_global_models']} models)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
