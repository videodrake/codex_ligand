# -*- coding: utf-8 -*-
"""Multi-scoring-function consensus re-scoring module.

Re-scores existing docked PDB structures using multiple PyRosetta energy
functions to reduce single-function systematic bias.  Models that rank
consistently high across all functions are flagged as "consensus hits",
providing stronger confidence than any single score.

Scoring functions used:
  - ref2015       (default, already used in production)
  - beta_nov16    (beta energy function, different parametrisation)
  - ref2015_cart  (Cartesian-space scoring)

Usage (HPC only — never run in dev environment):
  python -m egfr_pipeline.pyrosetta_docking.consensus_scoring \\
      --input_dir output/workflow_a/phase2_ppi_docking/3GT8_raw/prod_seed0 \\
      --output consensus_scores.csv

Design:
  - Reads filter_passed/*.pdb (or final_result/*.pdb)
  - Applies InterfaceAnalyzerMover with each score function
  - Outputs per-model dG from each function + consensus flag
  - "consensus_hit" = model ranks in top percentile across ALL functions
"""

import argparse
import csv
import logging
import math
import os
import sys
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger("pipeline.consensus")

# Score functions to evaluate — order matters for CSV column naming.
# Each entry: (display_name, rosetta_score_function_name)
SCORE_FUNCTIONS = [
    ("ref2015", "ref2015"),
    ("beta_nov16", "beta_nov16"),
    ("ref2015_cart", "ref2015_cart"),
]

# Default percentile threshold for "consensus hit" (top N%)
DEFAULT_CONSENSUS_PERCENTILE = 20

CSV_FIELDS = [
    "model_file",
    "dG_ref2015",
    "dG_beta_nov16",
    "dG_ref2015_cart",
    "dSASA_ref2015",
    "sc_ref2015",
    "rank_ref2015",
    "rank_beta_nov16",
    "rank_ref2015_cart",
    "max_rank",
    "consensus_hit",
]


# ---------------------------------------------------------------------------
# Core scoring
# ---------------------------------------------------------------------------

def _init_worker():
    """Initialise PyRosetta once per worker process."""
    try:
        from . import pyrosetta_init
        pyrosetta_init.safe_init()
    except Exception:
        # Standalone execution fallback
        import pyrosetta
        pyrosetta.init("-mute all", silent=True)


def _score_single_model(args: Tuple[str, List[Tuple[str, str]]]) -> Optional[Dict[str, Any]]:
    """Score one PDB with all score functions.

    Args:
        args: (pdb_path, score_function_specs)

    Returns:
        dict with model_file + dG per function, or None on failure.
    """
    pdb_path, sf_specs = args
    try:
        from pyrosetta import pose_from_pdb, create_score_function
        from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
        from pyrosetta.rosetta.core.pose import getPoseExtraScore
    except ImportError:
        logger.error("PyRosetta not available — cannot score")
        return None

    result: Dict[str, Any] = {"model_file": os.path.basename(pdb_path)}

    try:
        pose = pose_from_pdb(pdb_path)
    except Exception as e:
        logger.warning(f"Failed to load {pdb_path}: {e}")
        return None

    if pose.num_chains() < 2:
        logger.warning(f"Single-chain pose, skipping: {pdb_path}")
        return None

    # Detect interface definition
    conf = pose.conformation()
    c1 = pose.pdb_info().chain(conf.chain_begin(1))
    c2 = pose.pdb_info().chain(conf.chain_begin(2))
    interface_def = f"{c1}_{c2}"

    for display_name, sf_name in sf_specs:
        try:
            sfxn = create_score_function(sf_name)
            sfxn(pose)  # score the pose

            iam = InterfaceAnalyzerMover(interface_def)
            iam.set_scorefunction(sfxn)
            iam.set_pack_separated(True)

            work_pose = pose.clone()
            iam.apply(work_pose)

            dG = float(getPoseExtraScore(work_pose, "dG_separated"))
            result[f"dG_{display_name}"] = round(dG, 3)

            # Extra metrics only for primary function
            if display_name == "ref2015":
                try:
                    dSASA = float(getPoseExtraScore(work_pose, "dSASA_int"))
                except Exception:
                    dSASA = 0.0
                try:
                    sc = float(getPoseExtraScore(work_pose, "sc_value"))
                except Exception:
                    sc = 0.0
                result["dSASA_ref2015"] = round(dSASA, 1)
                result["sc_ref2015"] = round(sc, 4)

        except Exception as e:
            logger.warning(f"Scoring failed for {pdb_path} with {sf_name}: {e}")
            result[f"dG_{display_name}"] = None

    return result


# ---------------------------------------------------------------------------
# Consensus analysis
# ---------------------------------------------------------------------------

def compute_consensus(
    results: List[Dict[str, Any]],
    percentile: int = DEFAULT_CONSENSUS_PERCENTILE,
) -> List[Dict[str, Any]]:
    """Add rank and consensus_hit flag to scored results.

    A model is a "consensus hit" if it ranks in the top `percentile`%
    across ALL scoring functions.
    """
    sf_names = [name for name, _ in SCORE_FUNCTIONS]
    dg_keys = [f"dG_{name}" for name in sf_names]

    # Filter to models with all scores available
    complete = [r for r in results if all(r.get(k) is not None for k in dg_keys)]
    if not complete:
        logger.warning("No models with complete scores across all functions")
        for r in results:
            for name in sf_names:
                r[f"rank_{name}"] = ""
            r["max_rank"] = ""
            r["consensus_hit"] = "incomplete"
        return results

    n = len(complete)
    threshold_rank = max(1, math.ceil(n * percentile / 100.0))

    # Rank by each score function (lower dG = better = rank 1)
    for name in sf_names:
        key = f"dG_{name}"
        sorted_models = sorted(complete, key=lambda r: r[key])
        for rank, model in enumerate(sorted_models, 1):
            model[f"rank_{name}"] = rank

    # Max rank across functions (worst ranking)
    for model in complete:
        ranks = [model[f"rank_{name}"] for name in sf_names]
        model["max_rank"] = max(ranks)
        model["consensus_hit"] = "yes" if model["max_rank"] <= threshold_rank else "no"

    # Mark incomplete models
    incomplete_set = set(id(r) for r in results) - set(id(r) for r in complete)
    for r in results:
        if id(r) in incomplete_set:
            for name in sf_names:
                r.setdefault(f"rank_{name}", "")
            r.setdefault("max_rank", "")
            r.setdefault("consensus_hit", "incomplete")

    return results


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def find_model_pdbs(input_dir: str) -> List[str]:
    """Find PDB files to re-score.

    Priority: filter_passed/ > final_result/ > *.pdb in root.
    """
    root = Path(input_dir)
    candidates = []

    for subdir in ["filter_passed", "final_result"]:
        d = root / subdir
        if d.is_dir():
            pdbs = sorted(d.glob("*.pdb"))
            if pdbs:
                candidates = [str(p) for p in pdbs]
                logger.info(f"Found {len(candidates)} PDBs in {subdir}/")
                return candidates

    # Fallback: root-level PDBs
    pdbs = sorted(root.glob("*.pdb"))
    if pdbs:
        candidates = [str(p) for p in pdbs]
        logger.info(f"Found {len(candidates)} PDBs in root directory")

    return candidates


def write_consensus_csv(
    results: List[Dict[str, Any]],
    output_path: str,
) -> Path:
    """Write consensus scoring results to CSV."""
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    with open(out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    logger.info(f"Consensus scores written to {out} ({len(results)} models)")
    return out


def run_consensus_scoring(
    input_dir: str,
    output_path: str,
    n_workers: int = 0,
    percentile: int = DEFAULT_CONSENSUS_PERCENTILE,
) -> List[Dict[str, Any]]:
    """Main entry point: find PDBs, score with all functions, compute consensus.

    Args:
        input_dir: Directory containing docked PDB files.
        output_path: Path for output CSV.
        n_workers: Number of parallel workers (0 = auto).
        percentile: Top N% threshold for consensus hit.

    Returns:
        List of result dicts.
    """
    pdbs = find_model_pdbs(input_dir)
    if not pdbs:
        logger.error(f"No PDB files found in {input_dir}")
        return []

    sf_specs = [(name, sf) for name, sf in SCORE_FUNCTIONS]
    tasks = [(pdb, sf_specs) for pdb in pdbs]

    if n_workers == 0:
        n_workers = min(len(tasks), max(1, cpu_count() - 1))

    logger.info(f"Scoring {len(pdbs)} models with {len(SCORE_FUNCTIONS)} "
                f"functions using {n_workers} workers")

    if n_workers == 1:
        _init_worker()
        results = [_score_single_model(t) for t in tasks]
    else:
        with Pool(n_workers, initializer=_init_worker) as pool:
            results = list(pool.imap(_score_single_model, tasks))

    # Filter None results
    results = [r for r in results if r is not None]
    logger.info(f"Successfully scored {len(results)}/{len(pdbs)} models")

    results = compute_consensus(results, percentile)
    write_consensus_csv(results, output_path)

    # Summary
    n_consensus = sum(1 for r in results if r.get("consensus_hit") == "yes")
    logger.info(f"Consensus hits (top {percentile}% in all functions): "
                f"{n_consensus}/{len(results)}")

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Re-score docked PDBs with multiple scoring functions"
    )
    parser.add_argument(
        "--input_dir", required=True,
        help="Directory containing docked PDB files (e.g. prod_seed0/)",
    )
    parser.add_argument(
        "--output", default="consensus_scores.csv",
        help="Output CSV path (default: consensus_scores.csv)",
    )
    parser.add_argument(
        "--workers", type=int, default=0,
        help="Number of parallel workers (0 = auto)",
    )
    parser.add_argument(
        "--percentile", type=int, default=DEFAULT_CONSENSUS_PERCENTILE,
        help=f"Top N%% for consensus hit (default: {DEFAULT_CONSENSUS_PERCENTILE})",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    run_consensus_scoring(args.input_dir, args.output, args.workers, args.percentile)


if __name__ == "__main__":
    main()
