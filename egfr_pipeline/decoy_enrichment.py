# -*- coding: utf-8 -*-
"""Decoy enrichment validation framework.

Provides statistical evidence that docking results are non-random by comparing
real docking scores against negative controls (scrambled partner interfaces).

Three validation approaches:
  1. Scrambled Interface Control:
     Mutate MYO1D surface residues to random amino acids → dock → compare.
     If real scores are significantly better → sequence-specific binding.

  2. Z-score Classification:
     Z = (dG_observed - mean_control) / std_control
     Z < -2.0 → statistically significant binding (p < 0.023, one-tailed).

  3. Mann-Whitney U Test:
     Non-parametric test comparing real vs control score distributions.
     p < 0.05 → distributions are significantly different.

Usage (HPC only):
  # Step 1: Generate scrambled PDBs (dev or HPC)
  python -m egfr_pipeline.decoy_enrichment generate \\
      --partner_pdb input/PPI/phase1/MYO1D_TH1.pdb \\
      --output_dir output/decoy_controls/ \\
      --n_scrambles 5

  # Step 2: Run docking with scrambled partners (HPC via qsub)
  qsub -v CONTROL_DIR=output/decoy_controls config/run_decoy_docking.pbs

  # Step 3: Analyse results (dev or HPC)
  python -m egfr_pipeline.decoy_enrichment analyse \\
      --real_scores output/workflow_a/phase2_ppi_docking/3GT8_raw/prod_seed0/scored_all_models.csv \\
      --control_dir output/decoy_controls/ \\
      --output output/decoy_controls/enrichment_report.csv
"""

import argparse
import csv
import logging
import math
import os
import random
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger("pipeline.decoy")

# Standard amino acids (single-letter codes for display, 3-letter for PDB)
AA_THREE = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
]

# Default surface SASA threshold (Angstrom^2) for identifying surface residues
DEFAULT_SURFACE_SASA_THRESHOLD = 40.0

# Fraction of surface residues to mutate in each scramble
DEFAULT_SCRAMBLE_FRACTION = 0.5

ENRICHMENT_FIELDS = [
    "state",
    "seed",
    "metric",
    "real_mean",
    "real_std",
    "real_n",
    "control_mean",
    "control_std",
    "control_n",
    "z_score",
    "mann_whitney_u",
    "mann_whitney_p",
    "effect_size_r",
    "enrichment_verdict",
]


# ---------------------------------------------------------------------------
# Part 1: Scrambled partner generation
# ---------------------------------------------------------------------------

def identify_surface_residues_from_pdb(pdb_path: str) -> List[int]:
    """Identify surface residue numbers from PDB file using CA B-factors as proxy.

    For a proper SASA calculation PyRosetta is needed — this fallback uses
    a simple neighbour-count heuristic that works without PyRosetta.

    Returns list of PDB residue numbers (1-indexed).
    """
    ca_coords: Dict[int, Tuple[float, float, float]] = {}
    with open(pdb_path, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            try:
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ca_coords[resnum] = (x, y, z)
            except (ValueError, IndexError):
                continue

    if not ca_coords:
        logger.warning(f"No CA atoms found in {pdb_path}")
        return []

    # Neighbour count heuristic: fewer neighbours (< 12 within 10A) = surface
    surface = []
    residues = sorted(ca_coords.keys())
    for rn in residues:
        x1, y1, z1 = ca_coords[rn]
        n_neighbours = 0
        for rn2 in residues:
            if rn2 == rn:
                continue
            x2, y2, z2 = ca_coords[rn2]
            dist_sq = (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
            if dist_sq < 100.0:  # 10A radius
                n_neighbours += 1
        if n_neighbours < 12:
            surface.append(rn)

    logger.info(f"Identified {len(surface)}/{len(residues)} surface residues "
                f"in {os.path.basename(pdb_path)}")
    return surface


def generate_scrambled_pdb(
    pdb_path: str,
    output_path: str,
    surface_residues: List[int],
    scramble_fraction: float = DEFAULT_SCRAMBLE_FRACTION,
    seed: Optional[int] = None,
) -> Dict[int, Tuple[str, str]]:
    """Generate a scrambled PDB by mutating surface residues.

    This performs a simple residue-name substitution in the PDB file.
    The actual side-chain remodelling must be done by PyRosetta on the
    server (relax after mutation). This function only changes the residue
    name fields so that PyRosetta's ``MutateResidue`` mover can be applied
    in the docking pipeline.

    Returns dict of {resnum: (original_aa, new_aa)} mutations applied.
    """
    if seed is not None:
        rng = random.Random(seed)
    else:
        rng = random.Random()

    n_to_mutate = max(1, int(len(surface_residues) * scramble_fraction))
    selected = rng.sample(surface_residues, min(n_to_mutate, len(surface_residues)))

    # Read original residue names
    original_names: Dict[int, str] = {}
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    resnum = int(line[22:26])
                    resname = line[17:20].strip()
                    if resnum not in original_names:
                        original_names[resnum] = resname
                except (ValueError, IndexError):
                    continue

    # Assign random mutations (different from original)
    mutations: Dict[int, Tuple[str, str]] = {}
    for resnum in selected:
        orig = original_names.get(resnum, "ALA")
        candidates = [aa for aa in AA_THREE if aa != orig]
        new_aa = rng.choice(candidates)
        mutations[resnum] = (orig, new_aa)

    # Write modified PDB (only change residue name field)
    output_lines = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    resnum = int(line[22:26])
                    if resnum in mutations:
                        _, new_aa = mutations[resnum]
                        # PDB format: columns 17-20 = residue name (right-justified in 3 chars)
                        line = line[:17] + f"{new_aa:>3s}" + line[20:]
                except (ValueError, IndexError):
                    pass
            output_lines.append(line)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.writelines(output_lines)

    logger.info(f"Generated scrambled PDB: {output_path} "
                f"({len(mutations)} mutations)")
    return mutations


def generate_scrambled_set(
    partner_pdb: str,
    output_dir: str,
    n_scrambles: int = 5,
    scramble_fraction: float = DEFAULT_SCRAMBLE_FRACTION,
) -> List[str]:
    """Generate a set of scrambled partner PDBs for negative control docking.

    Returns list of output PDB paths.
    """
    surface = identify_surface_residues_from_pdb(partner_pdb)
    if not surface:
        logger.error("No surface residues found — cannot generate controls")
        return []

    stem = Path(partner_pdb).stem
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Write mutation log
    mutation_log_path = out_dir / "scramble_mutations.csv"
    mutation_rows = []

    paths = []
    for i in range(n_scrambles):
        out_path = str(out_dir / f"{stem}_scramble{i:02d}.pdb")
        mutations = generate_scrambled_pdb(
            partner_pdb, out_path, surface,
            scramble_fraction=scramble_fraction,
            seed=42 + i,
        )
        paths.append(out_path)

        for resnum, (orig, new) in sorted(mutations.items()):
            mutation_rows.append({
                "scramble_id": i,
                "residue_num": resnum,
                "original": orig,
                "mutated": new,
            })

    # Write mutation log
    with open(mutation_log_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["scramble_id", "residue_num", "original", "mutated"])
        writer.writeheader()
        writer.writerows(mutation_rows)
    logger.info(f"Mutation log: {mutation_log_path}")

    return paths


# ---------------------------------------------------------------------------
# Part 2: Statistical analysis
# ---------------------------------------------------------------------------

def load_scores_from_csv(
    csv_path: str,
    score_column: str = "dG_separated",
) -> List[float]:
    """Load dG scores from a scored_all_models.csv or similar CSV."""
    scores = []
    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            val = row.get(score_column, "")
            try:
                scores.append(float(val))
            except (ValueError, TypeError):
                continue
    return scores


def mann_whitney_u(x: List[float], y: List[float]) -> Tuple[float, float, float]:
    """Manual Mann-Whitney U test (no scipy dependency).

    Returns (U_statistic, p_value_approx, effect_size_r).
    Uses normal approximation for p-value (valid for n > 20).
    """
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return 0.0, 1.0, 0.0

    # Combine and rank
    combined = [(v, 0) for v in x] + [(v, 1) for v in y]
    combined.sort(key=lambda t: t[0])

    # Assign ranks (handle ties by averaging)
    n = len(combined)
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j < n and combined[j][0] == combined[i][0]:
            j += 1
        avg_rank = (i + j + 1) / 2.0  # 1-indexed average
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j

    # Sum of ranks for group x (group 0)
    r1 = sum(ranks[k] for k in range(n) if combined[k][1] == 0)

    u1 = r1 - nx * (nx + 1) / 2.0
    u2 = nx * ny - u1
    u_stat = min(u1, u2)

    # Normal approximation
    mu = nx * ny / 2.0
    sigma = math.sqrt(nx * ny * (nx + ny + 1) / 12.0)

    if sigma == 0:
        return u_stat, 1.0, 0.0

    z = (u_stat - mu) / sigma
    # One-tailed p-value approximation using error function
    p_value = 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))

    # Effect size r = Z / sqrt(N)
    effect_r = abs(z) / math.sqrt(nx + ny) if (nx + ny) > 0 else 0.0

    return u_stat, p_value, effect_r


def compute_enrichment_stats(
    real_scores: List[float],
    control_scores: List[float],
    label: str = "",
) -> Dict[str, Any]:
    """Compute enrichment statistics comparing real vs control scores."""
    if not real_scores or not control_scores:
        return {
            "metric": label,
            "real_mean": 0, "real_std": 0, "real_n": 0,
            "control_mean": 0, "control_std": 0, "control_n": 0,
            "z_score": 0, "mann_whitney_u": 0, "mann_whitney_p": 1.0,
            "effect_size_r": 0, "enrichment_verdict": "insufficient_data",
        }

    def _mean(xs):
        return sum(xs) / len(xs)

    def _std(xs, m):
        if len(xs) < 2:
            return 0.0
        return math.sqrt(sum((x - m)**2 for x in xs) / (len(xs) - 1))

    real_mean = _mean(real_scores)
    real_std = _std(real_scores, real_mean)
    ctrl_mean = _mean(control_scores)
    ctrl_std = _std(control_scores, ctrl_mean)

    # Z-score: how many control SDs is the real mean below the control mean?
    # For dG: more negative = better, so real_mean < ctrl_mean is good.
    z = (real_mean - ctrl_mean) / ctrl_std if ctrl_std > 0 else 0.0

    u_stat, p_val, effect_r = mann_whitney_u(real_scores, control_scores)

    # Verdict
    if z < -2.0 and p_val < 0.05:
        verdict = "significant"
    elif z < -1.5 and p_val < 0.10:
        verdict = "suggestive"
    elif z > 0:
        verdict = "no_enrichment"
    else:
        verdict = "not_significant"

    return {
        "metric": label,
        "real_mean": round(real_mean, 3),
        "real_std": round(real_std, 3),
        "real_n": len(real_scores),
        "control_mean": round(ctrl_mean, 3),
        "control_std": round(ctrl_std, 3),
        "control_n": len(control_scores),
        "z_score": round(z, 3),
        "mann_whitney_u": round(u_stat, 1),
        "mann_whitney_p": round(p_val, 6),
        "effect_size_r": round(effect_r, 3),
        "enrichment_verdict": verdict,
    }


def analyse_enrichment(
    real_csv: str,
    control_dir: str,
    output_path: str,
    state: str = "",
    seed: str = "",
) -> List[Dict[str, Any]]:
    """Compare real docking scores against control scores.

    Looks for scored_all_models.csv in each scramble subdirectory
    within control_dir.

    Returns list of enrichment stat dicts.
    """
    real_dg = load_scores_from_csv(real_csv, "dG_separated")
    if not real_dg:
        logger.error(f"No dG scores found in {real_csv}")
        return []

    # Collect control scores from all scramble runs
    control_dg: List[float] = []
    ctrl_dir = Path(control_dir)
    for csv_path in sorted(ctrl_dir.rglob("scored_all_models.csv")):
        scores = load_scores_from_csv(str(csv_path), "dG_separated")
        control_dg.extend(scores)
        logger.info(f"  Control: {csv_path.parent.name} → {len(scores)} scores")

    if not control_dg:
        # Also try consensus_scores.csv / dG_ref2015
        for csv_path in sorted(ctrl_dir.rglob("consensus_scores.csv")):
            scores = load_scores_from_csv(str(csv_path), "dG_ref2015")
            control_dg.extend(scores)

    if not control_dg:
        logger.error(f"No control scores found in {control_dir}")
        return []

    logger.info(f"Real: {len(real_dg)} scores, Control: {len(control_dg)} scores")

    stats = compute_enrichment_stats(real_dg, control_dg, label="dG_separated")
    stats["state"] = state
    stats["seed"] = seed

    # Write report
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=ENRICHMENT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerow(stats)

    logger.info(f"Enrichment report: {out}")
    _print_enrichment_summary(stats)

    return [stats]


def _print_enrichment_summary(stats: Dict[str, Any]) -> None:
    """Print human-readable enrichment summary."""
    print(f"\n{'='*60}")
    print(f"  Decoy Enrichment Validation")
    print(f"{'='*60}")
    print(f"  Real docking:    mean dG = {stats['real_mean']:.2f} ± {stats['real_std']:.2f} "
          f"(n={stats['real_n']})")
    print(f"  Control docking: mean dG = {stats['control_mean']:.2f} ± {stats['control_std']:.2f} "
          f"(n={stats['control_n']})")
    print(f"  Z-score:         {stats['z_score']:.2f}")
    print(f"  Mann-Whitney U:  {stats['mann_whitney_u']:.0f}  (p = {stats['mann_whitney_p']:.4f})")
    print(f"  Effect size (r): {stats['effect_size_r']:.3f}")
    print(f"  Verdict:         {stats['enrichment_verdict'].upper()}")
    if stats["enrichment_verdict"] == "significant":
        print(f"  → Real docking shows STATISTICALLY SIGNIFICANT enrichment")
        print(f"    over scrambled controls (Z < -2, p < 0.05).")
    elif stats["enrichment_verdict"] == "suggestive":
        print(f"  → Suggestive but not conclusive (Z < -1.5, p < 0.10).")
        print(f"    Consider adding more control runs.")
    elif stats["enrichment_verdict"] == "no_enrichment":
        print(f"  → WARNING: Real scores are NOT better than controls.")
        print(f"    Results may reflect non-specific binding.")
    else:
        print(f"  → Trend observed but not statistically significant.")
    print(f"{'='*60}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Decoy enrichment validation for PPI docking"
    )
    subparsers = parser.add_subparsers(dest="command")

    # Generate scrambled PDBs
    gen_parser = subparsers.add_parser("generate", help="Generate scrambled partner PDBs")
    gen_parser.add_argument("--partner_pdb", required=True, help="Partner PDB file")
    gen_parser.add_argument("--output_dir", required=True, help="Output directory")
    gen_parser.add_argument("--n_scrambles", type=int, default=5, help="Number of scrambles")
    gen_parser.add_argument("--scramble_fraction", type=float, default=DEFAULT_SCRAMBLE_FRACTION)

    # Analyse enrichment
    ana_parser = subparsers.add_parser("analyse", help="Analyse enrichment statistics")
    ana_parser.add_argument("--real_scores", required=True, help="Real scored_all_models.csv")
    ana_parser.add_argument("--control_dir", required=True, help="Directory with control results")
    ana_parser.add_argument("--output", required=True, help="Output enrichment CSV")
    ana_parser.add_argument("--state", default="", help="State label")
    ana_parser.add_argument("--seed", default="", help="Seed label")

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if args.command == "generate":
        paths = generate_scrambled_set(
            args.partner_pdb, args.output_dir,
            n_scrambles=args.n_scrambles,
            scramble_fraction=args.scramble_fraction,
        )
        print(f"Generated {len(paths)} scrambled PDBs in {args.output_dir}")
    elif args.command == "analyse":
        analyse_enrichment(
            args.real_scores, args.control_dir, args.output,
            state=args.state, seed=args.seed,
        )
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
