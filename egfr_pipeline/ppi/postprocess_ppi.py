"""PPI post-processing automation.

After PyRosetta docking completes, this module automates:
  1. Chain restoration (merged dimer → original chain IDs)
  2. PPI residue extraction (into project-level CSVs)
  3. Report regeneration (with updated PPI evidence)

Usage:
  python -m egfr_pipeline.ppi.postprocess_ppi \\
      --config config/example-project.yaml \\
      --docking-dir EGFR_dimer_beta_meander/ \\
      --mapping input/PPI/prepared/EGFR_dimer_beta_meander_mapping.csv \\
      --partner-name beta_meander \\
      --receptor-id 3GT8_raw
"""

import argparse
import shutil
import sys
from pathlib import Path
from typing import Optional

from egfr_pipeline.ppi.prepare_dimer_pdb import restore_chains, restore_csv


# ---------------------------------------------------------------------------
# Step 1: Restore chains in all docking outputs
# ---------------------------------------------------------------------------

def restore_all_results(
    docking_dir: Path,
    mapping_csv: Path,
) -> Path:
    """Restore original chain IDs in all docking output files.

    Creates a ``restored/`` subdirectory with restored PDBs and CSVs.
    Returns the path to the restored directory.
    """
    restored_dir = docking_dir / "restored"
    restored_dir.mkdir(parents=True, exist_ok=True)

    n_pdb = 0
    n_csv = 0

    # 1. final_result/*.pdb
    final_result_dir = docking_dir / "final_result"
    if final_result_dir.exists():
        out_final = restored_dir / "final_result"
        out_final.mkdir(parents=True, exist_ok=True)
        for pdb in sorted(final_result_dir.glob("Rank*.pdb")):
            restore_chains(pdb, mapping_csv, out_final / pdb.name)
            n_pdb += 1
        # Copy non-PDB files (reports, PyMOL scripts, etc.)
        for f in final_result_dir.iterdir():
            if f.suffix != ".pdb" and not (out_final / f.name).exists():
                dest = out_final / f.name
                if f.is_file():
                    shutil.copy2(f, dest)

    # 2. final_ranking.csv
    ranking_csv = docking_dir / "final_ranking.csv"
    if ranking_csv.exists():
        restore_csv(
            ranking_csv, mapping_csv,
            restored_dir / "final_ranking.csv",
            residue_columns=["Binding_Residues_A"],
        )
        n_csv += 1

    # 3. cluster_results/cluster_summary.csv
    cluster_csv = docking_dir / "cluster_results" / "cluster_summary.csv"
    if cluster_csv.exists():
        out_cluster = restored_dir / "cluster_results"
        out_cluster.mkdir(parents=True, exist_ok=True)
        restore_csv(
            cluster_csv, mapping_csv,
            out_cluster / "cluster_summary.csv",
            residue_columns=["Binding_Residues_A"],
        )
        n_csv += 1
        # Copy cluster PDBs
        for pdb in sorted((docking_dir / "cluster_results").glob("C*.pdb")):
            restore_chains(pdb, mapping_csv, out_cluster / pdb.name)
            n_pdb += 1

    # 4. Validate: check no offset residues remain in restored CSVs
    n_offset = _check_restored_offsets(restored_dir)
    if n_offset > 0:
        print(f"[postprocess] WARNING: {n_offset} offset residues (>1700) remain after restoration!")

    print(f"\n[postprocess] Restored {n_pdb} PDB files, {n_csv} CSV files")
    print(f"[postprocess] Output: {restored_dir}")
    return restored_dir


def _check_restored_offsets(restored_dir: Path) -> int:
    """Check restored CSVs for unrestored offset residues (resnum > 1700)."""
    import re
    count = 0
    for csv_path in restored_dir.rglob("*.csv"):
        try:
            with open(csv_path, encoding="utf-8") as f:
                content = f.read()
            # Find residue references like ALA1750, MET1971
            for m in re.finditer(r'[A-Z]{3}(\d{4,})', content):
                resnum = int(m.group(1))
                if resnum > 1700:
                    count += 1
        except (OSError, IOError):
            pass
    return count


# ---------------------------------------------------------------------------
# Step 2: Register restored results in config for PPI extraction
# ---------------------------------------------------------------------------

def _update_config_ppi_dir(
    config_path: str,
    receptor_id: str,
    restored_dir: Path,
) -> None:
    """Update config's ppi.pyrosetta_result_dirs with the restored path.

    This modifies the config in-memory only (does not write to disk)
    unless the caller saves it.
    """
    from egfr_pipeline.config import load_config, save_config

    config = load_config(config_path)
    ppi = config.setdefault("ppi", {})
    dirs = ppi.setdefault("pyrosetta_result_dirs", {})
    dirs[receptor_id] = str(restored_dir)
    save_config(config, config_path)
    print(f"[postprocess] Updated config ppi.pyrosetta_result_dirs[{receptor_id}] = {restored_dir}")


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def postprocess_ppi_results(
    config_path: str,
    docking_dir: str,
    mapping_csv: str,
    receptor_id: str,
    partner_name: str = "",
    skip_extract: bool = False,
    skip_report: bool = False,
) -> Path:
    """Full PPI post-processing pipeline.

    Args:
        config_path: Project YAML/JSON config path.
        docking_dir: PyRosetta docking output directory.
        mapping_csv: Chain mapping CSV from prepare_dimer_pdb.
        receptor_id: Receptor ID to register results under.
        partner_name: Partner name for labeling.
        skip_extract: If True, skip PPI residue extraction.
        skip_report: If True, skip report regeneration.

    Returns:
        Path to the restored directory.
    """
    docking_path = Path(docking_dir)
    mapping_path = Path(mapping_csv)

    if not docking_path.exists():
        print(f"[ERROR] Docking directory not found: {docking_path}", file=sys.stderr)
        return Path()

    if not mapping_path.exists():
        print(f"[ERROR] Mapping CSV not found: {mapping_path}", file=sys.stderr)
        return Path()

    print(f"\n{'='*60}")
    print(f"  PPI Post-Processing: {partner_name or docking_path.name}")
    print(f"{'='*60}")

    # Step 1: Restore chains
    print("\n--- Step 1: Chain restoration ---")
    restored_dir = restore_all_results(docking_path, mapping_path)

    # Step 2: Update config with restored dir
    print("\n--- Step 2: Register in config ---")
    _update_config_ppi_dir(config_path, receptor_id, restored_dir)

    # Step 3: PPI residue extraction
    if not skip_extract:
        print("\n--- Step 3: PPI residue extraction ---")
        try:
            from egfr_pipeline.ppi.pyrosetta_extract import extract_pyrosetta_batch
            res_csv, sum_csv = extract_pyrosetta_batch(config_path)
            print(f"  Residues: {res_csv}")
            print(f"  Summary:  {sum_csv}")
        except Exception as e:
            print(f"[WARN] PPI extraction failed: {e}")
    else:
        print("\n--- Step 3: Skipped (--skip-extract) ---")

    # Step 4: Report regeneration
    if not skip_report:
        print("\n--- Step 4: Report regeneration ---")
        try:
            from egfr_pipeline.report import generate_report
            report_path, combined_csv = generate_report(config_path)
            print(f"  Report:   {report_path}")
            print(f"  Combined: {combined_csv}")
        except Exception as e:
            print(f"[WARN] Report generation failed: {e}")
    else:
        print("\n--- Step 4: Skipped (--skip-report) ---")

    print(f"\n{'='*60}")
    print(f"  PPI Post-Processing Complete")
    print(f"  Restored: {restored_dir}")
    print(f"{'='*60}\n")

    return restored_dir


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Automate PPI post-docking chain restoration and extraction",
    )
    parser.add_argument("--config", required=True, help="Project config (YAML/JSON)")
    parser.add_argument("--docking-dir", required=True, help="PyRosetta docking output dir")
    parser.add_argument("--mapping", required=True, help="Chain mapping CSV")
    parser.add_argument("--receptor-id", required=True, help="Receptor ID for this result")
    parser.add_argument("--partner-name", default="", help="Partner name label")
    parser.add_argument("--skip-extract", action="store_true", help="Skip residue extraction")
    parser.add_argument("--skip-report", action="store_true", help="Skip report regeneration")

    args = parser.parse_args()
    postprocess_ppi_results(
        config_path=args.config,
        docking_dir=args.docking_dir,
        mapping_csv=args.mapping,
        receptor_id=args.receptor_id,
        partner_name=args.partner_name,
        skip_extract=args.skip_extract,
        skip_report=args.skip_report,
    )


if __name__ == "__main__":
    main()
