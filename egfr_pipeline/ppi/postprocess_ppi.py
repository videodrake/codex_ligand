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
import json
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from egfr_pipeline import paths
from egfr_pipeline.config import load_config
from egfr_pipeline.ppi.prepare_dimer_pdb import restore_chains, restore_csv


def _slugify_partner_name(value: str) -> str:
    token = (value or "").strip() or "unknown_partner"
    slug = "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in token).strip("_")
    return slug or "unknown_partner"


def _write_legacy_migration_marker(legacy_dir: Path, restored_dir: Path) -> None:
    marker_path = legacy_dir / "MOVED_TO.txt"
    if marker_path.exists():
        return
    marker_path.write_text(
        "This legacy restored directory is deprecated.\n"
        f"Current restored outputs are written to: {restored_dir}\n",
        encoding="utf-8",
    )


def _parse_run_dir_context(docking_dir: Path) -> tuple[str, str]:
    match = re.match(r"(?P<run_type>.+)_seed(?P<seed>\d+)$", docking_dir.name)
    if not match:
        return "", ""
    return match.group("run_type"), match.group("seed")


def _write_restored_manifest(
    restored_dir: Path,
    *,
    source_docking_dir: Path,
    receptor_id: str,
    partner_name: str,
) -> Path:
    run_type, seed = _parse_run_dir_context(source_docking_dir)
    state = source_docking_dir.parent.name if source_docking_dir.parent != source_docking_dir else ""
    payload = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "source_docking_dir": str(source_docking_dir),
        "source_state": state,
        "source_run_type": run_type,
        "source_seed": seed,
        "receptor_id": receptor_id,
        "partner_name": partner_name,
        "restored_dir": str(restored_dir),
    }
    path = restored_dir / "restored_manifest.json"
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Step 1: Restore chains in all docking outputs
# ---------------------------------------------------------------------------

def restore_all_results(
    docking_dir: Path,
    mapping_csv: Path,
    restored_root: Optional[Path] = None,
) -> Path:
    """Restore original chain IDs in all docking output files.

    Creates a restored output directory with restored PDBs and CSVs.
    Returns the path to the restored directory.
    """
    restored_dir = restored_root if restored_root is not None else docking_dir / "restored"
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
    partner_name: str = "",
    construct_type: str = "full_kinase_domain",
    orientation_validation_status: str = "not_available",
    source_manifest_path: str = "",
    source_state: str = "",
    source_run_type: str = "",
    source_seed: str = "",
) -> None:
    """Append a PPI result entry to config's pyrosetta_result_dirs list.

    Supports the multi-partner list format:
      ppi:
        pyrosetta_result_dirs:
          3GT8_raw:
            - path: .../restored
              partner: beta_meander
              construct_type: full_kinase_domain
              orientation_validation_status: not_available

    Migrates legacy string values to list format automatically.
    Upserts if the same path+partner entry already exists.
    """
    from egfr_pipeline.config import save_config

    config = load_config(config_path)
    ppi = config.get("ppi") or {}
    config["ppi"] = ppi
    dirs = ppi.get("pyrosetta_result_dirs") or {}
    ppi["pyrosetta_result_dirs"] = dirs

    existing = dirs.get(receptor_id, [])
    # Migrate legacy string format to list
    if isinstance(existing, str):
        existing = [{
            "path": existing,
            "partner": "legacy",
            "construct_type": "unknown_construct_type",
            "orientation_validation_status": "not_available",
        }] if existing else []

    new_entry = {
        "path": str(restored_dir),
        "partner": partner_name,
        "construct_type": construct_type,
        "orientation_validation_status": orientation_validation_status,
        "source_manifest_path": source_manifest_path,
        "source_state": source_state,
        "source_run_type": source_run_type,
        "source_seed": source_seed,
    }
    updated = False
    for entry in existing:
        if entry.get("path") == new_entry["path"] and entry.get("partner") == new_entry["partner"]:
            entry.update(new_entry)
            updated = True
            break
    if not updated:
        existing.append(new_entry)

    dirs[receptor_id] = existing
    save_config(config, config_path)
    label = f" ({partner_name})" if partner_name else ""
    print(
        f"[postprocess] Updated config ppi.pyrosetta_result_dirs[{receptor_id}]"
        f"{label} = {restored_dir}"
    )


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def postprocess_ppi_results(
    config_path: str,
    docking_dir: str,
    mapping_csv: str,
    receptor_id: str,
    partner_name: str = "",
    construct_type: str = "full_kinase_domain",
    orientation_validation_status: str = "not_available",
    restored_root: Optional[str] = None,
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
        construct_type: Receptor/partner system label for downstream handoff.
        orientation_validation_status: Current orientation filter status.
        restored_root: Optional explicit destination root for restored files.
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

    config = load_config(config_path)

    partner_slug = _slugify_partner_name(partner_name or docking_path.name)

    if restored_root:
        restored_output_dir = Path(restored_root)
    else:
        restored_output_dir = (
            paths.wa_phase3_ppi_postprocess(config)
            / "restored_runs"
            / receptor_id
            / partner_slug
        )
        print(f"[postprocess] Default restored output path: {restored_output_dir}")

    print(f"\n{'='*60}")
    print(f"  PPI Post-Processing: {partner_name or docking_path.name}")
    print(f"{'='*60}")

    # Step 1: Restore chains
    print("\n--- Step 1: Chain restoration ---")
    restored_dir = restore_all_results(
        docking_path,
        mapping_path,
        restored_root=restored_output_dir,
    )
    manifest_path = _write_restored_manifest(
        restored_dir,
        source_docking_dir=docking_path,
        receptor_id=receptor_id,
        partner_name=partner_name,
    )
    source_run_type, source_seed = _parse_run_dir_context(docking_path)
    source_state = docking_path.parent.name if docking_path.parent != docking_path else ""
    legacy_restored_dir = docking_path / "restored"
    if restored_dir != legacy_restored_dir and legacy_restored_dir.exists():
        print(
            "[postprocess] Legacy restored directory already exists at "
            f"{legacy_restored_dir}; new outputs are written to {restored_dir}."
        )
        _write_legacy_migration_marker(legacy_restored_dir, restored_dir)

    # Step 2: Update config with restored dir
    print("\n--- Step 2: Register in config ---")
    _update_config_ppi_dir(
        config_path,
        receptor_id,
        restored_dir,
        partner_name,
        construct_type=construct_type,
        orientation_validation_status=orientation_validation_status,
        source_manifest_path=str(manifest_path),
        source_state=source_state,
        source_run_type=source_run_type,
        source_seed=source_seed,
    )

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
    parser.add_argument(
        "--construct-type",
        default="full_kinase_domain",
        help="Construct type label to persist in config",
    )
    parser.add_argument(
        "--orientation-validation-status",
        default="not_available",
        help="Orientation validation status to persist in config",
    )
    parser.add_argument(
        "--restored-root",
        default=None,
        help="Optional destination root for restored outputs",
    )
    parser.add_argument("--skip-extract", action="store_true", help="Skip residue extraction")
    parser.add_argument("--skip-report", action="store_true", help="Skip report regeneration")

    args = parser.parse_args()
    postprocess_ppi_results(
        config_path=args.config,
        docking_dir=args.docking_dir,
        mapping_csv=args.mapping,
        receptor_id=args.receptor_id,
        partner_name=args.partner_name,
        construct_type=args.construct_type,
        orientation_validation_status=args.orientation_validation_status,
        restored_root=args.restored_root,
        skip_extract=args.skip_extract,
        skip_report=args.skip_report,
    )


if __name__ == "__main__":
    main()
