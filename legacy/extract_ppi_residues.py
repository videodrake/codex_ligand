#!/usr/bin/env python3
"""Task Group 6: Standardize PPI residue-level outputs.

Reads PyRosetta pipeline outputs (final_ranking.csv, cluster_summary.csv,
InterfaceEnergies CSVs) and produces a standardized receptor-side residue
summary that can be compared alongside Vina pocket residues.

Also provides a stub for AlphaFold-Multimer output parsing.

PPI outputs are explicitly **auxiliary evidence** — they inform receptor
surface analysis but do not override Vina-derived pocket definitions.
"""
import argparse
import csv
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


# ---------------------------------------------------------------------------
# Residue normalization (shared with compare_pockets.py)
# ---------------------------------------------------------------------------

_RESNAME_NORMALIZE = {
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "CYX": "CYS",
    "HIE": "HIS", "HID": "HIS", "HIP": "HIS",
}


def normalize_residue_id(residue_id: str) -> str:
    """Normalize 'A:MET971' -> 'MET971' (strip chain, normalize resname)."""
    if ":" in residue_id:
        residue_id = residue_id.split(":", 1)[1]
    i = 0
    while i < len(residue_id) and not residue_id[i].isdigit() and residue_id[i] != '-':
        i += 1
    resname = residue_id[:i]
    resnum = residue_id[i:]
    resname = _RESNAME_NORMALIZE.get(resname, resname)
    return f"{resname}{resnum}"


def extract_resnum(normalized_id: str) -> Optional[int]:
    """Extract residue number from normalized id like 'MET971' -> 971."""
    m = re.search(r'(-?\d+)$', normalized_id)
    return int(m.group(1)) if m else None


# ---------------------------------------------------------------------------
# PyRosetta output parsers
# ---------------------------------------------------------------------------

def load_csv(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def parse_binding_residues(raw: str) -> List[str]:
    """Parse 'A:LEU819,A:ALA822' into normalized list ['LEU819', 'ALA822']."""
    if not raw or raw in ("None", "No_Chain_2", "Analysis_Failed"):
        return []
    return [normalize_residue_id(r.strip()) for r in raw.split(",") if r.strip()]


def extract_pyrosetta_interface_residues(
    result_dir: Path,
    receptor_id: str,
) -> Dict[str, object]:
    """Extract receptor-side interface residues from a PyRosetta result directory.

    Reads final_ranking.csv and cluster_summary.csv to build:
    - Union of all receptor-side (Chain A) contact residues
    - Per-residue frequency (how many models contact this residue)
    - Per-residue best energy (from InterfaceEnergies CSVs if available)

    Returns a dict with all extracted data.
    """
    final_ranking = load_csv(result_dir / "final_ranking.csv")
    if not final_ranking:
        # Try nested path
        final_ranking = load_csv(result_dir / "final_result" / "final_ranking.csv")

    cluster_summary = load_csv(result_dir / "cluster_results" / "cluster_summary.csv")

    # --- Collect receptor-side binding residues across all models ---
    residue_counter: Counter = Counter()
    model_count = 0
    model_energies: Dict[str, List[float]] = defaultdict(list)

    for row in final_ranking:
        binding_a = row.get("Binding_Residues_A", "") or row.get("Binding_Residues", "")
        residues = parse_binding_residues(binding_a)
        if residues:
            model_count += 1
            residue_counter.update(residues)

        # Try to read per-residue interface energies
        csv_path = row.get("File_CSV", "")
        if csv_path:
            iface_csv = Path(csv_path).parent / Path(csv_path).name.replace("_Energies.csv", "_InterfaceEnergies.csv")
            if iface_csv.exists():
                for erow in load_csv(iface_csv):
                    resid = erow.get("Residue_ID", "")
                    chain = erow.get("Chain", "")
                    delta_e = erow.get("DeltaE_total", "")
                    if chain == "A" and delta_e:
                        norm = normalize_residue_id(f"{chain}:{erow.get('Residue_Name', '')}{resid}")
                        try:
                            model_energies[norm].append(float(delta_e))
                        except ValueError:
                            pass

    # --- Also collect from cluster_summary for broader coverage ---
    cluster_residue_counter: Counter = Counter()
    for row in cluster_summary:
        binding_a = row.get("Binding_Residues_A", "") or row.get("Binding_Residues", "")
        residues = parse_binding_residues(binding_a)
        cluster_residue_counter.update(residues)

    # Merge: union from both sources
    all_residues = set(residue_counter.keys()) | set(cluster_residue_counter.keys())

    # --- Build per-residue summary rows ---
    residue_rows = []
    for res in sorted(all_residues, key=lambda r: (extract_resnum(r) or 0, r)):
        freq_final = residue_counter.get(res, 0)
        freq_cluster = cluster_residue_counter.get(res, 0)
        occupancy = freq_final / model_count if model_count > 0 else 0.0
        energies = model_energies.get(res, [])
        mean_delta_e = sum(energies) / len(energies) if energies else None
        min_delta_e = min(energies) if energies else None
        resnum = extract_resnum(res)

        residue_rows.append({
            "receptor_id": receptor_id,
            "source": "pyrosetta_ppi",
            "residue_id": res,
            "residue_num": resnum if resnum is not None else "",
            "frequency_final_ranking": freq_final,
            "frequency_cluster_summary": freq_cluster,
            "n_models_final_ranking": model_count,
            "occupancy": round(occupancy, 4),
            "mean_interface_delta_e": round(mean_delta_e, 3) if mean_delta_e is not None else "",
            "best_interface_delta_e": round(min_delta_e, 3) if min_delta_e is not None else "",
        })

    # --- Summary metrics ---
    dg_values = []
    dsasa_values = []
    for row in final_ranking:
        try:
            dg_values.append(float(row.get("dG_separated", "nan")))
        except (ValueError, TypeError):
            pass
        try:
            dsasa_values.append(float(row.get("dSASA", "nan")))
        except (ValueError, TypeError):
            pass

    summary = {
        "receptor_id": receptor_id,
        "source": "pyrosetta_ppi",
        "n_final_models": len(final_ranking),
        "n_clusters": len(cluster_summary),
        "n_interface_residues": len(all_residues),
        "top_residues": ";".join(
            r for r, _ in residue_counter.most_common(10)
        ),
        "best_dg": round(min(dg_values), 3) if dg_values else "",
        "mean_dg": round(sum(dg_values) / len(dg_values), 3) if dg_values else "",
        "best_dsasa": round(max(dsasa_values), 1) if dsasa_values else "",
    }

    return {
        "residue_rows": residue_rows,
        "summary": summary,
    }


# ---------------------------------------------------------------------------
# AlphaFold-Multimer stub
# ---------------------------------------------------------------------------

def extract_afm_interface_residues(
    model_pdb: Path,
    receptor_chain: str = "A",
    partner_chain: str = "B",
    contact_cutoff: float = 8.0,
    receptor_id: str = "",
) -> Dict[str, object]:
    """Extract receptor-side contact residues from an AlphaFold-Multimer model PDB.

    Uses simple CA-CA distance (no PyRosetta dependency) for portability.
    For production use, consider adding PAE (predicted aligned error) weighting.

    Args:
        model_pdb: Path to AFM model PDB file
        receptor_chain: Chain ID for receptor
        partner_chain: Chain ID for binding partner
        contact_cutoff: CA-CA distance cutoff for contact (default 8.0 A)
        receptor_id: Receptor identifier for output labeling
    """
    if not model_pdb.exists():
        return {"residue_rows": [], "summary": {"receptor_id": receptor_id, "source": "alphafold_multimer", "error": "file_not_found"}}

    # Parse CA atoms per chain
    chain_atoms: Dict[str, List[Tuple[str, int, float, float, float]]] = defaultdict(list)
    with open(model_pdb, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            chain = line[21].strip()
            resname = line[17:20].strip()
            try:
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            chain_atoms[chain].append((resname, resnum, x, y, z))

    receptor_atoms = chain_atoms.get(receptor_chain, [])
    partner_atoms = chain_atoms.get(partner_chain, [])

    if not receptor_atoms or not partner_atoms:
        return {
            "residue_rows": [],
            "summary": {
                "receptor_id": receptor_id,
                "source": "alphafold_multimer",
                "error": f"missing_chain_{receptor_chain}_or_{partner_chain}",
            },
        }

    cutoff_sq = contact_cutoff ** 2
    contact_residues: Dict[str, float] = {}  # normalized_id -> min_distance

    for rname, rnum, rx, ry, rz in receptor_atoms:
        for _, _, px, py, pz in partner_atoms:
            dist_sq = (rx - px) ** 2 + (ry - py) ** 2 + (rz - pz) ** 2
            if dist_sq <= cutoff_sq:
                norm = normalize_residue_id(f"{receptor_chain}:{rname}{rnum}")
                dist = dist_sq ** 0.5
                if norm not in contact_residues or dist < contact_residues[norm]:
                    contact_residues[norm] = dist
                break  # one partner contact is enough for this residue

    residue_rows = []
    for res in sorted(contact_residues, key=lambda r: (extract_resnum(r) or 0, r)):
        residue_rows.append({
            "receptor_id": receptor_id,
            "source": "alphafold_multimer",
            "residue_id": res,
            "residue_num": extract_resnum(res) or "",
            "min_ca_distance": round(contact_residues[res], 2),
        })

    summary = {
        "receptor_id": receptor_id,
        "source": "alphafold_multimer",
        "n_interface_residues": len(residue_rows),
        "top_residues": ";".join(r["residue_id"] for r in residue_rows[:10]),
        "model_file": str(model_pdb),
    }

    return {"residue_rows": residue_rows, "summary": summary}


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

PPI_RESIDUE_FIELDS = [
    "receptor_id",
    "source",
    "residue_id",
    "residue_num",
    "frequency_final_ranking",
    "frequency_cluster_summary",
    "n_models_final_ranking",
    "occupancy",
    "mean_interface_delta_e",
    "best_interface_delta_e",
]

PPI_SUMMARY_FIELDS = [
    "receptor_id",
    "source",
    "n_final_models",
    "n_clusters",
    "n_interface_residues",
    "top_residues",
    "best_dg",
    "mean_dg",
    "best_dsasa",
]

AFM_RESIDUE_FIELDS = [
    "receptor_id",
    "source",
    "residue_id",
    "residue_num",
    "min_ca_distance",
]


def write_csv_rows(path: Path, rows: List[dict], fieldnames: List[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


# ---------------------------------------------------------------------------
# Batch extraction from project config
# ---------------------------------------------------------------------------

def load_config(path: str) -> dict:
    p = Path(path)
    if p.suffix in (".yaml", ".yml"):
        if not HAS_YAML:
            raise RuntimeError("pyyaml required for YAML config")
        with open(p, encoding="utf-8") as f:
            return yaml.safe_load(f)
    import json
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def extract_pyrosetta_batch(
    config_path: str,
    output_dir: Optional[str] = None,
) -> Tuple[Path, Path]:
    """Extract PyRosetta PPI residues for all receptors in config.

    Looks for PyRosetta result directories matching receptor IDs.
    Expected layout: <some_dir>/<receptor_id>/final_ranking.csv
    """
    config = load_config(config_path)
    out_root = Path(output_dir) if output_dir else Path(config.get("output_root", "./output"))
    project_name = config.get("project_name", "")
    if project_name:
        out_root = out_root / project_name

    ppi_config = config.get("ppi", {})
    pyrosetta_dirs = ppi_config.get("pyrosetta_result_dirs", {})

    all_residue_rows: List[dict] = []
    all_summaries: List[dict] = []

    for receptor in config.get("receptors", []):
        receptor_id = receptor["id"]
        result_dir_str = pyrosetta_dirs.get(receptor_id, "")
        if not result_dir_str:
            continue
        result_dir = Path(result_dir_str)
        if not result_dir.exists():
            print(f"[WARN] PyRosetta result dir not found for {receptor_id}: {result_dir}")
            continue

        data = extract_pyrosetta_interface_residues(result_dir, receptor_id)
        all_residue_rows.extend(data["residue_rows"])
        all_summaries.append(data["summary"])

    residue_csv = write_csv_rows(
        out_root / "ppi_pyrosetta_residues.csv",
        all_residue_rows,
        PPI_RESIDUE_FIELDS,
    )
    summary_csv = write_csv_rows(
        out_root / "ppi_pyrosetta_summary.csv",
        all_summaries,
        PPI_SUMMARY_FIELDS,
    )
    return residue_csv, summary_csv


def extract_afm_batch(
    config_path: str,
    output_dir: Optional[str] = None,
) -> Tuple[Path, Path]:
    """Extract AFM interface residues for all receptors in config.

    Expected config section:
        ppi:
          afm_models:
            3GT8_raw: path/to/model.pdb
            3GT8_cl38_48: path/to/model.pdb
    """
    config = load_config(config_path)
    out_root = Path(output_dir) if output_dir else Path(config.get("output_root", "./output"))
    project_name = config.get("project_name", "")
    if project_name:
        out_root = out_root / project_name

    ppi_config = config.get("ppi", {})
    afm_models = ppi_config.get("afm_models", {})
    afm_settings = ppi_config.get("afm_settings", {})
    receptor_chain = afm_settings.get("receptor_chain", "A")
    partner_chain = afm_settings.get("partner_chain", "B")
    contact_cutoff = float(afm_settings.get("contact_cutoff", 8.0))

    all_residue_rows: List[dict] = []
    all_summaries: List[dict] = []

    for receptor in config.get("receptors", []):
        receptor_id = receptor["id"]
        model_path_str = afm_models.get(receptor_id, "")
        if not model_path_str:
            continue
        model_pdb = Path(model_path_str)
        if not model_pdb.exists():
            print(f"[WARN] AFM model not found for {receptor_id}: {model_pdb}")
            continue

        data = extract_afm_interface_residues(
            model_pdb, receptor_chain, partner_chain, contact_cutoff, receptor_id,
        )
        all_residue_rows.extend(data["residue_rows"])
        all_summaries.append(data["summary"])

    residue_csv = write_csv_rows(
        out_root / "ppi_afm_residues.csv",
        all_residue_rows,
        AFM_RESIDUE_FIELDS,
    )
    summary_csv = write_csv_rows(
        out_root / "ppi_afm_summary.csv",
        all_summaries,
        ["receptor_id", "source", "n_interface_residues", "top_residues", "model_file"],
    )
    return residue_csv, summary_csv


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract standardized PPI residue summaries (Task Group 6).",
    )
    parser.add_argument("--config", default=None, help="Project YAML/JSON config path (required for batch mode)")
    parser.add_argument("--output-dir", default=None, help="Output directory override")
    parser.add_argument(
        "--source", choices=["pyrosetta", "afm", "all"], default="all",
        help="Which PPI source to extract (default: all)",
    )
    # Single-directory mode for quick extraction without full config
    parser.add_argument(
        "--pyrosetta-dir", default=None,
        help="Single PyRosetta result directory (bypasses config receptor lookup)",
    )
    parser.add_argument(
        "--receptor-id", default="unknown",
        help="Receptor ID label (used with --pyrosetta-dir)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Quick single-directory mode
    if args.pyrosetta_dir:
        result_dir = Path(args.pyrosetta_dir)
        data = extract_pyrosetta_interface_residues(result_dir, args.receptor_id)
        out_dir = Path(args.output_dir) if args.output_dir else result_dir
        res_csv = write_csv_rows(
            out_dir / "ppi_pyrosetta_residues.csv",
            data["residue_rows"],
            PPI_RESIDUE_FIELDS,
        )
        print(f"[OK] Wrote {len(data['residue_rows'])} residues: {res_csv}")
        return

    # Batch mode from config
    if not args.config:
        print("[ERROR] --config is required for batch mode")
        return

    if args.source in ("pyrosetta", "all"):
        try:
            r, s = extract_pyrosetta_batch(args.config, args.output_dir)
            print(f"[OK] PyRosetta residues: {r}")
            print(f"[OK] PyRosetta summary: {s}")
        except Exception as e:
            print(f"[WARN] PyRosetta extraction: {e}")

    if args.source in ("afm", "all"):
        try:
            r, s = extract_afm_batch(args.config, args.output_dir)
            print(f"[OK] AFM residues: {r}")
            print(f"[OK] AFM summary: {s}")
        except Exception as e:
            print(f"[WARN] AFM extraction: {e}")


if __name__ == "__main__":
    main()
