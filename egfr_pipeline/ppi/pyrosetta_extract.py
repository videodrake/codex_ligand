"""PyRosetta PPI residue extraction.

Reads PyRosetta pipeline outputs (final_ranking.csv, cluster_summary.csv,
InterfaceEnergies CSVs) and produces a standardized receptor-side residue
summary that can be compared alongside Vina pocket residues.

PPI outputs are explicitly **auxiliary evidence** -- they inform receptor
surface analysis but do not override Vina-derived pocket definitions.
"""
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from egfr_pipeline.residue_utils import normalize_residue_id, extract_resnum
from egfr_pipeline.config import load_config


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------

def load_csv(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------------------
# Field definitions
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


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_csv_rows(path: Path, rows: List[dict], fieldnames: List[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


# ---------------------------------------------------------------------------
# PyRosetta output parsers
# ---------------------------------------------------------------------------

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
# Batch extraction from project config
# ---------------------------------------------------------------------------

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
