"""AlphaFold-Multimer PPI residue extraction.

Extracts receptor-side contact residues from AlphaFold-Multimer model PDBs
using simple CA-CA distance (no PyRosetta dependency) for portability.

PPI outputs are explicitly **auxiliary evidence** -- they inform receptor
surface analysis but do not override Vina-derived pocket definitions.
"""
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from egfr_pipeline import paths
from egfr_pipeline.residue_utils import normalize_residue_id, extract_resnum
from egfr_pipeline.config import load_config
from egfr_pipeline.schemas import PPI_AFM_RESIDUES


AFM_RESIDUE_FIELDS = PPI_AFM_RESIDUES


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
# AlphaFold-Multimer extraction
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
# Batch extraction from project config
# ---------------------------------------------------------------------------

def extract_afm_batch(
    config_path: str,
    output_dir: Optional[str] = None,
) -> Tuple[Path, Path]:
    """Extract AFM interface residues for all receptors in config.

    Expected config section:
        ppi:
          afm_models:
            3GT8_raw: path/to/model.pdb
            EGFR_160-185: path/to/model.pdb
    """
    config = load_config(config_path)
    out_root = Path(output_dir) if output_dir else paths.wa_phase3_ppi_postprocess(config)

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
