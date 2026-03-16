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
    "partner_id",
    "source",
    "chain",
    "residue_id",
    "residue_num",
    "residue_name",
    "lobe_label",
    "construct_type",
    "orientation_validation_status",
    "frequency_final_ranking",
    "frequency_cluster_summary",
    "n_models_final_ranking",
    "occupancy",
    "mean_interface_delta_e",
    "best_interface_delta_e",
]

PPI_SUMMARY_FIELDS = [
    "receptor_id",
    "partner_id",
    "source",
    "construct_type",
    "orientation_validation_status",
    "n_final_models",
    "n_clusters",
    "n_interface_residues",
    "n_nlobe_interface_residues",
    "n_clobe_interface_residues",
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

NLOBE_CLOBE_BOUNDARY = 838


def _parse_receptor_residue(raw: str) -> dict:
    """Parse a receptor-side residue reference while preserving chain identity.

    Restored PyRosetta outputs may contain receptor residues on both chain A and
    chain B after EGFR dimer chain restoration. These should remain receptor-side
    evidence, not be treated as partner residues.
    """
    normalized = normalize_residue_id(raw.strip(), keep_chain=True)
    chain = ""
    residue_id = normalized
    if ":" in normalized:
        chain, residue_id = normalized.split(":", 1)

    resnum = extract_resnum(residue_id)
    residue_name = ""
    for idx, ch in enumerate(residue_id):
        if ch.isdigit() or ch == "-":
            residue_name = residue_id[:idx]
            break
    if not residue_name:
        residue_name = residue_id

    if resnum is None:
        lobe_label = "unknown"
    elif resnum < NLOBE_CLOBE_BOUNDARY:
        lobe_label = "N-lobe"
    else:
        lobe_label = "C-lobe"

    return {
        "chain": chain,
        "residue_id": residue_id,
        "residue_num": resnum,
        "residue_name": residue_name,
        "lobe_label": lobe_label,
    }


def parse_binding_residues(raw: str) -> List[str]:
    """Parse receptor-side residue refs while preserving restored chain IDs.

    Example:
      'A:LEU819,B:ASP855' -> ['A:LEU819', 'B:ASP855']
    """
    if not raw or raw in ("None", "No_Chain_2", "Analysis_Failed"):
        return []
    return [
        normalize_residue_id(r.strip(), keep_chain=True)
        for r in raw.split(",")
        if r.strip()
    ]


def _interface_csv_name(csv_path: str) -> str:
    name = Path(csv_path).name
    if name.endswith("_InterfaceEnergies.csv"):
        return name
    if name.endswith("_Energies.csv"):
        return name.replace("_Energies.csv", "_InterfaceEnergies.csv")
    return name


def _resolve_interface_csv_path(result_dir: Path, csv_path: str) -> Optional[Path]:
    """Resolve InterfaceEnergies CSV from a File_CSV value.

    PyRosetta final_ranking.csv may store File_CSV as:
      - a bare filename
      - a relative path
      - an absolute path
    while InterfaceEnergies files often live under ``final_result/``.
    """
    if not csv_path:
        return None

    raw = Path(csv_path)
    iface_name = _interface_csv_name(csv_path)
    candidates = []

    if raw.is_absolute():
        candidates.append(raw.with_name(iface_name))
    else:
        candidates.append(result_dir / raw.with_name(iface_name))
        candidates.append(result_dir / raw)

    candidates.append(result_dir / iface_name)
    candidates.append(result_dir / "final_result" / iface_name)

    seen = set()
    for candidate in candidates:
        normalized = str(candidate)
        if normalized in seen:
            continue
        seen.add(normalized)
        if candidate.exists():
            return candidate
    return None


def extract_pyrosetta_interface_residues(
    result_dir: Path,
    receptor_id: str,
    partner_id: str = "",
    construct_type: str = "full_kinase_domain",
    orientation_validation_status: str = "not_available",
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
            iface_csv = _resolve_interface_csv_path(result_dir, csv_path)
            if iface_csv is not None:
                for erow in load_csv(iface_csv):
                    resid = erow.get("Residue_ID", "")
                    chain = erow.get("Chain", "")
                    delta_e = erow.get("DeltaE_total", "")
                    if chain == "A" and delta_e:
                        residue_num = extract_resnum(resid)
                        residue_name = normalize_residue_id(
                            erow.get("Residue_Name", "")
                        )
                        if residue_num is None or not residue_name:
                            continue
                        norm = f"{chain}:{residue_name}{residue_num}"
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
    n_nlobe = 0
    n_clobe = 0
    for res in sorted(all_residues, key=lambda r: (extract_resnum(r) or 0, r)):
        residue_info = _parse_receptor_residue(res)
        freq_final = residue_counter.get(res, 0)
        freq_cluster = cluster_residue_counter.get(res, 0)
        occupancy = freq_final / model_count if model_count > 0 else 0.0
        energies = model_energies.get(res, [])
        mean_delta_e = sum(energies) / len(energies) if energies else None
        min_delta_e = min(energies) if energies else None
        if residue_info["lobe_label"] == "N-lobe":
            n_nlobe += 1
        elif residue_info["lobe_label"] == "C-lobe":
            n_clobe += 1

        residue_rows.append({
            "receptor_id": receptor_id,
            "partner_id": partner_id,
            "source": "pyrosetta_ppi",
            "chain": residue_info["chain"],
            "residue_id": residue_info["residue_id"],
            "residue_num": residue_info["residue_num"] if residue_info["residue_num"] is not None else "",
            "residue_name": residue_info["residue_name"],
            "lobe_label": residue_info["lobe_label"],
            "construct_type": construct_type,
            "orientation_validation_status": orientation_validation_status,
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
        "partner_id": partner_id,
        "source": "pyrosetta_ppi",
        "construct_type": construct_type,
        "orientation_validation_status": orientation_validation_status,
        "n_final_models": len(final_ranking),
        "n_clusters": len(cluster_summary),
        "n_interface_residues": len(all_residues),
        "n_nlobe_interface_residues": n_nlobe,
        "n_clobe_interface_residues": n_clobe,
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

def _normalize_result_dirs(raw) -> List[dict]:
    """Normalize pyrosetta_result_dirs value to list-of-dicts.

    Supports three formats:
      - str: legacy single path  -> [{"path": "...", "partner": ""}]
      - list of dicts: [{path, partner}, ...]
      - list of str: ["path1", "path2"]  -> [{"path": "...", "partner": ""}, ...]
    """
    if not raw:
        return []
    if isinstance(raw, str):
        return [{"path": raw, "partner": ""}]
    if isinstance(raw, list):
        entries = []
        for item in raw:
            if isinstance(item, dict):
                entries.append(item)
            elif isinstance(item, str):
                entries.append({"path": item, "partner": ""})
        return entries
    return []


def extract_pyrosetta_batch(
    config_path: str,
    output_dir: Optional[str] = None,
    pyrosetta_result_dirs: Optional[Dict[str, object]] = None,
) -> Tuple[Path, Path]:
    """Extract PyRosetta PPI residues for all receptors in config.

    Supports multiple PPI result directories per receptor:
      ppi:
        pyrosetta_result_dirs:
          3GT8_raw:
            - path: EGFR_dimer_beta_meander/restored
              partner: beta_meander
            - path: EGFR_dimer_TH1/restored
              partner: TH1

    Also supports legacy single-path format for backward compatibility.
    """
    config = load_config(config_path)
    out_root = Path(output_dir) if output_dir else Path(config.get("output_root", "./output"))
    project_name = config.get("project_name", "")
    if project_name:
        out_root = out_root / project_name

    ppi_config = config.get("ppi", {})
    pyrosetta_dirs = pyrosetta_result_dirs or ppi_config.get("pyrosetta_result_dirs", {})

    all_residue_rows: List[dict] = []
    all_summaries: List[dict] = []

    for receptor in config.get("receptors", []):
        receptor_id = receptor["id"]
        entries = _normalize_result_dirs(pyrosetta_dirs.get(receptor_id))
        if not entries:
            continue

        for entry in entries:
            result_dir = Path(entry["path"])
            partner = entry.get("partner", "")
            construct_type = entry.get("construct_type", "full_kinase_domain")
            orientation_validation_status = entry.get(
                "orientation_validation_status",
                "not_available",
            )
            if not result_dir.exists():
                print(f"[WARN] PyRosetta result dir not found for {receptor_id}"
                      f"{f' ({partner})' if partner else ''}: {result_dir}")
                continue

            data = extract_pyrosetta_interface_residues(
                result_dir,
                receptor_id,
                partner_id=partner,
                construct_type=construct_type,
                orientation_validation_status=orientation_validation_status,
            )

            # Tag source with partner name for multi-PPI distinction
            source_tag = f"pyrosetta_ppi:{partner}" if partner else "pyrosetta_ppi"
            for row in data["residue_rows"]:
                row["source"] = source_tag
            data["summary"]["source"] = source_tag

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
