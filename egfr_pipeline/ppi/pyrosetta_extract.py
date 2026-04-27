"""PyRosetta PPI residue extraction.

Reads PyRosetta pipeline outputs and produces:

- project-level receptor-side residue summaries for downstream verdict/report
- project-level partner-aware run summaries
- long-form per-model residue tables for provenance and reproducibility checks
- per-model score tables for run/seed auditing

PPI outputs are auxiliary evidence. They inform receptor surface analysis but do
not override Vina-derived pocket definitions.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from egfr_pipeline import paths
from egfr_pipeline.config import load_config
from egfr_pipeline.parsing_utils import safe_float
from egfr_pipeline.phase1.extract_interface import extract_run
from egfr_pipeline.residue_utils import extract_resnum, normalize_residue_id
from egfr_pipeline.schemas import (
    PPI_PYROSETTA_MODEL_TABLE,
    PPI_PYROSETTA_RESIDUE_LONG,
    PPI_PYROSETTA_RESIDUES,
    PPI_PYROSETTA_SUMMARY,
)


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------


def load_csv(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv_rows(path: Path, rows: List[dict], fieldnames: List[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


# ---------------------------------------------------------------------------
# Field definitions
# ---------------------------------------------------------------------------


PPI_RESIDUE_FIELDS = PPI_PYROSETTA_RESIDUES
PPI_SUMMARY_FIELDS = PPI_PYROSETTA_SUMMARY
PPI_RESIDUE_LONG_FIELDS = PPI_PYROSETTA_RESIDUE_LONG
PPI_MODEL_FIELDS = PPI_PYROSETTA_MODEL_TABLE


# ---------------------------------------------------------------------------
# PyRosetta output helpers
# ---------------------------------------------------------------------------


NLOBE_CLOBE_BOUNDARY = 838


def _parse_receptor_residue(raw: str) -> dict:
    """Parse a receptor-side residue reference while preserving chain identity."""
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
    """Parse receptor-side residue refs while preserving restored chain IDs."""
    if not raw or raw in ("None", "No_Chain_2", "Analysis_Failed", ""):
        return []
    return [
        normalize_residue_id(r.strip(), keep_chain=True)
        for r in raw.split(",")
        if r.strip()
    ]


def _normalize_result_dirs(raw) -> List[dict]:
    """Normalize pyrosetta_result_dirs value to list-of-dicts."""
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



def _normalize_seed_index(value: object) -> str:
    if value in ("", None):
        return ""
    return str(value).strip()


def _join_sorted(values: Iterable[object]) -> str:
    cleaned = []
    for value in values:
        if value in ("", None):
            continue
        cleaned.append(str(value).strip())
    unique = sorted(set(cleaned), key=lambda v: (v.isdigit(), int(v) if v.isdigit() else v))
    return ";".join(unique)


def _source_tag(partner_id: str) -> str:
    return f"pyrosetta_ppi:{partner_id}" if partner_id else "pyrosetta_ppi"


def _cluster_summary_path(result_dir: Path) -> Path:
    return result_dir / "cluster_results" / "cluster_summary.csv"


def _load_cluster_residue_counter(result_dir: Path) -> Tuple[Counter, int]:
    cluster_rows = load_csv(_cluster_summary_path(result_dir))
    counter: Counter = Counter()
    for row in cluster_rows:
        binding_a = row.get("Binding_Residues_A", "") or row.get("Binding_Residues", "")
        counter.update(parse_binding_residues(binding_a))
    return counter, len(cluster_rows)


def _load_run_metadata(
    result_dir: Path,
    *,
    receptor_id: str,
    partner_id: str,
    construct_type: str,
    orientation_validation_status: str,
) -> dict:
    meta = {}
    meta_path = result_dir / "pyrosetta_run_metadata.json"
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            meta = {}

    resolved_partner = partner_id or str(meta.get("partner_id", "")).strip()
    resolved_construct = construct_type or str(meta.get("construct_type", "")).strip() or "full_kinase_domain"

    meta["receptor_id"] = receptor_id or str(meta.get("receptor_id", "")).strip()
    meta["partner_id"] = resolved_partner
    meta["construct_type"] = resolved_construct
    meta["orientation_validation_status"] = orientation_validation_status or "not_available"
    meta["seed_index"] = meta.get("seed_index", 0)
    meta["run_dir"] = str(result_dir)
    meta["run_label"] = result_dir.name
    return meta


def _build_run_payload(
    result_dir: Path,
    *,
    receptor_id: str,
    partner_id: str,
    construct_type: str,
    orientation_validation_status: str,
) -> Optional[dict]:
    metadata = _load_run_metadata(
        result_dir,
        receptor_id=receptor_id,
        partner_id=partner_id,
        construct_type=construct_type,
        orientation_validation_status=orientation_validation_status,
    )
    residue_rows, model_rows = extract_run(result_dir, metadata)
    cluster_counter, cluster_row_count = _load_cluster_residue_counter(result_dir)

    if not residue_rows and not model_rows and cluster_row_count == 0:
        return None

    source = _source_tag(metadata.get("partner_id", ""))
    run_label = metadata.get("run_label", result_dir.name)
    run_dir_str = metadata.get("run_dir", str(result_dir))
    seed_index = _normalize_seed_index(metadata.get("seed_index", ""))
    orientation = metadata.get("orientation_validation_status", "not_available")

    enriched_residue_rows = []
    for row in residue_rows:
        enriched = dict(row)
        enriched["source"] = source
        enriched["orientation_validation_status"] = orientation
        enriched["run_label"] = run_label
        enriched["run_dir"] = run_dir_str
        enriched["seed_index"] = _normalize_seed_index(enriched.get("seed_index", seed_index))
        enriched_residue_rows.append(enriched)

    enriched_model_rows = []
    for row in model_rows:
        enriched = dict(row)
        enriched["source"] = source
        enriched["orientation_validation_status"] = orientation
        enriched["run_label"] = run_label
        enriched["run_dir"] = run_dir_str
        enriched["seed_index"] = _normalize_seed_index(enriched.get("seed_index", seed_index))
        enriched_model_rows.append(enriched)

    return {
        "metadata": metadata,
        "source": source,
        "run_label": run_label,
        "run_dir": run_dir_str,
        "seed_index": seed_index,
        "residue_rows": enriched_residue_rows,
        "model_rows": enriched_model_rows,
        "cluster_counter": cluster_counter,
        "cluster_row_count": cluster_row_count,
    }


def _payload_group_key(payload: dict) -> Tuple[str, str, str, str, str]:
    meta = payload["metadata"]
    return (
        str(meta.get("receptor_id", "")),
        str(meta.get("partner_id", "")),
        payload["source"],
        str(meta.get("construct_type", "full_kinase_domain")),
        str(meta.get("orientation_validation_status", "not_available")),
    )


def _sort_residue_rows(rows: List[dict]) -> List[dict]:
    return sorted(
        rows,
        key=lambda row: (
            row.get("receptor_id", ""),
            row.get("partner_id", ""),
            safe_float(row.get("residue_num")) or -999999,
            row.get("chain", ""),
            row.get("residue_id", ""),
        ),
    )


def _sort_long_rows(rows: List[dict]) -> List[dict]:
    return sorted(
        rows,
        key=lambda row: (
            row.get("receptor_id", ""),
            row.get("partner_id", ""),
            row.get("seed_index", ""),
            safe_float(row.get("rank")) or 999999,
            row.get("model_id", ""),
            row.get("chain", ""),
            safe_float(row.get("residue_num")) or -999999,
            row.get("residue_id", ""),
        ),
    )


def _sort_summary_rows(rows: List[dict]) -> List[dict]:
    return sorted(
        rows,
        key=lambda row: (
            row.get("receptor_id", ""),
            row.get("partner_id", ""),
            row.get("construct_type", ""),
        ),
    )


def _aggregate_run_payloads(run_payloads: List[dict]) -> Tuple[List[dict], List[dict]]:
    groups: Dict[Tuple[str, str, str, str, str], dict] = {}

    for payload in run_payloads:
        key = _payload_group_key(payload)
        receptor_id, partner_id, source, construct_type, orientation = key
        state = groups.setdefault(
            key,
            {
                "receptor_id": receptor_id,
                "partner_id": partner_id,
                "source": source,
                "construct_type": construct_type,
                "orientation_validation_status": orientation,
                "run_labels": [],
                "seed_indices": set(),
                "n_models_total": 0,
                "n_clusters_total": 0,
                "dg_values": [],
                "dsasa_values": [],
                "residue_aggs": {},
            },
        )

        state["run_labels"].append(payload["run_label"])
        if payload["seed_index"]:
            state["seed_indices"].add(payload["seed_index"])
        state["n_models_total"] += len(payload["model_rows"])
        state["n_clusters_total"] += payload["cluster_row_count"]

        for model_row in payload["model_rows"]:
            dg = safe_float(model_row.get("dG_separated"))
            if dg is not None:
                state["dg_values"].append(dg)
            dsasa = safe_float(model_row.get("dSASA"))
            if dsasa is not None:
                state["dsasa_values"].append(dsasa)

        for row in payload["residue_rows"]:
            res_key = (
                row.get("chain", ""),
                row.get("residue_id", ""),
                row.get("residue_num", ""),
                row.get("residue_name", ""),
                row.get("lobe_label", ""),
            )
            agg = state["residue_aggs"].setdefault(
                res_key,
                {
                    "chain": row.get("chain", ""),
                    "residue_id": row.get("residue_id", ""),
                    "residue_num": row.get("residue_num", ""),
                    "residue_name": row.get("residue_name", ""),
                    "lobe_label": row.get("lobe_label", ""),
                    "frequency_final_ranking": 0,
                    "frequency_cluster_summary": 0,
                    "supporting_runs": set(),
                    "supporting_seed_indices": set(),
                    "energies": [],
                },
            )
            agg["frequency_final_ranking"] += 1
            agg["supporting_runs"].add(payload["run_label"])
            if payload["seed_index"]:
                agg["supporting_seed_indices"].add(payload["seed_index"])
            delta_e = safe_float(row.get("delta_e_total"))
            if delta_e is not None:
                agg["energies"].append(delta_e)

        for residue_ref, count in payload["cluster_counter"].items():
            residue_info = _parse_receptor_residue(residue_ref)
            res_key = (
                residue_info["chain"],
                residue_info["residue_id"],
                residue_info["residue_num"] if residue_info["residue_num"] is not None else "",
                residue_info["residue_name"],
                residue_info["lobe_label"],
            )
            agg = state["residue_aggs"].setdefault(
                res_key,
                {
                    "chain": residue_info["chain"],
                    "residue_id": residue_info["residue_id"],
                    "residue_num": residue_info["residue_num"] if residue_info["residue_num"] is not None else "",
                    "residue_name": residue_info["residue_name"],
                    "lobe_label": residue_info["lobe_label"],
                    "frequency_final_ranking": 0,
                    "frequency_cluster_summary": 0,
                    "supporting_runs": set(),
                    "supporting_seed_indices": set(),
                    "energies": [],
                },
            )
            agg["frequency_cluster_summary"] += count
            agg["supporting_runs"].add(payload["run_label"])
            if payload["seed_index"]:
                agg["supporting_seed_indices"].add(payload["seed_index"])

    residue_rows: List[dict] = []
    summary_rows: List[dict] = []

    for state in groups.values():
        run_count = len(state["run_labels"])
        n_models_total = state["n_models_total"]
        residue_aggs = state["residue_aggs"]

        top_residue_pairs = sorted(
            residue_aggs.values(),
            key=lambda agg: (-agg["frequency_final_ranking"], agg["residue_id"]),
        )[:10]

        n_nlobe = sum(1 for agg in residue_aggs.values() if agg["lobe_label"] == "N-lobe")
        n_clobe = sum(1 for agg in residue_aggs.values() if agg["lobe_label"] == "C-lobe")

        for agg in residue_aggs.values():
            n_support = len(agg["supporting_runs"])
            occupancy = (
                agg["frequency_final_ranking"] / n_models_total if n_models_total > 0 else 0.0
            )
            mean_delta_e = (
                sum(agg["energies"]) / len(agg["energies"]) if agg["energies"] else None
            )
            best_delta_e = min(agg["energies"]) if agg["energies"] else None

            residue_rows.append(
                {
                    "receptor_id": state["receptor_id"],
                    "partner_id": state["partner_id"],
                    "source": state["source"],
                    "chain": agg["chain"],
                    "residue_id": agg["residue_id"],
                    "residue_num": agg["residue_num"],
                    "residue_name": agg["residue_name"],
                    "lobe_label": agg["lobe_label"],
                    "construct_type": state["construct_type"],
                    "orientation_validation_status": state["orientation_validation_status"],
                    "n_runs_total": run_count,
                    "n_runs_supporting": n_support,
                    "frac_runs_supporting": round(n_support / run_count, 4) if run_count else 0.0,
                    "supporting_seed_indices": _join_sorted(agg["supporting_seed_indices"]),
                    "frequency_final_ranking": agg["frequency_final_ranking"],
                    "frequency_cluster_summary": agg["frequency_cluster_summary"],
                    "n_models_final_ranking": n_models_total,
                    "occupancy": round(occupancy, 4),
                    "mean_interface_delta_e": round(mean_delta_e, 3) if mean_delta_e is not None else "",
                    "best_interface_delta_e": round(best_delta_e, 3) if best_delta_e is not None else "",
                }
            )

        dg_values = state["dg_values"]
        dsasa_values = state["dsasa_values"]
        summary_rows.append(
            {
                "receptor_id": state["receptor_id"],
                "partner_id": state["partner_id"],
                "source": state["source"],
                "construct_type": state["construct_type"],
                "orientation_validation_status": state["orientation_validation_status"],
                "n_runs_total": run_count,
                "n_runs_completed": run_count,
                "seed_indices": _join_sorted(state["seed_indices"]),
                "n_final_models": n_models_total,
                "n_clusters": state["n_clusters_total"],
                "n_interface_residues": len(residue_aggs),
                "n_nlobe_interface_residues": n_nlobe,
                "n_clobe_interface_residues": n_clobe,
                "top_residues": ";".join(agg["residue_id"] for agg in top_residue_pairs),
                "best_dg": round(min(dg_values), 3) if dg_values else "",
                "mean_dg": round(sum(dg_values) / len(dg_values), 3) if dg_values else "",
                "best_dsasa": round(max(dsasa_values), 1) if dsasa_values else "",
            }
        )

    return _sort_residue_rows(residue_rows), _sort_summary_rows(summary_rows)


def extract_pyrosetta_interface_residues(
    result_dir: Path,
    receptor_id: str,
    partner_id: str = "",
    construct_type: str = "full_kinase_domain",
    orientation_validation_status: str = "not_available",
) -> Dict[str, object]:
    """Extract receptor-side interface residues from a single PyRosetta run dir."""
    payload = _build_run_payload(
        result_dir,
        receptor_id=receptor_id,
        partner_id=partner_id,
        construct_type=construct_type,
        orientation_validation_status=orientation_validation_status,
    )
    if payload is None:
        return {"residue_rows": [], "summary": {}}

    residue_rows, summary_rows = _aggregate_run_payloads([payload])
    summary = summary_rows[0] if summary_rows else {}
    return {
        "residue_rows": residue_rows,
        "summary": summary,
        "residue_long_rows": _sort_long_rows(payload["residue_rows"]),
        "model_rows": _sort_long_rows(payload["model_rows"]),
    }


# ---------------------------------------------------------------------------
# Batch extraction from project config
# ---------------------------------------------------------------------------


def extract_pyrosetta_batch(
    config_path: str,
    output_dir: Optional[str] = None,
    pyrosetta_result_dirs: Optional[Dict[str, object]] = None,
) -> Tuple[Path, Path]:
    """Extract PyRosetta PPI evidence for all receptors in config."""
    config = load_config(config_path)
    out_root = Path(output_dir) if output_dir else paths.wa_phase3_ppi_postprocess(config)

    ppi_config = config.get("ppi", {})
    pyrosetta_dirs = pyrosetta_result_dirs or ppi_config.get("pyrosetta_result_dirs", {})

    run_payloads: List[dict] = []

    for receptor in config.get("receptors", []):
        receptor_id = receptor["id"]
        entries = _normalize_result_dirs(pyrosetta_dirs.get(receptor_id))
        if not entries:
            continue

        for entry in entries:
            result_dir = Path(entry["path"])
            partner = str(entry.get("partner", "")).strip()
            construct_type = str(entry.get("construct_type", "full_kinase_domain")).strip() or "full_kinase_domain"
            orientation_validation_status = str(
                entry.get("orientation_validation_status", "not_available")
            ).strip() or "not_available"

            if not result_dir.exists():
                print(
                    f"[WARN] PyRosetta result dir not found for {receptor_id}"
                    f"{f' ({partner})' if partner else ''}: {result_dir}"
                )
                continue

            payload = _build_run_payload(
                result_dir,
                receptor_id=receptor_id,
                partner_id=partner,
                construct_type=construct_type,
                orientation_validation_status=orientation_validation_status,
            )
            if payload is None:
                print(
                    f"[WARN] PyRosetta result dir has no extractable final_ranking.csv for "
                    f"{receptor_id}{f' ({partner})' if partner else ''}: {result_dir}"
                )
                continue
            run_payloads.append(payload)

    residue_rows, summary_rows = _aggregate_run_payloads(run_payloads)
    residue_long_rows = _sort_long_rows(
        [row for payload in run_payloads for row in payload["residue_rows"]]
    )
    model_rows = _sort_long_rows(
        [row for payload in run_payloads for row in payload["model_rows"]]
    )

    residue_csv = write_csv_rows(
        out_root / "ppi_pyrosetta_residues.csv",
        residue_rows,
        PPI_RESIDUE_FIELDS,
    )
    summary_csv = write_csv_rows(
        out_root / "ppi_pyrosetta_summary.csv",
        summary_rows,
        PPI_SUMMARY_FIELDS,
    )
    write_csv_rows(
        out_root / "ppi_pyrosetta_residue_long.csv",
        residue_long_rows,
        PPI_RESIDUE_LONG_FIELDS,
    )
    write_csv_rows(
        out_root / "ppi_pyrosetta_model_table.csv",
        model_rows,
        PPI_MODEL_FIELDS,
    )
    return residue_csv, summary_csv
