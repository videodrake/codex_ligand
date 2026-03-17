"""Helpers for PyRosetta run traceability and score exports."""

from __future__ import annotations

import csv
import configparser
import json
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Set


RECEPTOR_PATTERNS = (
    "3GT8_raw",
    "EGFR_160-185",
    "EGFR_170-200",
)


def _get_metadata_value(
    config: configparser.ConfigParser, key: str, fallback: Optional[str] = None
) -> Optional[str]:
    if not config.has_section("Metadata"):
        return fallback
    value = config.get("Metadata", key, fallback=fallback)
    if value is None:
        return None
    value = value.strip()
    return value or fallback


def infer_receptor_id(config: configparser.ConfigParser, input_pdb: str) -> str:
    explicit = _get_metadata_value(config, "receptor_id")
    if explicit:
        return explicit

    stem = Path(input_pdb).stem
    for receptor_id in RECEPTOR_PATTERNS:
        if receptor_id.lower() in stem.lower():
            return receptor_id
    return "unknown_receptor"


def infer_partner_id(config: configparser.ConfigParser, input_pdb: str) -> str:
    explicit = _get_metadata_value(config, "partner_id")
    if explicit:
        return explicit

    stem = Path(input_pdb).stem.lower()
    if "beta_meander" in stem:
        return "MYO1D_beta_meander"
    if "th1" in stem:
        return "MYO1D_TH1"
    return Path(input_pdb).stem


def infer_construct_type(config: configparser.ConfigParser, input_pdb: str) -> str:
    explicit = _get_metadata_value(config, "construct_type")
    if explicit:
        return explicit

    stem = Path(input_pdb).stem.lower()
    if "clobe" in stem or "c_lobe" in stem or "fragment" in stem:
        return "legacy_clobe_fragment"
    if "full_kinase" in stem or "kinase_domain" in stem:
        return "full_kinase_domain"
    return "unknown_construct_type"


def infer_partner_construct(config: configparser.ConfigParser, input_pdb: str) -> str:
    explicit = _get_metadata_value(config, "partner_construct")
    if explicit:
        return explicit

    stem = Path(input_pdb).stem.lower()
    if "beta_meander" in stem:
        return "beta_meander"
    if "th1" in stem:
        return "TH1_domain"
    return "unknown_partner_construct"


def infer_numbering_system(config: configparser.ConfigParser) -> str:
    return _get_metadata_value(config, "numbering_system", "unspecified")


def _parse_chain_ids(value: object) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(item).strip() for item in value if str(item).strip()]
    return [part.strip() for part in str(value).split(",") if part.strip()]


def _sanitize_output_label(value: object) -> str:
    normalized = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip()).strip("_")
    return normalized or "unknown"


def build_output_root_name(
    config: configparser.ConfigParser, input_pdb: str, filename_stem: str
) -> str:
    mode = config.get("Output", "root_dir_mode", fallback="legacy_stem")
    if mode != "metadata_tagged":
        return filename_stem

    run_label = config.get("Output", "run_label", fallback="run")
    parts = [
        filename_stem,
        infer_receptor_id(config, input_pdb),
        infer_partner_construct(config, input_pdb),
        infer_construct_type(config, input_pdb),
        run_label,
    ]
    return "__".join(_sanitize_output_label(part) for part in parts)


def _summarize_pdb_chains(input_pdb: str) -> Dict[str, Dict[str, int]]:
    chain_residues: Dict[str, Set[int]] = {}
    with open(input_pdb, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if len(line) < 26:
                continue
            chain_id = line[21].strip() or "?"
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                continue
            chain_residues.setdefault(chain_id, set()).add(resnum)

    summary: Dict[str, Dict[str, int]] = {}
    for chain_id, residues in sorted(chain_residues.items()):
        summary[chain_id] = {
            "min_resnum": min(residues),
            "max_resnum": max(residues),
            "n_unique_residues": len(residues),
        }
    return summary


def _load_sidecar_info(input_pdb: str) -> Dict[str, object]:
    pdb_path = Path(input_pdb)
    info_path = pdb_path.with_name(f"{pdb_path.stem}_info.json")
    mapping_path = pdb_path.with_name(f"{pdb_path.stem}_mapping.csv")

    sidecar: Dict[str, object] = {
        "info_json_path": str(info_path.resolve()),
        "info_json_exists": info_path.exists(),
        "mapping_csv_path": str(mapping_path.resolve()),
        "mapping_csv_exists": mapping_path.exists(),
    }

    if info_path.exists():
        with open(info_path, "r", encoding="utf-8") as handle:
            sidecar["info_json"] = json.load(handle)

    if mapping_path.exists():
        mapping_summary: Dict[str, Dict[str, object]] = {}
        with open(mapping_path, "r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                new_chain = row.get("new_chain", "").strip() or "?"
                original_chain = row.get("original_chain", "").strip() or "?"
                try:
                    new_resnum = int(row.get("new_resnum", ""))
                    original_resnum = int(row.get("original_resnum", ""))
                except ValueError:
                    continue

                chain_summary = mapping_summary.setdefault(
                    new_chain,
                    {
                        "new_min_resnum": new_resnum,
                        "new_max_resnum": new_resnum,
                        "original_min_resnum": original_resnum,
                        "original_max_resnum": original_resnum,
                        "original_chain_ids": set(),
                    },
                )
                chain_summary["new_min_resnum"] = min(chain_summary["new_min_resnum"], new_resnum)
                chain_summary["new_max_resnum"] = max(chain_summary["new_max_resnum"], new_resnum)
                chain_summary["original_min_resnum"] = min(
                    chain_summary["original_min_resnum"], original_resnum
                )
                chain_summary["original_max_resnum"] = max(
                    chain_summary["original_max_resnum"], original_resnum
                )
                chain_summary["original_chain_ids"].add(original_chain)

        normalized_summary: Dict[str, Dict[str, object]] = {}
        for chain_id, chain_summary in sorted(mapping_summary.items()):
            normalized_summary[chain_id] = {
                "new_min_resnum": chain_summary["new_min_resnum"],
                "new_max_resnum": chain_summary["new_max_resnum"],
                "original_min_resnum": chain_summary["original_min_resnum"],
                "original_max_resnum": chain_summary["original_max_resnum"],
                "original_chain_ids": sorted(chain_summary["original_chain_ids"]),
            }
        sidecar["mapping_summary"] = normalized_summary

    return sidecar


def build_run_metadata(
    *,
    config: configparser.ConfigParser,
    config_file: str,
    input_pdb: str,
    root_dir: str,
    filename: str,
    requested_cpus: int,
    used_cpus: int,
    master_seed: Optional[int],
    dir_filter: str,
    dir_cluster: str,
    dir_final: str,
) -> Dict[str, object]:
    construct_type = infer_construct_type(config, input_pdb)
    return {
        "config_file": os.path.abspath(config_file),
        "config_basename": os.path.basename(config_file),
        "input_pdb": os.path.abspath(input_pdb),
        "input_pdb_basename": filename,
        "root_dir": os.path.abspath(root_dir),
        "receptor_id": infer_receptor_id(config, input_pdb),
        "partner_id": infer_partner_id(config, input_pdb),
        "construct_type": construct_type,
        "receptor_construct": _get_metadata_value(
            config, "receptor_construct", construct_type
        ),
        "partner_construct": infer_partner_construct(config, input_pdb),
        "receptor_chain_ids": _get_metadata_value(config, "receptor_chain_ids", "A"),
        "partner_chain_ids": _get_metadata_value(config, "partner_chain_ids", "B"),
        "numbering_system": infer_numbering_system(config),
        "n_cpus_requested": requested_cpus,
        "n_cpus_used": used_cpus,
        "random_seed": master_seed,
        "root_dir_mode": config.get("Output", "root_dir_mode", fallback="legacy_stem"),
        "run_label": config.get("Output", "run_label", fallback="run"),
        "total_global_models": config.getint("Docking", "total_global_models"),
        "save_top_n": config.getint("Output", "save_top_n"),
        "save_filter_max": config.getint("Output", "save_filter_max", fallback=1000),
        "cluster_threshold": config.get("Cluster", "threshold", fallback="auto"),
        "cluster_top_n": config.get("Cluster", "cluster_top_n", fallback="auto"),
        "output_dirs": {
            "filter_passed": os.path.abspath(dir_filter),
            "cluster_results": os.path.abspath(dir_cluster),
            "final_result": os.path.abspath(dir_final),
        },
        "run_status": "initialized",
    }


def build_input_validation_report(
    *,
    config: configparser.ConfigParser,
    input_pdb: str,
    run_metadata: Dict[str, object],
) -> Dict[str, object]:
    input_path = Path(input_pdb)
    receptor_chain_ids = _parse_chain_ids(run_metadata.get("receptor_chain_ids"))
    partner_chain_ids = _parse_chain_ids(run_metadata.get("partner_chain_ids"))
    warnings: List[str] = []

    report: Dict[str, object] = {
        "input_pdb": str(input_path.resolve()),
        "input_pdb_exists": input_path.exists(),
        "receptor_id": run_metadata.get("receptor_id", "unknown_receptor"),
        "partner_id": run_metadata.get("partner_id", "unknown_partner"),
        "construct_type": run_metadata.get("construct_type", "unknown_construct_type"),
        "receptor_construct": run_metadata.get("receptor_construct", ""),
        "partner_construct": run_metadata.get("partner_construct", ""),
        "numbering_system": run_metadata.get("numbering_system", "unspecified"),
        "expected_receptor_chain_ids": receptor_chain_ids,
        "expected_partner_chain_ids": partner_chain_ids,
        "checks": {},
        "warnings": warnings,
    }

    if not input_path.exists():
        report["status"] = "fail"
        report["checks"] = {"input_pdb_exists": False}
        warnings.append("Configured input PDB does not exist.")
        return report

    observed_chains = _summarize_pdb_chains(input_pdb)
    observed_chain_ids = sorted(observed_chains.keys())
    report["observed_chains"] = observed_chains

    receptor_chains_present = all(chain_id in observed_chains for chain_id in receptor_chain_ids)
    partner_chains_present = all(chain_id in observed_chains for chain_id in partner_chain_ids)

    sidecar = _load_sidecar_info(input_pdb)
    report["sidecar"] = sidecar

    checks = {
        "input_pdb_exists": True,
        "observed_chain_count": len(observed_chain_ids),
        "expected_receptor_chains_present": receptor_chains_present,
        "expected_partner_chains_present": partner_chains_present,
        "info_json_exists": sidecar["info_json_exists"],
        "mapping_csv_exists": sidecar["mapping_csv_exists"],
        "numbering_system_declared": report["numbering_system"] != "unspecified",
    }
    report["checks"] = checks

    if not receptor_chains_present:
        warnings.append(
            f"Configured receptor chain IDs {receptor_chain_ids} were not fully found in input PDB chains {observed_chain_ids}."
        )
    if not partner_chains_present:
        warnings.append(
            f"Configured partner chain IDs {partner_chain_ids} were not fully found in input PDB chains {observed_chain_ids}."
        )
    if not sidecar["info_json_exists"]:
        warnings.append("Prepared input sidecar info JSON is missing.")
    if not sidecar["mapping_csv_exists"]:
        warnings.append("Prepared input mapping CSV is missing.")
    if report["numbering_system"] == "unspecified":
        warnings.append("Numbering system is not explicitly declared in config metadata.")

    partner_construct = str(report["partner_construct"])
    if "legacy_beta_meander" in partner_construct:
        partner_start = None
        for chain_id in partner_chain_ids:
            if chain_id in observed_chains:
                partner_start = observed_chains[chain_id]["min_resnum"]
                break
        warnings.append(
            "Partner input is still a legacy beta-meander construct"
            + (f" starting at residue {partner_start}" if partner_start is not None else "")
            + "; extended ~955-start input has not been wired in yet."
        )

    report["status"] = "warn" if warnings else "pass"
    return report


def build_input_validation_markdown(report: Dict[str, object]) -> str:
    """Render a small human-readable summary for pre-qsub inspection."""
    lines = [
        "# Phase 1 Input Validation Summary",
        "",
        f"- Status: {report.get('status', 'unknown')}",
        f"- Input PDB: {report.get('input_pdb', '')}",
        f"- Receptor ID: {report.get('receptor_id', '')}",
        f"- Partner ID: {report.get('partner_id', '')}",
        f"- Receptor Construct: {report.get('receptor_construct', '')}",
        f"- Partner Construct: {report.get('partner_construct', '')}",
        f"- Numbering System: {report.get('numbering_system', '')}",
        "",
        "## Checks",
    ]

    checks = report.get("checks", {})
    if isinstance(checks, dict) and checks:
        for key in sorted(checks.keys()):
            lines.append(f"- {key}: {checks[key]}")
    else:
        lines.append("- No checks recorded")

    observed_chains = report.get("observed_chains", {})
    lines.extend(["", "## Observed Chains"])
    if isinstance(observed_chains, dict) and observed_chains:
        for chain_id in sorted(observed_chains.keys()):
            chain_summary = observed_chains[chain_id]
            if isinstance(chain_summary, dict):
                lines.append(
                    "- Chain "
                    f"{chain_id}: "
                    f"{chain_summary.get('min_resnum', '')}-"
                    f"{chain_summary.get('max_resnum', '')} "
                    f"({chain_summary.get('n_unique_residues', '')} residues)"
                )
    else:
        lines.append("- No chain summary recorded")

    sidecar = report.get("sidecar", {})
    lines.extend(["", "## Sidecar Files"])
    if isinstance(sidecar, dict):
        lines.append(f"- info_json_exists: {sidecar.get('info_json_exists', False)}")
        lines.append(f"- mapping_csv_exists: {sidecar.get('mapping_csv_exists', False)}")
        mapping_summary = sidecar.get("mapping_summary", {})
        if isinstance(mapping_summary, dict) and mapping_summary:
            lines.append("")
            lines.append("### Mapping Summary")
            for chain_id in sorted(mapping_summary.keys()):
                chain_summary = mapping_summary[chain_id]
                if isinstance(chain_summary, dict):
                    lines.append(
                        "- Chain "
                        f"{chain_id}: new "
                        f"{chain_summary.get('new_min_resnum', '')}-"
                        f"{chain_summary.get('new_max_resnum', '')}, original "
                        f"{chain_summary.get('original_min_resnum', '')}-"
                        f"{chain_summary.get('original_max_resnum', '')}, "
                        f"source chains {chain_summary.get('original_chain_ids', [])}"
                    )
    else:
        lines.append("- No sidecar summary recorded")

    warnings = report.get("warnings", [])
    lines.extend(["", "## Warnings"])
    if isinstance(warnings, list) and warnings:
        lines.extend(f"- {warning}" for warning in warnings)
    else:
        lines.append("- None")

    return "\n".join(lines) + "\n"


def build_decoy_score_rows(
    fully_scored: List[Dict[str, object]], run_metadata: Dict[str, object]
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for item in fully_scored:
        rows.append(
            {
                "decoy_id": item.get("id", ""),
                "source_file": item.get("File", ""),
                "receptor_id": run_metadata.get("receptor_id", "unknown_receptor"),
                "partner_id": run_metadata.get("partner_id", "unknown_partner"),
                "construct_type": run_metadata.get("construct_type", "unknown_construct_type"),
                "receptor_construct": run_metadata.get("receptor_construct", ""),
                "partner_construct": run_metadata.get("partner_construct", ""),
                "receptor_chain_ids": run_metadata.get("receptor_chain_ids", ""),
                "partner_chain_ids": run_metadata.get("partner_chain_ids", ""),
                "total_score": item.get("Total_Score", item.get("total_score", "")),
                "I_sc": item.get("dG_separated", ""),
                "dG_separated": item.get("dG_separated", ""),
                "dSASA": item.get("dSASA", ""),
                "sc_value": item.get("sc_value", ""),
                "packstat": item.get("packstat", ""),
                "delta_unsatHbonds": item.get("delta_unsatHbonds", ""),
                "nres_int": item.get("nres_int", ""),
                "hbonds_int": item.get("hbonds_int", ""),
                "L_RMSD": item.get("L_RMSD", ""),
                "center_x": item.get("center_x", ""),
                "center_y": item.get("center_y", ""),
                "center_z": item.get("center_z", ""),
                "binding_residues_A": item.get("Binding_Residues_A", ""),
                "binding_residues_B": item.get("Binding_Residues_B", ""),
            }
        )
    return rows
