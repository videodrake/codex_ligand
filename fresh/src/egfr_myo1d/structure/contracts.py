"""Structural input contracts and EGFR receptor QC."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from egfr_myo1d.io.hashing import describe_file
from egfr_myo1d.structure.pdb_parser import PDBParseError, PDBStructure, parse_pdb


class ContractError(ValueError):
    """Raised for malformed structural contract files."""


CAP_OR_MODIFIED_RESNAMES = {"ACE", "NME", "MSE", "SEP", "TPO", "PTR"}


def default_contract() -> dict[str, Any]:
    return {
        "receptors": [
            {
                "receptor_id": "EGFR_160-185",
                "path": "EGFR_160-185.pdb",
                "expected_assembly": "dimer",
                "expected_protomers": [
                    {"protomer_id": "A", "chain_ids": ["A"], "required_ranges": [[699, 700]]},
                    {"protomer_id": "B", "chain_ids": ["B"], "required_ranges": [[699, 700]]},
                ],
                "numbering_system": "project_uniprot_or_documented_mapping",
                "allow_single_chain_dimer_with_mapping": True,
                "known_mutation_checks": [
                    {
                        "label": "3GT8_V924R",
                        "residue_number": 924,
                        "expected_wt": "VAL",
                        "warn_if_resname": "ARG",
                    }
                ],
            }
        ],
        "membrane_frame": "membrane_frame.json",
        "myo1d_partner": "AF-O94832-F1-model_v6.pdb",
    }


def load_contract(input_root: Path) -> tuple[dict[str, Any], Path | None]:
    contract_path = Path(input_root) / "contract.json"
    if not contract_path.exists():
        return default_contract(), None
    try:
        data = json.loads(contract_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ContractError(f"malformed JSON contract: {contract_path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ContractError("contract root must be a JSON object")
    if "receptors" in data and not isinstance(data["receptors"], list):
        raise ContractError("contract field 'receptors' must be a list")
    return data, contract_path


def resolve_contract_path(input_root: Path, relative_path: str) -> Path:
    candidate = (Path(input_root) / relative_path).resolve()
    root = Path(input_root).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise ContractError(f"contract path escapes input root: {relative_path}") from exc
    return candidate


def _missing_input_status(mode: str, profile: str) -> str:
    if mode == "smoke_env" and profile == "codex_dev":
        return "WARN"
    return "FAIL"


def _required_range_missing_status(mode: str, profile: str) -> str:
    if mode == "smoke_env" and profile == "codex_dev":
        return "WARN"
    return "FAIL"


def _residue_range_present(numbers: set[int], start: int, end: int) -> tuple[bool, list[int]]:
    required = set(range(start, end + 1))
    missing = sorted(required - numbers)
    return not missing, missing


def _chain_residue_summary(structure: PDBStructure, structure_id: str, path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for chain_id in structure.chain_ids:
        numbers = sorted(structure.residue_numbers_for_chain(chain_id))
        jumps = [
            f"{left}->{right}"
            for left, right in zip(numbers, numbers[1:])
            if right - left > 50
        ]
        reset_risk = bool(numbers and min(numbers) <= 5 and max(numbers) < 200)
        rows.append(
            {
                "structure_id": structure_id,
                "path": str(path),
                "chain_id": chain_id,
                "min_residue": min(numbers) if numbers else "",
                "max_residue": max(numbers) if numbers else "",
                "residue_count": len(numbers),
                "large_jumps": ";".join(jumps),
                "residue_reset_risk": reset_risk,
                "status": "WARN" if jumps or reset_risk else "PASS",
                "message": "large residue-number jumps or reset risk detected" if jumps or reset_risk else "",
            }
        )
    return rows


def validate_receptor_contract(
    receptor: dict[str, Any],
    input_root: Path,
    mode: str,
    profile: str,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    receptor_id = receptor.get("receptor_id", "unnamed_receptor")
    rel_path = receptor.get("path")
    if not isinstance(rel_path, str):
        raise ContractError(f"receptor {receptor_id} is missing string field 'path'")
    path = resolve_contract_path(input_root, rel_path)
    file_desc = describe_file(path)
    manifest_record = {
        "input_id": receptor_id,
        "input_type": "egfr_receptor",
        "path": rel_path,
        "absolute_path": str(path),
        **file_desc,
    }

    chain_rows: list[dict[str, Any]] = []
    residue_rows: list[dict[str, Any]] = []
    messages: list[dict[str, Any]] = []
    report: dict[str, Any] = {
        "receptor_id": receptor_id,
        "path": rel_path,
        "status": "PASS",
        "observed_chain_ids": [],
        "atom_count": 0,
        "residue_count": 0,
        "messages": messages,
    }

    if not path.exists():
        status = _missing_input_status(mode, profile)
        messages.append({"severity": status, "code": "missing_receptor", "message": f"missing receptor input {rel_path}"})
        report["status"] = status
        manifest_record["parse_status"] = "missing"
        manifest_record["status"] = status
        return report, chain_rows, residue_rows, manifest_record

    try:
        structure = parse_pdb(path)
    except PDBParseError as exc:
        messages.append({"severity": "FAIL", "code": "pdb_parse_error", "message": str(exc)})
        report["status"] = "FAIL"
        manifest_record["parse_status"] = "failed"
        manifest_record["status"] = "FAIL"
        return report, chain_rows, residue_rows, manifest_record

    observed_chains = list(structure.chain_ids)
    report["observed_chain_ids"] = observed_chains
    report["atom_count"] = len(structure.atoms)
    report["residue_count"] = len(structure.residues)
    report["duplicate_atom_identity_count"] = structure.duplicate_atom_identity_count()
    manifest_record["parse_status"] = "parsed"

    expected_protomers = receptor.get("expected_protomers", [])
    allow_single_chain = bool(receptor.get("allow_single_chain_dimer_with_mapping", False))
    expected_assembly = receptor.get("expected_assembly", "")

    if expected_assembly == "dimer" and len(observed_chains) < 2:
        severity = "WARN" if allow_single_chain else _required_range_missing_status(mode, profile)
        messages.append(
            {
                "severity": severity,
                "code": "single_chain_dimer_ambiguity",
                "message": "dimer contract observed fewer than two chain IDs; mapping/normalization required later",
            }
        )

    if structure.duplicate_atom_identity_count() > 0:
        messages.append(
            {
                "severity": "WARN",
                "code": "duplicate_atom_identities",
                "message": "duplicate atom identities detected; possible same-chain duplicate dimer input",
            }
        )

    for protomer in expected_protomers:
        protomer_id = protomer.get("protomer_id", "")
        expected_chain_ids = protomer.get("chain_ids", [])
        present_chains = [chain_id for chain_id in expected_chain_ids if chain_id in observed_chains]
        row_status = "PASS" if present_chains else "WARN"
        row_message = ""
        if not present_chains:
            row_message = "expected chain IDs not observed"
            severity = "WARN" if allow_single_chain and len(observed_chains) == 1 else _required_range_missing_status(mode, profile)
            messages.append(
                {
                    "severity": severity,
                    "code": "expected_protomer_chain_missing",
                    "message": f"protomer {protomer_id} expected chains {expected_chain_ids} but observed {observed_chains}",
                }
            )

        chain_rows.append(
            {
                "receptor_id": receptor_id,
                "path": rel_path,
                "expected_assembly": expected_assembly,
                "protomer_id": protomer_id,
                "expected_chain_ids": ";".join(expected_chain_ids),
                "observed_chain_ids": ";".join(observed_chains),
                "status": row_status,
                "message": row_message,
            }
        )

        chains_for_range = present_chains or (observed_chains if allow_single_chain and len(observed_chains) == 1 else [])
        for chain_id in chains_for_range:
            numbers = structure.residue_numbers_for_chain(chain_id)
            for start, end in protomer.get("required_ranges", []):
                present, missing = _residue_range_present(numbers, int(start), int(end))
                if not present:
                    severity = _required_range_missing_status(mode, profile)
                    messages.append(
                        {
                            "severity": severity,
                            "code": "required_residue_range_missing",
                            "message": f"{receptor_id} chain {chain_id} missing residues in range {start}-{end}: {missing[:10]}",
                        }
                    )
                if numbers and min(numbers) <= 5 and max(numbers) < int(start):
                    severity = _required_range_missing_status(mode, profile)
                    messages.append(
                        {
                            "severity": severity,
                            "code": "residue_number_reset_risk",
                            "message": f"{receptor_id} chain {chain_id} looks renumbered {min(numbers)}-{max(numbers)} while expected {start}-{end}",
                        }
                    )

    for mutation_check in receptor.get("known_mutation_checks", []):
        residue_number = int(mutation_check["residue_number"])
        warn_if_resname = mutation_check.get("warn_if_resname")
        observed_resnames = structure.residue_resnames(residue_number)
        if warn_if_resname in observed_resnames:
            messages.append(
                {
                    "severity": "WARN",
                    "code": mutation_check.get("label", "known_mutation"),
                    "message": f"residue {residue_number} contains {warn_if_resname}; do not treat as normal WT evidence",
                }
            )

    cap_residues = [
        residue
        for residue in structure.residues
        if "HETATM" in residue.record_types_present and residue.resname in CAP_OR_MODIFIED_RESNAMES
    ]
    report["hetatm_residue_count"] = sum(
        1 for residue in structure.residues if "HETATM" in residue.record_types_present
    )
    report["cap_or_modified_residue_count"] = len(cap_residues)
    if cap_residues:
        messages.append(
            {
                "severity": "WARN",
                "code": "hetatm_cap_or_modified_residue_present",
                "message": "HETATM cap or modified residues are present and included in QC counts",
            }
        )

    residue_rows.extend(_chain_residue_summary(structure, receptor_id, path))
    for row in residue_rows:
        if row["structure_id"] == receptor_id and row["status"] == "WARN":
            messages.append(
                {
                    "severity": "WARN",
                    "code": "residue_mapping_audit_warning",
                    "message": f"chain {row['chain_id']}: {row['message']}",
                }
            )

    if any(message["severity"] == "FAIL" for message in messages):
        report["status"] = "FAIL"
    elif messages:
        report["status"] = "WARN"
    manifest_record["status"] = report["status"]
    return report, chain_rows, residue_rows, manifest_record

