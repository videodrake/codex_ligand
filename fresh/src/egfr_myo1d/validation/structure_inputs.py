"""Task 3 structural input contract validation."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_phase_status
from egfr_myo1d.core.manifest import initialize_manifests, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.io.hashing import describe_file
from egfr_myo1d.structure.contracts import (
    ContractError,
    load_contract,
    resolve_contract_path,
    validate_receptor_contract,
)
from egfr_myo1d.structure.geometry import validate_membrane_frame
from egfr_myo1d.structure.myo1d_annotation import validate_myo1d_file
from egfr_myo1d.structure.pdb_parser import PDBParseError


def combined_status(statuses: list[str]) -> str:
    if "FAIL" in statuses:
        return "FAIL"
    if "WARN" in statuses:
        return "WARN"
    return "PASS"


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str], ctx: RunContext) -> None:
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    with safe_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _missing_status(mode: str, profile: str) -> str:
    return "WARN" if mode == "smoke_env" and profile == "codex_dev" else "FAIL"


def validate_myo1d_contract(
    contract: dict[str, Any],
    input_root: Path,
    mode: str,
    profile: str,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    myo1d_entry = contract.get("myo1d_partner")
    if isinstance(myo1d_entry, dict):
        rel_path = myo1d_entry.get("path")
    else:
        rel_path = myo1d_entry
    if not isinstance(rel_path, str):
        return (
            {
                "status": _missing_status(mode, profile),
                "messages": [
                    {
                        "severity": _missing_status(mode, profile),
                        "code": "myo1d_partner_missing_from_contract",
                        "message": "contract does not define myo1d_partner",
                    }
                ],
            },
            None,
        )

    path = resolve_contract_path(input_root, rel_path)
    manifest = {
        "input_id": "MYO1D_partner",
        "input_type": "myo1d_partner",
        "path": rel_path,
        "absolute_path": str(path),
        **describe_file(path),
    }
    if not path.exists():
        status = _missing_status(mode, profile)
        manifest["status"] = status
        manifest["parse_status"] = "missing"
        return (
            {
                "input_path": rel_path,
                "status": status,
                "messages": [
                    {
                        "severity": status,
                        "code": "missing_myo1d_partner",
                        "message": f"missing MYO1D partner input {rel_path}",
                    }
                ],
            },
            manifest,
        )
    try:
        report = validate_myo1d_file(path, mode, profile)
    except PDBParseError as exc:
        report = {
            "input_path": rel_path,
            "status": "FAIL",
            "messages": [{"severity": "FAIL", "code": "pdb_parse_error", "message": str(exc)}],
        }
        manifest["parse_status"] = "failed"
        manifest["status"] = "FAIL"
        return report, manifest

    report["input_path"] = rel_path
    manifest["parse_status"] = "parsed"
    manifest["status"] = report["status"]
    return report, manifest


def validate_membrane_contract(
    contract: dict[str, Any],
    input_root: Path,
    mode: str,
    profile: str,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    rel_path = contract.get("membrane_frame")
    if not isinstance(rel_path, str):
        status = _missing_status(mode, profile)
        return (
            {
                "status": status,
                "messages": [
                    {
                        "severity": status,
                        "code": "membrane_frame_missing_from_contract",
                        "message": "contract does not define membrane_frame",
                    }
                ],
            },
            None,
        )

    path = resolve_contract_path(input_root, rel_path)
    manifest = {
        "input_id": "membrane_frame",
        "input_type": "membrane_frame",
        "path": rel_path,
        "absolute_path": str(path),
        **describe_file(path),
    }
    if not path.exists():
        status = _missing_status(mode, profile)
        manifest["status"] = status
        manifest["parse_status"] = "missing"
        return (
            {
                "input_path": rel_path,
                "status": status,
                "messages": [
                    {
                        "severity": status,
                        "code": "missing_membrane_frame",
                        "message": f"missing membrane-frame JSON {rel_path}",
                    }
                ],
            },
            manifest,
        )

    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        manifest["status"] = "FAIL"
        manifest["parse_status"] = "failed"
        return (
            {
                "input_path": rel_path,
                "status": "FAIL",
                "messages": [
                    {
                        "severity": "FAIL",
                        "code": "malformed_membrane_frame_json",
                        "message": str(exc),
                    }
                ],
            },
            manifest,
        )

    report = validate_membrane_frame(data)
    report["input_path"] = rel_path
    manifest["status"] = report["status"]
    manifest["parse_status"] = "parsed"
    return report, manifest


def _failure_outputs(
    ctx: RunContext,
    mode: str,
    profile: str,
    input_root: Path,
    message: str,
) -> dict[str, Any]:
    manifest = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "mode": mode,
        "profile": profile,
        "input_root": str(input_root),
        "status": "FAIL",
        "inputs": [],
        "messages": [{"severity": "FAIL", "code": "contract_error", "message": message}],
    }
    report = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "status": "FAIL",
        "messages": manifest["messages"],
        "receptors": [],
        "myo1d": None,
        "membrane_frame": None,
    }
    write_json(ctx.manifest_dir / "structure_input_manifest.json", manifest, ctx)
    write_json(ctx.manifest_dir / "structure_qc_report.json", report, ctx)
    write_json(ctx.qc_dir / "membrane_frame_report.json", {"status": "FAIL", "messages": manifest["messages"]}, ctx)
    write_json(ctx.qc_dir / "myo1d_annotation_report.json", {"status": "FAIL", "messages": manifest["messages"]}, ctx)
    write_csv(ctx.qc_dir / "chain_protomer_audit.csv", [], [
        "receptor_id", "path", "expected_assembly", "protomer_id",
        "expected_chain_ids", "observed_chain_ids", "status", "message",
    ], ctx)
    write_csv(ctx.qc_dir / "residue_mapping_audit.csv", [], [
        "structure_id", "path", "chain_id", "min_residue", "max_residue",
        "residue_count", "large_jumps", "residue_reset_risk", "status", "message",
    ], ctx)
    append_phase_status(ctx, "structure_inputs", "FAIL", message, {"profile": profile, "mode": mode})
    return report


def validate_structure_inputs(
    ctx: RunContext,
    mode: str,
    profile: str,
    input_root: Path,
) -> dict[str, Any]:
    initialize_manifests(ctx, mode)
    root = Path(input_root)
    if not root.is_absolute():
        root = (ctx.repo_root / root).resolve()

    if not root.exists():
        status = _missing_status(mode, profile)
        if status == "FAIL":
            return _failure_outputs(ctx, mode, profile, root, f"missing input root: {root}")

    try:
        contract, contract_path = load_contract(root)
    except ContractError as exc:
        return _failure_outputs(ctx, mode, profile, root, str(exc))

    receptor_reports: list[dict[str, Any]] = []
    chain_rows: list[dict[str, Any]] = []
    residue_rows: list[dict[str, Any]] = []
    input_records: list[dict[str, Any]] = []
    messages: list[dict[str, Any]] = []

    if contract_path is not None:
        contract_desc = describe_file(contract_path)
        input_records.append(
            {
                "input_id": "contract",
                "input_type": "contract",
                "path": contract_path.name,
                "absolute_path": str(contract_path),
                "status": "PASS",
                "parse_status": "parsed",
                **contract_desc,
            }
        )

    for receptor in contract.get("receptors", []):
        report, receptor_chain_rows, receptor_residue_rows, manifest_record = validate_receptor_contract(
            receptor, root, mode, profile
        )
        receptor_reports.append(report)
        chain_rows.extend(receptor_chain_rows)
        residue_rows.extend(receptor_residue_rows)
        input_records.append(manifest_record)
        messages.extend(report.get("messages", []))

    myo1d_report, myo1d_manifest = validate_myo1d_contract(contract, root, mode, profile)
    if myo1d_manifest:
        input_records.append(myo1d_manifest)
    messages.extend(myo1d_report.get("messages", []))

    membrane_report, membrane_manifest = validate_membrane_contract(contract, root, mode, profile)
    if membrane_manifest:
        input_records.append(membrane_manifest)
    messages.extend(membrane_report.get("messages", []))

    status = combined_status(
        [report.get("status", "PASS") for report in receptor_reports]
        + [myo1d_report.get("status", "PASS"), membrane_report.get("status", "PASS")]
    )
    manifest = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "mode": mode,
        "profile": profile,
        "input_root": str(root),
        "contract_path": str(contract_path) if contract_path else None,
        "status": status,
        "inputs": input_records,
        "messages": messages,
    }
    qc_report = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "mode": mode,
        "profile": profile,
        "status": status,
        "summary": {
            "receptor_count": len(receptor_reports),
            "message_counts": {
                "PASS": sum(1 for item in input_records if item.get("status") == "PASS"),
                "WARN": sum(1 for item in input_records if item.get("status") == "WARN"),
                "FAIL": sum(1 for item in input_records if item.get("status") == "FAIL"),
            },
        },
        "receptors": receptor_reports,
        "myo1d": myo1d_report,
        "membrane_frame": membrane_report,
        "messages": messages,
        "policy": {
            "docking_enabled": False,
            "normalization_enabled": False,
            "scoring_enabled": False,
        },
    }

    write_json(ctx.manifest_dir / "structure_input_manifest.json", manifest, ctx)
    write_json(ctx.manifest_dir / "structure_qc_report.json", qc_report, ctx)
    write_json(ctx.qc_dir / "membrane_frame_report.json", membrane_report, ctx)
    write_json(ctx.qc_dir / "myo1d_annotation_report.json", myo1d_report, ctx)
    write_csv(ctx.qc_dir / "chain_protomer_audit.csv", chain_rows, [
        "receptor_id", "path", "expected_assembly", "protomer_id",
        "expected_chain_ids", "observed_chain_ids", "status", "message",
    ], ctx)
    write_csv(ctx.qc_dir / "residue_mapping_audit.csv", residue_rows, [
        "structure_id", "path", "chain_id", "min_residue", "max_residue",
        "residue_count", "large_jumps", "residue_reset_risk", "status", "message",
    ], ctx)
    append_phase_status(
        ctx,
        "structure_inputs",
        status,
        f"structure input validation completed with {status}",
        {"mode": mode, "profile": profile, "input_root": str(root)},
    )
    return qc_report

