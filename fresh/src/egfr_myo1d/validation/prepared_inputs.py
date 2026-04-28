"""Task 4 PPI prepared-input generation and QC."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import initialize_manifests, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError, ensure_within
from egfr_myo1d.io.hashing import describe_file
from egfr_myo1d.preparation.constructs import (
    CAP_RESNAMES,
    STANDARD_AA,
    detect_terminal_artifact,
    prepare_myo1d_construct,
    residue_rows_for_construct,
    validate_active_face_presence,
)
from egfr_myo1d.preparation.masks import build_egfr_masks, write_json as write_contract_json
from egfr_myo1d.preparation.pdb_writer import write_pdb_atoms
from egfr_myo1d.preparation.restraints import write_active_face_contract
from egfr_myo1d.structure.pdb_parser import PDBParseError, PDBStructure, parse_pdb


TASK3_WARNING_CLASSES = {
    "single_chain_dimer_ambiguity",
    "duplicate_atom_identities",
    "expected_chains_AB_missing",
    "v924r_like_mutation",
    "myo1d_terminal_artifact",
    "myo1d_tail_masked_comparator",
}


def read_contract(input_root: Path, contract_path: Path | None) -> tuple[dict[str, Any], Path]:
    if contract_path is None:
        path = input_root / "ppi_input_contract.json"
    elif contract_path.is_absolute():
        path = contract_path
    elif contract_path.exists():
        path = contract_path
    else:
        path = input_root / contract_path
    path = path.resolve()
    root = input_root.resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise RunContextError(f"contract path escapes input root: {path}") from exc
    try:
        return json.loads(path.read_text(encoding="utf-8")), path
    except json.JSONDecodeError as exc:
        raise ValueError(f"malformed PPI input contract: {path}: {exc}") from exc


def resolve_input(root: Path, rel_path: str) -> Path:
    path = (root / rel_path).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError as exc:
        raise RunContextError(f"input path escapes input root: {rel_path}") from exc
    return path


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str], ctx: RunContext) -> None:
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    with safe_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def classify_warning(
    warning_class: str,
    fixture_id: str,
    contract: dict[str, Any],
    mode: str,
    profile: str,
    strict: bool = False,
) -> tuple[str, str]:
    whitelist = set(contract.get("task3_closeout", {}).get("fixture_warning_whitelist", []))
    fixture_allowed = fixture_id in whitelist and mode == "smoke_env" and profile == "codex_dev" and not strict
    if fixture_allowed:
        return "fixture_warning", "WARN"
    return "production_blocker", "FAIL"


def audit_receptor(
    structure: PDBStructure,
    receptor_contract: dict[str, Any],
    contract: dict[str, Any],
    mode: str,
    profile: str,
    strict: bool = False,
    fixture_id: str | None = None,
) -> tuple[str, list[dict[str, Any]], list[dict[str, Any]]]:
    receptor_id = receptor_contract["id"]
    fixture_role = fixture_id or receptor_id
    expected_chains = receptor_contract["expected_chains"]
    warn_mutations = receptor_contract.get("warn_mutations", [])
    messages: list[dict[str, Any]] = []
    observed_chains = set(structure.chain_ids)
    duplicate_count = structure.duplicate_atom_identity_count()

    if len(observed_chains) < receptor_contract.get("required_protomer_count", 2):
        classification, severity = classify_warning("single_chain_dimer_ambiguity", fixture_role, contract, mode, profile, strict)
        messages.append({"severity": severity, "classification": classification, "code": "single_chain_dimer_ambiguity", "message": "fewer than two receptor chains observed"})
    if duplicate_count:
        classification, severity = classify_warning("duplicate_atom_identities", fixture_role, contract, mode, profile, strict)
        messages.append({"severity": severity, "classification": classification, "code": "duplicate_atom_identities", "message": f"{duplicate_count} duplicate atom identities observed"})
    if not set(expected_chains).issubset(observed_chains):
        classification, severity = classify_warning("expected_chains_AB_missing", fixture_role, contract, mode, profile, strict)
        messages.append({"severity": severity, "classification": classification, "code": "expected_chains_AB_missing", "message": f"expected chains {expected_chains}; observed {sorted(observed_chains)}"})

    mutation_residues: set[int] = set()
    for check in warn_mutations:
        residue_number = int(check["residue_number"])
        warn_if = check["warn_if"]
        observed = {residue.resname for residue in structure.residues if residue.residue_number == residue_number}
        if warn_if in observed:
            classification, severity = classify_warning("v924r_like_mutation", fixture_role, contract, mode, profile, strict)
            mutation_residues.add(residue_number)
            messages.append({"severity": severity, "classification": classification, "code": "v924r_like_mutation", "message": check.get("reason", f"{warn_if}{residue_number} detected")})

    rows: list[dict[str, Any]] = []
    for residue in structure.residues:
        protomer_id = residue.chain_id if residue.chain_id in expected_chains else "unknown"
        biological = residue.resname in STANDARD_AA
        is_cap = residue.resname in CAP_RESNAMES
        warn_mut = residue.residue_number in mutation_residues
        rows.append(
            {
                "receptor_id": receptor_id,
                "fixture_role": fixture_role,
                "chain_id": residue.chain_id,
                "protomer_id": protomer_id,
                "residue_number": residue.residue_number,
                "insertion_code": residue.insertion_code,
                "residue_name": residue.resname,
                "record_type": "/".join(residue.record_types_present),
                "biological": biological,
                "is_cap": is_cap,
                "is_warn_mutation": warn_mut,
                "warning_classification": ";".join(sorted({m["classification"] for m in messages})) if messages else "",
                "production_policy": "fail_or_quarantine" if messages else "accepted",
                "notes": "preserve; do not mutate" if warn_mut else "",
            }
        )

    status = "FAIL" if any(message["severity"] == "FAIL" for message in messages) else ("WARN" if messages else "PASS")
    return status, rows, messages


def output_record(ctx: RunContext, path: Path, role: str) -> dict[str, Any]:
    desc = describe_file(path)
    return {
        "role": role,
        "path": ctx.relative_to_repo(path),
        "sha256": desc["sha256"],
        "size_bytes": desc["size_bytes"],
        "present": desc["present"],
    }


def summarize_qc(messages: list[dict[str, Any]], extra_statuses: list[str]) -> dict[str, Any]:
    pass_count = sum(1 for status in extra_statuses if status == "PASS")
    warn_count = sum(1 for status in extra_statuses if status == "WARN") + sum(1 for m in messages if m.get("severity") == "WARN")
    fail_count = sum(1 for status in extra_statuses if status == "FAIL") + sum(1 for m in messages if m.get("severity") == "FAIL")
    fixture_warning = sum(1 for m in messages if m.get("classification") == "fixture_warning")
    production_blocker = sum(1 for m in messages if m.get("classification") == "production_blocker")
    comparator = sum(1 for m in messages if m.get("classification") == "explicitly_whitelisted_comparator")
    verdict = "FAIL" if fail_count else ("WARN" if warn_count else "PASS")
    return {
        "PASS": pass_count,
        "WARN": warn_count,
        "FAIL": fail_count,
        "fixture_warning": fixture_warning,
        "production_blocker": production_blocker,
        "explicitly_whitelisted_comparator": comparator,
        "verdict": verdict,
    }


def validate_negative_controls(contract: dict[str, Any], input_root: Path, mode: str, profile: str, strict: bool) -> list[dict[str, Any]]:
    messages: list[dict[str, Any]] = []
    for item in contract.get("negative_controls", []):
        control_id = item["id"]
        path = resolve_input(input_root, item["path"])
        if not path.exists():
            continue
        if control_id.startswith("myo1d_"):
            structure = parse_pdb(path)
            terminal = detect_terminal_artifact(structure, control_id, path, allow_fixture_warning=not strict)
            if terminal["verdict"] in {"WARN", "FAIL"}:
                classification, severity = classify_warning("myo1d_terminal_artifact", control_id, contract, mode, profile, strict)
                messages.append({"severity": severity, "classification": classification, "code": "myo1d_terminal_artifact", "message": terminal["notes"]})
        else:
            structure = parse_pdb(path)
            status, _rows, receptor_messages = audit_receptor(
                structure,
                contract["receptor"],
                contract,
                mode,
                profile,
                strict=strict,
                fixture_id=control_id,
            )
            messages.extend(receptor_messages)
    return messages


def prepare_ppi_inputs(
    ctx: RunContext,
    mode: str,
    profile: str,
    input_root: Path,
    contract_path: Path | None = None,
    strict: bool = False,
) -> dict[str, Any]:
    initialize_manifests(ctx, mode)
    root = input_root if input_root.is_absolute() else (ctx.repo_root / input_root)
    root = root.resolve()
    contract, resolved_contract_path = read_contract(root, contract_path)
    prepared_root = ctx.run_dir / "prepared"
    egfr_dir = prepared_root / "egfr"
    myo_dir = prepared_root / "myo1d"
    restraints_dir = prepared_root / "restraints"
    for path in [egfr_dir, myo_dir, restraints_dir]:
        ctx.require_within_run_dir(path)
        path.mkdir(parents=True, exist_ok=True)

    receptor_contract = contract["receptor"]
    receptor_path = resolve_input(root, receptor_contract["path"])
    receptor_structure = parse_pdb(receptor_path)
    receptor_status, egfr_audit_rows, receptor_messages = audit_receptor(
        receptor_structure, receptor_contract, contract, mode, profile, strict=strict
    )
    if receptor_status == "FAIL":
        raise ValueError("primary receptor failed Task 4 precondition checks")
    prepared_receptor = egfr_dir / "egfr_receptor_normalized.pdb"
    write_pdb_atoms(prepared_receptor, receptor_structure.atoms)

    myo_contract = contract["myo1d"]
    primary = myo_contract["primary_construct"]
    comparator = myo_contract["comparator_construct"]
    primary_path = resolve_input(root, primary["path"])
    comparator_path = resolve_input(root, comparator["path"])
    primary_out = myo_dir / f"{primary['id']}.pdb"
    comparator_out = myo_dir / f"{comparator['id']}.pdb"
    primary_record, primary_rows, primary_terminal = prepare_myo1d_construct(
        primary_path,
        primary_out,
        primary["id"],
        tuple(primary["biological_residue_range"]),
        myo_contract,
        "production_primary",
    )
    comparator_record, comparator_rows, comparator_terminal = prepare_myo1d_construct(
        comparator_path,
        comparator_out,
        comparator["id"],
        tuple(comparator["biological_residue_range"]),
        myo_contract,
        "tail_masked_comparator",
    )
    comparator_message = {
        "severity": "WARN",
        "classification": "explicitly_whitelisted_comparator",
        "code": "myo1d_tail_masked_comparator",
        "message": "955-1006 comparator emitted only as tail-masked/noise-monitor construct",
    }

    bad_messages = validate_negative_controls(contract, root, mode, profile, strict)

    active_face_path = restraints_dir / "myo1d_active_face_contract.json"
    active_face_contract = write_active_face_contract(active_face_path, myo_contract)
    active_status, face_rows = validate_active_face_presence(parse_pdb(primary_out), myo_contract)

    membrane_frame_path = resolve_input(root, contract["membrane_frame"]["path"])
    membrane_mask, non_target_mask, membrane_rows = build_egfr_masks(contract, membrane_frame_path)
    membrane_mask_path = restraints_dir / "egfr_membrane_exclusion_mask.json"
    non_target_mask_path = restraints_dir / "egfr_non_target_region_mask.json"
    write_contract_json(membrane_mask_path, membrane_mask)
    write_contract_json(non_target_mask_path, non_target_mask)

    terminal_rows = [primary_terminal, comparator_terminal]
    for item in contract.get("negative_controls", []):
        if item["id"].startswith("myo1d_"):
            path = resolve_input(root, item["path"])
            if path.exists():
                terminal_rows.append(detect_terminal_artifact(parse_pdb(path), item["id"], path, allow_fixture_warning=True))

    write_csv(ctx.qc_dir / "egfr_receptor_normalization_audit.csv", egfr_audit_rows, [
        "receptor_id", "fixture_role", "chain_id", "protomer_id", "residue_number",
        "insertion_code", "residue_name", "record_type", "biological", "is_cap",
        "is_warn_mutation", "warning_classification", "production_policy", "notes",
    ], ctx)
    write_csv(ctx.qc_dir / "myo1d_construct_audit.csv", primary_rows + comparator_rows, [
        "construct_id", "source_path", "chain_id", "residue_number", "insertion_code",
        "residue_name", "record_type", "biological", "is_cap", "annotation_class",
        "output_role", "score_bonus_allowed", "notes",
    ], ctx)
    write_csv(ctx.qc_dir / "terminal_artifact_audit.csv", terminal_rows, [
        "construct_id", "n_terminal_residue", "c_terminal_residue", "has_cap_records",
        "has_free_artifact_terminal", "terminal_contact_region_policy", "verdict", "notes",
    ], ctx)
    write_csv(ctx.qc_dir / "face_orientation_audit.csv", face_rows, [
        "construct_id", "residue_number", "annotation_group", "present",
        "score_bonus_allowed", "status", "notes",
    ], ctx)
    write_csv(ctx.qc_dir / "membrane_exclusion_audit.csv", membrane_rows, [
        "receptor_id", "mask_name", "source", "policy", "threshold",
        "residue_count", "status", "notes",
    ], ctx)

    messages = receptor_messages + [comparator_message] + bad_messages
    statuses = [receptor_status, active_status, membrane_mask["membrane_frame_qc"]["status"]]
    qc_summary = summarize_qc(messages, statuses)

    outputs = {
        "egfr_receptor": output_record(ctx, prepared_receptor, "prepared_receptor"),
        "myo1d_primary": output_record(ctx, primary_out, "myo1d_primary"),
        "myo1d_comparator": output_record(ctx, comparator_out, "myo1d_tail_masked_comparator"),
        "active_face_contract": output_record(ctx, active_face_path, "active_face_contract"),
        "membrane_exclusion_mask": output_record(ctx, membrane_mask_path, "egfr_membrane_exclusion_mask"),
        "non_target_region_mask": output_record(ctx, non_target_mask_path, "egfr_non_target_region_mask"),
    }
    inputs = {
        "contract": {"path": str(resolved_contract_path), **describe_file(resolved_contract_path)},
        "receptor": {"path": str(receptor_path), **describe_file(receptor_path)},
        "myo1d_primary": {"path": str(primary_path), **describe_file(primary_path)},
        "myo1d_comparator": {"path": str(comparator_path), **describe_file(comparator_path)},
        "membrane_frame": {"path": str(membrane_frame_path), **describe_file(membrane_frame_path)},
    }
    manifest = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": "prepare-ppi-inputs",
        "task3_closeout": {
            "warning_whitelist_applied": True,
            "fixture_warning_whitelist": contract["task3_closeout"]["fixture_warning_whitelist"],
            "production_same_warning_policy": contract["task3_closeout"]["production_same_warning_policy"],
        },
        "inputs": inputs,
        "outputs": outputs,
        "egfr": {
            "receptor_id": receptor_contract["id"],
            "chains": list(receptor_structure.chain_ids),
            "protomer_count": len(receptor_structure.chain_ids),
            "chain_protomer_identity_preserved": set(receptor_contract["expected_chains"]).issubset(set(receptor_structure.chain_ids)),
            "source_file_mapping_preserved": True,
            "residue_numbering_preserved": True,
            "dimer_context_preserved": len(receptor_structure.chain_ids) >= 2,
        },
        "myo1d": {
            "primary_construct": primary["id"],
            "comparator_construct": comparator["id"],
            "active_face_annotation_present": active_status == "PASS",
            "score_bonus_allowed": False,
        },
        "qc_summary": qc_summary,
    }
    report = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "status": qc_summary["verdict"],
        "messages": messages,
        "task3_closeout": manifest["task3_closeout"],
        "egfr_audit_rows": len(egfr_audit_rows),
        "myo1d_audit_rows": len(primary_rows + comparator_rows),
        "terminal_audit_rows": len(terminal_rows),
        "active_face_contract": active_face_contract,
        "membrane_mask": membrane_mask,
        "qc_summary": qc_summary,
        "not_implemented": [
            "docking",
            "pocket_discovery",
            "scoring",
            "candidate_nomination",
            "qsub_pbs_generation",
            "cleanup_deletion",
        ],
    }
    write_json(ctx.manifest_dir / "prepared_input_manifest.json", manifest, ctx)
    write_json(ctx.manifest_dir / "preparation_qc_report.json", report, ctx)
    append_phase_status(
        ctx,
        "prepare-ppi-inputs",
        qc_summary["verdict"],
        f"PPI input preparation completed with {qc_summary['verdict']}",
        {"mode": mode, "profile": profile, "strict": strict},
    )
    append_job_status(ctx, "prepare_ppi_inputs_local", qc_summary["verdict"], details={"message": "local pure-Python preparation"})
    return report
