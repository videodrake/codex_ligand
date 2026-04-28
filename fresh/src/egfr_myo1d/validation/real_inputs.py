"""Task 5 real-input readiness bridge.

This module validates real EGFR/MYO1D/membrane-frame inputs against the Task 4
PPI input contract style. It is a readiness/audit layer only: it does not write
prepared PDBs, normalize structures, run docking, or invoke external tools.
"""

import csv
import json
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import initialize_manifests, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError
from egfr_myo1d.io.hashing import describe_file
from egfr_myo1d.preparation.constructs import CAP_RESNAMES, STANDARD_AA, residue_annotation
from egfr_myo1d.preparation.masks import build_egfr_masks
from egfr_myo1d.preparation.restraints import build_active_face_contract
from egfr_myo1d.structure.geometry import validate_membrane_frame
from egfr_myo1d.structure.pdb_parser import PDBParseError, PDBStructure, parse_pdb
from egfr_myo1d.validation.prepared_inputs import read_contract, resolve_input


NOT_IMPLEMENTED = [
    "ppi_docking",
    "pyrosetta_docking_or_relaxation",
    "lightdock",
    "vina",
    "pocket_discovery",
    "compound_docking",
    "scoring",
    "candidate_nomination",
    "qsub_pbs_generation",
    "cleanup_deletion",
    "structure_mutation_repair",
    "silent_renumbering_or_normalization",
]


def write_csv(path, rows, fieldnames, ctx):
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    with safe_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def real_warning_policy(
    code,
    input_id,
    contract,
    mode,
    profile,
    strict,
):
    """Return severity, classification, production_blocker, quarantine."""
    whitelist = set(contract.get("task3_closeout", {}).get("fixture_warning_whitelist", []))
    fixture_allowed = input_id in whitelist and profile == "codex_dev" and not strict
    if fixture_allowed:
        return "WARN", "fixture_warning", False, False
    return "FAIL", "quarantine", True, True


def message(
    input_id,
    severity,
    classification,
    code,
    text,
    production_blocker=False,
    quarantine=False,
):
    return {
        "input_id": input_id,
        "severity": severity,
        "classification": classification,
        "production_blocker": production_blocker,
        "quarantine": quarantine,
        "code": code,
        "message": text,
    }


def status_from_messages(messages):
    if any(item["severity"] == "FAIL" for item in messages):
        return "FAIL"
    if any(item["severity"] == "WARN" for item in messages):
        return "WARN"
    return "PASS"


def count_summary(messages, pass_count=0):
    fail_count = sum(1 for item in messages if item["severity"] == "FAIL")
    warn_count = sum(1 for item in messages if item["severity"] == "WARN")
    return {
        "PASS": pass_count,
        "WARN": warn_count,
        "FAIL": fail_count,
        "fixture_warning": sum(1 for item in messages if item["classification"] == "fixture_warning"),
        "production_blocker": sum(1 for item in messages if item.get("production_blocker")),
        "quarantine": sum(1 for item in messages if item.get("quarantine")),
        "verdict": "FAIL" if fail_count else ("WARN" if warn_count else "PASS"),
    }


def residue_reset_risk(structure, expected_range):
    expected_start = int(expected_range[0])
    for chain_id in structure.chain_ids:
        numbers = structure.residue_numbers_for_chain(chain_id)
        if numbers and min(numbers) <= 5 and max(numbers) < expected_start:
            return True
    return False


def validate_real_receptor(
    contract,
    input_root,
    mode,
    profile,
    strict,
):
    receptor = contract["receptor"]
    input_id = receptor["id"]
    path = resolve_input(input_root, receptor["path"])
    desc = describe_file(path)
    manifest_record = {"input_id": input_id, "role": "receptor", "path": str(path)}
    manifest_record.update(desc)
    rows = []
    messages = []
    if not path.exists():
        messages.append(message(input_id, "FAIL", "missing_input", "missing_receptor", "missing receptor input: {0}".format(path), True, True))
        return {"status": "FAIL", "chains": [], "protomer_count": 0, "messages": messages}, rows, messages, manifest_record

    try:
        structure = parse_pdb(path)
    except PDBParseError as exc:
        messages.append(message(input_id, "FAIL", "malformed_input", "pdb_parse_error", str(exc), True, True))
        return {"status": "FAIL", "chains": [], "protomer_count": 0, "messages": messages}, rows, messages, manifest_record

    expected_chains = set(receptor.get("expected_chains", []))
    observed_chains = set(structure.chain_ids)
    if len(observed_chains) < int(receptor.get("required_protomer_count", 2)):
        sev, cls, blocker, quarantine = real_warning_policy("single_chain_dimer_ambiguity", input_id, contract, mode, profile, strict)
        messages.append(message(input_id, sev, cls, "single_chain_dimer_ambiguity", "real receptor has fewer than two explicit chains", blocker, quarantine))
    if not expected_chains.issubset(observed_chains):
        sev, cls, blocker, quarantine = real_warning_policy("expected_chains_AB_missing", input_id, contract, mode, profile, strict)
        messages.append(message(input_id, sev, cls, "expected_chains_AB_missing", "expected chains {0}; observed {1}".format(sorted(expected_chains), sorted(observed_chains)), blocker, quarantine))
    duplicate_count = structure.duplicate_atom_identity_count()
    if duplicate_count:
        sev, cls, blocker, quarantine = real_warning_policy("duplicate_atom_identities", input_id, contract, mode, profile, strict)
        messages.append(message(input_id, sev, cls, "duplicate_atom_identities", "{0} duplicate atom identities observed".format(duplicate_count), blocker, quarantine))
    if residue_reset_risk(structure, receptor["biological_residue_range"]):
        messages.append(message(input_id, "FAIL", "quarantine", "residue_number_reset_risk", "residue numbering appears reset rather than project numbering", True, True))

    warn_mutations = receptor.get("warn_mutations", [])
    mutation_residues = set()
    for check in warn_mutations:
        residue_number = int(check["residue_number"])
        warn_if = check["warn_if"]
        observed = {residue.resname for residue in structure.residues if residue.residue_number == residue_number}
        if warn_if in observed:
            mutation_residues.add(residue_number)
            sev, cls, blocker, quarantine = real_warning_policy("v924r_like_mutation", input_id, contract, mode, profile, strict)
            messages.append(message(input_id, sev, cls, "v924r_like_mutation", check.get("reason", "{0}{1} detected; do not repair silently".format(warn_if, residue_number)), blocker, quarantine))

    for residue in structure.residues:
        protomer_id = residue.chain_id if residue.chain_id in expected_chains else "unknown"
        rows.append(
            {
                "receptor_id": input_id,
                "path": str(path),
                "chain_id": residue.chain_id,
                "protomer_id": protomer_id,
                "residue_number": residue.residue_number,
                "insertion_code": residue.insertion_code,
                "residue_name": residue.resname,
                "record_type": "/".join(residue.record_types_present),
                "biological": residue.resname in STANDARD_AA,
                "is_warn_mutation": residue.residue_number in mutation_residues,
                "chain_identity_ok": expected_chains.issubset(observed_chains),
                "dimer_context_ok": len(observed_chains) >= 2,
                "residue_numbering_reset_risk": residue_reset_risk(structure, receptor["biological_residue_range"]),
                "verdict": status_from_messages(messages),
                "notes": "do not mutate or renumber" if residue.residue_number in mutation_residues else "",
            }
        )

    report = {
        "status": status_from_messages(messages),
        "receptor_id": input_id,
        "path": str(path),
        "chains": list(structure.chain_ids),
        "protomer_count": len(structure.chain_ids),
        "chain_protomer_identity_preserved": expected_chains.issubset(observed_chains),
        "dimer_context_preserved": len(observed_chains) >= 2,
        "residue_numbering_preserved": not residue_reset_risk(structure, receptor["biological_residue_range"]),
        "mutations_repaired": False,
        "messages": messages,
    }
    return report, rows, messages, manifest_record


def validate_real_myo1d(
    contract,
    input_root,
    mode,
    profile,
    strict,
):
    myo = contract["myo1d"]
    primary = myo["primary_construct"]
    comparator = myo.get("comparator_construct")
    messages = []
    rows = []
    inputs = {}

    def check_construct(construct, role):
        input_id = construct["id"]
        path = resolve_input(input_root, construct["path"])
        inputs[role] = {"input_id": input_id, "role": role, "path": str(path)}
        inputs[role].update(describe_file(path))
        if not path.exists():
            messages.append(message(input_id, "FAIL", "missing_input", "missing_myo1d_construct", "missing MYO1D {0}: {1}".format(role, path), True, True))
            return {"status": "FAIL", "path": str(path), "residues": []}
        try:
            structure = parse_pdb(path)
        except PDBParseError as exc:
            messages.append(message(input_id, "FAIL", "malformed_input", "pdb_parse_error", str(exc), True, True))
            return {"status": "FAIL", "path": str(path), "residues": []}
        biological_numbers = sorted({residue.residue_number for residue in structure.residues if residue.resname not in CAP_RESNAMES})
        start, end = construct["biological_residue_range"]
        if not biological_numbers or min(biological_numbers) != int(start) or max(biological_numbers) != int(end):
            code = "myo1d_terminal_artifact" if biological_numbers and min(biological_numbers) == 962 else "myo1d_range_mismatch"
            sev, cls, blocker, quarantine = real_warning_policy(code, input_id, contract, mode, profile, strict)
            if role == "primary" and code != "myo1d_terminal_artifact":
                sev, cls, blocker, quarantine = "FAIL", "quarantine", True, True
            messages.append(message(input_id, sev, cls, code, "{0} construct observed range {1}-{2}; expected {3}-{4}".format(role, biological_numbers[:1], biological_numbers[-1:], start, end), blocker, quarantine))
        if role == "primary" and biological_numbers and min(biological_numbers) == 962:
            sev, cls, blocker, quarantine = real_warning_policy("myo1d_terminal_artifact", input_id, contract, mode, profile, strict)
            messages.append(message(input_id, sev, cls, "myo1d_terminal_artifact", "primary MYO1D starts at 962 and lacks 955-960 structural buffer", blocker, quarantine))
        if role == "comparator":
            messages.append(message(input_id, "WARN", "comparator_noise_monitor", "myo1d_955_1006_comparator_only", "955-1006 is comparator/noise-monitor only, not production primary"))
        residue_numbers = set(biological_numbers)
        for required in myo["active_face"] + myo["sheet12_support"]:
            if required not in residue_numbers:
                messages.append(message(input_id, "FAIL", "quarantine", "missing_myo1d_annotation_residue", "missing annotation residue {0}".format(required), True, True))
        for residue in structure.residues:
            rows.append(
                {
                    "construct_id": input_id,
                    "role": role,
                    "path": str(path),
                    "chain_id": residue.chain_id,
                    "residue_number": residue.residue_number,
                    "insertion_code": residue.insertion_code,
                    "residue_name": residue.resname,
                    "record_type": "/".join(residue.record_types_present),
                    "annotation_class": "cap" if residue.resname in CAP_RESNAMES else residue_annotation(residue.residue_number, myo),
                    "score_bonus_allowed": False,
                    "verdict": "PASS",
                    "notes": "comparator/noise monitor only" if role == "comparator" and residue.residue_number >= 998 else "annotation only; no scoring bonus",
                }
            )
        return {
            "status": "PASS",
            "path": str(path),
            "residue_range": [min(biological_numbers), max(biological_numbers)] if biological_numbers else [],
            "residue_count": len(biological_numbers),
        }

    primary_report = check_construct(primary, "primary")
    comparator_report = check_construct(comparator, "comparator") if comparator else None
    report = {
        "status": status_from_messages([item for item in messages if item["severity"] == "FAIL"]),
        "primary_construct": primary["id"],
        "primary": primary_report,
        "comparator_construct": comparator["id"] if comparator else None,
        "comparator": comparator_report,
        "active_face_annotation_present": not any(item["code"] == "missing_myo1d_annotation_residue" for item in messages),
        "score_bonus_allowed": False,
        "messages": messages,
    }
    return report, rows, messages, inputs


def validate_real_membrane_frame(
    contract,
    input_root,
):
    frame_path = resolve_input(input_root, contract["membrane_frame"]["path"])
    desc = describe_file(frame_path)
    manifest_record = {"input_id": "membrane_frame", "role": "membrane_frame", "path": str(frame_path)}
    manifest_record.update(desc)
    rows = []
    messages = []
    if not frame_path.exists():
        messages.append(message("membrane_frame", "FAIL", "missing_input", "missing_membrane_frame", "missing membrane frame: {0}".format(frame_path), True, True))
        return {"status": "FAIL", "messages": messages}, rows, messages, manifest_record
    try:
        frame = json.loads(frame_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        messages.append(message("membrane_frame", "FAIL", "malformed_input", "malformed_membrane_frame", str(exc), True, True))
        return {"status": "FAIL", "messages": messages}, rows, messages, manifest_record
    frame_report = validate_membrane_frame(frame)
    if frame_report["status"] == "FAIL":
        messages.append(message("membrane_frame", "FAIL", "quarantine", "invalid_membrane_frame", "membrane-frame vector validation failed", True, True))
    elif frame_report["status"] == "WARN":
        messages.append(message("membrane_frame", "WARN", "metadata_warning", "membrane_frame_warning", "membrane-frame metadata has warnings"))
    for field in ["membrane_normal", "dimer_axis"]:
        rows.append(
            {
                "field": field,
                "path": str(frame_path),
                "status": frame_report["status"],
                "norm": frame_report.get("{0}_norm".format(field)),
                "unit_vector": json.dumps(frame_report.get("{0}_unit".format(field))),
                "coordinate_convention": frame_report.get("coordinate_convention", ""),
                "notes": "; ".join(item["message"] for item in frame_report.get("messages", [])),
            }
        )
    report = {"status": status_from_messages(messages), "path": str(frame_path), "frame_qc": frame_report, "messages": messages}
    return report, rows, messages, manifest_record


def failure_outputs(
    ctx,
    mode,
    profile,
    input_root,
    contract_path,
    text,
):
    msg = message("contract", "FAIL", "missing_or_malformed_contract", "contract_error", text, True, True)
    summary = count_summary([msg])
    manifest = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "mode": mode,
        "profile": profile,
        "input_root": str(input_root),
        "contract": {"path": str(contract_path) if contract_path else None, "present": False, "sha256": None},
        "inputs": {},
        "qc_summary": summary,
        "blockers": [msg],
        "not_implemented": NOT_IMPLEMENTED,
    }
    report = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "mode": mode,
        "profile": profile,
        "input_root": str(input_root),
        "status": "FAIL",
        "qc_summary": summary,
        "blockers": [msg],
        "not_implemented": NOT_IMPLEMENTED,
    }
    write_json(ctx.manifest_dir / "real_input_manifest.json", manifest, ctx)
    write_json(ctx.manifest_dir / "real_input_readiness_report.json", report, ctx)
    write_csv(ctx.qc_dir / "real_egfr_readiness_audit.csv", [], REAL_EGFR_FIELDS, ctx)
    write_csv(ctx.qc_dir / "real_myo1d_readiness_audit.csv", [], REAL_MYO_FIELDS, ctx)
    write_csv(ctx.qc_dir / "real_membrane_frame_audit.csv", [], REAL_MEMBRANE_FIELDS, ctx)
    write_csv(ctx.qc_dir / "real_input_blockers.csv", [msg], REAL_BLOCKER_FIELDS, ctx)
    append_phase_status(ctx, "validate-real-inputs", "FAIL", text, {"mode": mode, "profile": profile})
    return report


REAL_EGFR_FIELDS = [
    "receptor_id", "path", "chain_id", "protomer_id", "residue_number", "insertion_code",
    "residue_name", "record_type", "biological", "is_warn_mutation", "chain_identity_ok",
    "dimer_context_ok", "residue_numbering_reset_risk", "verdict", "notes",
]
REAL_MYO_FIELDS = [
    "construct_id", "role", "path", "chain_id", "residue_number", "insertion_code",
    "residue_name", "record_type", "annotation_class", "score_bonus_allowed", "verdict", "notes",
]
REAL_MEMBRANE_FIELDS = ["field", "path", "status", "norm", "unit_vector", "coordinate_convention", "notes"]
REAL_BLOCKER_FIELDS = ["input_id", "severity", "classification", "production_blocker", "quarantine", "code", "message"]


def validate_real_inputs(
    ctx,
    mode,
    profile,
    input_root,
    contract_path=None,
    strict=False,
):
    initialize_manifests(ctx, mode)
    root = input_root if input_root.is_absolute() else (ctx.repo_root / input_root)
    root = root.resolve()
    try:
        contract, resolved_contract_path = read_contract(root, contract_path)
    except (FileNotFoundError, RunContextError, ValueError) as exc:
        return failure_outputs(ctx, mode, profile, root, contract_path, str(exc))

    messages = []
    egfr_report, egfr_rows, egfr_messages, egfr_input = validate_real_receptor(contract, root, mode, profile, strict or profile == "hpc_strict")
    messages.extend(egfr_messages)
    myo_report, myo_rows, myo_messages, myo_inputs = validate_real_myo1d(contract, root, mode, profile, strict or profile == "hpc_strict")
    messages.extend(myo_messages)
    membrane_report, membrane_rows, membrane_messages, membrane_input = validate_real_membrane_frame(contract, root)
    messages.extend(membrane_messages)

    mask_preview = None
    if membrane_report["status"] != "FAIL":
        try:
            membrane_mask, non_target_mask, _rows = build_egfr_masks(contract, Path(membrane_input["path"]))
            mask_preview = {
                "membrane_exclusion_mask_serializable": True,
                "non_target_region_mask_serializable": True,
                "membrane_exclusion_mask": membrane_mask,
                "non_target_region_mask": non_target_mask,
            }
        except Exception as exc:  # keep readiness deterministic and explicit
            messages.append(message("mask_contracts", "FAIL", "quarantine", "mask_serialization_failed", str(exc), True, True))
            mask_preview = {"membrane_exclusion_mask_serializable": False, "non_target_region_mask_serializable": False}

    active_face_contract = build_active_face_contract(contract["myo1d"])
    if active_face_contract["score_bonus_allowed"] is not False:
        messages.append(message("active_face", "FAIL", "quarantine", "score_bonus_not_allowed", "active-face residues must not create score bonuses", True, True))

    contract_input = {"path": str(resolved_contract_path)}
    contract_input.update(describe_file(resolved_contract_path))
    present_inputs = {
        "contract": contract_input,
        "receptor": egfr_input,
        "membrane_frame": membrane_input,
    }
    present_inputs.update(myo_inputs)
    pass_count = sum(1 for item in present_inputs.values() if item.get("present"))
    summary = count_summary(messages, pass_count=pass_count)
    blockers = [item for item in messages if item["severity"] == "FAIL"]
    report = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "mode": mode,
        "profile": profile,
        "input_root": str(root),
        "contract": contract_input,
        "inputs": present_inputs,
        "status": summary["verdict"],
        "qc_summary": summary,
        "blockers": blockers,
        "messages": messages,
        "egfr": egfr_report,
        "myo1d": myo_report,
        "membrane_frame": membrane_report,
        "active_face_contract": active_face_contract,
        "mask_contract_preview": mask_preview,
        "not_implemented": NOT_IMPLEMENTED,
    }
    manifest = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "phase": "validate-real-inputs",
        "mode": mode,
        "profile": profile,
        "input_root": str(root),
        "contract": report["contract"],
        "inputs": present_inputs,
        "qc_summary": summary,
        "final_verdict": summary["verdict"],
        "blocker_list": blockers,
        "not_implemented": NOT_IMPLEMENTED,
    }

    write_json(ctx.manifest_dir / "real_input_manifest.json", manifest, ctx)
    write_json(ctx.manifest_dir / "real_input_readiness_report.json", report, ctx)
    write_csv(ctx.qc_dir / "real_egfr_readiness_audit.csv", egfr_rows, REAL_EGFR_FIELDS, ctx)
    write_csv(ctx.qc_dir / "real_myo1d_readiness_audit.csv", myo_rows, REAL_MYO_FIELDS, ctx)
    write_csv(ctx.qc_dir / "real_membrane_frame_audit.csv", membrane_rows, REAL_MEMBRANE_FIELDS, ctx)
    write_csv(ctx.qc_dir / "real_input_blockers.csv", blockers, REAL_BLOCKER_FIELDS, ctx)
    append_phase_status(
        ctx,
        "validate-real-inputs",
        summary["verdict"],
        "real input readiness completed with {0}".format(summary['verdict']),
        {"mode": mode, "profile": profile, "strict": strict},
    )
    append_job_status(ctx, "validate_real_inputs_local", summary["verdict"], details={"message": "local pure-Python readiness validation"})
    return report

