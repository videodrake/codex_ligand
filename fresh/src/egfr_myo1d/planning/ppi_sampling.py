"""Deterministic spec-only PPI sampling plan builder."""

from pathlib import Path
from typing import Any

from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.io.hashing import describe_file
from egfr_myo1d.validation.prepared_inputs import resolve_input


PLANNER_VERSION = "task6_ppi_sampling_plan_v0_1"
METHODS = ["pyrosetta_global_ppi", "lightdock_gso_ppi"]
PRIMARY_SEEDS = [0, 1]
COMPARATOR_SEEDS = [0]


def _method_tag(method):
    if method == "pyrosetta_global_ppi":
        return "pyrosetta"
    if method == "lightdock_gso_ppi":
        return "lightdock"
    return method.replace("_", "-")


def _as_repo_path(ctx, path):
    return ctx.relative_to_repo(path)


def _file_entry(ctx, path, source_path):
    entry = {"path": _as_repo_path(ctx, path)}
    entry.update(describe_file(source_path))
    return entry


def source_file_entries(ctx, input_root, contract, contract_path):
    entries = {"contract": _file_entry(ctx, contract_path, contract_path)}
    receptor_path = resolve_input(input_root, contract["receptor"]["path"])
    primary_path = resolve_input(input_root, contract["myo1d"]["primary_construct"]["path"])
    comparator_path = resolve_input(input_root, contract["myo1d"]["comparator_construct"]["path"])
    membrane_path = resolve_input(input_root, contract["membrane_frame"]["path"])
    entries.update(
        {
            "receptor": _file_entry(ctx, receptor_path, receptor_path),
            "myo1d_primary": _file_entry(ctx, primary_path, primary_path),
            "myo1d_comparator": _file_entry(ctx, comparator_path, comparator_path),
            "membrane_frame": _file_entry(ctx, membrane_path, membrane_path),
        }
    )
    for item in contract.get("negative_controls", []):
        control_path = resolve_input(input_root, item["path"])
        entries["negative_control:{0}".format(item['id'])] = _file_entry(ctx, control_path, control_path)
    return entries


def negative_control_quarantine_rows(contract):
    rows = []
    for item in contract.get("negative_controls", []):
        control_id = item.get("id", "")
        if "962_1006" in control_id or "terminal_artifact" in control_id:
            rows.append(
                {
                    "input_id": control_id,
                    "severity": "WARN",
                    "classification": "quarantined_regression_fixture",
                    "code": "myo1d_terminal_artifact_no_jobs",
                    "message": "962-start MYO1D terminal-artifact control receives zero PPI jobs",
                    "planned_jobs": 0,
                    "production_primary_allowed": False,
                    "notes": "intentional regression input only",
                }
            )
        elif "single_chain" in control_id or "ambiguous" in control_id:
            rows.append(
                {
                    "input_id": control_id,
                    "severity": "WARN",
                    "classification": "quarantined_regression_fixture",
                    "code": "ambiguous_receptor_no_jobs",
                    "message": "ambiguous single-chain receptor control receives zero PPI jobs",
                    "planned_jobs": 0,
                    "production_primary_allowed": False,
                    "notes": "intentional regression input only",
                }
            )
        elif "v924r" in control_id.lower():
            rows.append(
                {
                    "input_id": control_id,
                    "severity": "WARN",
                    "classification": "quarantined_regression_fixture",
                    "code": "v924r_receptor_no_jobs",
                    "message": "V924R-like receptor control is not planned as production primary",
                    "planned_jobs": 0,
                    "production_primary_allowed": False,
                    "notes": "do not mutate or silently treat as WT",
                }
            )
    return rows


def readiness_blocker_rows(readiness_report):
    rows = []
    for item in readiness_report.get("blockers", []):
        rows.append(
            {
                "input_id": item.get("input_id", ""),
                "severity": item.get("severity", "FAIL"),
                "classification": item.get("classification", "readiness_blocker"),
                "code": item.get("code", ""),
                "message": item.get("message", ""),
                "planned_jobs": 0,
                "production_primary_allowed": False,
                "notes": "blocked by Task 5 real-input readiness",
            }
        )
    return rows


def accepted_protomers(contract, readiness_report):
    expected = list(contract["receptor"].get("expected_chains", []))
    observed = set(readiness_report.get("egfr", {}).get("chains", []))
    return [
        {"protomer_id": chain_id, "chain_ids": [chain_id]}
        for chain_id in expected
        if chain_id in observed
    ]


def build_job_specs(
    ctx,
    mode,
    profile,
    input_root,
    contract,
    readiness_report,
):
    receptor = contract["receptor"]
    myo = contract["myo1d"]
    primary = myo["primary_construct"]
    comparator = myo["comparator_construct"]
    receptor_path = resolve_input(input_root, receptor["path"])
    primary_path = resolve_input(input_root, primary["path"])
    comparator_path = resolve_input(input_root, comparator["path"])
    membrane_id = Path(contract["membrane_frame"]["path"]).stem
    jobs = []
    for method in METHODS:
        method_tag = _method_tag(method)
        for protomer in accepted_protomers(contract, readiness_report):
            for seed in PRIMARY_SEEDS:
                job_id = "ppi_smoke_primary_{0}_seed{1:03d}_{2}_spec".format(protomer['protomer_id'], seed, method_tag)
                jobs.append(
                    {
                        "job_id": job_id,
                        "execution_mode": "spec_only",
                        "execution_allowed": False,
                        "method": method,
                        "run_mode": mode,
                        "profile": profile,
                        "receptor_id": receptor["id"],
                        "receptor_path": _as_repo_path(ctx, receptor_path),
                        "receptor_assembly": "dimer",
                        "receptor_protomer_id": protomer["protomer_id"],
                        "receptor_chain_ids": protomer["chain_ids"],
                        "partner_construct_id": primary["id"],
                        "partner_role": "production_primary",
                        "partner_path": _as_repo_path(ctx, primary_path),
                        "partner_chain_id": "M",
                        "seed": seed,
                        "active_face_residues": list(myo["active_face"]),
                        "sheet12_support_residues": list(myo["sheet12_support"]),
                        "short_c_terminal_cap_residues": [number for number in myo["contact_monitoring_cap"] if number <= 1001],
                        "tail_noise_residues": [],
                        "score_bonus_allowed": False,
                        "key_residue_bonus_weight": 0.0,
                        "membrane_frame_id": membrane_id,
                        "future_output_dir": _as_repo_path(ctx, ctx.scratch_dir / "ppi" / job_id),
                        "blocked": False,
                        "block_reasons": [],
                    }
                )
            for seed in COMPARATOR_SEEDS:
                job_id = "ppi_smoke_noise_{0}_seed{1:03d}_{2}_spec".format(protomer['protomer_id'], seed, method_tag)
                jobs.append(
                    {
                        "job_id": job_id,
                        "execution_mode": "spec_only",
                        "execution_allowed": False,
                        "method": method,
                        "run_mode": mode,
                        "profile": profile,
                        "receptor_id": receptor["id"],
                        "receptor_path": _as_repo_path(ctx, receptor_path),
                        "receptor_assembly": "dimer",
                        "receptor_protomer_id": protomer["protomer_id"],
                        "receptor_chain_ids": protomer["chain_ids"],
                        "partner_construct_id": comparator["id"],
                        "partner_role": "noise_monitor",
                        "partner_path": _as_repo_path(ctx, comparator_path),
                        "partner_chain_id": "M",
                        "seed": seed,
                        "active_face_residues": list(myo["active_face"]),
                        "sheet12_support_residues": list(myo["sheet12_support"]),
                        "short_c_terminal_cap_residues": [number for number in myo["contact_monitoring_cap"] if number <= 1001],
                        "tail_noise_residues": [number for number in myo["contact_monitoring_cap"] if number >= 998],
                        "score_bonus_allowed": False,
                        "key_residue_bonus_weight": 0.0,
                        "membrane_frame_id": membrane_id,
                        "future_output_dir": _as_repo_path(ctx, ctx.scratch_dir / "ppi" / job_id),
                        "blocked": False,
                        "block_reasons": [],
                    }
                )
    return jobs


def build_sampling_plan(
    ctx,
    mode,
    profile,
    input_root,
    contract,
    readiness_report,
    jobs,
    blockers,
):
    primary_count = sum(1 for job in jobs if job["partner_role"] == "production_primary")
    comparator_count = sum(1 for job in jobs if job["partner_role"] == "noise_monitor")
    return {
        "planner_version": PLANNER_VERSION,
        "execution_mode": "spec_only",
        "execution_allowed": False,
        "no_docking_executed": True,
        "run_id": ctx.run_id,
        "mode": mode,
        "profile": profile,
        "input_root": str(input_root),
        "methods": METHODS,
        "smoke_budget": {
            "protomers": contract["receptor"].get("expected_chains", []),
            "primary_seeds": PRIMARY_SEEDS,
            "comparator_noise_monitor_seeds": COMPARATOR_SEEDS,
        },
        "receptor": {
            "receptor_id": contract["receptor"]["id"],
            "assembly": "dimer",
            "accepted_protomers": accepted_protomers(contract, readiness_report),
            "readiness_status": readiness_report.get("egfr", {}).get("status"),
        },
        "myo1d": {
            "primary_construct": contract["myo1d"]["primary_construct"]["id"],
            "comparator_construct": contract["myo1d"]["comparator_construct"]["id"],
            "comparator_role": "noise_monitor",
            "active_face_residues": list(contract["myo1d"]["active_face"]),
            "sheet12_support_residues": list(contract["myo1d"]["sheet12_support"]),
            "short_c_terminal_cap": [number for number in contract["myo1d"]["contact_monitoring_cap"] if number <= 1001],
            "extended_tail_noise_zone": [number for number in contract["myo1d"]["contact_monitoring_cap"] if number >= 998],
            "score_bonus_allowed": False,
            "key_residue_bonus_weight": 0.0,
        },
        "job_count_total": len(jobs),
        "job_count_primary": primary_count,
        "job_count_comparator_noise_monitor": comparator_count,
        "job_count_blocked": len(blockers),
        "jobs": jobs,
        "blocked_or_quarantined_inputs": blockers,
    }


def job_audit_rows(jobs, blockers):
    rows = []
    for job in jobs:
        rows.append(
            {
                "record_type": "job_spec",
                "job_id": job["job_id"],
                "method": job["method"],
                "receptor_id": job["receptor_id"],
                "receptor_protomer_id": job["receptor_protomer_id"],
                "partner_construct_id": job["partner_construct_id"],
                "partner_role": job["partner_role"],
                "seed": job["seed"],
                "status": "PLANNED_SPEC_ONLY",
                "execution_allowed": job["execution_allowed"],
                "notes": "future job spec only; no engine execution",
            }
        )
    for item in blockers:
        rows.append(
            {
                "record_type": "blocked_or_quarantined_input",
                "job_id": "",
                "method": "",
                "receptor_id": item.get("input_id", ""),
                "receptor_protomer_id": "",
                "partner_construct_id": item.get("input_id", ""),
                "partner_role": "blocked",
                "seed": "",
                "status": item.get("severity", "WARN"),
                "execution_allowed": False,
                "notes": item.get("message", ""),
            }
        )
    return rows
