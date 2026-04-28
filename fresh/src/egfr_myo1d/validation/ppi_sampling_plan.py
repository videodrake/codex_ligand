"""Task 6 PPI sampling-plan validation and serialization."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from egfr_myo1d.core.logging_utils import append_job_status, append_phase_status
from egfr_myo1d.core.manifest import initialize_manifests, now_iso, write_json
from egfr_myo1d.core.run_context import RunContext, RunContextError
from egfr_myo1d.io.hashing import describe_file
from egfr_myo1d.planning.pose_qc_policy import audit_pose_policy, build_pose_acceptance_policy
from egfr_myo1d.planning.ppi_sampling import (
    PLANNER_VERSION,
    build_job_specs,
    build_sampling_plan,
    job_audit_rows,
    negative_control_quarantine_rows,
    readiness_blocker_rows,
    source_file_entries,
)
from egfr_myo1d.validation.prepared_inputs import read_contract
from egfr_myo1d.validation.real_inputs import NOT_IMPLEMENTED as TASK5_NOT_IMPLEMENTED
from egfr_myo1d.validation.real_inputs import validate_real_inputs


NON_GOALS_CONFIRMED = [
    "no_ppi_docking_executed",
    "no_pyrosetta_import_or_execution",
    "no_lightdock_import_or_execution",
    "no_vina_execution",
    "no_fpocket_or_p2rank_execution",
    "no_compound_docking_or_scoring",
    "no_qsub_pbs_generation",
    "no_cleanup_deletion",
    "no_mutation_repair_or_renumbering",
]


JOB_CSV_FIELDS = [
    "job_id", "execution_mode", "execution_allowed", "method", "run_mode", "profile",
    "receptor_id", "receptor_path", "receptor_assembly", "receptor_protomer_id",
    "receptor_chain_ids", "partner_construct_id", "partner_role", "partner_path",
    "partner_chain_id", "seed", "active_face_residues", "sheet12_support_residues",
    "short_c_terminal_cap_residues", "tail_noise_residues", "score_bonus_allowed",
    "key_residue_bonus_weight", "membrane_frame_id", "future_output_dir", "blocked",
    "block_reasons",
]
AUDIT_FIELDS = [
    "record_type", "job_id", "method", "receptor_id", "receptor_protomer_id",
    "partner_construct_id", "partner_role", "seed", "status", "execution_allowed", "notes",
]
POSE_AUDIT_FIELDS = ["pose_class", "decision", "threshold_summary", "status", "notes"]
BLOCKER_FIELDS = [
    "input_id", "severity", "classification", "code", "message",
    "planned_jobs", "production_primary_allowed", "notes",
]


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str], ctx: RunContext) -> None:
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    with safe_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            encoded = {}
            for key, value in row.items():
                if isinstance(value, (list, dict)):
                    encoded[key] = json.dumps(value, sort_keys=True)
                else:
                    encoded[key] = value
            writer.writerow(encoded)


def write_jsonl(path: Path, rows: list[dict[str, Any]], ctx: RunContext) -> None:
    safe_path = ctx.require_within_run_dir(path)
    safe_path.parent.mkdir(parents=True, exist_ok=True)
    with safe_path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")


def _resolve_root(ctx: RunContext, input_root: Path) -> Path:
    return input_root.resolve() if input_root.is_absolute() else (ctx.repo_root / input_root).resolve()


def _failure_outputs(
    ctx: RunContext,
    mode: str,
    profile: str,
    input_root: Path,
    contract_path: Path | None,
    message: str,
) -> dict[str, Any]:
    ppi_dir = ctx.require_within_run_dir(ctx.run_dir / "prepared" / "ppi")
    ppi_dir.mkdir(parents=True, exist_ok=True)
    blocker = {
        "input_id": "contract",
        "severity": "FAIL",
        "classification": "planning_blocker",
        "code": "contract_error",
        "message": message,
        "planned_jobs": 0,
        "production_primary_allowed": False,
        "notes": "Task 6 requires a readable PPI input contract",
    }
    policy = {
        "policy_version": "task6_spec_only_v0_1",
        "current_task_executes_pose_classification": False,
        "future_pose_classes": [],
        "score_bonus_allowed": False,
        "non_scoring_notice": "contract unavailable; no pose policy could be specialized",
    }
    manifest = {
        "run_id": ctx.run_id,
        "mode": mode,
        "profile": profile,
        "input_root": str(input_root),
        "contract_path": str(contract_path) if contract_path else None,
        "created_at": now_iso(),
        "command_line": None,
        "planner_version": PLANNER_VERSION,
        "source_contract_sha256": None,
        "source_file_sha256_entries": {},
        "job_count_total": 0,
        "job_count_primary": 0,
        "job_count_comparator_noise_monitor": 0,
        "job_count_blocked": 1,
        "execution_mode": "spec_only",
        "no_docking_executed": True,
    }
    outputs = {
        "ppi_sampling_plan": ctx.relative_to_repo(ppi_dir / "ppi_sampling_plan.json"),
        "ppi_job_specs_jsonl": ctx.relative_to_repo(ppi_dir / "ppi_job_specs.jsonl"),
        "ppi_job_specs_csv": ctx.relative_to_repo(ppi_dir / "ppi_job_specs.csv"),
        "pose_acceptance_policy": ctx.relative_to_repo(ppi_dir / "pose_acceptance_policy.json"),
        "ppi_sampling_plan_manifest": ctx.relative_to_repo(ctx.manifest_dir / "ppi_sampling_plan_manifest.json"),
        "ppi_sampling_plan_report": ctx.relative_to_repo(ctx.manifest_dir / "ppi_sampling_plan_report.json"),
        "ppi_sampling_plan_audit": ctx.relative_to_repo(ctx.qc_dir / "ppi_sampling_plan_audit.csv"),
        "ppi_pose_qc_policy_audit": ctx.relative_to_repo(ctx.qc_dir / "ppi_pose_qc_policy_audit.csv"),
        "ppi_sampling_blockers": ctx.relative_to_repo(ctx.qc_dir / "ppi_sampling_blockers.csv"),
    }
    report = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "status": "FAIL",
        "warnings": [],
        "blockers": [blocker],
        "accepted_receptor_count": 0,
        "accepted_partner_count": 0,
        "planned_job_count": 0,
        "blocked_job_count": 1,
        "outputs": outputs,
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
        "not_implemented": TASK5_NOT_IMPLEMENTED,
    }
    write_json(ppi_dir / "ppi_sampling_plan.json", {"jobs": [], "blocked_or_quarantined_inputs": [blocker]}, ctx)
    write_jsonl(ppi_dir / "ppi_job_specs.jsonl", [], ctx)
    write_csv(ppi_dir / "ppi_job_specs.csv", [], JOB_CSV_FIELDS, ctx)
    write_json(ppi_dir / "pose_acceptance_policy.json", policy, ctx)
    write_json(ctx.manifest_dir / "ppi_sampling_plan_manifest.json", manifest, ctx)
    write_json(ctx.manifest_dir / "ppi_sampling_plan_report.json", report, ctx)
    write_csv(ctx.qc_dir / "ppi_sampling_plan_audit.csv", [], AUDIT_FIELDS, ctx)
    write_csv(ctx.qc_dir / "ppi_pose_qc_policy_audit.csv", [], POSE_AUDIT_FIELDS, ctx)
    write_csv(ctx.qc_dir / "ppi_sampling_blockers.csv", [blocker], BLOCKER_FIELDS, ctx)
    append_phase_status(ctx, "plan-ppi-sampling", "FAIL", message, {"mode": mode, "profile": profile})
    append_job_status(ctx, "plan_ppi_sampling_local", "FAIL", details={"message": message})
    return report


def plan_ppi_sampling(
    ctx: RunContext,
    mode: str,
    profile: str,
    input_root: Path,
    contract_path: Path | None = None,
    strict: bool = False,
    command_line: str | None = None,
) -> dict[str, Any]:
    initialize_manifests(ctx, mode)
    root = _resolve_root(ctx, input_root)
    ppi_dir = ctx.require_within_run_dir(ctx.run_dir / "prepared" / "ppi")
    ppi_dir.mkdir(parents=True, exist_ok=True)
    try:
        contract, resolved_contract_path = read_contract(root, contract_path)
    except (FileNotFoundError, RunContextError, ValueError) as exc:
        return _failure_outputs(ctx, mode, profile, root, contract_path, str(exc))

    readiness = validate_real_inputs(
        ctx,
        mode,
        profile,
        root,
        resolved_contract_path,
        strict=strict,
    )
    blocker_rows = readiness_blocker_rows(readiness)
    fixture_rows = negative_control_quarantine_rows(contract)
    blockers_for_csv = blocker_rows + fixture_rows
    if readiness["status"] == "FAIL":
        jobs: list[dict[str, Any]] = []
    else:
        jobs = build_job_specs(ctx, mode, profile, root, contract, readiness)

    policy = build_pose_acceptance_policy(contract)
    plan = build_sampling_plan(ctx, mode, profile, root, contract, readiness, jobs, blockers_for_csv)
    source_entries = source_file_entries(ctx, root, contract, resolved_contract_path)
    source_contract_sha = source_entries["contract"]["sha256"]
    primary_count = sum(1 for job in jobs if job["partner_role"] == "production_primary")
    comparator_count = sum(1 for job in jobs if job["partner_role"] == "noise_monitor")
    hard_blockers = [item for item in blocker_rows if item.get("severity") == "FAIL"]
    warnings = [item for item in blockers_for_csv if item.get("severity") == "WARN"]
    status = "FAIL" if hard_blockers else ("PASS_WITH_WARNINGS" if warnings or readiness["status"] == "WARN" else "PASS")
    phase_status = "FAIL" if status == "FAIL" else ("WARN" if status == "PASS_WITH_WARNINGS" else "PASS")

    outputs = {
        "ppi_sampling_plan": ctx.relative_to_repo(ppi_dir / "ppi_sampling_plan.json"),
        "ppi_job_specs_jsonl": ctx.relative_to_repo(ppi_dir / "ppi_job_specs.jsonl"),
        "ppi_job_specs_csv": ctx.relative_to_repo(ppi_dir / "ppi_job_specs.csv"),
        "pose_acceptance_policy": ctx.relative_to_repo(ppi_dir / "pose_acceptance_policy.json"),
        "ppi_sampling_plan_manifest": ctx.relative_to_repo(ctx.manifest_dir / "ppi_sampling_plan_manifest.json"),
        "ppi_sampling_plan_report": ctx.relative_to_repo(ctx.manifest_dir / "ppi_sampling_plan_report.json"),
        "ppi_sampling_plan_audit": ctx.relative_to_repo(ctx.qc_dir / "ppi_sampling_plan_audit.csv"),
        "ppi_pose_qc_policy_audit": ctx.relative_to_repo(ctx.qc_dir / "ppi_pose_qc_policy_audit.csv"),
        "ppi_sampling_blockers": ctx.relative_to_repo(ctx.qc_dir / "ppi_sampling_blockers.csv"),
    }
    manifest = {
        "run_id": ctx.run_id,
        "mode": mode,
        "profile": profile,
        "input_root": str(root),
        "contract_path": str(resolved_contract_path),
        "created_at": now_iso(),
        "command_line": command_line,
        "planner_version": PLANNER_VERSION,
        "source_contract_sha256": source_contract_sha,
        "source_file_sha256_entries": source_entries,
        "job_count_total": len(jobs),
        "job_count_primary": primary_count,
        "job_count_comparator_noise_monitor": comparator_count,
        "job_count_blocked": len(blockers_for_csv),
        "execution_mode": "spec_only",
        "no_docking_executed": True,
    }
    report = {
        "run_id": ctx.run_id,
        "created_at": now_iso(),
        "status": status,
        "warnings": warnings + [item for item in readiness.get("messages", []) if item.get("severity") == "WARN"],
        "blockers": hard_blockers,
        "accepted_receptor_count": 1 if jobs else 0,
        "accepted_partner_count": 2 if jobs else 0,
        "planned_job_count": len(jobs),
        "blocked_job_count": len(blockers_for_csv),
        "outputs": outputs,
        "readiness_status": readiness["status"],
        "non_goals_confirmed": NON_GOALS_CONFIRMED,
        "not_implemented": TASK5_NOT_IMPLEMENTED,
    }
    summary = [
        "# Task 6 PPI Sampling Plan Summary",
        "",
        f"- Status: {status}",
        f"- Planned spec-only jobs: {len(jobs)}",
        f"- Blocked/quarantined input records: {len(blockers_for_csv)}",
        "- Docking executed: no",
        "- Score bonuses allowed: no",
        "",
    ]

    write_json(ppi_dir / "ppi_sampling_plan.json", plan, ctx)
    write_jsonl(ppi_dir / "ppi_job_specs.jsonl", jobs, ctx)
    write_csv(ppi_dir / "ppi_job_specs.csv", jobs, JOB_CSV_FIELDS, ctx)
    write_json(ppi_dir / "pose_acceptance_policy.json", policy, ctx)
    write_json(ctx.manifest_dir / "ppi_sampling_plan_manifest.json", manifest, ctx)
    write_json(ctx.manifest_dir / "ppi_sampling_plan_report.json", report, ctx)
    write_csv(ctx.qc_dir / "ppi_sampling_plan_audit.csv", job_audit_rows(jobs, blockers_for_csv), AUDIT_FIELDS, ctx)
    write_csv(ctx.qc_dir / "ppi_pose_qc_policy_audit.csv", audit_pose_policy(policy), POSE_AUDIT_FIELDS, ctx)
    write_csv(ctx.qc_dir / "ppi_sampling_blockers.csv", blockers_for_csv, BLOCKER_FIELDS, ctx)
    safe_summary = ctx.require_within_run_dir(ctx.reports_dir / "task6_ppi_sampling_plan_summary.md")
    safe_summary.write_text("\n".join(summary), encoding="utf-8")
    append_phase_status(
        ctx,
        "plan-ppi-sampling",
        phase_status,
        f"PPI sampling plan completed with {status}",
        {"mode": mode, "profile": profile, "strict": strict, "job_count": len(jobs)},
    )
    append_job_status(
        ctx,
        "plan_ppi_sampling_local",
        phase_status,
        details={"message": "local pure-Python spec-only PPI sampling plan"},
    )
    return report
