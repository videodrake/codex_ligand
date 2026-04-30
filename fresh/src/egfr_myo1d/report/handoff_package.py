"""M3 handoff manifest generation."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from egfr_myo1d.core.run_context import RunContext


HANDOFF_FIELDS = ["artifact_role", "relative_path", "exists", "sha256", "rows_or_size", "required_for_handoff", "confidentiality_scanned", "notes"]

REQUIRED_ARTIFACTS = {
    "final_candidate_hypotheses": "phase3_compounds/tables/final_m3_candidate_hypotheses.csv",
    "pocket_compound_evidence_table": "phase3_compounds/tables/pocket_compound_evidence_table.csv",
    "rejected_candidate_reasons": "phase3_compounds/tables/rejected_candidate_reasons.csv",
    "compound_anchor_convergence": "phase3_compounds/tables/compound_anchor_convergence.csv",
    "compound_pose_clusters": "phase3_compounds/tables/compound_pose_clusters.csv",
    "compound_pose_attribution": "phase3_compounds/tables/compound_pose_attribution.csv",
    "final_candidate_gate_qc_json": "phase3_compounds/qc/final_candidate_gate_qc.json",
    "evidence_integration_qc": "phase3_compounds/qc/evidence_integration_qc.csv",
    "candidate_tiering_qc": "phase3_compounds/qc/candidate_tiering_qc.csv",
    "compound_manifest_public": "phase3_compounds/manifests/compound_manifest_public.csv",
    "vina_job_manifest": "phase3_compounds/manifests/vina_job_manifest.csv",
    "docking_box_manifest": "phase3_compounds/manifests/docking_box_manifest.csv",
    "milestone3_summary": "report/milestone3_summary.md",
    "final_candidate_summary": "report/final_candidate_summary.md",
    "candidate_hypothesis_summary": "report/candidate_hypothesis_summary.md",
    "qc_summary": "report/qc_summary.md",
    "methods_traceability": "report/methods_traceability.md",
    "reviewer_risk_notes": "report/reviewer_risk_notes.md",
    "cleanup_report": "report/cleanup_report.json",
    "phase_status": "logs/phase_status.jsonl",
    "job_status": "logs/job_status.jsonl",
    "master_log": "logs/master.log",
    "phase3_log": "logs/phase3_compounds.log",
    "error_summary": "logs/errors/error_summary.txt",
    "failed_jobs": "logs/errors/failed_jobs.csv",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rows_or_size(path: Path) -> str:
    if not path.exists():
        return "0"
    if path.suffix.lower() == ".csv":
        try:
            with path.open("r", encoding="utf-8", newline="") as handle:
                return str(max(0, sum(1 for _ in handle) - 1))
        except OSError:
            return str(path.stat().st_size)
    return str(path.stat().st_size)


def build_handoff_manifest(
    ctx: RunContext,
    *,
    m2_run_id: str,
    profile: str | None,
    overall_status: str,
    post_m3_ready: bool,
    m4_ready: bool,
    confidentiality_pass: bool,
    claim_safety_pass: bool,
    cleanup_safety_pass: bool,
    blockers: list[str],
    warnings: list[str],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for role, rel in REQUIRED_ARTIFACTS.items():
        path = ctx.run_dir / rel
        exists = path.is_file()
        required = role in {
            "final_candidate_hypotheses",
            "pocket_compound_evidence_table",
            "final_candidate_gate_qc_json",
            "milestone3_summary",
            "cleanup_report",
            "phase_status",
            "job_status",
        }
        rows.append(
            {
                "artifact_role": role,
                "relative_path": rel,
                "exists": "true" if exists else "false",
                "sha256": sha256(path) if exists else "",
                "rows_or_size": rows_or_size(path) if exists else "0",
                "required_for_handoff": "true" if required else "false",
                "confidentiality_scanned": "true",
                "notes": "present" if exists else ("missing required" if required else "missing optional"),
            }
        )
    required_present = all(row["exists"] == "true" for row in rows if row["required_for_handoff"] == "true")
    summary = {
        "schema_version": "m3_handoff_manifest_v1",
        "run_id": ctx.run_id,
        "m2_run_id": m2_run_id,
        "profile": profile,
        "generated_at": None,
        "overall_status": overall_status,
        "post_m3_ready": post_m3_ready,
        "m4_ready": m4_ready,
        "required_artifacts_present": required_present,
        "confidentiality_pass": confidentiality_pass,
        "claim_safety_pass": claim_safety_pass,
        "cleanup_safety_pass": cleanup_safety_pass,
        "artifacts": rows,
        "blockers": blockers,
        "warnings": warnings,
        "next_stage": "Milestone 4 - fresh compound library or focused analogue expansion",
    }
    return rows, summary


def write_handoff_manifests(ctx: RunContext, rows: list[dict[str, Any]], summary: dict[str, Any], generated_at: str) -> tuple[Path, Path, Path, Path]:
    summary["generated_at"] = generated_at
    root_csv = ctx.run_dir / "report" / "milestone3_handoff_manifest.csv"
    root_json = ctx.run_dir / "report" / "milestone3_handoff_manifest.json"
    phase_csv = ctx.run_dir / "phase3_compounds" / "reports" / "milestone3_handoff_manifest.csv"
    phase_json = ctx.run_dir / "phase3_compounds" / "reports" / "milestone3_handoff_manifest.json"
    for path in [root_csv, root_json, phase_csv, phase_json]:
        ctx.require_within_run_dir(path).parent.mkdir(parents=True, exist_ok=True)
    for path in [root_csv, phase_csv]:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=HANDOFF_FIELDS, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    for path in [root_json, phase_json]:
        path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return root_csv, root_json, phase_csv, phase_json
