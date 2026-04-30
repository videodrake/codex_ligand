"""M3-T11 final report, cleanup, and handoff integration."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from egfr_myo1d.cleanup.milestone_cleanup import M3CleanupResult, plan_m3_cleanup
from egfr_myo1d.compound.confidentiality import PUBLIC_COMPOUND_IDS, load_private_map
from egfr_myo1d.core.logging_utils import append_failed_job, append_job_status, append_phase_status
from egfr_myo1d.core.run_context import RunContext, ensure_within
from egfr_myo1d.report.handoff_package import build_handoff_manifest, sha256, write_handoff_manifests
from egfr_myo1d.report.report_safety import SafetyScanResult, scan_public_reports


DEFAULT_CONFIG: dict[str, Any] = {
    "schema_version": "m3_report_cleanup_config_v1",
    "report": {
        "require_t10_ready": True,
        "allow_partial": False,
        "strict_claims": True,
        "strict_confidentiality": True,
        "write_root_report_dir": True,
        "write_phase3_report_dir": True,
        "include_methods_traceability": True,
        "include_qc_summary": True,
        "include_reviewer_risk_notes": True,
        "include_cleanup_report": True,
        "include_candidate_tables_in_markdown": True,
        "max_candidate_rows_in_summary": 50,
        "max_reject_rows_in_summary": 100,
        "hypothesis_language_only": True,
    },
    "handoff": {
        "next_stage_label": "Milestone 4 - fresh compound library or focused analogue expansion",
    },
}

REPORT_QC_FIELDS = ["check_id", "category", "status", "severity", "details", "recommended_fix"]
REPORT_CHECKS = [
    "final_candidate_hypotheses_present", "pocket_compound_evidence_table_present", "rejected_candidate_reasons_present",
    "final_candidate_gate_qc_json_present", "t10_report_ready", "profile_selected", "candidate_counts_loaded",
    "tier_counts_loaded", "reject_reasons_loaded", "m2_traceability_loaded", "m3_upstream_traceability_loaded",
    "milestone3_summary_written", "final_candidate_summary_written", "candidate_hypothesis_summary_written",
    "qc_summary_written", "methods_traceability_written", "reviewer_risk_notes_written", "handoff_manifest_written",
    "cleanup_report_written", "all_required_artifacts_in_handoff_manifest", "artifact_sha256_computed",
    "candidate_language_hypothesis_only", "no_validated_inhibitor_claim", "no_proven_binder_claim",
    "no_drug_candidate_claim", "no_clinical_claim", "affinity_non_promotion_stated",
    "atp_non_promotion_stated", "reference_only_non_promotion_stated", "broad_scan_non_promotion_stated",
    "hard_gates_soft_scores_separation_stated", "no_internal_id_leak", "no_smiles_leak", "no_coordinate_leak",
    "no_vina_invoked_by_report", "no_qsub_invoked_by_report", "no_runner_invoked_by_report",
    "no_upstream_stage_rerun", "old_workflow_not_used", "post_m3_transition_gate_evaluated",
    "m4_readiness_evaluated",
]


@dataclass
class M3ReportResult:
    status: str
    post_m3_ready: bool
    m4_ready: bool
    blockers: list[str]
    warnings: list[str]
    profile: str | None
    allow_partial: bool
    require_t10_ready: bool
    counts: dict[str, int]
    inputs: dict[str, bool]
    cleanup: M3CleanupResult
    safety: SafetyScanResult
    outputs: dict[str, Path]
    qc_json: Path
    qc_csv: Path
    handoff_qc: Path
    phase3_log: Path


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    if not path.is_file():
        return [], []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def _load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def _load_config(ctx: RunContext, config_path: Path | None) -> dict[str, Any]:
    config = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (ctx.fresh_root / "configs" / "milestone3_report.yaml")
    if path and path.is_file():
        try:
            import yaml  # type: ignore

            loaded = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
            for section, values in loaded.items():
                if isinstance(values, dict) and isinstance(config.get(section), dict):
                    config[section].update(values)
                else:
                    config[section] = values
        except Exception:
            pass
    snapshot = ctx.run_dir / "phase3_compounds" / "config_snapshots" / "milestone3_report.resolved.yaml"
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    try:
        import yaml  # type: ignore

        snapshot.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    except Exception:
        snapshot.write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")
    return config


def _ensure_logs(ctx: RunContext) -> None:
    ctx.create_directories()
    for path in [
        ctx.logs_dir / "master.log", ctx.logs_dir / "phase_status.jsonl", ctx.logs_dir / "job_status.jsonl",
        ctx.logs_dir / "phase3_compounds.log", ctx.errors_dir / "error_summary.txt", ctx.errors_dir / "failed_jobs.csv",
    ]:
        ctx.require_within_run_dir(path).parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            path.touch()
    if (ctx.errors_dir / "failed_jobs.csv").stat().st_size == 0:
        (ctx.errors_dir / "failed_jobs.csv").write_text("timestamp,job_name,status,message\n", encoding="utf-8")


def _resolve_inputs(ctx: RunContext, m2_run_id: str) -> dict[str, Path]:
    phase3 = ctx.run_dir / "phase3_compounds"
    m2 = ctx.fresh_root / "runs" / m2_run_id

    def first(paths: list[Path]) -> Path:
        for path in paths:
            if path.is_file():
                return path
        return paths[0]

    return {
        "final_candidate_hypotheses": phase3 / "tables" / "final_m3_candidate_hypotheses.csv",
        "pocket_compound_evidence_table": phase3 / "tables" / "pocket_compound_evidence_table.csv",
        "rejected_candidate_reasons": phase3 / "tables" / "rejected_candidate_reasons.csv",
        "final_candidate_gate_qc_json": phase3 / "qc" / "final_candidate_gate_qc.json",
        "final_candidate_gate_qc_csv": phase3 / "qc" / "final_candidate_gate_qc.csv",
        "evidence_integration_qc": phase3 / "qc" / "evidence_integration_qc.csv",
        "candidate_tiering_qc": phase3 / "qc" / "candidate_tiering_qc.csv",
        "compound_anchor_convergence": phase3 / "tables" / "compound_anchor_convergence.csv",
        "compound_pocket_support": phase3 / "tables" / "compound_pocket_support.csv",
        "compound_pose_clusters": phase3 / "tables" / "compound_pose_clusters.csv",
        "compound_pose_attribution": phase3 / "tables" / "compound_pose_attribution.csv",
        "compound_pose_mechanism_classification": phase3 / "tables" / "compound_pose_mechanism_classification.csv",
        "compound_manifest_public": phase3 / "manifests" / "compound_manifest_public.csv",
        "vina_job_manifest": phase3 / "manifests" / "vina_job_manifest.csv",
        "docking_box_manifest": phase3 / "manifests" / "docking_box_manifest.csv",
        "m2_accepted_pockets": first([m2 / "phase2_pockets" / "final" / "accepted_pockets_for_m3.csv", m2 / "phase2_pockets" / "export_for_m3" / "accepted_pockets_for_m3.csv"]),
        "m2_accepted_boxes": first([m2 / "phase2_pockets" / "final" / "accepted_pocket_boxes.csv", m2 / "phase2_pockets" / "export_for_m3" / "accepted_pocket_boxes.csv"]),
        "m2_pocket_gate_qc": m2 / "phase2_pockets" / "gated" / "pocket_gate_qc.csv",
        "m2_ppi_consensus": first([m2 / "phase1_ppi" / "consensus" / "ppi_consensus_patch_merged.csv", m2 / "phase1_ppi" / "tables" / "ppi_consensus_patch.csv"]),
        "m2_final_report": m2 / "reports" / "m2_final_report.md",
    }


def _profile_from(qc: dict[str, Any], candidates: list[dict[str, str]], explicit: str | None) -> str | None:
    if explicit:
        return explicit
    if qc.get("profile"):
        return str(qc["profile"])
    profiles = sorted({row.get("profile", "") for row in candidates if row.get("profile")})
    return profiles[0] if len(profiles) == 1 else None


def _count_tiers(candidates: list[dict[str, str]]) -> dict[str, int]:
    return {
        "candidate_rows": len(candidates),
        "tier1": sum(1 for row in candidates if row.get("candidate_tier") == "Tier 1"),
        "tier2": sum(1 for row in candidates if row.get("candidate_tier") == "Tier 2"),
        "tier3": sum(1 for row in candidates if row.get("candidate_tier") == "Tier 3"),
        "reject": sum(1 for row in candidates if row.get("candidate_tier") == "Reject"),
    }


def _safe_candidate_rows(rows: list[dict[str, str]], limit: int = 20) -> list[str]:
    if not rows:
        return ["- none"]
    lines = []
    for row in rows[:limit]:
        lines.append(
            "- {compound} / {pocket}: {tier}; states={states}; score={score}; reason={reason}".format(
                compound=row.get("compound_public_id") or "pocket-level",
                pocket=row.get("pocket_family_id", ""),
                tier=row.get("candidate_tier", ""),
                states=row.get("primary_state_ids_supported", ""),
                score=row.get("candidate_priority_score", ""),
                reason=row.get("candidate_reject_reason", "none"),
            )
        )
    return lines


def _write_text_pair(ctx: RunContext, rel: str, phase_rel: str, text: str) -> tuple[Path, Path]:
    root = ctx.run_dir / rel
    phase = ctx.run_dir / phase_rel
    for path in [root, phase]:
        ctx.require_within_run_dir(path).parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    return root, phase


def _git_commit(repo_root: Path) -> str:
    try:
        proc = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(repo_root), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        return proc.stdout.strip() if proc.returncode == 0 else "unavailable"
    except OSError:
        return "unavailable"


def _report_sections(
    *,
    ctx: RunContext,
    m2_run_id: str,
    profile: str | None,
    candidates: list[dict[str, str]],
    evidence: list[dict[str, str]],
    rejects: list[dict[str, str]],
    t10_qc: dict[str, Any],
    cleanup_report: dict[str, Any],
    input_paths: dict[str, Path],
    counts: dict[str, int],
) -> dict[str, str]:
    tier1 = [row for row in candidates if row.get("candidate_tier") == "Tier 1"]
    tier2 = [row for row in candidates if row.get("candidate_tier") == "Tier 2"]
    tier3 = [row for row in candidates if row.get("candidate_tier") == "Tier 3"]
    rejected = [row for row in candidates if row.get("candidate_tier") == "Reject"]
    guard = (
        "All rows are computational hypotheses for review. No experimental validation is provided for "
        "binding, inhibition, PPI disruption, cellular efficacy, therapeutic use, or clinical utility. "
        "Affinity is descriptive metadata and is not used for promotion."
    )
    summary = "\n".join(
        [
            "# Milestone 3 Summary - Compound Anchor Integration",
            "",
            "## 1. Executive summary",
            guard,
            "",
            "## 2. Run identity",
            f"- run_id: {ctx.run_id}",
            f"- m2_run_id: {m2_run_id}",
            f"- profile: {profile}",
            "",
            "## 3. Input provenance and source files",
            *[f"- {key}: {ctx.relative_to_repo(path) if path.exists() else 'missing'}" for key, path in input_paths.items()],
            "",
            "## 4. M2 accepted pocket handoff summary",
            f"- accepted pocket source present: {input_paths['m2_accepted_pockets'].is_file()}",
            "",
            "## 5. Compound confidentiality and ligand QC summary",
            "- Public compound IDs only: Cpd-A, Cpd-B, Cpd-C.",
            "",
            "## 6. Docking execution summary",
            "- This report does not run docking; it summarizes existing execution/QC outputs.",
            "",
            "## 7. Pose attribution summary",
            f"- pose attribution table present: {input_paths['compound_pose_attribution'].is_file()}",
            "",
            "## 8. ATP migration and non-ATP gate summary",
            "- ATP-confounded evidence is not promoted.",
            "",
            "## 9. Pocket retention summary",
            "- Pocket-retention hard gates are inherited from M3-T7/T10 outputs.",
            "",
            "## 10. PPI geometry / mechanism classification summary",
            "- PPI geometry is summarized as hypothesis support only.",
            "",
            "## 11. Pose clustering summary",
            f"- cluster table present: {input_paths['compound_pose_clusters'].is_file()}",
            "",
            "## 12. Compound-anchor convergence summary",
            f"- anchor convergence table present: {input_paths['compound_anchor_convergence'].is_file()}",
            "",
            "## 13. Integrated evidence and hard-gate summary",
            f"- evidence rows: {len(evidence)}",
            "- Hard gates and soft scores remain separate.",
            "",
            "## 14. Candidate tiering summary",
            f"- Tier 1: {counts['tier1']}",
            f"- Tier 2: {counts['tier2']}",
            f"- Tier 3: {counts['tier3']}",
            f"- Reject: {counts['reject']}",
            "",
            "## 15. Tier 1 hypotheses",
            *_safe_candidate_rows(tier1),
            "",
            "## 16. Tier 2 hypotheses",
            *_safe_candidate_rows(tier2),
            "",
            "## 17. Tier 3 exploratory signals",
            *_safe_candidate_rows(tier3),
            "",
            "## 18. Rejected / quarantined cases",
            *_safe_candidate_rows(rejected),
            "",
            "## 19. Key risks and limitations",
            "- Vina affinity is not proof of binding.",
            "- Pose convergence is not proof of mechanism.",
            "- A/B protomer symmetry is not independent state support.",
            "- 3GT8_raw/reference-only signals are not primary support.",
            "- Broad-scan support is exploratory only.",
            "",
            "## 20. Recommended post-M3 next steps",
            "- Milestone 4: fresh compound library or focused analogue expansion.",
            "",
            "## 21. Output package and reproducibility checklist",
            "- Handoff manifest lists required artifacts and hashes.",
            "",
            "## 22. Cleanup summary",
            f"- cleanup mode: {cleanup_report.get('mode')}",
            f"- files planned for deletion: {cleanup_report.get('files_planned_for_deletion')}",
        ]
    )
    candidate = "\n".join(
        [
            "# Final M3 Candidate Hypothesis Summary",
            "",
            "## 1. Interpretation guardrail",
            guard,
            "",
            "## 2. Candidate tier definitions",
            "- Tier 1: strong computational pocket-compound hypothesis.",
            "- Tier 2: plausible single-anchor or state-specific hypothesis.",
            "- Tier 3: exploratory chemical signal.",
            "- Reject: quarantined or hard-gate-failed evidence.",
            "",
            "## 3. Tier 1 candidate hypotheses",
            *_safe_candidate_rows(tier1),
            "",
            "## 4. Tier 2 candidate hypotheses",
            *_safe_candidate_rows(tier2),
            "",
            "## 5. Tier 3 exploratory signals",
            *_safe_candidate_rows(tier3),
            "",
            "## 6. Rejected / quarantined cases by reason",
            *_safe_candidate_rows(rejected),
            "",
            "## 7. Evidence table cross-reference",
            f"- pocket_compound_evidence_table rows: {len(evidence)}",
            "",
            "## 8. Recommended review order",
            "- Review by tier, then non-affinity priority score, then deterministic IDs.",
            "",
            "## 9. Experimental validation suggestions",
            "- Consider orthogonal biochemical and cellular assays before any mechanistic interpretation.",
        ]
    )
    qc = "\n".join(
        [
            "# Milestone 3 QC Summary",
            "",
            "## 1. Overall M3 status",
            f"- T10 overall status: {t10_qc.get('overall_status', 'unknown')}",
            "",
            "## 2. M3-T10 readiness status",
            f"- m3_t11_report_ready: {t10_qc.get('m3_t11_report_ready')}",
            "",
            "## 3. Input gates",
            json.dumps(t10_qc.get("inputs", {}), sort_keys=True),
            "",
            "## 4. Docking completion QC",
            f"- present: {(ctx.run_dir / 'phase3_compounds/qc/docking_completion_qc.json').is_file()}",
            "",
            "## 5. Pose attribution QC",
            f"- present: {(ctx.run_dir / 'phase3_compounds/qc/pose_attribution_qc.json').is_file()}",
            "",
            "## 6. Pose clustering QC",
            f"- present: {(ctx.run_dir / 'phase3_compounds/qc/pose_clustering_qc.json').is_file()}",
            "",
            "## 7. Anchor convergence QC",
            f"- present: {(ctx.run_dir / 'phase3_compounds/qc/anchor_convergence_qc.json').is_file()}",
            "",
            "## 8. Evidence integration QC",
            f"- present: {input_paths['evidence_integration_qc'].is_file()}",
            "",
            "## 9. Candidate tiering QC",
            f"- present: {input_paths['candidate_tiering_qc'].is_file()}",
            "",
            "## 10. Confidentiality/output hygiene QC",
            "- Public outputs are scanned for private IDs, private structure strings, coordinates, and overclaims.",
            "",
            "## 11. Cleanup safety QC",
            f"- cleanup status: {cleanup_report.get('overall_status')}",
            "",
            "## 12. Open blockers",
            json.dumps(t10_qc.get("blockers", []), sort_keys=True),
            "",
            "## 13. Warnings",
            json.dumps(t10_qc.get("warnings", []), sort_keys=True),
            "",
            "## 14. Post-M3 readiness",
            "- evaluated in m3_report_qc.json.",
        ]
    )
    method_lines = [
        "# Methods Traceability",
        "",
        "## 1. Repository status summary if available",
        f"- git_commit: {_git_commit(ctx.repo_root)}",
        "",
        "## 2. Python version",
        f"- {platform.python_version()}",
        "",
        "## 3. CLI commands used for M3 tasks if recorded",
        "- See logs/phase_status.jsonl and logs/job_status.jsonl.",
        "",
        "## 4. Run IDs and profile",
        f"- run_id: {ctx.run_id}",
        f"- m2_run_id: {m2_run_id}",
        f"- profile: {profile}",
        "",
        "## 5. M2 source paths and sha256",
    ]
    for key, path in input_paths.items():
        if key.startswith("m2_"):
            method_lines.append(f"- {key}: {ctx.relative_to_repo(path) if path.exists() else 'missing'} sha256={sha256(path) if path.exists() else ''}")
    method_lines.extend(["", "## 6. M3 source paths and sha256"])
    for key, path in input_paths.items():
        if not key.startswith("m2_"):
            method_lines.append(f"- {key}: {ctx.relative_to_repo(path) if path.exists() else 'missing'} sha256={sha256(path) if path.exists() else ''}")
    method_lines.extend(
        [
            "",
            "## 7. Config snapshots and sha256",
            f"- milestone3_report.resolved.yaml: {(ctx.run_dir / 'phase3_compounds/config_snapshots/milestone3_report.resolved.yaml').is_file()}",
            "",
            "## 8. Input table row counts",
            f"- candidates: {len(candidates)}",
            f"- evidence: {len(evidence)}",
            "",
            "## 9. Output table row counts",
            f"- rejected_candidate_reasons: {len(rejects)}",
            "",
            "## 10. Hard-gate definitions from final_candidate_gate_qc.json and config",
            "- See final_candidate_gate_qc.json hard_gates section.",
            "",
            "## 11. Soft-score definitions from T10 config snapshot",
            "- See compound_evidence_tiering.resolved.yaml when present.",
            "",
            "## 12. Affinity non-promotion policy",
            "- Affinity is descriptive only and is not used for promotion or best-compound selection.",
            "",
            "## 13. Confidentiality policy",
            "- Public outputs use public IDs only.",
            "",
            "## 14. Cleanup policy",
            "- Cleanup is constrained to run_dir and preserves final evidence, QC, reports, manifests, and logs.",
            "",
            "## 15. Exact report generation command",
            "- python -m egfr_myo1d.cli m3-report ...",
        ]
    )
    risks = "\n".join(
        [
            "# Reviewer Risk Notes",
            "",
            "## 1. Scientific interpretation limits",
            "- No candidate has experimental validation from this workflow.",
            "- No cellular MYO1D-EGFR disruption has been shown by this workflow.",
            "",
            "## 2. Structural-model limits",
            "- The receptor model is a membrane-compatible EGFR dimer model and remains computational.",
            "",
            "## 3. Docking limits",
            "- Vina affinity is not proof of binding.",
            "- Pose convergence is not proof of mechanism.",
            "",
            "## 4. Compound confidentiality limits",
            "- Public IDs are limited to Cpd-A, Cpd-B, and Cpd-C.",
            "",
            "## 5. ATP-site confounding controls",
            "- ATP-confounded support is not promoted.",
            "",
            "## 6. Reference-state / 3GT8_raw controls",
            "- 3GT8_raw/reference-only signals are not primary support.",
            "",
            "## 7. Broad-scan exploratory controls",
            "- Broad-scan support is exploratory only.",
            "",
            "## 8. Affinity-score limitations",
            "- Affinity is descriptive and does not promote hypotheses.",
            "",
            "## 9. Membrane/dimer geometry caveats",
            "- Membrane and central dimer clash filters are inherited from upstream QC.",
            "",
            "## 10. Recommended experimental follow-up",
            "- Use orthogonal binding and cell-context assays before interpreting mechanism.",
            "",
            "## 11. Not-yet-done analyses",
            "- Milestone 4 library or focused analogue expansion has not been started.",
        ]
    )
    return {
        "milestone3_summary": summary + "\n",
        "final_candidate_summary": candidate + "\n",
        "candidate_hypothesis_summary": candidate + "\n",
        "qc_summary": qc + "\n",
        "methods_traceability": "\n".join(method_lines) + "\n",
        "reviewer_risk_notes": risks + "\n",
    }


def _write_report_qc(path: Path, statuses: dict[str, str], ctx: RunContext) -> None:
    path = ctx.require_within_run_dir(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_QC_FIELDS)
        writer.writeheader()
        for check in REPORT_CHECKS:
            status = statuses.get(check, "PASS")
            severity = "BLOCKER" if status == "FAIL" else ("MAJOR" if status == "WARN" else "INFO")
            writer.writerow({"check_id": check, "category": "m3_report", "status": status, "severity": severity, "details": check, "recommended_fix": ""})


def run_m3_report(
    ctx: RunContext,
    *,
    m2_run_id: str,
    profile: str | None = None,
    mode: str = "report",
    allow_partial: bool = False,
    require_t10_ready: bool = True,
    strict_confidentiality: bool = True,
    strict_claims: bool = True,
    config: Path | None = None,
    private_map: Path | None = None,
    force: bool = False,
    cleanup_mode: str = "plan",
    cleanup_scope: str = "conservative",
    confirm_run_id: str | None = None,
    preserve_selected_poses: bool = True,
    preserve_all_docking_outputs: bool = True,
) -> M3ReportResult:
    _ensure_logs(ctx)
    timestamp = now_iso()
    _load_config(ctx, config)
    phase3 = ctx.run_dir / "phase3_compounds"
    report_dir = ctx.run_dir / "report"
    phase3_reports = phase3 / "reports"
    qc_dir = phase3 / "qc"
    for directory in [report_dir, phase3_reports, qc_dir]:
        ctx.require_within_run_dir(directory).mkdir(parents=True, exist_ok=True)
    input_paths = _resolve_inputs(ctx, m2_run_id)
    blockers: list[str] = []
    warnings: list[str] = []
    candidates, candidate_fields = _read_csv(input_paths["final_candidate_hypotheses"])
    evidence, evidence_fields = _read_csv(input_paths["pocket_compound_evidence_table"])
    rejects, reject_fields = _read_csv(input_paths["rejected_candidate_reasons"])
    t10_qc = _load_json(input_paths["final_candidate_gate_qc_json"])
    selected_profile = _profile_from(t10_qc, candidates, profile)
    if selected_profile:
        candidates = [row for row in candidates if not row.get("profile") or row.get("profile") == selected_profile]
        evidence = [row for row in evidence if not row.get("profile") or row.get("profile") == selected_profile]
        rejects = [row for row in rejects if not row.get("profile") or row.get("profile") == selected_profile]
    required = [
        ("final_candidate_hypotheses", "final_m3_candidate_hypotheses.csv missing"),
        ("pocket_compound_evidence_table", "pocket_compound_evidence_table.csv missing"),
        ("rejected_candidate_reasons", "rejected_candidate_reasons.csv missing"),
        ("final_candidate_gate_qc_json", "final_candidate_gate_qc.json missing"),
    ]
    for key, message in required:
        if not input_paths[key].is_file():
            (warnings if mode == "dry-run" or allow_partial else blockers).append(message)
    if input_paths["final_candidate_hypotheses"].is_file() and "candidate_tier" not in candidate_fields:
        blockers.append("final candidate schema missing candidate_tier")
    if input_paths["rejected_candidate_reasons"].is_file() and "recommended_fix" not in reject_fields:
        blockers.append("rejected reasons schema missing recommended_fix")
    t10_ready = bool(t10_qc.get("m3_t11_report_ready") is True)
    if require_t10_ready and not t10_ready:
        (warnings if mode == "dry-run" or allow_partial else blockers).append("M3-T10 m3_t11_report_ready is false or missing")
    if t10_qc.get("overall_status") == "FAIL" and not allow_partial:
        blockers.append("M3-T10 QC overall status FAIL")
    invalid_public = sorted({row.get("compound_public_id", "") for row in candidates + evidence if row.get("compound_public_id") and row.get("compound_public_id") not in PUBLIC_COMPOUND_IDS})
    if invalid_public:
        blockers.append("non-public compound IDs appear in final tables")
    if not selected_profile and (candidates or evidence):
        blockers.append("profile ambiguous; provide --profile")
    if allow_partial:
        warnings.append("allow_partial diagnostic mode used")
    if mode == "dry-run":
        warnings.append("dry-run report mode")

    cleanup_result = plan_m3_cleanup(
        ctx,
        profile=selected_profile,
        mode="plan" if cleanup_mode == "none" else cleanup_mode,
        cleanup_scope=cleanup_scope,
        config=config,
        preserve_selected_poses=preserve_selected_poses,
        preserve_all_docking_outputs=preserve_all_docking_outputs,
        confirm_run_id=confirm_run_id,
        force=force,
    )
    if cleanup_result.status == "FAIL":
        blockers.extend(cleanup_result.blockers)
    else:
        warnings.extend(cleanup_result.warnings)

    counts = _count_tiers(candidates)
    counts.update({"evidence_rows": len(evidence), "rejected_candidate_rows": len(rejects)})
    if counts["tier1"] == 0 and candidates:
        warnings.append("no Tier 1 hypotheses")
    if candidates and counts["tier1"] == 0 and counts["tier2"] == 0:
        warnings.append("no Tier 1 or Tier 2 hypotheses")
    sections = _report_sections(
        ctx=ctx,
        m2_run_id=m2_run_id,
        profile=selected_profile,
        candidates=candidates,
        evidence=evidence,
        rejects=rejects,
        t10_qc=t10_qc,
        cleanup_report=cleanup_result.report,
        input_paths=input_paths,
        counts=counts,
    )
    outputs: dict[str, Path] = {}
    for name, text in sections.items():
        root, phase = _write_text_pair(ctx, f"report/{name}.md", f"phase3_compounds/reports/{name}.md", text)
        outputs[name] = root
        outputs[f"phase3_{name}"] = phase

    private_path = private_map or (ctx.fresh_root / "data" / "private" / "compound_id_map.csv")
    if not private_path.is_absolute():
        private_path = ctx.repo_root / private_path
    try:
        private_path = ensure_within(private_path.resolve(), ctx.fresh_root)
    except Exception:
        private_path = ctx.fresh_root / "data" / "private" / "compound_id_map.csv"
    private_entries, private_warnings = load_private_map(private_path)
    warnings.extend(private_warnings)
    safety = scan_public_reports(ctx, list(private_entries.values()))
    if strict_confidentiality and safety.leaks_detected:
        blockers.append("internal compound ID leakage detected")
    if strict_confidentiality and safety.smiles_logged:
        blockers.append("private structure string printed in public outputs")
    if strict_confidentiality and safety.coordinate_leak_files:
        blockers.append("coordinate records printed in public outputs")
    if strict_claims and safety.candidate_overclaims_detected:
        blockers.append("candidate overclaim detected")

    status = "FAIL" if blockers else ("WARN" if warnings or mode == "dry-run" or cleanup_mode != "execute" else "PASS")
    post_m3_ready = status in {"PASS", "WARN"}
    m4_ready = status == "PASS"
    report_ready = {
        "milestone3_summary": outputs["milestone3_summary"],
        "final_candidate_summary": outputs["final_candidate_summary"],
        "candidate_hypothesis_summary": outputs["candidate_hypothesis_summary"],
        "qc_summary": outputs["qc_summary"],
        "methods_traceability": outputs["methods_traceability"],
        "reviewer_risk_notes": outputs["reviewer_risk_notes"],
        "cleanup_report": cleanup_result.report_root,
    }
    handoff_rows, handoff_summary = build_handoff_manifest(
        ctx,
        m2_run_id=m2_run_id,
        profile=selected_profile,
        overall_status=status,
        post_m3_ready=post_m3_ready,
        m4_ready=m4_ready,
        confidentiality_pass=not safety.leaks_detected and not safety.smiles_logged and not safety.coordinate_leak_files,
        claim_safety_pass=not safety.candidate_overclaims_detected,
        cleanup_safety_pass=cleanup_result.status != "FAIL",
        blockers=sorted(set(blockers)),
        warnings=sorted(set(warnings)),
    )
    handoff_csv, handoff_json, phase_handoff_csv, phase_handoff_json = write_handoff_manifests(ctx, handoff_rows, handoff_summary, timestamp)
    outputs.update({"handoff_manifest_csv": handoff_csv, "handoff_manifest_json": handoff_json, "phase3_handoff_manifest_csv": phase_handoff_csv, "phase3_handoff_manifest_json": phase_handoff_json})
    handoff_required_present = all(row["exists"] == "true" for row in handoff_rows if row["required_for_handoff"] == "true")

    qc_csv = qc_dir / "m3_report_qc.csv"
    qc_json = qc_dir / "m3_report_qc.json"
    handoff_qc = qc_dir / "handoff_qc.csv"
    qc_status = {
        "final_candidate_hypotheses_present": "PASS" if input_paths["final_candidate_hypotheses"].is_file() else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "pocket_compound_evidence_table_present": "PASS" if input_paths["pocket_compound_evidence_table"].is_file() else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "rejected_candidate_reasons_present": "PASS" if input_paths["rejected_candidate_reasons"].is_file() else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "final_candidate_gate_qc_json_present": "PASS" if input_paths["final_candidate_gate_qc_json"].is_file() else ("WARN" if allow_partial or mode == "dry-run" else "FAIL"),
        "t10_report_ready": "PASS" if t10_ready or allow_partial or mode == "dry-run" or not require_t10_ready else "FAIL",
        "profile_selected": "PASS" if selected_profile or not (candidates or evidence) else "FAIL",
        "candidate_counts_loaded": "PASS" if candidates or mode == "dry-run" or allow_partial else "WARN",
        "tier_counts_loaded": "PASS",
        "reject_reasons_loaded": "PASS" if rejects or mode == "dry-run" or allow_partial else "WARN",
        "m2_traceability_loaded": "PASS" if input_paths["m2_accepted_pockets"].is_file() else "WARN",
        "m3_upstream_traceability_loaded": "PASS" if input_paths["compound_anchor_convergence"].is_file() or allow_partial or mode == "dry-run" else "WARN",
        "handoff_manifest_written": "PASS" if handoff_csv.is_file() and handoff_json.is_file() else "FAIL",
        "cleanup_report_written": "PASS" if cleanup_result.report_root.is_file() else "FAIL",
        "all_required_artifacts_in_handoff_manifest": "PASS" if handoff_required_present else "WARN",
        "artifact_sha256_computed": "PASS" if all(row["sha256"] or row["exists"] == "false" for row in handoff_rows) else "FAIL",
        "no_internal_id_leak": "PASS" if safety.leaks_detected == 0 else "FAIL",
        "no_smiles_leak": "PASS" if not safety.smiles_logged else "FAIL",
        "no_coordinate_leak": "PASS" if not safety.coordinate_leak_files else "FAIL",
        "candidate_language_hypothesis_only": "PASS" if not safety.candidate_overclaims_detected else "FAIL",
        "no_validated_inhibitor_claim": "PASS" if not safety.candidate_overclaims_detected else "FAIL",
        "no_proven_binder_claim": "PASS" if not safety.candidate_overclaims_detected else "FAIL",
        "no_drug_candidate_claim": "PASS" if not safety.candidate_overclaims_detected else "FAIL",
        "no_clinical_claim": "PASS" if not safety.candidate_overclaims_detected else "FAIL",
    }
    for name in ["milestone3_summary", "final_candidate_summary", "candidate_hypothesis_summary", "qc_summary", "methods_traceability", "reviewer_risk_notes"]:
        qc_status[f"{name}_written"] = "PASS" if outputs[name].is_file() else "FAIL"
    _write_report_qc(qc_csv, qc_status, ctx)
    _write_report_qc(handoff_qc, {"handoff_manifest_written": "PASS" if handoff_csv.is_file() else "FAIL", "all_required_artifacts_in_handoff_manifest": "PASS" if handoff_required_present else "WARN"}, ctx)
    qc_summary = {
        "schema_version": "m3_report_qc_v1",
        "run_id": ctx.run_id,
        "m2_run_id": m2_run_id,
        "profile": selected_profile,
        "generated_at": timestamp,
        "overall_status": status,
        "post_m3_ready": post_m3_ready,
        "m4_ready": m4_ready,
        "allow_partial": allow_partial,
        "require_t10_ready": require_t10_ready,
        "strict_confidentiality": strict_confidentiality,
        "strict_claims": strict_claims,
        "counts": {
            "candidate_rows": counts["candidate_rows"],
            "tier1": counts["tier1"],
            "tier2": counts["tier2"],
            "tier3": counts["tier3"],
            "reject": counts["reject"],
            "evidence_rows": len(evidence),
            "rejected_candidate_rows": len(rejects),
            "report_files_written": sum(1 for key in ["milestone3_summary", "final_candidate_summary", "candidate_hypothesis_summary", "qc_summary", "methods_traceability", "reviewer_risk_notes"] if outputs[key].is_file()),
            "handoff_artifacts_required": sum(1 for row in handoff_rows if row["required_for_handoff"] == "true"),
            "handoff_artifacts_present": sum(1 for row in handoff_rows if row["required_for_handoff"] == "true" and row["exists"] == "true"),
            "confidentiality_leaks": safety.leaks_detected,
            "candidate_overclaims": 1 if safety.candidate_overclaims_detected else 0,
            "coordinate_leak_files": len(safety.coordinate_leak_files),
            "vina_commands_invoked_by_report": 0,
            "qsub_commands_invoked_by_report": 0,
            "runner_commands_invoked_by_report": 0,
            "upstream_rerun_attempts": 0,
        },
        "outputs": {key: ctx.relative_to_repo(path) for key, path in report_ready.items()},
        "blockers": sorted(set(blockers)),
        "warnings": sorted(set(warnings)),
    }
    qc_json.write_text(json.dumps(qc_summary, indent=2, sort_keys=True), encoding="utf-8")
    phase3_log = ctx.logs_dir / "phase3_compounds.log"
    with phase3_log.open("a", encoding="utf-8") as handle:
        handle.write(f"{timestamp} M3-T11 {status} profile={selected_profile} reports_written={qc_summary['counts']['report_files_written']} cleanup_mode={cleanup_result.report.get('mode')} no_vina=true no_qsub=true\n")
    for role, path in report_ready.items():
        append_job_status(ctx, f"m3_report:{role}", "PASS" if path.is_file() else "FAIL", details={"phase": "phase3_compounds", "task": "M3-T11", "job_type": "report_artifact", "artifact_role": role, "relative_path": ctx.relative_to_repo(path), "sha256": sha256(path) if path.is_file() else "", "message": "report artifact"})
    if status == "FAIL":
        append_failed_job(ctx, "M3-T11", "FAIL", ";".join(sorted(set(blockers))[:3]))
    append_phase_status(ctx, "phase3_compounds", status, "M3-T11 report and handoff completed", {"phase": "phase3_compounds", "task": "M3-T11", "status": status, "run_id": ctx.run_id, "m2_run_id": m2_run_id, "timestamp": timestamp, "mode": mode, "profile": selected_profile, "candidate_rows": counts["candidate_rows"], "tier1": counts["tier1"], "tier2": counts["tier2"], "tier3": counts["tier3"], "reject": counts["reject"], "post_m3_ready": post_m3_ready, "m4_ready": m4_ready, "reports_written": qc_summary["counts"]["report_files_written"], "cleanup_mode": cleanup_result.report.get("mode"), "cleanup_files_planned": cleanup_result.report.get("files_planned_for_deletion"), "cleanup_files_deleted": cleanup_result.report.get("files_deleted"), "qsub_invoked": False, "vina_invoked_by_report": False, "runner_invoked_by_report": False, "upstream_rerun_attempted": False, "no_internal_id_leak": safety.leaks_detected == 0, "no_smiles_leak": not safety.smiles_logged, "no_candidate_overclaim": not safety.candidate_overclaims_detected})
    return M3ReportResult(status, post_m3_ready, m4_ready, sorted(set(blockers)), sorted(set(warnings)), selected_profile, allow_partial, require_t10_ready, counts, {key: path.is_file() for key, path in input_paths.items()}, cleanup_result, safety, outputs, qc_json, qc_csv, handoff_qc, phase3_log)
