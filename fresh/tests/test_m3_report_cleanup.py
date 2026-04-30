import csv
import json
from pathlib import Path

from egfr_myo1d.cleanup.milestone_cleanup import plan_m3_cleanup
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.report.milestone3_report import run_m3_report


def make_ctx(tmp_path, run_id="m3_report"):
    repo = tmp_path / "repo"
    fresh = repo / "fresh"
    run_dir = fresh / "runs" / run_id
    ctx = RunContext(
        repo_root=repo,
        fresh_root=fresh,
        run_id=run_id,
        run_dir=run_dir,
        manifest_dir=run_dir / "manifest",
        logs_dir=run_dir / "logs",
        jobs_log_dir=run_dir / "logs" / "jobs",
        errors_dir=run_dir / "logs" / "errors",
        qc_dir=run_dir / "qc",
        reports_dir=run_dir / "reports",
        scratch_dir=run_dir / "scratch",
        tmp_dir=run_dir / "tmp",
    )
    initialize_logs(ctx)
    (repo / ".gitignore").write_text("", encoding="utf-8")
    private = fresh / "data" / "private"
    private.mkdir(parents=True, exist_ok=True)
    (private / "compound_id_map.csv").write_text("public_id,internal_id,notes\nCpd-A,SECRET_A,test\n", encoding="utf-8")
    return ctx


def write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def candidate_row(tier="Tier 1", compound="Cpd-A", pocket="fam_1", reject="none", score="0.900000"):
    code = {
        "Tier 1": "TIER1_STRONG_CHEMICAL_ANCHOR_POCKET_HYPOTHESIS",
        "Tier 2": "TIER2_PLAUSIBLE_SINGLE_ANCHOR_OR_STATE_SPECIFIC_HYPOTHESIS",
        "Tier 3": "TIER3_EXPLORATORY_CHEMICAL_SIGNAL",
        "Reject": "REJECTED_OR_QUARANTINED",
    }[tier]
    return {
        "run_id": "m3_report",
        "m2_run_id": "m2_run",
        "profile": "mini",
        "candidate_hypothesis_id": f"candidate_{compound}_{pocket}",
        "candidate_scope": "compound_pocket",
        "candidate_tier": tier,
        "candidate_tier_code": code,
        "tier_rank": {"Tier 1": 1, "Tier 2": 2, "Tier 3": 3, "Reject": 4}[tier],
        "compound_public_id": compound,
        "pocket_family_id": pocket,
        "candidate_name_public": f"{compound}::{pocket}",
        "candidate_hypothesis_statement": "Computational hypothesis prioritized for review; experimental follow-up is required.",
        "m2_pocket_accepted": "true",
        "non_atp_pass": "true",
        "ppi_relationship_pass": "true",
        "lower_lateral_pass": "true",
        "dimer_accessibility_pass": "true",
        "pose_retention_pass": "true",
        "atp_migration_absent": "true",
        "compound_convergence_pass": "true",
        "primary_state_support": "true",
        "reference_only_support": "false",
        "confidentiality_pass": "true",
        "no_membrane_or_dimer_conflict": "true",
        "all_tier1_hard_gates_pass": "true" if tier == "Tier 1" else "false",
        "all_tier2_hard_gates_pass": "true" if tier in {"Tier 1", "Tier 2"} else "false",
        "candidate_reject_reason": reject,
        "candidate_accept_reason": "hypothesis gates pass" if reject == "none" else "",
        "candidate_priority_score": score,
        "primary_state_ids_supported": "EGFR_160-185",
        "reference_state_ids_supported": "",
        "symmetry_support": "false",
        "protomer_support_summary": "EGFR_160-185:A",
        "compound_public_ids_supported_for_pocket": compound,
        "num_compounds_supported_for_pocket": "2" if tier == "Tier 1" else "1",
        "num_primary_states_supported": "1",
        "best_compound_public_id": compound,
        "best_compound_selection_rule": "deterministic_public_id_order",
        "affinity_used_for_best_compound_selection": "false",
        "affinity_used_for_tier_assignment": "false",
        "affinity_used_for_candidate_promotion": "false",
        "dominant_mechanism_class": "rim_blocker",
        "mechanism_class_counts": '{"rim_blocker": 2}',
        "cluster_ids_supporting": "cluster_1",
        "support_ids_supporting": "support_1",
        "evidence_ids_supporting": "evidence_1",
        "evaluated_at": "2026-04-30T00:00:00+00:00",
    }


def evidence_row(compound="Cpd-A", pocket="fam_1"):
    return {
        "run_id": "m3_report",
        "m2_run_id": "m2_run",
        "profile": "mini",
        "evidence_id": "evidence_1",
        "candidate_hypothesis_id": f"candidate_{compound}_{pocket}",
        "candidate_scope": "compound_pocket",
        "compound_public_id": compound,
        "pocket_family_id": pocket,
        "state_id": "EGFR_160-185",
        "state_role": "primary",
        "primary_state_flag": "true",
        "reference_state_flag": "false",
        "candidate_reject_reason": "none",
        "candidate_hard_gate_status": "PASS",
        "m2_pocket_accepted": "true",
        "compound_trace_present": "true",
        "anchor_trace_present": "true",
        "accepted_pocket_trace_present": "true",
        "evaluated_at": "2026-04-30T00:00:00+00:00",
    }


def write_complete_t10(ctx, m2_run_id="m2_run", tier="Tier 1"):
    phase3 = ctx.run_dir / "phase3_compounds"
    cand = candidate_row(tier=tier, reject="none" if tier != "Reject" else "ATP_migration")
    ev = evidence_row()
    rej = {
        "run_id": ctx.run_id,
        "m2_run_id": m2_run_id,
        "profile": "mini",
        "candidate_hypothesis_id": "candidate_reject",
        "candidate_scope": "compound_pocket",
        "compound_public_id": "Cpd-B",
        "pocket_family_id": "fam_1",
        "state_ids": "EGFR_160-185",
        "candidate_tier": "Reject",
        "candidate_reject_reason": "ATP_migration",
        "hard_gate_failures": "atp_migration_absent",
        "evidence_notes": "quarantined",
        "recommended_fix": "Review upstream gates.",
        "evaluated_at": "2026-04-30T00:00:00+00:00",
    }
    write_csv(phase3 / "tables" / "final_m3_candidate_hypotheses.csv", list(cand), [cand])
    write_csv(phase3 / "tables" / "pocket_compound_evidence_table.csv", list(ev), [ev])
    write_csv(phase3 / "tables" / "rejected_candidate_reasons.csv", list(rej), [rej])
    write_csv(phase3 / "tables" / "compound_anchor_convergence.csv", ["profile", "anchor_id"], [{"profile": "mini", "anchor_id": "anchor_1"}])
    write_csv(phase3 / "tables" / "compound_pose_clusters.csv", ["profile", "cluster_id"], [{"profile": "mini", "cluster_id": "cluster_1"}])
    write_csv(phase3 / "tables" / "compound_pose_attribution.csv", ["profile", "job_id"], [{"profile": "mini", "job_id": "job_1"}])
    write_csv(phase3 / "manifests" / "compound_manifest_public.csv", ["compound_public_id"], [{"compound_public_id": "Cpd-A"}])
    write_csv(phase3 / "manifests" / "vina_job_manifest.csv", ["job_id"], [{"job_id": "job_1"}])
    write_csv(phase3 / "manifests" / "docking_box_manifest.csv", ["box_id"], [{"box_id": "box_1"}])
    qc = {
        "overall_status": "PASS",
        "m3_t11_report_ready": True,
        "profile": "mini",
        "candidate_tiering": {"candidate_rows": 1, "tier1": 1 if tier == "Tier 1" else 0, "tier2": 1 if tier == "Tier 2" else 0, "tier3": 1 if tier == "Tier 3" else 0, "reject": 1 if tier == "Reject" else 0},
        "inputs": {"accepted_pockets_found": True},
        "blockers": [],
        "warnings": [],
    }
    qcd = phase3 / "qc"
    qcd.mkdir(parents=True, exist_ok=True)
    (qcd / "final_candidate_gate_qc.json").write_text(json.dumps(qc), encoding="utf-8")
    write_csv(qcd / "final_candidate_gate_qc.csv", ["check_id", "status"], [{"check_id": "ok", "status": "PASS"}])
    write_csv(qcd / "evidence_integration_qc.csv", ["check_id", "status"], [{"check_id": "ok", "status": "PASS"}])
    write_csv(qcd / "candidate_tiering_qc.csv", ["check_id", "status"], [{"check_id": "ok", "status": "PASS"}])
    m2 = ctx.fresh_root / "runs" / m2_run_id
    write_csv(m2 / "phase2_pockets" / "final" / "accepted_pockets_for_m3.csv", ["pocket_family_id", "state_id"], [{"pocket_family_id": "fam_1", "state_id": "EGFR_160-185"}])
    write_csv(m2 / "phase2_pockets" / "gated" / "pocket_gate_qc.csv", ["check_id", "status"], [{"check_id": "ok", "status": "PASS"}])
    write_csv(m2 / "phase1_ppi" / "consensus" / "ppi_consensus_patch_merged.csv", ["hotspot_id"], [{"hotspot_id": "h1"}])


def test_missing_inputs_dry_run_warns_and_report_fails_without_partial(tmp_path):
    ctx = make_ctx(tmp_path)
    dry = run_m3_report(ctx, m2_run_id="m2_run", profile="mini", mode="dry-run")
    fail = run_m3_report(ctx, m2_run_id="m2_run", profile="mini", mode="report")

    assert dry.status in {"WARN", "FAIL"}
    assert dry.qc_json.is_file()
    assert fail.status == "FAIL"
    assert any("final_m3_candidate" in item for item in fail.blockers)


def test_complete_t10_outputs_generate_all_reports_and_handoff(tmp_path):
    ctx = make_ctx(tmp_path)
    write_complete_t10(ctx)
    result = run_m3_report(ctx, m2_run_id="m2_run", profile="mini", cleanup_mode="execute", cleanup_scope="test", confirm_run_id=ctx.run_id)

    assert result.status == "PASS"
    assert result.m4_ready is True
    for name in ["milestone3_summary", "final_candidate_summary", "candidate_hypothesis_summary", "qc_summary", "methods_traceability", "reviewer_risk_notes", "handoff_manifest_csv", "handoff_manifest_json"]:
        assert result.outputs[name].is_file()
    summary = result.outputs["milestone3_summary"].read_text(encoding="utf-8")
    assert "# Milestone 3 Summary" in summary
    assert "## 15. Tier 1 hypotheses" in summary
    assert "Affinity is descriptive" in summary
    assert "ATP-confounded evidence is not promoted" in summary
    handoff = json.loads(result.outputs["handoff_manifest_json"].read_text(encoding="utf-8"))
    assert handoff["required_artifacts_present"] is True


def test_t10_not_ready_fails_unless_partial_or_no_require(tmp_path):
    ctx = make_ctx(tmp_path)
    write_complete_t10(ctx)
    qc_path = ctx.run_dir / "phase3_compounds" / "qc" / "final_candidate_gate_qc.json"
    data = json.loads(qc_path.read_text(encoding="utf-8"))
    data["m3_t11_report_ready"] = False
    qc_path.write_text(json.dumps(data), encoding="utf-8")

    assert run_m3_report(ctx, m2_run_id="m2_run", profile="mini").status == "FAIL"
    assert run_m3_report(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True).status in {"WARN", "PASS"}
    assert run_m3_report(ctx, m2_run_id="m2_run", profile="mini", require_t10_ready=False).status in {"WARN", "PASS"}


def test_safety_scan_detects_private_id_claim_and_coordinates(tmp_path):
    ctx = make_ctx(tmp_path)
    write_complete_t10(ctx)
    report_dir = ctx.run_dir / "phase3_compounds" / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "bad.md").write_text("SECRET_A\nvalidated inhibitor\nATOM      1  C   LIG A   1       1.000   2.000   3.000\n", encoding="utf-8")
    result = run_m3_report(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True)

    assert result.status == "FAIL"
    assert result.safety.leaks_detected >= 1
    assert result.safety.candidate_overclaims_detected is True
    assert result.safety.coordinate_leak_files


def test_cleanup_plan_preserves_essential_outputs_and_finds_tmp(tmp_path):
    ctx = make_ctx(tmp_path)
    write_complete_t10(ctx)
    tmp_file = ctx.run_dir / "phase3_compounds" / "scratch" / "x.tmp"
    tmp_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_file.write_text("temporary", encoding="utf-8")

    result = plan_m3_cleanup(ctx, profile="mini", mode="plan", cleanup_scope="test")
    planned = {item["relative_path"] for item in result.report["planned_deletions"]}
    assert "phase3_compounds/scratch/x.tmp" in planned
    assert not any(path.startswith("phase3_compounds/tables/") for path in planned)
    assert tmp_file.exists()


def test_cleanup_execute_requires_confirm_and_deletes_only_planned_tmp(tmp_path):
    ctx = make_ctx(tmp_path)
    write_complete_t10(ctx)
    tmp_file = ctx.run_dir / "phase3_compounds" / "tmp" / "delete.tmp"
    tmp_file.parent.mkdir(parents=True, exist_ok=True)
    tmp_file.write_text("temporary", encoding="utf-8")

    fail = plan_m3_cleanup(ctx, profile="mini", mode="execute", cleanup_scope="test")
    assert fail.status == "FAIL"
    assert tmp_file.exists()
    ok = plan_m3_cleanup(ctx, profile="mini", mode="execute", cleanup_scope="test", confirm_run_id=ctx.run_id)
    assert ok.status == "PASS"
    assert not tmp_file.exists()
    assert ok.report["files_deleted"] == 1


def test_cleanup_blocks_symlink_escape(tmp_path):
    ctx = make_ctx(tmp_path)
    write_complete_t10(ctx)
    outside = tmp_path / "outside.txt"
    outside.write_text("nope", encoding="utf-8")
    link = ctx.run_dir / "phase3_compounds" / "tmp" / "escape.tmp"
    link.parent.mkdir(parents=True, exist_ok=True)
    try:
        link.symlink_to(outside)
    except (OSError, NotImplementedError):
        return
    result = plan_m3_cleanup(ctx, profile="mini", mode="plan", cleanup_scope="test")
    assert result.status == "FAIL"
    assert outside.exists()


def test_report_outputs_are_deterministic_in_plan_mode(tmp_path):
    ctx = make_ctx(tmp_path)
    write_complete_t10(ctx, tier="Tier 3")
    first = run_m3_report(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True)
    text1 = first.outputs["milestone3_summary"].read_text(encoding="utf-8")
    second = run_m3_report(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True)
    text2 = second.outputs["milestone3_summary"].read_text(encoding="utf-8")
    assert text1 == text2
    assert second.counts["tier3"] == 1
