import csv
import json
from pathlib import Path

from egfr_myo1d.compound.anchor_convergence import ANCHOR_FIELDS, SUPPORT_FIELDS
from egfr_myo1d.compound.evidence_integration import run_m3_evidence_tiering
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext


def make_ctx(tmp_path, run_id="m3_t10"):
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


def write_m2(ctx, m2_run_id="m2_run", states=("EGFR_160-185",)):
    m2 = ctx.fresh_root / "runs" / m2_run_id
    pockets = [
        {
            "pocket_family_id": "fam_1",
            "state_id": state,
            "non_atp_pass": "true",
            "lower_lateral_pass": "true",
            "dimer_accessibility_pass": "true",
            "ppi_relationship_pass": "true",
        }
        for state in states
    ]
    write_csv(m2 / "phase2_pockets" / "final" / "accepted_pockets_for_m3.csv", list(pockets[0]), pockets)
    write_csv(m2 / "phase2_pockets" / "gated" / "pocket_gate_qc.csv", ["check_id", "status"], [{"check_id": "ok", "status": "PASS"}])
    write_csv(m2 / "phase1_ppi" / "consensus" / "ppi_consensus_patch_merged.csv", ["hotspot_id"], [{"hotspot_id": "ppi_1"}])


def support_row(cid="Cpd-A", state="EGFR_160-185", support_id="support_1", allowed="true", reject="none", atp="0.000000", mem="0.000000", dimer="0.000000"):
    row = {field: "" for field in SUPPORT_FIELDS}
    row.update(
        {
            "run_id": "m3_t10",
            "m2_run_id": "m2_run",
            "profile": "mini",
            "support_id": support_id,
            "compound_public_id": cid,
            "pocket_family_id": "fam_1",
            "state_id": state,
            "state_role": "reference" if state == "3GT8_raw" else "primary",
            "primary_state_flag": "false" if state == "3GT8_raw" else "true",
            "reference_state_flag": "true" if state == "3GT8_raw" else "false",
            "protomer_ids_with_support": "A",
            "protomer_symmetry_support": "false",
            "box_ids_with_support": "box_1",
            "cluster_ids_supporting": f"cluster_{support_id}",
            "num_anchor_allowed_clusters": "1" if allowed == "true" else "0",
            "num_valid_pose_clusters": "1",
            "max_cluster_size": "2",
            "max_cluster_fraction": "0.750000",
            "median_cluster_fraction": "0.750000",
            "within_pocket_fraction_min": "1.000000",
            "atp_migration_fraction_max": atp,
            "ppi_relationship_confirmed_fraction_max": "1.000000",
            "ppi_relationship_confirmed_fraction_median": "1.000000",
            "ppi_hotspot_contact_median_max": "1.000000",
            "ppi_hotspot_contact_mean_max": "1.000000",
            "nearest_ppi_hotspot_distance_median_min": "3.000000",
            "ppi_rim_distance_median_min": "3.000000",
            "membrane_penetration_fraction_max": mem,
            "dimer_interface_clash_fraction_max": dimer,
            "mechanism_class_counts": json.dumps({"rim_blocker": 2}),
            "mechanism_class_majority": "rim_blocker",
            "vina_affinity_median_kcal_mol": "-7.000000",
            "vina_affinity_mean_kcal_mol": "-7.000000",
            "vina_affinity_best_kcal_mol": "-7.000000",
            "vina_affinity_worst_kcal_mol": "-7.000000",
            "affinity_used_for_ranking": "false",
            "support_status": "PASS" if allowed == "true" else "WARN",
            "support_reject_reason": reject,
            "allowed_for_anchor_convergence": allowed,
            "allowed_for_evidence_integration": allowed,
        }
    )
    return row


def anchor_row(anchor_class="multi_compound_primary_state_supported", status="PASS", allowed="true", compounds="Cpd-A;Cpd-B", support_ids="support_1;support_2", primary="1"):
    row = {field: "" for field in ANCHOR_FIELDS}
    row.update(
        {
            "run_id": "m3_t10",
            "m2_run_id": "m2_run",
            "profile": "mini",
            "anchor_id": "anchor_1",
            "anchor_scope": "primary_cross_state_pocket",
            "pocket_family_id": "fam_1",
            "state_ids_supported": "EGFR_160-185",
            "primary_state_support": primary,
            "primary_state_ids_supported": "EGFR_160-185" if primary != "0" else "",
            "reference_state_support": "0" if primary != "0" else "1",
            "reference_state_ids_supported": "" if primary != "0" else "3GT8_raw",
            "symmetry_support": "false",
            "num_compounds_tested": "3",
            "num_compounds_with_valid_pose": "2",
            "num_compounds_with_converged_pose": "2",
            "num_compounds_allowed_for_anchor_convergence": str(len([c for c in compounds.split(";") if c])),
            "compound_public_ids_tested": "Cpd-A;Cpd-B;Cpd-C",
            "compound_public_ids_supported": compounds,
            "best_compound_public_id": "Cpd-A",
            "best_compound_selection_rule": "deterministic_public_id_order",
            "affinity_used_for_best_compound_selection": "false",
            "median_affinity_across_supported": "-7.000000",
            "mean_affinity_across_supported": "-7.000000",
            "best_affinity_across_supported": "-7.000000",
            "worst_affinity_across_supported": "-7.000000",
            "affinity_used_for_ranking": "false",
            "cluster_ids_supporting": "cluster_support_1;cluster_support_2",
            "support_ids_supporting": support_ids,
            "mechanism_class_counts": json.dumps({"rim_blocker": 2}),
            "dominant_mechanism_class": "rim_blocker",
            "max_cluster_fraction_across_supported": "0.750000",
            "median_cluster_fraction_across_supported": "0.750000",
            "min_within_pocket_fraction_across_supported": "1.000000",
            "max_atp_migration_fraction_across_supported": "0.000000",
            "max_membrane_penetration_fraction_across_supported": "0.000000",
            "max_dimer_interface_clash_fraction_across_supported": "0.000000",
            "max_ppi_relationship_confirmed_fraction": "1.000000",
            "anchor_converged": allowed,
            "anchor_convergence_class": anchor_class,
            "anchor_convergence_status": status,
            "anchor_reject_reason": "none" if allowed == "true" else "weak_or_scattered_support",
            "allowed_for_evidence_integration": allowed,
        }
    )
    return row


def write_m3(ctx, supports, anchors, ready=True):
    tables = ctx.run_dir / "phase3_compounds" / "tables"
    qc = ctx.run_dir / "phase3_compounds" / "qc"
    write_csv(tables / "compound_pocket_support.csv", SUPPORT_FIELDS, supports)
    write_csv(tables / "compound_anchor_convergence.csv", ANCHOR_FIELDS, anchors)
    write_csv(tables / "compound_pose_clusters.csv", ["profile", "cluster_id"], [{"profile": "mini", "cluster_id": "cluster_1"}])
    qc.mkdir(parents=True, exist_ok=True)
    (qc / "anchor_convergence_qc.json").write_text(json.dumps({"overall_status": "PASS", "m3_t10_evidence_integration_ready": ready}), encoding="utf-8")
    write_csv(ctx.run_dir / "phase3_compounds" / "manifests" / "compound_manifest_public.csv", ["compound_public_id"], [{"compound_public_id": "Cpd-A"}, {"compound_public_id": "Cpd-B"}, {"compound_public_id": "Cpd-C"}])


def test_dry_run_missing_inputs_warns_and_integrate_missing_fails(tmp_path):
    ctx = make_ctx(tmp_path)
    dry = run_m3_evidence_tiering(ctx, m2_run_id="m2_run", profile="mini", mode="dry-run", write_empty_evidence=True)
    fail = run_m3_evidence_tiering(ctx, m2_run_id="m2_run", profile="mini", mode="integrate")

    assert dry.status == "WARN"
    assert dry.evidence_csv.is_file()
    assert dry.qc_json.is_file()
    assert fail.status == "FAIL"
    assert any("accepted_pockets" in item for item in fail.blockers)


def test_synthetic_multi_compound_primary_support_becomes_tier1_without_affinity_promotion(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m2(ctx)
    write_m3(ctx, [support_row("Cpd-A", support_id="support_1"), support_row("Cpd-B", support_id="support_2")], [anchor_row()])

    result = run_m3_evidence_tiering(ctx, m2_run_id="m2_run", profile="mini")
    candidates = read_csv(result.candidates_csv)
    evidence = read_csv(result.evidence_csv)

    assert result.status == "PASS"
    assert result.m3_t11_report_ready is True
    assert any(row["candidate_tier"] == "Tier 1" for row in candidates)
    assert all(row["affinity_used_for_candidate_promotion"] == "false" for row in candidates)
    assert all(row["affinity_used_for_best_compound_selection"] == "false" for row in candidates)
    assert all("validated" not in row["candidate_hypothesis_statement"].lower() for row in candidates)
    assert len(evidence) >= 2


def test_prior_smiles_diagnostics_and_neutral_candidate_terms_do_not_block(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m2(ctx)
    write_m3(ctx, [support_row("Cpd-A", support_id="support_1"), support_row("Cpd-B", support_id="support_2")], [anchor_row()])
    (ctx.errors_dir / "error_summary.txt").write_text("previous blocker: SMILES printed in public outputs/logs\n", encoding="utf-8")
    report = ctx.run_dir / "phase3_compounds" / "reports" / "old_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("final candidate hypotheses table was generated for computational review only\n", encoding="utf-8")

    result = run_m3_evidence_tiering(ctx, m2_run_id="m2_run", profile="mini")

    assert result.status == "PASS"
    assert not any("SMILES printed" in item for item in result.blockers)
    assert not any("candidate overclaim" in item for item in result.blockers)


def test_missing_anchor_qc_and_not_ready_fail_unless_partial(tmp_path):
    ctx = make_ctx(tmp_path / "missing_qc")
    write_m2(ctx)
    write_m3(ctx, [support_row()], [anchor_row()])
    (ctx.run_dir / "phase3_compounds" / "qc" / "anchor_convergence_qc.json").unlink()
    assert run_m3_evidence_tiering(ctx, m2_run_id="m2_run", profile="mini").status == "FAIL"
    assert run_m3_evidence_tiering(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True).status in {"PASS", "WARN"}

    ctx2 = make_ctx(tmp_path / "not_ready")
    write_m2(ctx2)
    write_m3(ctx2, [support_row()], [anchor_row()], ready=False)
    assert run_m3_evidence_tiering(ctx2, m2_run_id="m2_run", profile="mini").status == "FAIL"
    assert run_m3_evidence_tiering(ctx2, m2_run_id="m2_run", profile="mini", allow_partial=True).status in {"PASS", "WARN"}


def test_nonpublic_and_affinity_flags_fail(tmp_path):
    ctx = make_ctx(tmp_path / "badid")
    write_m2(ctx)
    write_m3(ctx, [support_row("SECRET")], [anchor_row(compounds="SECRET")])
    assert run_m3_evidence_tiering(ctx, m2_run_id="m2_run", profile="mini").status == "FAIL"

    ctx2 = make_ctx(tmp_path / "aff")
    write_m2(ctx2)
    bad_support = support_row()
    bad_support["affinity_used_for_ranking"] = "true"
    bad_anchor = anchor_row()
    bad_anchor["affinity_used_for_best_compound_selection"] = "true"
    write_m3(ctx2, [bad_support], [bad_anchor])
    result = run_m3_evidence_tiering(ctx2, m2_run_id="m2_run", profile="mini")
    assert result.status == "FAIL"
    assert any("affinity" in item for item in result.blockers)


def test_atp_membrane_reference_and_missing_trace_are_not_promoted(tmp_path):
    ctx = make_ctx(tmp_path / "atp")
    write_m2(ctx)
    write_m3(ctx, [support_row(atp="0.500000", reject="ATP_migration")], [anchor_row()])
    candidates = read_csv(run_m3_evidence_tiering(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True).candidates_csv)
    assert all(row["candidate_tier"] not in {"Tier 1", "Tier 2"} for row in candidates)

    ctx2 = make_ctx(tmp_path / "mem")
    write_m2(ctx2)
    write_m3(ctx2, [support_row(mem="0.500000")], [anchor_row()])
    candidates2 = read_csv(run_m3_evidence_tiering(ctx2, m2_run_id="m2_run", profile="mini", allow_partial=True).candidates_csv)
    assert all(row["candidate_tier"] not in {"Tier 1", "Tier 2"} for row in candidates2)

    ctx3 = make_ctx(tmp_path / "ref")
    write_m2(ctx3, states=("3GT8_raw",))
    write_m3(ctx3, [support_row("Cpd-A", state="3GT8_raw")], [anchor_row(anchor_class="reference_only", status="WARN", compounds="Cpd-A", support_ids="support_1", primary="0")])
    candidates3 = read_csv(run_m3_evidence_tiering(ctx3, m2_run_id="m2_run", profile="mini", allow_partial=True).candidates_csv)
    assert all(row["candidate_tier"] not in {"Tier 1", "Tier 2"} for row in candidates3)

    ctx4 = make_ctx(tmp_path / "missing_trace")
    write_m3(ctx4, [support_row()], [anchor_row()])
    candidates4 = read_csv(run_m3_evidence_tiering(ctx4, m2_run_id="m2_run", profile="mini", allow_partial=True).candidates_csv)
    assert all(row["candidate_tier"] == "Reject" for row in candidates4)


def test_tier2_patterns_and_weak_support(tmp_path):
    ctx = make_ctx(tmp_path / "robust")
    write_m2(ctx, states=("EGFR_160-185", "EGFR_170-200"))
    write_m3(ctx, [support_row("Cpd-A", "EGFR_160-185", "support_1"), support_row("Cpd-A", "EGFR_170-200", "support_2")], [anchor_row(anchor_class="single_compound_state_robust", compounds="Cpd-A")])
    candidates = read_csv(run_m3_evidence_tiering(ctx, m2_run_id="m2_run", profile="mini").candidates_csv)
    assert any(row["candidate_tier"] == "Tier 2" for row in candidates)

    ctx2 = make_ctx(tmp_path / "weak")
    write_m2(ctx2)
    write_m3(ctx2, [support_row(allowed="false", reject="cluster_not_allowed_for_anchor_convergence")], [anchor_row(anchor_class="weak_or_scattered", status="WARN", allowed="false", compounds="Cpd-A", support_ids="support_1")])
    candidates2 = read_csv(run_m3_evidence_tiering(ctx2, m2_run_id="m2_run", profile="mini", allow_partial=True).candidates_csv)
    assert all(row["candidate_tier"] in {"Tier 3", "Reject"} for row in candidates2)


def test_broad_scan_ignored_by_default_and_never_tier1_when_allowed(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m2(ctx)
    write_m3(ctx, [], [], ready=True)
    broad = ctx.run_dir / "phase3_compounds" / "tables" / "broad_anchor_scan_support.csv"
    write_csv(broad, ["compound_public_id", "pocket_family_id", "state_id"], [{"compound_public_id": "Cpd-A", "pocket_family_id": "fam_1", "state_id": "EGFR_160-185"}])
    result = run_m3_evidence_tiering(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True, broad_anchor_scan_support=broad)
    assert result.inputs["broad_anchor_scan_support_found"] is False
    result2 = run_m3_evidence_tiering(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True, allow_broad_scan_support=True, broad_anchor_scan_support=broad)
    candidates = read_csv(result2.candidates_csv)
    assert any(row["candidate_scope"] == "exploratory_broad_scan" for row in candidates)
    assert all(row["candidate_tier"] != "Tier 1" for row in candidates)


def test_public_outputs_have_no_private_ids_coordinates_or_overclaims_and_logs_written(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m2(ctx)
    write_m3(ctx, [support_row()], [anchor_row(compounds="Cpd-A", support_ids="support_1")])
    result = run_m3_evidence_tiering(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True)

    public_text = "\n".join(
        path.read_text(encoding="utf-8", errors="ignore")
        for root in [ctx.run_dir / "phase3_compounds" / "tables", ctx.run_dir / "phase3_compounds" / "qc", ctx.run_dir / "phase3_compounds" / "reports", ctx.logs_dir]
        if root.exists()
        for path in root.rglob("*")
        if path.is_file()
    )
    assert "SECRET_A" not in public_text
    assert "ATOM      1" not in public_text
    assert "validated inhibitor" not in public_text.lower()
    assert result.qc_json.is_file()
    assert (ctx.logs_dir / "phase_status.jsonl").read_text(encoding="utf-8")
    assert (ctx.logs_dir / "job_status.jsonl").read_text(encoding="utf-8")
    assert (ctx.run_dir / "phase3_compounds" / "config_snapshots" / "compound_evidence_tiering.resolved.yaml").is_file()
