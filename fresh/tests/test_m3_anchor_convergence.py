import csv
import hashlib
import json
from pathlib import Path

from egfr_myo1d.compound.anchor_convergence import run_m3_anchor_convergence
from egfr_myo1d.compound.pose_clustering import CLUSTER_FIELDS, MEMBERSHIP_FIELDS
from egfr_myo1d.compound.pose_attribution import ATTRIBUTION_FIELDS
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext


def make_ctx(tmp_path, run_id="m3_anchor"):
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


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_refs(ctx, states=("EGFR_160-185", "EGFR_170-200"), m2_run_id="m2_run"):
    m2 = ctx.fresh_root / "runs" / m2_run_id
    pockets = []
    boxes = []
    for state in states:
        pockets.append({"pocket_family_id": "fam_1", "state_id": state, "non_atp_pass": "true", "ppi_relationship_class": "PPI-adjacent"})
        boxes.append({"pocket_family_id": "fam_1", "state_id": state, "protomer_id": "A", "box_id": "box_1"})
    write_csv(m2 / "phase2_pockets" / "final" / "accepted_pockets_for_m3.csv", ["pocket_family_id", "state_id", "non_atp_pass", "ppi_relationship_class"], pockets)
    write_csv(m2 / "phase2_pockets" / "final" / "accepted_pocket_boxes.csv", ["pocket_family_id", "state_id", "protomer_id", "box_id"], boxes)


def cluster_row(cid="Cpd-A", state="EGFR_160-185", cluster_id="cluster_1", protomer="A", allowed="true", converged="true", status="PASS_CONVERGED", reject="none", affinity="-7.0", reference=False, atp="0.000000", mem="0.000000", dimer="0.000000", fraction="0.750000"):
    row = {field: "" for field in CLUSTER_FIELDS}
    row.update(
        {
            "run_id": "m3_anchor",
            "m2_run_id": "m2_run",
            "profile": "mini",
            "cluster_id": cluster_id,
            "cluster_algorithm": "graph_connected_components",
            "cluster_method": "hybrid",
            "compound_public_id": cid,
            "state_id": state,
            "state_role": "reference" if reference or state == "3GT8_raw" else "primary",
            "pocket_family_id": "fam_1",
            "protomer_id": protomer,
            "box_ids_in_cluster": "box_1",
            "eligible_pose_count_in_group": "2",
            "cluster_size": "2",
            "cluster_fraction": fraction,
            "vina_affinity_median_kcal_mol": affinity,
            "vina_affinity_mean_kcal_mol": affinity,
            "vina_affinity_best_kcal_mol": affinity,
            "vina_affinity_worst_kcal_mol": affinity,
            "affinity_used_for_ranking": "false",
            "within_pocket_fraction": "1.000000",
            "atp_migration_fraction": atp,
            "ppi_relationship_confirmed_fraction": "1.000000",
            "ppi_hotspot_contact_median": "1",
            "ppi_hotspot_contact_mean": "1",
            "nearest_ppi_hotspot_distance_median": "3.0",
            "ppi_rim_distance_median": "3.0",
            "membrane_penetration_fraction": mem,
            "dimer_interface_clash_fraction": dimer,
            "mechanism_class_majority": "rim_blocker",
            "mechanism_class_counts": json.dumps({"rim_blocker": 2}),
            "converged_pose_cluster": converged,
            "convergence_status": status,
            "cluster_reject_reason": reject,
            "allowed_for_anchor_convergence": allowed,
        }
    )
    return row


def write_inputs(ctx, clusters, ready=True, manifest=True, qc_status="PASS"):
    membership = ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_cluster_membership.csv"
    write_csv(membership, MEMBERSHIP_FIELDS, [])
    attr = ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_attribution.csv"
    write_csv(attr, ATTRIBUTION_FIELDS, [])
    attr_sha = sha256(attr)
    prepared_clusters = []
    for row in clusters:
        prepared = dict(row)
        if not prepared.get("source_compound_pose_attribution_sha256"):
            prepared["source_compound_pose_attribution_sha256"] = attr_sha
        prepared_clusters.append(prepared)
    cluster_path = ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_clusters.csv"
    write_csv(cluster_path, CLUSTER_FIELDS, prepared_clusters)
    qc = ctx.run_dir / "phase3_compounds" / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    (qc / "pose_clustering_qc.json").write_text(json.dumps({"overall_status": qc_status, "m3_t9_anchor_convergence_ready": ready}), encoding="utf-8")
    if manifest:
        write_csv(ctx.run_dir / "phase3_compounds" / "manifests" / "compound_manifest_public.csv", ["compound_public_id"], [{"compound_public_id": "Cpd-A"}, {"compound_public_id": "Cpd-B"}, {"compound_public_id": "Cpd-C"}])
    return cluster_path


def test_dry_run_missing_clusters_warns_and_converge_missing_fails(tmp_path):
    ctx = make_ctx(tmp_path)
    dry = run_m3_anchor_convergence(ctx, m2_run_id="m2_run", profile="mini", mode="dry-run", write_empty_convergence=True)
    fail = run_m3_anchor_convergence(ctx, m2_run_id="m2_run", profile="mini", mode="converge")

    assert dry.status == "WARN"
    assert dry.support_csv.is_file()
    assert fail.status == "FAIL"
    assert any("compound_pose_clusters" in item for item in fail.blockers)


def test_two_compounds_same_primary_pocket_passes_without_affinity_ranking(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx, states=("EGFR_160-185",))
    write_inputs(ctx, [cluster_row("Cpd-A", cluster_id="c1", affinity="-1.0"), cluster_row("Cpd-B", cluster_id="c2", affinity="-12.0")])

    result = run_m3_anchor_convergence(ctx, m2_run_id="m2_run", profile="mini")
    anchors = read_csv(result.anchor_csv)
    support = read_csv(result.support_csv)

    assert result.status == "PASS"
    assert result.m3_t10_evidence_integration_ready is True
    assert any(row["anchor_convergence_class"] in {"multi_compound_primary_state_supported", "single_state_multi_compound"} for row in anchors)
    assert all(row["affinity_used_for_ranking"] == "false" for row in anchors)
    assert all(row["affinity_used_for_best_compound_selection"] == "false" for row in anchors)
    assert any(row["best_compound_public_id"] == "Cpd-A" for row in anchors if row["allowed_for_evidence_integration"] == "true")
    assert len(support) == 2


def test_schema_duplicate_nonpublic_and_affinity_ranking_fail(tmp_path):
    ctx = make_ctx(tmp_path / "schema")
    bad = ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_clusters.csv"
    write_csv(bad, ["cluster_id"], [{"cluster_id": "x"}])
    assert run_m3_anchor_convergence(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True).status == "FAIL"

    ctx2 = make_ctx(tmp_path / "dupe")
    write_refs(ctx2)
    row = cluster_row(cluster_id="same")
    write_inputs(ctx2, [row, dict(row)])
    assert run_m3_anchor_convergence(ctx2, m2_run_id="m2_run", profile="mini").status == "FAIL"

    ctx3 = make_ctx(tmp_path / "badid")
    write_refs(ctx3)
    write_inputs(ctx3, [cluster_row(cid="SECRET")])
    assert run_m3_anchor_convergence(ctx3, m2_run_id="m2_run", profile="mini").status == "FAIL"

    ctx4 = make_ctx(tmp_path / "aff")
    write_refs(ctx4)
    bad_aff = cluster_row()
    bad_aff["affinity_used_for_ranking"] = "true"
    write_inputs(ctx4, [bad_aff])
    assert run_m3_anchor_convergence(ctx4, m2_run_id="m2_run", profile="mini").status == "FAIL"


def test_qc_gate_and_missing_accepted_pockets(tmp_path):
    ctx = make_ctx(tmp_path / "qc")
    write_refs(ctx)
    write_inputs(ctx, [cluster_row()], ready=False)
    assert run_m3_anchor_convergence(ctx, m2_run_id="m2_run", profile="mini").status == "FAIL"
    assert run_m3_anchor_convergence(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True).status in {"PASS", "WARN"}

    ctx2 = make_ctx(tmp_path / "nopocket")
    write_inputs(ctx2, [cluster_row()])
    assert run_m3_anchor_convergence(ctx2, m2_run_id="m2_run", profile="mini").status == "FAIL"


def test_t8_qc_must_pass_and_cluster_source_hash_must_match(tmp_path):
    ctx = make_ctx(tmp_path / "qc_status")
    write_refs(ctx)
    write_inputs(ctx, [cluster_row()], ready=True, qc_status="WARN")
    blocked = run_m3_anchor_convergence(ctx, m2_run_id="m2_run", profile="mini")
    assert blocked.status == "FAIL"
    assert any("pose_clustering_qc.json has m3_t9_anchor_convergence_ready=false" in item for item in blocked.blockers)

    ctx2 = make_ctx(tmp_path / "hash")
    write_refs(ctx2)
    write_inputs(ctx2, [cluster_row()])
    attr = ctx2.run_dir / "phase3_compounds" / "tables" / "compound_pose_attribution.csv"
    attr.write_text(attr.read_text(encoding="utf-8") + "# changed after clustering\n", encoding="utf-8")

    mismatch = run_m3_anchor_convergence(ctx2, m2_run_id="m2_run", profile="mini")
    assert mismatch.status == "FAIL"
    assert any("cluster source attribution hash mismatch" in item for item in mismatch.blockers)


def test_support_filters_and_rejected_clusters_are_represented(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx)
    rows = [
        cluster_row(cluster_id="ok"),
        cluster_row(cluster_id="not_allowed", allowed="false"),
        cluster_row(cluster_id="atp", atp="0.500000"),
        cluster_row(cluster_id="mem", mem="0.500000"),
        cluster_row(cluster_id="dimer", dimer="0.500000"),
    ]
    write_inputs(ctx, rows)
    result = run_m3_anchor_convergence(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True)
    support = read_csv(result.support_csv)

    assert len(support) == 1
    assert support[0]["num_clusters_total"] == "5"
    assert support[0]["num_anchor_allowed_clusters"] == "1"
    assert support[0]["allowed_for_anchor_convergence"] == "true"


def test_state_symmetry_reference_and_single_compound_robust_classes(tmp_path):
    ctx = make_ctx(tmp_path / "robust")
    write_refs(ctx, states=("EGFR_160-185", "EGFR_170-200"))
    write_inputs(ctx, [cluster_row("Cpd-A", "EGFR_160-185", "c1"), cluster_row("Cpd-A", "EGFR_170-200", "c2")])
    anchors = read_csv(run_m3_anchor_convergence(ctx, m2_run_id="m2_run", profile="mini").anchor_csv)
    assert any(row["anchor_convergence_class"] == "single_compound_state_robust" for row in anchors)

    ctx2 = make_ctx(tmp_path / "sym")
    write_refs(ctx2, states=("EGFR_160-185",))
    write_inputs(ctx2, [cluster_row("Cpd-A", "EGFR_160-185", "a", protomer="A"), cluster_row("Cpd-A", "EGFR_160-185", "b", protomer="B")])
    anchor = [row for row in read_csv(run_m3_anchor_convergence(ctx2, m2_run_id="m2_run", profile="mini", allow_partial=True).anchor_csv) if row["anchor_scope"] == "primary_cross_state_pocket"][0]
    assert anchor["primary_state_support"] == "1"
    assert anchor["symmetry_support"] == "true"

    ctx3 = make_ctx(tmp_path / "ref")
    write_refs(ctx3, states=("3GT8_raw",))
    write_inputs(ctx3, [cluster_row("Cpd-A", "3GT8_raw", "ref", reference=True)])
    result = run_m3_anchor_convergence(ctx3, m2_run_id="m2_run", profile="mini", allow_partial=True)
    anchors_ref = read_csv(result.anchor_csv)
    assert any(row["anchor_convergence_class"] == "reference_only" for row in anchors_ref)
    assert all(row["allowed_for_evidence_integration"] == "false" for row in anchors_ref if row["anchor_convergence_class"] == "reference_only")


def test_no_support_yields_warn_and_no_final_tables_or_leaks(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx)
    write_inputs(ctx, [cluster_row(cluster_id="scattered", allowed="false", converged="false", status="WARN_SINGLETON_NOT_CONVERGED", reject="singleton_cluster")])
    result = run_m3_anchor_convergence(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True)

    assert result.status == "WARN"
    assert result.m3_t10_evidence_integration_ready is False
    assert not (ctx.run_dir / "phase3_compounds" / "tables" / "pocket_compound_evidence_table.csv").exists()
    assert not (ctx.run_dir / "phase3_compounds" / "tables" / "final_m3_candidate_hypotheses.csv").exists()
    public_text = "\n".join(
        path.read_text(encoding="utf-8", errors="ignore")
        for root in [ctx.run_dir / "phase3_compounds" / "tables", ctx.run_dir / "phase3_compounds" / "qc", ctx.run_dir / "phase3_compounds" / "reports", ctx.logs_dir]
        if root.exists()
        for path in root.rglob("*")
        if path.is_file()
    )
    assert "SECRET_A" not in public_text
    assert "ATOM      1" not in public_text
    assert (ctx.run_dir / "phase3_compounds" / "config_snapshots" / "compound_anchor_convergence.resolved.yaml").is_file()
