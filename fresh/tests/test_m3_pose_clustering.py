import csv
import hashlib
import json
from pathlib import Path

from egfr_myo1d.compound.pose_clustering import run_m3_pose_clustering
from egfr_myo1d.compound.pose_attribution import ATTRIBUTION_FIELDS
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext


POSE_MULTI = """MODEL 1
REMARK VINA RESULT: -7.500 0.000 0.000
ATOM      1  C   LIG A   1       0.000   0.000   8.000  0.00  0.00     0.000 C
ATOM      2  N   LIG A   1       1.000   0.000   8.000  0.00  0.00     0.000 N
ENDMDL
MODEL 2
REMARK VINA RESULT: -7.200 0.000 0.000
ATOM      1  C   LIG A   1       0.200   0.100   8.000  0.00  0.00     0.000 C
ATOM      2  N   LIG A   1       1.200   0.100   8.000  0.00  0.00     0.000 N
ENDMDL
MODEL 3
REMARK VINA RESULT: -6.100 0.000 0.000
ATOM      1  C   LIG A   1      12.000  12.000   8.000  0.00  0.00     0.000 C
ATOM      2  N   LIG A   1      13.000  12.000   8.000  0.00  0.00     0.000 N
ENDMDL
"""


def make_ctx(tmp_path, run_id="m3_pose_cluster"):
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
    (repo / ".gitignore").write_text("fresh/runs/*/phase3_compounds/docking_outputs/**\n", encoding="utf-8")
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


def write_refs(ctx, m2_run_id="m2_run"):
    m2 = ctx.fresh_root / "runs" / m2_run_id
    write_csv(
        m2 / "phase2_pockets" / "final" / "accepted_pocket_boxes.csv",
        ["pocket_family_id", "state_id", "protomer_id", "box_id", "box_center_x", "box_center_y", "box_center_z", "box_size_x", "box_size_y", "box_size_z"],
        [{"pocket_family_id": "fam_1", "state_id": "EGFR_160-185", "protomer_id": "A", "box_id": "box_1", "box_center_x": 0, "box_center_y": 0, "box_center_z": 8, "box_size_x": 8, "box_size_y": 8, "box_size_z": 8}],
    )


def write_pose(ctx, text=POSE_MULTI):
    pose = ctx.run_dir / "phase3_compounds" / "docking_outputs" / "focused_pocket_first" / "mini" / "job_1" / "pose_out.pdbqt"
    pose.parent.mkdir(parents=True, exist_ok=True)
    pose.write_text(text, encoding="utf-8")
    pose.with_name("vina.log").write_text("log\n", encoding="utf-8")
    return pose


def attr_row(ctx, pose, job_id="job_1", rank="1", model="1", cid="Cpd-A", allowed="true", hard="true", status="PASS", reject="none", pocket="true", atp="false", membrane="false", dimer="false", state="EGFR_160-185"):
    row = {field: "" for field in ATTRIBUTION_FIELDS}
    row.update(
        {
            "run_id": ctx.run_id,
            "m2_run_id": "m2_run",
            "profile": "mini",
            "job_id": job_id,
            "compound_public_id": cid,
            "state_id": state,
            "state_role": "reference" if state == "3GT8_raw" else "primary",
            "pocket_family_id": "fam_1",
            "box_id": "box_1",
            "protomer_id": "A",
            "vina_repeat_id": rank,
            "pose_rank": rank,
            "pose_model_index": model,
            "vina_affinity_kcal_mol": str(-7.5 + float(rank) / 10.0),
            "vina_log_affinity_kcal_mol": str(-7.5 + float(rank) / 10.0),
            "pose_file": ctx.relative_to_repo(pose),
            "output_pdbqt_file": ctx.relative_to_repo(pose),
            "output_pdbqt_sha256": sha256(pose) if pose.is_file() else "",
            "vina_log_file": ctx.relative_to_repo(pose.with_name("vina.log")),
            "vina_log_sha256": sha256(pose.with_name("vina.log")) if pose.with_name("vina.log").is_file() else "",
            "pose_center_x": "0.5",
            "pose_center_y": "0.0",
            "pose_center_z": "8.0",
            "pose_hard_gate_pass": hard,
            "pose_reject_reason": reject,
            "attribution_status": status,
            "allowed_for_clustering": allowed,
            "pocket_retention_pass": pocket,
            "atp_migration_flag": atp,
            "ppi_relationship_confirmed": "true",
            "ppi_hotspot_contact_count": "1",
            "ppi_contact_residue_count": "1",
            "nearest_ppi_hotspot_distance": "3.0",
            "ppi_rim_distance": "3.0",
            "membrane_penetration_flag": membrane,
            "dimer_interface_clash_flag": dimer,
            "pose_mechanism_class": "rim_blocker",
        }
    )
    return row


def write_attr(ctx, rows, ready=True, status="PASS"):
    attr = ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_attribution.csv"
    write_csv(attr, ATTRIBUTION_FIELDS, rows)
    qc = ctx.run_dir / "phase3_compounds" / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    (qc / "pose_attribution_qc.json").write_text(json.dumps({"overall_status": status, "m3_t8_clustering_ready": ready}), encoding="utf-8")
    return attr


def test_dry_run_missing_attribution_warns_and_cluster_missing_fails(tmp_path):
    ctx = make_ctx(tmp_path)
    dry = run_m3_pose_clustering(ctx, m2_run_id="m2_run", profile="mini", mode="dry-run", write_empty_clusters=True)
    fail = run_m3_pose_clustering(ctx, m2_run_id="m2_run", profile="mini", mode="cluster")

    assert dry.status == "WARN"
    assert dry.clusters_csv.is_file()
    assert fail.status == "FAIL"
    assert any("compound_pose_attribution" in item for item in fail.blockers)


def test_cluster_synthetic_nearby_poses_converges_and_writes_membership(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx)
    pose = write_pose(ctx)
    attr = write_attr(ctx, [attr_row(ctx, pose, job_id="job_1", rank="1", model="1"), attr_row(ctx, pose, job_id="job_2", rank="2", model="2")])

    result = run_m3_pose_clustering(ctx, m2_run_id="m2_run", pose_attribution_table=attr, profile="mini")
    clusters = read_csv(result.clusters_csv)
    membership = read_csv(result.membership_csv)

    assert result.status == "PASS"
    assert result.m3_t9_anchor_convergence_ready is True
    assert len(clusters) == 1
    assert clusters[0]["converged_pose_cluster"] == "true"
    assert clusters[0]["affinity_used_for_ranking"] == "false"
    assert len(membership) == 2
    assert all(row["cluster_id"] for row in membership)


def test_empty_downstream_tables_and_prior_smiles_diagnostics_do_not_block(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx)
    pose = write_pose(ctx)
    attr = write_attr(ctx, [attr_row(ctx, pose, job_id="job_1", rank="1", model="1"), attr_row(ctx, pose, job_id="job_2", rank="2", model="2")])
    write_csv(ctx.run_dir / "phase3_compounds" / "tables" / "compound_anchor_convergence.csv", ["anchor_id"], [])
    (ctx.errors_dir / "error_summary.txt").write_text("previous blocker: SMILES-like fields printed in public outputs\n", encoding="utf-8")

    result = run_m3_pose_clustering(ctx, m2_run_id="m2_run", pose_attribution_table=attr, profile="mini")

    assert result.status == "PASS"
    assert not any("forbidden downstream" in item for item in result.blockers)
    assert not any("SMILES-like fields" in item for item in result.blockers)


def test_schema_duplicate_nonpublic_and_path_validation_fail(tmp_path):
    ctx = make_ctx(tmp_path / "schema")
    bad = ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_attribution.csv"
    write_csv(bad, ["job_id"], [{"job_id": "x"}])
    assert run_m3_pose_clustering(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True).status == "FAIL"

    ctx2 = make_ctx(tmp_path / "dupe")
    write_refs(ctx2)
    pose = write_pose(ctx2)
    row = attr_row(ctx2, pose)
    write_attr(ctx2, [row, dict(row)])
    assert run_m3_pose_clustering(ctx2, m2_run_id="m2_run", profile="mini").status == "FAIL"

    ctx3 = make_ctx(tmp_path / "badid")
    write_refs(ctx3)
    pose3 = write_pose(ctx3)
    write_attr(ctx3, [attr_row(ctx3, pose3, cid="SECRET")])
    assert run_m3_pose_clustering(ctx3, m2_run_id="m2_run", profile="mini").status == "FAIL"

    ctx4 = make_ctx(tmp_path / "path")
    write_refs(ctx4)
    outside = ctx4.repo_root / "outside.pdbqt"
    outside.write_text(POSE_MULTI, encoding="utf-8")
    write_attr(ctx4, [attr_row(ctx4, outside)])
    assert run_m3_pose_clustering(ctx4, m2_run_id="m2_run", profile="mini").status == "FAIL"


def test_t7_qc_gate_and_missing_accepted_boxes(tmp_path):
    ctx = make_ctx(tmp_path / "qc")
    write_refs(ctx)
    pose = write_pose(ctx)
    write_attr(ctx, [attr_row(ctx, pose)], ready=False)

    blocked = run_m3_pose_clustering(ctx, m2_run_id="m2_run", profile="mini")
    partial = run_m3_pose_clustering(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True)
    assert blocked.status == "FAIL"
    assert partial.status in {"WARN", "PASS"}

    ctx2 = make_ctx(tmp_path / "nobox")
    pose2 = write_pose(ctx2)
    write_attr(ctx2, [attr_row(ctx2, pose2)])
    assert run_m3_pose_clustering(ctx2, m2_run_id="m2_run", profile="mini").status == "FAIL"


def test_t7_qc_must_pass_and_pose_hash_must_match(tmp_path):
    ctx = make_ctx(tmp_path / "qc_fail")
    write_refs(ctx)
    pose = write_pose(ctx)
    write_attr(ctx, [attr_row(ctx, pose)], ready=True, status="WARN")

    result = run_m3_pose_clustering(ctx, m2_run_id="m2_run", profile="mini")
    assert result.status == "FAIL"
    assert any("pose_attribution_qc.json has m3_t8_clustering_ready=false" in item for item in result.blockers)

    ctx2 = make_ctx(tmp_path / "hash")
    write_refs(ctx2)
    pose2 = write_pose(ctx2)
    row = attr_row(ctx2, pose2)
    pose2.write_text(POSE_MULTI.replace("0.200", "5.200"), encoding="utf-8")
    write_attr(ctx2, [row])

    mismatch = run_m3_pose_clustering(ctx2, m2_run_id="m2_run", profile="mini")
    membership = read_csv(mismatch.membership_csv)[0]
    assert mismatch.status == "FAIL"
    assert membership["not_clustered_reason"] == "pose_hash_mismatch"


def test_ineligible_rows_are_represented_with_explicit_reasons(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx)
    pose = write_pose(ctx)
    rows = [
        attr_row(ctx, pose, job_id="ok", rank="1", model="1"),
        attr_row(ctx, pose, job_id="reject", rank="2", model="2", hard="false", reject="outside_accepted_pocket", pocket="false"),
        attr_row(ctx, pose, job_id="atp", rank="3", model="3", atp="true"),
        attr_row(ctx, pose, job_id="noallow", rank="4", model="1", allowed="false"),
    ]
    write_attr(ctx, rows)
    result = run_m3_pose_clustering(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True)
    membership = read_csv(result.membership_csv)

    assert len(membership) == 4
    reasons = {row["job_id"]: row["not_clustered_reason"] for row in membership if row["not_clustered_reason"] != "none"}
    assert reasons["reject"] == "outside_accepted_pocket"
    assert reasons["atp"] == "ATP_migration"
    assert reasons["noallow"] == "not_allowed_for_clustering"


def test_centroid_scattered_and_singleton_clusters_warn_not_promoted(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx)
    pose = write_pose(ctx)
    attr = write_attr(ctx, [attr_row(ctx, pose, job_id="near", rank="1", model="1"), attr_row(ctx, pose, job_id="far", rank="3", model="3")])
    result = run_m3_pose_clustering(ctx, m2_run_id="m2_run", pose_attribution_table=attr, profile="mini", cluster_method="centroid")
    clusters = read_csv(result.clusters_csv)

    assert result.status == "WARN"
    assert len(clusters) == 2
    assert all(row["convergence_status"] == "WARN_SINGLETON_NOT_CONVERGED" for row in clusters)
    assert result.clustering["scattered_groups"] == 1


def test_missing_eligible_pose_file_fails_but_membership_is_written(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx)
    missing = ctx.run_dir / "phase3_compounds" / "docking_outputs" / "focused_pocket_first" / "mini" / "job_1" / "missing.pdbqt"
    write_attr(ctx, [attr_row(ctx, missing)])

    result = run_m3_pose_clustering(ctx, m2_run_id="m2_run", profile="mini", write_empty_clusters=True)
    membership = read_csv(result.membership_csv)

    assert result.status == "FAIL"
    assert membership[0]["not_clustered_reason"] == "missing_pose_file"
    assert membership[0]["cluster_membership_status"] == "FAIL"


def test_no_coordinate_lines_or_downstream_tables_created_and_snapshot_written(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx)
    pose = write_pose(ctx)
    write_attr(ctx, [attr_row(ctx, pose, job_id="job_1", rank="1", model="1"), attr_row(ctx, pose, job_id="job_2", rank="2", model="2")])
    result = run_m3_pose_clustering(ctx, m2_run_id="m2_run", profile="mini")

    public_text = "\n".join(
        path.read_text(encoding="utf-8", errors="ignore")
        for root in [ctx.run_dir / "phase3_compounds" / "tables", ctx.run_dir / "phase3_compounds" / "qc", ctx.run_dir / "phase3_compounds" / "reports", ctx.logs_dir]
        if root.exists()
        for path in root.rglob("*")
        if path.is_file() and "docking_outputs" not in path.parts
    )
    assert "ATOM      1" not in public_text
    assert "SECRET_A" not in public_text
    assert not (ctx.run_dir / "phase3_compounds" / "tables" / "compound_anchor_convergence.csv").exists()
    assert not (ctx.run_dir / "phase3_compounds" / "tables" / "final_m3_candidate_hypotheses.csv").exists()
    assert (ctx.run_dir / "phase3_compounds" / "config_snapshots" / "compound_pose_clustering.resolved.yaml").is_file()
    assert result.qc_json.is_file()
