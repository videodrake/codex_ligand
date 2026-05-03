import csv
import hashlib
import json
from pathlib import Path

from egfr_myo1d.compound.pdbqt_parse import extract_pose_atoms
from egfr_myo1d.compound.pose_attribution import run_m3_pose_attribution
from egfr_myo1d.compound.vina_collect import RAW_FIELDS
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext


POSE_MULTI = """MODEL 1
REMARK VINA RESULT: -7.500 0.000 0.000
ATOM      1  C   LIG A   1       0.000   0.000   8.000  0.00  0.00     0.000 C
ATOM      2  N   LIG A   1       1.000   0.000   8.000  0.00  0.00     0.000 N
ENDMDL
MODEL 2
REMARK VINA RESULT: -6.250 0.000 0.000
ATOM      1  C   LIG A   1       9.000   9.000   8.000  0.00  0.00     0.000 C
ATOM      2  O   LIG A   1      10.000   9.000   8.000  0.00  0.00     0.000 O
ENDMDL
"""

POSE_SINGLE = """REMARK VINA RESULT: -5.000 0.000 0.000
ATOM      1  C   LIG A   1       0.000   0.000   8.000  0.00  0.00     0.000 C
ATOM      2  O   LIG A   1       1.000   0.000   8.000  0.00  0.00     0.000 O
END
"""


def make_ctx(tmp_path, run_id="m3_pose_attr"):
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


def write_receptor_pdb(path, state_residue=900, xyz=(50, 50, 50), chain="A", resname="GLY"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "ATOM  {serial:5d}  CA  {resname:>3s} {chain:1s}{resseq:4d}    "
        "{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\nEND\n".format(
            serial=1,
            resname=resname,
            chain=chain,
            resseq=state_residue,
            x=xyz[0],
            y=xyz[1],
            z=xyz[2],
        ),
        encoding="utf-8",
    )


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_refs(ctx, m2_run_id="m2_run", atp=(40, 40, 40), ppi=(3, 0, 8), membrane_bounds=(-2, 2), mapping_point=(50, 50, 50)):
    m2 = ctx.fresh_root / "runs" / m2_run_id
    write_csv(m2 / "phase2_pockets" / "final" / "accepted_pockets_for_m3.csv", ["pocket_family_id", "state_id"], [{"pocket_family_id": "fam_1", "state_id": "EGFR_160-185"}])
    write_csv(
        m2 / "phase2_pockets" / "final" / "accepted_pocket_boxes.csv",
        ["pocket_family_id", "state_id", "protomer_id", "box_id", "box_center_x", "box_center_y", "box_center_z", "box_size_x", "box_size_y", "box_size_z"],
        [{"pocket_family_id": "fam_1", "state_id": "EGFR_160-185", "protomer_id": "A", "box_id": "box_1", "box_center_x": 0, "box_center_y": 0, "box_center_z": 8, "box_size_x": 8, "box_size_y": 8, "box_size_z": 8}],
    )
    write_csv(m2 / "phase2_pockets" / "atp_reference" / "atp_site_reference.csv", ["atp_id", "atp_center_x", "atp_center_y", "atp_center_z"], [{"atp_id": "ATP_REF", "atp_center_x": atp[0], "atp_center_y": atp[1], "atp_center_z": atp[2]}])
    write_csv(m2 / "phase1_ppi" / "consensus" / "ppi_consensus_patch_merged.csv", ["hotspot_id", "residue_public", "x", "y", "z"], [{"hotspot_id": "HS1", "residue_public": "A:900", "x": ppi[0], "y": ppi[1], "z": ppi[2]}])
    membrane = m2 / "membrane" / "membrane_frame.json"
    membrane.parent.mkdir(parents=True, exist_ok=True)
    membrane.write_text(json.dumps({"origin": [0, 0, 0], "normal": [0, 0, 1], "core_z_min": membrane_bounds[0], "core_z_max": membrane_bounds[1]}), encoding="utf-8")
    write_csv(m2 / "phase1_receptor" / "receptor_mapping.csv", ["protomer_id", "role", "x", "y", "z"], [{"protomer_id": "A", "role": "central_interface", "x": mapping_point[0], "y": mapping_point[1], "z": mapping_point[2]}])


def raw_row(ctx, pose, job_id="job_1", pose_rank="1", model="1", cid="Cpd-A", parse_status="PASS", allowed="true", state="EGFR_160-185"):
    row = {field: "" for field in RAW_FIELDS}
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
            "vina_repeat_id": "1",
            "pose_rank": pose_rank,
            "pose_model_index": model,
            "vina_affinity_kcal_mol": "-7.5",
            "vina_log_affinity_kcal_mol": "-7.5",
            "pose_file": ctx.relative_to_repo(pose),
            "output_pdbqt_file": ctx.relative_to_repo(pose),
            "output_pdbqt_sha256": sha256(pose) if pose.is_file() else "",
            "vina_log_file": ctx.relative_to_repo(pose.with_name("vina.log")),
            "vina_log_sha256": sha256(pose.with_name("vina.log")) if pose.with_name("vina.log").is_file() else "",
            "pose_center_x": "0.5",
            "pose_center_y": "0.0",
            "pose_center_z": "8.0",
            "pose_atom_count": "2",
            "pose_heavy_atom_count": "2",
            "source_vina_job_manifest": "manifest",
            "parse_status": parse_status,
            "collected_at": "2026-04-30T00:00:00+00:00",
            "allowed_for_attribution": allowed,
        }
    )
    return row


def write_raw(ctx, rows):
    raw = ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_raw.csv"
    write_csv(raw, RAW_FIELDS, rows)
    qc = ctx.run_dir / "phase3_compounds" / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    (qc / "docking_completion_qc.json").write_text(json.dumps({"overall_status": "PASS", "m3_t7_attribution_ready": True}), encoding="utf-8")
    return raw


def write_pose(ctx, text=POSE_MULTI, name="pose_out.pdbqt"):
    pose = ctx.run_dir / "phase3_compounds" / "docking_outputs" / "focused_pocket_first" / "mini" / "job_1" / name
    pose.parent.mkdir(parents=True, exist_ok=True)
    pose.write_text(text, encoding="utf-8")
    pose.with_name("vina.log").write_text("log\n", encoding="utf-8")
    return pose


def test_dry_run_missing_raw_warns_and_attribute_missing_raw_fails(tmp_path):
    ctx = make_ctx(tmp_path)
    dry = run_m3_pose_attribution(ctx, m2_run_id="m2_run", profile="mini", mode="dry-run", write_empty_attribution=True)
    fail = run_m3_pose_attribution(ctx, m2_run_id="m2_run", profile="mini", mode="attribute")

    assert dry.status == "WARN"
    assert dry.attribution_csv.is_file()
    assert fail.status == "FAIL"
    assert any("compound_pose_raw" in item for item in fail.blockers)


def test_extract_pose_atoms_multimodel_and_single_model(tmp_path):
    multi = tmp_path / "multi.pdbqt"
    single = tmp_path / "single.pdbqt"
    multi.write_text(POSE_MULTI, encoding="utf-8")
    single.write_text(POSE_SINGLE, encoding="utf-8")

    assert len(extract_pose_atoms(multi, pose_model_index=2, pose_rank=2)) == 2
    assert extract_pose_atoms(multi, pose_model_index=2, pose_rank=2)[0].x == 9.0
    assert len(extract_pose_atoms(single, pose_model_index=None, pose_rank=1)) == 2


def test_attribute_synthetic_pose_writes_all_qc_tables(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx)
    pose = write_pose(ctx)
    raw = write_raw(ctx, [raw_row(ctx, pose)])

    result = run_m3_pose_attribution(ctx, m2_run_id="m2_run", pose_raw_table=raw, profile="mini")
    attribution = read_csv(result.attribution_csv)

    assert result.status == "PASS"
    assert result.m3_t8_clustering_ready is True
    assert len(attribution) == 1
    assert attribution[0]["pose_hard_gate_pass"] == "true"
    assert attribution[0]["pose_mechanism_class"] in {"orthosteric_or_direct_PPI_patch_blocker", "rim_blocker"}
    assert len(read_csv(result.atp_qc_csv)) == 1
    assert len(read_csv(result.pocket_qc_csv)) == 1
    assert len(read_csv(result.ppi_qc_csv)) == 1
    assert len(read_csv(result.membrane_qc_csv)) == 1


def test_rejected_pose_rows_do_not_block_clustering_when_pass_rows_exist(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx)
    pose = write_pose(ctx, text=POSE_MULTI)
    raw = write_raw(
        ctx,
        [
            raw_row(ctx, pose, job_id="job_pass", pose_rank="1", model="1"),
            raw_row(ctx, pose, job_id="job_reject", pose_rank="2", model="2"),
        ],
    )

    result = run_m3_pose_attribution(ctx, m2_run_id="m2_run", pose_raw_table=raw, profile="mini")
    attribution = read_csv(result.attribution_csv)

    assert result.status == "PASS"
    assert result.m3_t8_clustering_ready is True
    assert sum(1 for row in attribution if row["pose_hard_gate_pass"] == "true") == 1
    assert sum(1 for row in attribution if row["attribution_status"] == "REJECT") == 1


def test_actual_m2_export_paths_and_centroid_columns_resolve(tmp_path):
    ctx = make_ctx(tmp_path)
    m2 = ctx.fresh_root / "runs" / "m2_run"
    write_csv(
        m2 / "phase2_pockets" / "export_for_m3" / "accepted_pockets_for_m3.csv",
        ["pocket_family_id", "state_id"],
        [{"pocket_family_id": "fam_1", "state_id": "EGFR_160-185"}, {"pocket_family_id": "fam_1", "state_id": "EGFR_170-200"}],
    )
    receptor_160 = m2 / "normalized" / "receptors" / "EGFR_160-185_dockable_669_1014_explicit_AB.pdb"
    receptor_170 = m2 / "normalized" / "receptors" / "EGFR_170-200_dockable_669_1014_explicit_AB.pdb"
    write_receptor_pdb(receptor_160)
    write_receptor_pdb(receptor_170)
    write_csv(
        m2 / "phase2_pockets" / "export_for_m3" / "accepted_pocket_boxes.csv",
        ["pocket_family_id", "state_id", "protomer_id", "box_id", "box_center_x", "box_center_y", "box_center_z", "box_size_x", "box_size_y", "box_size_z", "receptor_pdb", "receptor_mapping_csv"],
        [
            {"pocket_family_id": "fam_1", "state_id": "EGFR_160-185", "protomer_id": "A", "box_id": "box_1", "box_center_x": 0, "box_center_y": 0, "box_center_z": 8, "box_size_x": 8, "box_size_y": 8, "box_size_z": 8, "receptor_pdb": ctx.relative_to_repo(receptor_160), "receptor_mapping_csv": "fresh/runs/m2_run/qc/EGFR_160-185_receptor_mapping.csv"},
            {"pocket_family_id": "fam_1", "state_id": "EGFR_170-200", "protomer_id": "A", "box_id": "box_1", "box_center_x": 0, "box_center_y": 0, "box_center_z": 8, "box_size_x": 8, "box_size_y": 8, "box_size_z": 8, "receptor_pdb": ctx.relative_to_repo(receptor_170), "receptor_mapping_csv": "fresh/runs/m2_run/qc/EGFR_170-200_receptor_mapping.csv"},
        ],
    )
    write_csv(
        m2 / "phase2_pockets" / "atp_reference" / "atp_site_reference.csv",
        ["uniprot_residue_number", "residue_name"],
        [{"uniprot_residue_number": "900", "residue_name": "ATP_SITE"}],
    )
    write_csv(
        m2 / "phase2_pockets" / "atp_reference" / "atp_site_centroid_by_state.csv",
        ["atp_id", "state_id", "egfr_protomer_id", "atp_centroid_x", "atp_centroid_y", "atp_centroid_z"],
        [
            {"atp_id": "ATP_160", "state_id": "EGFR_160-185", "egfr_protomer_id": "A", "atp_centroid_x": 40, "atp_centroid_y": 40, "atp_centroid_z": 40},
            {"atp_id": "ATP_170", "state_id": "EGFR_170-200", "egfr_protomer_id": "A", "atp_centroid_x": 40, "atp_centroid_y": 40, "atp_centroid_z": 40},
        ],
    )
    write_csv(
        m2 / "phase2_pockets" / "export_for_m3" / "ppi_consensus_patch.csv",
        ["ppi_patch_id", "state_id", "egfr_protomer_id", "egfr_residue_number", "egfr_contact_centroid_x", "egfr_contact_centroid_y", "egfr_contact_centroid_z"],
        [
            {"ppi_patch_id": "ppi_160", "state_id": "EGFR_160-185", "egfr_protomer_id": "A", "egfr_residue_number": "900", "egfr_contact_centroid_x": 3, "egfr_contact_centroid_y": 0, "egfr_contact_centroid_z": 8},
            {"ppi_patch_id": "ppi_170", "state_id": "EGFR_170-200", "egfr_protomer_id": "A", "egfr_residue_number": "900", "egfr_contact_centroid_x": 3, "egfr_contact_centroid_y": 0, "egfr_contact_centroid_z": 8},
        ],
    )
    membrane = m2 / "manifest" / "membrane_frame.json"
    membrane.parent.mkdir(parents=True, exist_ok=True)
    membrane.write_text(json.dumps({"origin": [0, 0, 0], "normal": [0, 0, 1], "core_z_min": -2, "core_z_max": 2}), encoding="utf-8")
    mapping_fields = ["state", "source_file", "protomer_id", "source_chain", "source_resseq", "source_icode", "source_resname", "runtime_chain", "runtime_resseq", "atom_count", "role"]
    write_csv(m2 / "qc" / "EGFR_160-185_receptor_mapping.csv", mapping_fields, [{"state": "EGFR_160-185", "source_file": "src.pdb", "protomer_id": "A", "source_chain": "A", "source_resseq": "900", "source_icode": "_", "source_resname": "GLY", "runtime_chain": "A", "runtime_resseq": "900", "atom_count": "1", "role": "primary"}])
    write_csv(m2 / "qc" / "EGFR_170-200_receptor_mapping.csv", mapping_fields, [{"state": "EGFR_170-200", "source_file": "src.pdb", "protomer_id": "A", "source_chain": "A", "source_resseq": "900", "source_icode": "_", "source_resname": "GLY", "runtime_chain": "A", "runtime_resseq": "900", "atom_count": "1", "role": "primary"}])
    write_csv(ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_clusters.csv", ["cluster_id"], [])
    broad = ctx.run_dir / "phase3_compounds" / "docking_outputs" / "broad_anchor_scan_optional"
    broad.mkdir(parents=True, exist_ok=True)
    (broad / ".gitkeep").write_text("", encoding="utf-8")
    (ctx.errors_dir / "error_summary.txt").write_text("previous blocker: SMILES-like fields printed in public outputs\n", encoding="utf-8")
    pose = write_pose(ctx)
    raw = write_raw(ctx, [raw_row(ctx, pose, job_id="job_160"), raw_row(ctx, pose, job_id="job_170", state="EGFR_170-200")])

    result = run_m3_pose_attribution(ctx, m2_run_id="m2_run", pose_raw_table=raw, profile="mini")
    attribution = read_csv(result.attribution_csv)

    assert result.status == "PASS"
    assert result.m3_t8_clustering_ready is True
    assert {row["state_id"] for row in attribution} == {"EGFR_160-185", "EGFR_170-200"}
    assert {row["pose_hard_gate_pass"] for row in attribution} == {"true"}
    assert not any("SMILES" in item for item in result.blockers)


def test_schema_duplicate_public_id_and_pose_path_gates(tmp_path):
    ctx = make_ctx(tmp_path / "schema")
    pose = write_pose(ctx)
    bad_raw = ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_raw.csv"
    write_csv(bad_raw, ["job_id"], [{"job_id": "x"}])
    assert run_m3_pose_attribution(ctx, m2_run_id="m2_run", profile="mini").status == "FAIL"

    ctx2 = make_ctx(tmp_path / "dupe")
    write_refs(ctx2)
    pose2 = write_pose(ctx2)
    row = raw_row(ctx2, pose2)
    write_raw(ctx2, [row, dict(row)])
    assert run_m3_pose_attribution(ctx2, m2_run_id="m2_run", profile="mini").status == "FAIL"

    ctx3 = make_ctx(tmp_path / "badid")
    write_refs(ctx3)
    pose3 = write_pose(ctx3)
    write_raw(ctx3, [raw_row(ctx3, pose3, cid="SECRET")])
    assert run_m3_pose_attribution(ctx3, m2_run_id="m2_run", profile="mini").status == "FAIL"

    ctx4 = make_ctx(tmp_path / "path")
    write_refs(ctx4)
    outside = ctx4.repo_root / "outside.pdbqt"
    outside.write_text(POSE_SINGLE, encoding="utf-8")
    write_raw(ctx4, [raw_row(ctx4, outside)])
    assert run_m3_pose_attribution(ctx4, m2_run_id="m2_run", profile="mini").status == "FAIL"


def test_requires_m3_t6_ready_and_pose_hash_match(tmp_path):
    ctx = make_ctx(tmp_path / "not_ready")
    write_refs(ctx)
    pose = write_pose(ctx)
    write_raw(ctx, [raw_row(ctx, pose)])
    qc = ctx.run_dir / "phase3_compounds" / "qc" / "docking_completion_qc.json"
    qc.write_text(json.dumps({"overall_status": "WARN", "m3_t7_attribution_ready": False}), encoding="utf-8")

    not_ready = run_m3_pose_attribution(ctx, m2_run_id="m2_run", profile="mini")
    assert not_ready.status == "FAIL"
    assert any("M3-T6 docking completion QC did not PASS" in item for item in not_ready.blockers)

    ctx2 = make_ctx(tmp_path / "hash")
    write_refs(ctx2)
    pose2 = write_pose(ctx2)
    row = raw_row(ctx2, pose2)
    pose2.write_text(POSE_SINGLE, encoding="utf-8")
    write_raw(ctx2, [row])

    mismatch = run_m3_pose_attribution(ctx2, m2_run_id="m2_run", profile="mini")
    attr = read_csv(mismatch.attribution_csv)[0]
    assert mismatch.status == "FAIL"
    assert attr["pose_reject_reason"] == "pose_hash_mismatch"


def test_missing_references_fail_attribute_mode(tmp_path):
    ctx = make_ctx(tmp_path)
    pose = write_pose(ctx)
    write_raw(ctx, [raw_row(ctx, pose)])

    result = run_m3_pose_attribution(ctx, m2_run_id="m2_run", profile="mini")

    assert result.status == "FAIL"
    assert any("atp_reference" in item or "accepted_boxes" in item for item in result.blockers)


def test_pose_references_are_state_and_protomer_matched(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx)
    m2 = ctx.fresh_root / "runs" / "m2_run"
    write_csv(
        m2 / "phase1_ppi" / "consensus" / "ppi_consensus_patch_merged.csv",
        ["hotspot_id", "residue_public", "state_id", "protomer_id", "x", "y", "z"],
        [{"hotspot_id": "HS_OTHER", "residue_public": "B:900", "state_id": "EGFR_other", "protomer_id": "B", "x": 0, "y": 0, "z": 8}],
    )
    pose = write_pose(ctx)
    write_raw(ctx, [raw_row(ctx, pose)])

    result = run_m3_pose_attribution(ctx, m2_run_id="m2_run", profile="mini")
    attr = read_csv(result.attribution_csv)[0]

    assert result.status == "WARN"
    assert result.m3_t8_clustering_ready is False
    assert attr["pose_reject_reason"] == "missing_mapping"
    assert attr["ppi_relationship_confirmed"] == "false"


def test_outside_pocket_atp_membrane_and_dimer_rejections(tmp_path):
    ctx = make_ctx(tmp_path / "outside")
    write_refs(ctx)
    pose = write_pose(ctx, text=POSE_MULTI)
    write_raw(ctx, [raw_row(ctx, pose, pose_rank="2", model="2")])
    out = run_m3_pose_attribution(ctx, m2_run_id="m2_run", profile="mini")
    attr = read_csv(out.attribution_csv)[0]
    assert attr["pose_reject_reason"] == "outside_accepted_pocket"

    ctx2 = make_ctx(tmp_path / "atp")
    write_refs(ctx2, atp=(0, 0, 8))
    pose2 = write_pose(ctx2)
    write_raw(ctx2, [raw_row(ctx2, pose2)])
    attr2 = read_csv(run_m3_pose_attribution(ctx2, m2_run_id="m2_run", profile="mini").attribution_csv)[0]
    assert attr2["atp_migration_flag"] == "true"
    assert attr2["pose_mechanism_class"] == "ATP_like_reject"

    ctx3 = make_ctx(tmp_path / "mem")
    write_refs(ctx3, membrane_bounds=(7, 9))
    pose3 = write_pose(ctx3)
    write_raw(ctx3, [raw_row(ctx3, pose3)])
    attr3 = read_csv(run_m3_pose_attribution(ctx3, m2_run_id="m2_run", profile="mini").attribution_csv)[0]
    assert attr3["pose_reject_reason"] == "membrane_penetration"

    ctx4 = make_ctx(tmp_path / "clash")
    write_refs(ctx4, mapping_point=(0, 0, 8))
    pose4 = write_pose(ctx4)
    write_raw(ctx4, [raw_row(ctx4, pose4)])
    attr4 = read_csv(run_m3_pose_attribution(ctx4, m2_run_id="m2_run", profile="mini").attribution_csv)[0]
    assert attr4["pose_reject_reason"] == "central_dimer_clash"


def test_parse_fail_and_not_allowed_rows_are_not_dropped(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx)
    pose = write_pose(ctx)
    write_raw(ctx, [raw_row(ctx, pose, job_id="parse_fail", parse_status="FAIL"), raw_row(ctx, pose, job_id="not_allowed", allowed="false")])

    result = run_m3_pose_attribution(ctx, m2_run_id="m2_run", profile="mini", allow_partial=True, require_all_poses=False)
    attribution = read_csv(result.attribution_csv)

    assert len(attribution) == 2
    assert {row["pose_reject_reason"] for row in attribution} == {"parse_failure", "not_allowed_for_attribution"}
    assert len(read_csv(result.atp_qc_csv)) == 2


def test_no_coordinate_lines_or_downstream_tables_are_created(tmp_path):
    ctx = make_ctx(tmp_path)
    write_refs(ctx)
    pose = write_pose(ctx)
    write_raw(ctx, [raw_row(ctx, pose)])
    result = run_m3_pose_attribution(ctx, m2_run_id="m2_run", profile="mini")

    public_text = "\n".join(
        path.read_text(encoding="utf-8", errors="ignore")
        for root in [ctx.run_dir / "phase3_compounds" / "tables", ctx.run_dir / "phase3_compounds" / "qc", ctx.run_dir / "phase3_compounds" / "reports", ctx.logs_dir]
        if root.exists()
        for path in root.rglob("*")
        if path.is_file() and "docking_outputs" not in path.parts
    )
    assert "ATOM      1" not in public_text
    assert "SECRET_A" not in public_text
    assert not (ctx.run_dir / "phase3_compounds" / "tables" / "compound_pose_clusters.csv").exists()
    assert not (ctx.run_dir / "phase3_compounds" / "tables" / "final_m3_candidate_hypotheses.csv").exists()
    assert result.mechanism_classes["orthosteric_or_direct_PPI_patch_blocker"] + result.mechanism_classes["rim_blocker"] >= 1
