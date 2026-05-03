import csv
import json
from pathlib import Path

from egfr_myo1d.compound.m3_readiness import run_m3_readiness
from egfr_myo1d.core.run_context import RunContext


def make_ctx(tmp_path, run_id="m3_run"):
    repo = tmp_path / "repo"
    fresh = repo / "fresh"
    run_dir = fresh / "runs" / run_id
    return RunContext(
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


def write_csv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_m1_inputs(fresh_root):
    for state in ["EGFR_160-185", "EGFR_170-200", "3GT8_raw"]:
        state_dir = fresh_root / "data" / "normalized" / "receptors" / state
        state_dir.mkdir(parents=True, exist_ok=True)
        dockable = "dockable_reference_explicit_AB.pdb" if state == "3GT8_raw" else "dockable_669_1014_explicit_AB.pdb"
        (state_dir / dockable).write_text("ATOM      1  CA  GLY A   1       0.0 0.0 0.0\n", encoding="utf-8")
        (state_dir / "runtime_offset_receptor_only.pdb").write_text("ATOM      1  CA  GLY A   1       0.0 0.0 0.0\n", encoding="utf-8")
        write_csv(
            state_dir / "receptor_mapping.csv",
            ["state_id", "egfr_chain_id", "egfr_protomer_id", "receptor_residue_number", "uniprot_residue_number"],
            [
                {
                    "state_id": state,
                    "egfr_chain_id": "A",
                    "egfr_protomer_id": "A",
                    "receptor_residue_number": "900",
                    "uniprot_residue_number": "900",
                }
            ],
        )
    geometry = fresh_root / "data" / "normalized" / "geometry"
    geometry.mkdir(parents=True, exist_ok=True)
    (geometry / "membrane_frame.json").write_text('{"frame_id":"test","normal":[0,0,1]}\n', encoding="utf-8")
    write_csv(geometry / "frame_qc.csv", ["status"], [{"status": "PASS"}])
    write_csv(geometry / "dimer_geometry_qc.csv", ["status"], [{"status": "PASS"}])


def write_ligands(fresh_root):
    ligand_dir = fresh_root / "data" / "raw" / "ligands"
    ligand_dir.mkdir(parents=True, exist_ok=True)
    for compound_id in ["Cpd-A", "Cpd-B", "Cpd-C"]:
        (ligand_dir / f"{compound_id}.smi").write_text("C\n", encoding="utf-8")


def write_m2_package(
    fresh_root,
    m2_run_id="m2_run",
    accepted_rows=None,
    box_size="12.0",
    m2_review_go=True,
):
    m2_dir = fresh_root / "runs" / m2_run_id
    export = m2_dir / "phase2_pockets" / "export_for_m3"
    export.mkdir(parents=True, exist_ok=True)
    if accepted_rows is None:
        accepted_rows = [
            {
                "pocket_family_id": "fam_primary",
                "export_tier": "primary",
                "gate_class": "accepted_primary_pocket",
                "representative_candidate_pocket_id": "cand1",
                "representative_state_id": "EGFR_160-185",
                "representative_state_role": "primary",
                "representative_protomer_id": "A",
                "member_candidate_pocket_ids": "cand1",
                "member_state_ids": "EGFR_160-185",
                "member_state_roles": "primary",
                "member_protomer_ids": "A",
                "primary_state_count": "1",
                "reference_state_count": "0",
                "hard_gate_pass": "True",
                "non_atp_pass": "True",
                "ppi_relationship_pass": "True",
                "lower_lateral_pass": "True",
                "dimer_accessibility_pass": "True",
                "primary_state_support_or_review_flag": "primary_supported",
                "ppi_relationship_class": "orthosteric",
                "atp_overlap_fraction": "0.0",
                "atp_centroid_distance": "20.0",
                "z_percentile_within_receptor_surface": "0.25",
                "lateral_score": "0.5",
                "interface_overlap_fraction": "0.0",
                "sasa_ratio_dimer_vs_isolated_protomer": "",
                "family_center_x": "1.0",
                "family_center_y": "2.0",
                "family_center_z": "3.0",
                "residue_union_uniprot": "900",
                "residue_union_state_residue_ids": "EGFR_160-185:A:900",
                "box_export_status": "PASS",
                "m3_primary_allowed": "True",
                "m3_review_allowed": "False",
                "reference_only": "False",
                "acceptance_reason": "synthetic accepted primary",
                "review_notes": "",
                "warnings": "",
            }
        ]
    accepted_fields = list(accepted_rows[0].keys()) if accepted_rows else [
        "pocket_family_id",
        "gate_class",
        "export_tier",
        "representative_state_id",
        "representative_protomer_id",
        "hard_gate_pass",
        "non_atp_pass",
        "lower_lateral_pass",
        "dimer_accessibility_pass",
        "reference_only",
    ]
    write_csv(export / "accepted_pockets_for_m3.csv", accepted_fields, accepted_rows)
    write_csv(
        export / "accepted_pocket_boxes.csv",
        [
            "pocket_family_id",
            "box_id",
            "export_tier",
            "state_id",
            "state_role",
            "protomer_id",
            "chain_id",
            "receptor_pdb",
            "receptor_mapping_csv",
            "box_center_x",
            "box_center_y",
            "box_center_z",
            "box_size_x",
            "box_size_y",
            "box_size_z",
            "box_margin_angstrom",
            "box_source",
            "box_qc_status",
            "box_qc_notes",
            "non_atp_pass",
            "ppi_relationship_class",
            "lower_lateral_pass",
            "dimer_accessibility_pass",
            "reference_only",
            "m3_primary_allowed",
            "warnings",
        ],
        [
            {
                "pocket_family_id": "fam_primary",
                "box_id": "box_fam_primary",
                "export_tier": "primary",
                "state_id": "EGFR_160-185",
                "state_role": "primary",
                "protomer_id": "A",
                "chain_id": "A",
                "receptor_pdb": "fresh/data/normalized/receptors/EGFR_160-185/dockable_669_1014_explicit_AB.pdb",
                "receptor_mapping_csv": "fresh/data/normalized/receptors/EGFR_160-185/receptor_mapping.csv",
                "box_center_x": "1.0",
                "box_center_y": "2.0",
                "box_center_z": "3.0",
                "box_size_x": box_size,
                "box_size_y": "12.0",
                "box_size_z": "12.0",
                "box_margin_angstrom": "4.0",
                "box_source": "synthetic",
                "box_qc_status": "PASS",
                "box_qc_notes": "",
                "non_atp_pass": "True",
                "ppi_relationship_class": "orthosteric",
                "lower_lateral_pass": "True",
                "dimer_accessibility_pass": "True",
                "reference_only": "False",
                "m3_primary_allowed": "True",
                "warnings": "",
            }
        ],
    )
    (export / "accepted_pocket_residue_sets.json").write_text('{"pockets":[]}\n', encoding="utf-8")
    write_csv(
        export / "ppi_consensus_patch.csv",
        ["state_id", "egfr_chain_id", "egfr_protomer_id", "egfr_residue_number", "uniprot_residue_number"],
        [{"state_id": "EGFR_160-185", "egfr_chain_id": "A", "egfr_protomer_id": "A", "egfr_residue_number": "900", "uniprot_residue_number": "900"}],
    )
    write_csv(export / "pocket_gate_qc.csv", ["pocket_family_id", "gate_class"], [{"pocket_family_id": "fam_primary", "gate_class": "accepted_primary_pocket"}])
    (export / "m2_final_report.md").write_text("# M2 final report\n", encoding="utf-8")
    (export / "m2_export_manifest.json").write_text('{"schema_version":"test"}\n', encoding="utf-8")
    atp_dir = m2_dir / "phase2_pockets" / "atp_reference"
    write_csv(
        atp_dir / "atp_site_reference.csv",
        ["uniprot_residue_number", "residue_name", "evidence_status"],
        [{"uniprot_residue_number": "800", "residue_name": "LYS", "evidence_status": "PASS"}],
    )
    review_dir = m2_dir / "reviews"
    review_dir.mkdir(parents=True, exist_ok=True)
    review_payload = {
        "m3_go_no_go": "GO" if m2_review_go else "NO_GO",
        "m3_docking_allowed": bool(m2_review_go),
        "blockers": [] if m2_review_go else ["synthetic review blocker"],
    }
    (review_dir / "m2_review_status.json").write_text(json.dumps(review_payload), encoding="utf-8")
    return m2_dir


def write_m2_run_local_receptor_inputs(m2_dir, state="EGFR_160-185"):
    receptor_dir = m2_dir / "prepared" / "m2_1_ppi_inputs" / state / "receptor"
    receptor_dir.mkdir(parents=True, exist_ok=True)
    (receptor_dir / f"{state}_runtime_offset_receptor_only.pdb").write_text(
        "ATOM      1  CA  GLY A 900       1.000   2.000   3.000  1.00 20.00           C\n"
        "ATOM      2  CA  GLY B 900       4.000   5.000   6.000  1.00 20.00           C\n",
        encoding="utf-8",
    )
    write_csv(
        m2_dir / "qc" / f"{state}_receptor_mapping.csv",
        ["state_id", "egfr_chain_id", "egfr_protomer_id", "receptor_residue_number", "uniprot_residue_number"],
        [
            {
                "state_id": state,
                "egfr_chain_id": "A",
                "egfr_protomer_id": "A",
                "receptor_residue_number": "900",
                "uniprot_residue_number": "900",
            }
        ],
    )


def write_m2_final_status_allows_m3(m2_dir):
    final = m2_dir / "phase2_pockets" / "final"
    final.mkdir(parents=True, exist_ok=True)
    (final / "m2_final_status.json").write_text(
        json.dumps(
            {
                "schema_version": "m2_8_final_status_v1",
                "m3_docking_allowed": True,
                "synthetic_fixture": False,
                "blockers": [],
            }
        ),
        encoding="utf-8",
    )


def test_m3_readiness_valid_minimal_inputs_pass_in_docking_mode(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m1_inputs(ctx.fresh_root)
    write_ligands(ctx.fresh_root)
    write_m2_package(ctx.fresh_root)

    result = run_m3_readiness(ctx, "m2_run", mode="docking")
    summary = read_json(result.summary_json)

    assert result.status == "PASS"
    assert summary["m3_docking_allowed"] is True
    assert summary["counts"]["docking_ready_pocket_rows"] == 1
    assert (ctx.run_dir / "phase3_compounds" / "prepared_ligands" / "Cpd-A").is_dir()
    assert not list((ctx.run_dir / "phase3_compounds").rglob("*.pdbqt"))


def test_m3_readiness_accepts_zone_composite_boxes_for_surface_zone_handoff(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m1_inputs(ctx.fresh_root)
    write_ligands(ctx.fresh_root)
    m2_dir = write_m2_package(ctx.fresh_root)
    export = m2_dir / "phase2_pockets" / "export_for_m3"
    box_rows = read_csv(export / "accepted_pocket_boxes.csv")
    write_csv(
        export / "zone_composite_boxes.csv",
        [
            "zone_id",
            "pocket_family_id",
            "box_id",
            "state_id",
            "protomer_id",
            "box_center_x",
            "box_center_y",
            "box_center_z",
            "box_size_x",
            "box_size_y",
            "box_size_z",
            "box_source",
            "source_box_id",
        ],
        [
            {
                "zone_id": "zone_fam_primary",
                "pocket_family_id": row["pocket_family_id"],
                "box_id": row["box_id"],
                "state_id": row["state_id"],
                "protomer_id": row["protomer_id"],
                "box_center_x": row["box_center_x"],
                "box_center_y": row["box_center_y"],
                "box_center_z": row["box_center_z"],
                "box_size_x": row["box_size_x"],
                "box_size_y": row["box_size_y"],
                "box_size_z": row["box_size_z"],
                "box_source": "zone_composite_boxes",
                "source_box_id": row["box_id"],
            }
            for row in box_rows
        ],
    )
    (export / "accepted_pocket_boxes.csv").unlink()

    result = run_m3_readiness(ctx, "m2_run", mode="docking")

    assert result.status == "PASS"
    path_rows = {row["input_id"]: row for row in read_csv(result.path_resolution_csv)}
    assert path_rows["accepted_pocket_boxes"]["observed_path"].endswith("zone_composite_boxes.csv")


def test_m3_readiness_missing_accepted_pockets_fails(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m1_inputs(ctx.fresh_root)
    write_ligands(ctx.fresh_root)
    m2_dir = write_m2_package(ctx.fresh_root)
    (m2_dir / "phase2_pockets" / "export_for_m3" / "accepted_pockets_for_m3.csv").unlink()

    result = run_m3_readiness(ctx, "m2_run", mode="docking")

    assert result.status == "FAIL"
    assert any("accepted_pockets_for_m3.csv" in item for item in result.blockers)


def test_m3_readiness_zero_accepted_pockets_fails(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m1_inputs(ctx.fresh_root)
    write_ligands(ctx.fresh_root)
    write_m2_package(ctx.fresh_root, accepted_rows=[])

    result = run_m3_readiness(ctx, "m2_run", mode="docking")

    assert result.status == "FAIL"
    assert any("zero accepted rows" in item for item in result.blockers)


def test_m3_readiness_missing_ligands_warn_in_dry_run_and_fail_in_docking(tmp_path):
    dry_ctx = make_ctx(tmp_path / "dry")
    write_m1_inputs(dry_ctx.fresh_root)
    write_m2_package(dry_ctx.fresh_root)
    dry = run_m3_readiness(dry_ctx, "m2_run", mode="dry-run")
    assert dry.status == "WARN"
    assert dry.m3_docking_allowed is False
    assert dry.counts["ligands_missing"] == 3

    dock_ctx = make_ctx(tmp_path / "dock")
    write_m1_inputs(dock_ctx.fresh_root)
    write_m2_package(dock_ctx.fresh_root)
    dock = run_m3_readiness(dock_ctx, "m2_run", mode="docking")
    assert dock.status == "FAIL"
    assert dock.counts["ligands_missing"] == 3


def test_m3_readiness_invalid_box_size_fails(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m1_inputs(ctx.fresh_root)
    write_ligands(ctx.fresh_root)
    write_m2_package(ctx.fresh_root, box_size="0")

    result = run_m3_readiness(ctx, "m2_run", mode="docking")

    assert result.status == "FAIL"
    assert any("malformed or missing accepted pocket boxes" in item for item in result.blockers)


def test_m3_readiness_3gt8_raw_only_rows_do_not_allow_docking(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m1_inputs(ctx.fresh_root)
    write_ligands(ctx.fresh_root)
    reference_row = {
        "pocket_family_id": "fam_ref",
        "export_tier": "reference_only",
        "gate_class": "reference_only",
        "representative_state_id": "3GT8_raw",
        "representative_protomer_id": "A",
        "primary_state_count": "0",
        "reference_state_count": "1",
        "hard_gate_pass": "True",
        "non_atp_pass": "True",
        "lower_lateral_pass": "True",
        "dimer_accessibility_pass": "True",
        "reference_only": "True",
    }
    write_m2_package(ctx.fresh_root, accepted_rows=[reference_row])

    result = run_m3_readiness(ctx, "m2_run", mode="docking")

    assert result.status == "FAIL"
    assert result.counts["reference_only_rows"] == 1
    assert any("no accepted primary or secondary" in item for item in result.blockers)


def test_m3_readiness_m2_review_no_go_blocks(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m1_inputs(ctx.fresh_root)
    write_ligands(ctx.fresh_root)
    write_m2_package(ctx.fresh_root, m2_review_go=False)

    result = run_m3_readiness(ctx, "m2_run", mode="docking")

    assert result.status == "FAIL"
    assert any("M2 review blocks M3 readiness" in item for item in result.blockers)


def test_m3_readiness_missing_m2_review_blocks_docking_mode(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m1_inputs(ctx.fresh_root)
    write_ligands(ctx.fresh_root)
    m2_dir = write_m2_package(ctx.fresh_root)
    (m2_dir / "reviews" / "m2_review_status.json").unlink()

    result = run_m3_readiness(ctx, "m2_run", mode="docking")

    assert result.status == "FAIL"
    assert any("M2 review status file not found" in item for item in result.blockers)


def test_m3_readiness_uses_m2_final_status_and_run_local_receptor_inputs(tmp_path):
    ctx = make_ctx(tmp_path)
    write_ligands(ctx.fresh_root)
    m2_dir = write_m2_package(ctx.fresh_root)
    (m2_dir / "reviews" / "m2_review_status.json").unlink()
    write_m2_final_status_allows_m3(m2_dir)
    write_m2_run_local_receptor_inputs(m2_dir)
    manifest = m2_dir / "manifest"
    manifest.mkdir(parents=True, exist_ok=True)
    (manifest / "membrane_frame.json").write_text('{"frame_id":"run-local","normal":[0,0,1]}\n', encoding="utf-8")

    result = run_m3_readiness(ctx, "m2_run", mode="docking")
    summary = read_json(result.summary_json)

    assert result.status == "PASS"
    assert result.m3_docking_allowed is True
    assert not summary["blockers"]
    checks = read_csv(result.report_csv)
    assert any(row["check_id"] == "M2_REVIEW_STATUS" and row["status"] == "PASS" for row in checks)
    assert summary["inputs"]["receptor_mappings"]["EGFR_160-185"]["runtime"]["path"].endswith(
        "prepared/m2_1_ppi_inputs/EGFR_160-185/receptor/EGFR_160-185_runtime_offset_receptor_only.pdb"
    )


def test_m3_readiness_honors_m3_primary_allow_flag(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m1_inputs(ctx.fresh_root)
    write_ligands(ctx.fresh_root)
    blocked_row = {
        "pocket_family_id": "fam_primary",
        "export_tier": "primary",
        "gate_class": "accepted_primary_pocket",
        "representative_state_id": "EGFR_160-185",
        "representative_protomer_id": "A",
        "primary_state_count": "2",
        "reference_state_count": "0",
        "hard_gate_pass": "True",
        "non_atp_pass": "True",
        "lower_lateral_pass": "True",
        "dimer_accessibility_pass": "True",
        "m3_primary_allowed": "False",
        "m3_review_allowed": "False",
        "reference_only": "False",
    }
    write_m2_package(ctx.fresh_root, accepted_rows=[blocked_row])

    result = run_m3_readiness(ctx, "m2_run", mode="docking")

    assert result.status == "FAIL"
    assert result.counts["docking_ready_pocket_rows"] == 0
    assert any("no accepted primary or secondary" in item for item in result.blockers)


def test_m3_readiness_missing_accepted_state_receptor_pdb_blocks(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m1_inputs(ctx.fresh_root)
    write_ligands(ctx.fresh_root)
    write_m2_package(ctx.fresh_root)
    dockable = ctx.fresh_root / "data" / "normalized" / "receptors" / "EGFR_160-185" / "dockable_669_1014_explicit_AB.pdb"
    dockable.unlink()

    result = run_m3_readiness(ctx, "m2_run", mode="docking")

    assert result.status == "FAIL"
    assert any("dockable receptor PDB" in item for item in result.blockers)


def test_m3_readiness_outputs_and_skeleton_stay_inside_run_dir(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m1_inputs(ctx.fresh_root)
    write_m2_package(ctx.fresh_root)

    result = run_m3_readiness(ctx, "m2_run", mode="dry-run")

    for path in [result.report_csv, result.summary_json, result.report_md, result.path_resolution_csv, result.skeleton_qc_csv]:
        path.resolve().relative_to(ctx.run_dir.resolve())
        assert path.is_file()
    skeleton_rows = read_csv(result.skeleton_qc_csv)
    assert any(row["path"].endswith("phase3_compounds/docking_outputs/focused_pocket_first") for row in skeleton_rows)


def test_m3_readiness_non_goal_flags_remain_true(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m1_inputs(ctx.fresh_root)
    write_ligands(ctx.fresh_root)
    write_m2_package(ctx.fresh_root)

    result = run_m3_readiness(ctx, "m2_run", mode="docking")
    summary = read_json(result.summary_json)

    assert all(summary["non_goals_preserved"].values())
    assert not (ctx.run_dir / "phase4_scoring").exists()
    generated_docking_files = [
        path
        for path in (ctx.run_dir / "phase3_compounds" / "docking_outputs").rglob("*")
        if path.is_file() and path.name != ".gitkeep"
    ]
    assert generated_docking_files == []
