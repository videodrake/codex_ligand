import csv
import json
from pathlib import Path

from egfr_myo1d.compound.receptor_prepare import ReceptorToolStatus, run_m3_receptor_boxes
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext


PDB_TEXT = (
    "ATOM      1  CA  GLY A 900       1.000   2.000   3.000  1.00 20.00           C\n"
    "ATOM      2  CA  GLY B 900       4.000   5.000   6.000  1.00 20.00           C\n"
    "END\n"
)
PDBQT_TEXT = (
    "ATOM      1  CA  GLY A 900       1.000   2.000   3.000  1.00 20.00     0.000 C\n"
    "ATOM      2  CA  GLY B 900       4.000   5.000   6.000  1.00 20.00     0.000 C\n"
    "END\n"
)


def make_ctx(tmp_path, run_id="m3_t3"):
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
    (ctx.repo_root / ".gitignore").write_text("", encoding="utf-8")
    private = ctx.fresh_root / "data" / "private"
    private.mkdir(parents=True, exist_ok=True)
    (private / "compound_id_map.csv").write_text("public_id,internal_id,notes\nCpd-A,SECRET_A,test\n", encoding="utf-8")
    return ctx


def write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_m3_t2_pass(ctx, status="PASS", allowed=True):
    qc = ctx.run_dir / "phase3_compounds" / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    (qc / "m3_ligand_prep_summary.json").write_text(
        json.dumps({"overall_status": status, "m3_t3_allowed": allowed}),
        encoding="utf-8",
    )


def write_m2_fixture(ctx, state="EGFR_160-185", pdb_text=PDB_TEXT, suffix=".pdb", box_updates=None, pocket_updates=None):
    m2 = ctx.fresh_root / "runs" / "m2_run"
    receptor_dir = m2 / "phase1_receptors" / "normalized"
    receptor_dir.mkdir(parents=True, exist_ok=True)
    source = receptor_dir / f"{state}_dockable_669_1014_explicit_AB{suffix}"
    source.write_text(pdb_text, encoding="utf-8")
    mapping = receptor_dir / "receptor_mapping.csv"
    write_csv(
        mapping,
        ["state", "source_file", "protomer_id", "source_chain", "source_resseq", "runtime_chain", "runtime_resseq", "role"],
        [
            {"state": state, "source_file": str(source), "protomer_id": "A", "source_chain": "A", "source_resseq": "900", "runtime_chain": "A", "runtime_resseq": "900", "role": "primary"},
            {"state": state, "source_file": str(source), "protomer_id": "B", "source_chain": "B", "source_resseq": "900", "runtime_chain": "B", "runtime_resseq": "1900", "role": "primary"},
        ],
    )
    phase0 = m2 / "phase0_inputs"
    phase0.mkdir(parents=True, exist_ok=True)
    (phase0 / "membrane_frame.json").write_text('{"frame":"synthetic"}\n', encoding="utf-8")
    final = m2 / "phase2_pockets" / "final"
    pocket = {
        "pocket_family_id": "fam_primary",
        "pocket_acceptance_status": "accepted",
        "non_atp_pass": "true",
        "lower_lateral_pass": "true",
        "dimer_accessibility_pass": "true",
    }
    if pocket_updates:
        pocket.update(pocket_updates)
    write_csv(final / "accepted_pockets_for_m3.csv", list(pocket.keys()), [pocket])
    box = {
        "pocket_family_id": "fam_primary",
        "box_id": "box_fam_primary",
        "state_id": state,
        "state_role": "primary",
        "protomer_id": "A",
        "chain_id": "A",
        "receptor_pdb": ctx.relative_to_repo(source),
        "receptor_mapping_csv": ctx.relative_to_repo(mapping),
        "box_center_x": "2.0",
        "box_center_y": "3.0",
        "box_center_z": "4.0",
        "box_size_x": "10.0",
        "box_size_y": "11.0",
        "box_size_z": "12.0",
        "ppi_relationship_class": "PPI-adjacent",
        "geometry_class": "lower_lateral",
        "non_atp_pass": "true",
        "lower_lateral_pass": "true",
        "dimer_accessibility_pass": "true",
        "pocket_acceptance_status": "accepted",
    }
    if box_updates:
        box.update(box_updates)
    write_csv(final / "accepted_pocket_boxes.csv", list(box.keys()), [box])
    return source, mapping


def fake_tools():
    return ReceptorToolStatus(
        prepare_receptor4_available=True,
        prepare_receptor4_version="test prepare_receptor4",
        prepare_receptor4_binary="prepare_receptor4.py",
        mk_prepare_receptor_available=False,
        mk_prepare_receptor_version=None,
        mk_prepare_receptor_binary=None,
        meeko_available=False,
        meeko_version=None,
        openbabel_available=False,
        openbabel_version=None,
        rdkit_available=False,
        rdkit_version=None,
        vina_available=True,
        vina_invoked=False,
    )


def patch_fake_receptor_conversion(monkeypatch):
    monkeypatch.setattr("egfr_myo1d.compound.receptor_prepare.detect_receptor_tools", fake_tools)

    def fake_prepare(source, output, state_id, ctx, tools):
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(PDBQT_TEXT, encoding="utf-8")
        stdout = ctx.jobs_log_dir / f"m3_t3_{state_id}_receptor_prep.stdout.txt"
        stderr = ctx.jobs_log_dir / f"m3_t3_{state_id}_receptor_prep.stderr.txt"
        stdout.write_text("synthetic receptor prep completed\n", encoding="utf-8")
        stderr.write_text("", encoding="utf-8")
        return True, "prepare_receptor4", "prepare_receptor4", "test prepare_receptor4", stdout, stderr

    monkeypatch.setattr("egfr_myo1d.compound.receptor_prepare._copy_or_prepare_pdbqt", fake_prepare)


def test_valid_synthetic_dimer_receptor_is_prepared_and_boxes_pass(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    source, _ = write_m2_fixture(ctx)
    before = source.read_text(encoding="utf-8")
    patch_fake_receptor_conversion(monkeypatch)

    result = run_m3_receptor_boxes(ctx, "m2_run", box_margin=1.5)
    rows = read_csv(result.docking_box_manifest_csv)

    assert result.status == "PASS"
    assert result.m3_t4_allowed is True
    assert result.counts["states_prepared"] == 1
    assert result.counts["boxes_pass"] == 1
    assert rows[0]["box_size_x"] == "13"
    assert source.read_text(encoding="utf-8") == before
    assert (ctx.run_dir / "phase3_compounds" / "receptor_pdbqt" / "EGFR_160-185" / "receptor.pdbqt").is_file()


def test_m2_1_prepared_receptor_source_is_resolved(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    source, _ = write_m2_fixture(ctx)
    source.unlink()
    prepared_dir = ctx.fresh_root / "runs" / "m2_run" / "prepared" / "m2_1_ppi_inputs" / "EGFR_160-185" / "receptor"
    prepared_dir.mkdir(parents=True, exist_ok=True)
    prepared_source = prepared_dir / "EGFR_160-185_runtime_offset_receptor_only.pdb"
    prepared_source.write_text(PDB_TEXT, encoding="utf-8")
    patch_fake_receptor_conversion(monkeypatch)

    result = run_m3_receptor_boxes(ctx, "m2_run")
    rows = read_csv(result.receptor_preparation_manifest_csv)

    assert result.status == "PASS"
    assert rows[0]["source_receptor_file"].endswith(
        "prepared/m2_1_ppi_inputs/EGFR_160-185/receptor/EGFR_160-185_runtime_offset_receptor_only.pdb"
    )


def test_m2_1_prepared_receptor_wins_over_multiple_candidates(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    write_m2_fixture(ctx)
    prepared_dir = ctx.fresh_root / "runs" / "m2_run" / "prepared" / "m2_1_ppi_inputs" / "EGFR_160-185" / "receptor"
    prepared_dir.mkdir(parents=True, exist_ok=True)
    prepared_source = prepared_dir / "EGFR_160-185_runtime_offset_receptor_only.pdb"
    prepared_source.write_text(PDB_TEXT, encoding="utf-8")
    patch_fake_receptor_conversion(monkeypatch)

    result = run_m3_receptor_boxes(ctx, "m2_run")
    rows = read_csv(result.receptor_preparation_manifest_csv)

    assert result.status == "PASS"
    assert rows[0]["source_receptor_file"].endswith(
        "prepared/m2_1_ppi_inputs/EGFR_160-185/receptor/EGFR_160-185_runtime_offset_receptor_only.pdb"
    )


def test_empty_m3_skeleton_dirs_do_not_block_receptor_boxes(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    write_m2_fixture(ctx)
    for rel in [
        "docking_inputs",
        "docking_inputs/focused_pocket_first",
        "docking_outputs",
        "docking_outputs/focused_pocket_first",
    ]:
        (ctx.run_dir / "phase3_compounds" / rel).mkdir(parents=True, exist_ok=True)
    patch_fake_receptor_conversion(monkeypatch)

    result = run_m3_receptor_boxes(ctx, "m2_run")

    assert result.status == "PASS"
    assert not any("docking/Vina/scoring/candidate outputs" in blocker for blocker in result.blockers)


def test_supplied_receptor_pdbqt_is_validated_and_copied(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    write_m2_fixture(ctx, pdb_text=PDBQT_TEXT, suffix=".pdbqt", box_updates={"supplied_receptor_pdbqt_manual_review": "true"})
    monkeypatch.setattr("egfr_myo1d.compound.receptor_prepare.detect_receptor_tools", fake_tools)

    result = run_m3_receptor_boxes(ctx, "m2_run")
    rows = read_csv(result.receptor_preparation_manifest_csv)

    assert result.status == "WARN"
    assert rows[0]["preparation_method"] == "supplied_receptor_pdbqt_validated"
    assert "supplied receptor PDBQT" in "; ".join(result.warnings)


def test_supplied_receptor_pdbqt_requires_manual_review(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    write_m2_fixture(ctx, pdb_text=PDBQT_TEXT, suffix=".pdbqt")
    monkeypatch.setattr("egfr_myo1d.compound.receptor_prepare.detect_receptor_tools", fake_tools)

    result = run_m3_receptor_boxes(ctx, "m2_run")

    assert result.status == "FAIL"
    assert any("manual review" in blocker for blocker in result.blockers)


def test_monomer_only_receptor_source_fails(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    write_m2_fixture(ctx, pdb_text="ATOM      1  CA  GLY A 900       1.000   2.000   3.000  1.00 20.00           C\nEND\n")
    patch_fake_receptor_conversion(monkeypatch)

    result = run_m3_receptor_boxes(ctx, "m2_run")

    assert result.status == "FAIL"
    assert result.counts["monomer_only_sources_rejected"] == 1


def test_missing_mapping_and_membrane_frame_fail_prepare(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    _, mapping = write_m2_fixture(ctx)
    mapping.unlink()
    (ctx.fresh_root / "runs" / "m2_run" / "phase0_inputs" / "membrane_frame.json").unlink()
    patch_fake_receptor_conversion(monkeypatch)

    result = run_m3_receptor_boxes(ctx, "m2_run")

    assert result.status == "FAIL"
    assert any("receptor mapping" in item for item in result.blockers)
    assert any("membrane_frame" in item for item in result.blockers)


def test_missing_accepted_boxes_fails_prepare_and_warns_dry_run(tmp_path):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    result = run_m3_receptor_boxes(ctx, "m2_run", mode="prepare")
    dry = run_m3_receptor_boxes(ctx, "m2_run", mode="dry-run")

    assert result.status == "FAIL"
    assert dry.status == "WARN"
    assert dry.counts["boxes_input"] == 0


def test_bad_box_dimensions_and_missing_metadata_fail(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    write_m2_fixture(ctx, box_updates={"box_size_x": "0", "state_id": "", "protomer_id": "", "pocket_family_id": ""})
    patch_fake_receptor_conversion(monkeypatch)

    result = run_m3_receptor_boxes(ctx, "m2_run")

    assert result.status == "FAIL"
    assert result.counts["boxes_fail"] == 1


def test_box_not_traceable_or_missing_receptor_fails(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    write_m2_fixture(ctx, box_updates={"pocket_family_id": "fam_missing", "state_id": "EGFR_170-200"})
    patch_fake_receptor_conversion(monkeypatch)

    result = run_m3_receptor_boxes(ctx, "m2_run")

    assert result.status == "FAIL"
    assert result.counts["boxes_fail"] == 1
    assert any("traceable" in item or "validation failed" in item for item in result.blockers)


def test_non_m3_authorized_pocket_row_is_not_traceable(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    write_m2_fixture(
        ctx,
        pocket_updates={
            "pocket_acceptance_status": "",
            "export_tier": "primary",
            "m3_primary_allowed": "false",
            "m3_review_allowed": "false",
            "reference_only": "false",
            "hard_gate_pass": "true",
        },
    )
    patch_fake_receptor_conversion(monkeypatch)

    result = run_m3_receptor_boxes(ctx, "m2_run")
    rows = read_csv(result.docking_box_manifest_csv)

    assert result.status == "FAIL"
    assert rows[0]["traceability_status"] == "FAIL"
    assert rows[0]["source_pocket_status"] == "unknown"


def test_reference_only_box_is_rejected_from_vina_smoke(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    write_m2_fixture(
        ctx,
        box_updates={
            "export_tier": "reference_only",
            "reference_only": "true",
            "m3_primary_allowed": "false",
            "m3_review_allowed": "false",
        },
        pocket_updates={
            "export_tier": "reference_only",
            "gate_class": "reference_only",
            "reference_only": "true",
            "m3_primary_allowed": "false",
            "m3_review_allowed": "false",
        },
    )
    patch_fake_receptor_conversion(monkeypatch)

    result = run_m3_receptor_boxes(ctx, "m2_run")
    rows = read_csv(result.docking_box_manifest_csv)

    assert result.status == "FAIL"
    assert rows[0]["allowed_for_vina_smoke"] == "false"
    assert "authorized" in rows[0]["box_qc_notes"] or rows[0]["traceability_status"] == "FAIL"


def test_tool_unavailable_blocks_pdb_conversion_but_not_supplied_pdbqt(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    write_m2_fixture(ctx)

    unavailable = ReceptorToolStatus(False, None, None, False, None, None, False, None, False, None, False, None, False)
    monkeypatch.setattr("egfr_myo1d.compound.receptor_prepare.detect_receptor_tools", lambda: unavailable)
    result = run_m3_receptor_boxes(ctx, "m2_run")
    assert result.status == "FAIL"

    ctx2 = make_ctx(tmp_path / "second")
    write_m3_t2_pass(ctx2)
    write_m2_fixture(ctx2, pdb_text=PDBQT_TEXT, suffix=".pdbqt", box_updates={"supplied_receptor_pdbqt_manual_review": "true"})
    result2 = run_m3_receptor_boxes(ctx2, "m2_run")
    assert result2.status == "WARN"
    assert result2.counts["states_prepared"] == 1


def test_m3_t2_fail_blocks_m3_t4_and_allow_independent_records_bypass(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx, status="FAIL", allowed=False)
    write_m2_fixture(ctx)
    patch_fake_receptor_conversion(monkeypatch)

    blocked = run_m3_receptor_boxes(ctx, "m2_run")
    bypass = run_m3_receptor_boxes(ctx, "m2_run", force=True, allow_independent=True, require_m3_t2_pass=False)

    assert blocked.status == "FAIL"
    assert blocked.m3_t4_allowed is False
    assert bypass.status == "WARN"
    assert "allow_independent" in json.loads(bypass.summary_json.read_text(encoding="utf-8"))["warnings"][0]


def test_no_coordinates_or_forbidden_outputs_created_and_vina_not_invoked(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_m3_t2_pass(ctx)
    write_m2_fixture(ctx)
    patch_fake_receptor_conversion(monkeypatch)

    result = run_m3_receptor_boxes(ctx, "m2_run")
    combined = "\n".join(
        path.read_text(encoding="utf-8", errors="ignore")
        for path in [
            result.receptor_preparation_manifest_csv,
            result.docking_box_manifest_csv,
            result.receptor_prep_qc_csv,
            result.docking_box_qc_csv,
            result.summary_json,
            result.report_md,
            result.phase3_log,
        ]
    )

    assert result.tools.vina_available is True
    assert result.tools.vina_invoked is False
    assert "ATOM      1" not in combined
    for forbidden in ["docking_inputs", "docking_outputs", "vina_raw", "pose_attribution"]:
        assert not (ctx.run_dir / "phase3_compounds" / forbidden).exists()
