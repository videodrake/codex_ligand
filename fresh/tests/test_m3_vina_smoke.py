import csv
import json
from pathlib import Path

from egfr_myo1d.compound.vina_adapter import VinaAvailability, VinaRunResult
from egfr_myo1d.compound.vina_smoke import run_m3_vina_smoke
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext


LIGAND_PDBQT = "ATOM      1  C   LIG A   1       0.000   0.000   0.000  0.00  0.00     0.000 C\nEND\n"
RECEPTOR_PDBQT = (
    "ATOM      1  CA  GLY A 900       1.000   2.000   3.000  1.00 20.00     0.000 C\n"
    "ATOM      2  CA  GLY B 900       4.000   5.000   6.000  1.00 20.00     0.000 C\nEND\n"
)
POSE_PDBQT = "MODEL 1\nATOM      1  C   LIG A   1       1.000   2.000   3.000  0.00  0.00     0.000 C\nENDMDL\n"


def make_ctx(tmp_path, run_id="m3_vina_smoke"):
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
        for row in rows:
            writer.writerow(row)


def read_csv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def sha(path):
    import hashlib

    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_inputs(ctx, state="EGFR_160-185", state_role="primary", t2_allowed=True, t3_allowed=True, ligand_status="PASS", box_updates=None, ligand_outside=False, receptor_outside=False):
    phase3 = ctx.run_dir / "phase3_compounds"
    manifest = phase3 / "manifests"
    qc = phase3 / "qc"
    ligand_dir = phase3 / "prepared_ligands" / "Cpd-A"
    receptor_dir = phase3 / "receptor_pdbqt" / state
    ligand_dir.mkdir(parents=True, exist_ok=True)
    receptor_dir.mkdir(parents=True, exist_ok=True)
    ligand = (ctx.repo_root / "outside_ligand.pdbqt") if ligand_outside else ligand_dir / "Cpd-A.pdbqt"
    receptor = (ctx.repo_root / "outside_receptor.pdbqt") if receptor_outside else receptor_dir / "receptor.pdbqt"
    ligand.write_text(LIGAND_PDBQT, encoding="utf-8")
    receptor.write_text(RECEPTOR_PDBQT, encoding="utf-8")
    qc.mkdir(parents=True, exist_ok=True)
    (qc / "m3_ligand_prep_summary.json").write_text(json.dumps({"overall_status": "PASS", "m3_t3_allowed": t2_allowed}), encoding="utf-8")
    (qc / "m3_receptor_box_summary.json").write_text(json.dumps({"overall_status": "PASS", "m3_t4_allowed": t3_allowed}), encoding="utf-8")
    write_csv(
        manifest / "ligand_preparation_manifest.csv",
        ["compound_public_id", "prepared_pdbqt_file", "prepared_pdbqt_sha256", "pdbqt_validation_success", "preparation_status"],
        [
            {"compound_public_id": "Cpd-B", "prepared_pdbqt_file": "", "prepared_pdbqt_sha256": "", "pdbqt_validation_success": "false", "preparation_status": "FAIL"},
            {"compound_public_id": "Cpd-A", "prepared_pdbqt_file": ctx.relative_to_repo(ligand), "prepared_pdbqt_sha256": sha(ligand), "pdbqt_validation_success": "true", "preparation_status": ligand_status},
        ],
    )
    write_csv(
        manifest / "receptor_preparation_manifest.csv",
        [
            "state_id",
            "state_role",
            "prepared_receptor_pdbqt_file",
            "prepared_receptor_pdbqt_sha256",
            "pdbqt_validation_success",
            "preparation_status",
            "allowed_for_vina_smoke",
            "source_receptor_is_old_workflow",
            "source_receptor_is_monomer_only",
        ],
        [
            {
                "state_id": state,
                "state_role": state_role,
                "prepared_receptor_pdbqt_file": ctx.relative_to_repo(receptor),
                "prepared_receptor_pdbqt_sha256": sha(receptor),
                "pdbqt_validation_success": "true",
                "preparation_status": "PASS",
                "allowed_for_vina_smoke": "true",
                "source_receptor_is_old_workflow": "false",
                "source_receptor_is_monomer_only": "false",
            }
        ],
    )
    box = {
        "pocket_family_id": "fam_primary",
        "state_id": state,
        "state_role": state_role,
        "protomer_id": "A",
        "box_id": "box_1",
        "box_center_x": "1.0",
        "box_center_y": "2.0",
        "box_center_z": "3.0",
        "box_size_x": "10.0",
        "box_size_y": "10.0",
        "box_size_z": "10.0",
        "non_atp_pass": "true",
        "lower_lateral_pass": "true",
        "dimer_accessibility_pass": "true",
        "receptor_pdbqt_file": ctx.relative_to_repo(receptor),
        "receptor_pdbqt_sha256": sha(receptor),
        "traceability_status": "PASS",
        "box_qc_status": "PASS",
        "allowed_for_vina_smoke": "true",
    }
    if box_updates:
        box.update(box_updates)
    write_csv(manifest / "docking_box_manifest.csv", list(box.keys()), [box])
    return ligand, receptor


def patch_vina_success(monkeypatch, ctx):
    monkeypatch.setattr(
        "egfr_myo1d.compound.vina_smoke.detect_vina",
        lambda exe="vina": VinaAvailability(executable=exe, available=True, version="AutoDock Vina test"),
    )

    def fake_run(argv, stdout_file, stderr_file, timeout_seconds, now_iso):
        assert (ctx.run_dir / "phase3_compounds" / "manifests" / "vina_smoke_manifest.csv").is_file()
        out = Path(argv[argv.index("--out") + 1])
        log = Path(argv[argv.index("--log") + 1])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(POSE_PDBQT, encoding="utf-8")
        log.write_text("mode | affinity | dist\n1 -7.0 0.0\n", encoding="utf-8")
        stdout_file.write_text("Vina completed\n", encoding="utf-8")
        stderr_file.write_text("", encoding="utf-8")
        return VinaRunResult(True, 0, False, now_iso(), now_iso(), 0.01)

    monkeypatch.setattr("egfr_myo1d.compound.vina_smoke.run_vina_once", fake_run)


def test_dry_run_writes_manifest_and_does_not_invoke_vina(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    monkeypatch.setattr("egfr_myo1d.compound.vina_smoke.detect_vina", lambda exe="vina": VinaAvailability(exe, False, None))

    result = run_m3_vina_smoke(ctx, "m2_run", mode="dry-run")
    rows = read_csv(result.manifest_csv)

    assert result.status == "WARN"
    assert result.vina["invoked"] is False
    assert rows[0]["command_manifest_written_before_execution"] == "true"
    assert rows[0]["vina_invocation_count"] == "0"


def test_run_mode_invokes_exactly_one_mocked_vina_and_validates_outputs(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    patch_vina_success(monkeypatch, ctx)

    result = run_m3_vina_smoke(ctx, "m2_run")
    rows = read_csv(result.manifest_csv)

    assert result.status == "PASS"
    assert result.m3_t5_allowed is True
    assert result.vina["invocation_count"] == 1
    assert rows[0]["output_pdbqt_exists"] == "true"
    assert rows[0]["minimal_pose_detected"] == "true"
    assert rows[0]["affinity_table_detected"] == "true"
    assert Path(ctx.repo_root / rows[0]["stdout_file"]).is_file()
    assert Path(ctx.repo_root / rows[0]["stderr_file"]).is_file()


def test_readiness_skeleton_gitkeep_does_not_count_as_forbidden_output(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    for rel in [
        "docking_inputs/production",
        "docking_outputs/focused_pocket_first/production",
        "docking_outputs/broad_anchor_scan_optional",
        "vina_raw/production",
    ]:
        directory = ctx.run_dir / "phase3_compounds" / rel
        directory.mkdir(parents=True, exist_ok=True)
        (directory / ".gitkeep").write_text("", encoding="utf-8")
    patch_vina_success(monkeypatch, ctx)

    result = run_m3_vina_smoke(ctx, "m2_run")

    assert result.status == "PASS"
    assert "production/broad/pose/candidate outputs created" not in result.blockers


def test_deterministic_and_explicit_selection(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    patch_vina_success(monkeypatch, ctx)

    result = run_m3_vina_smoke(ctx, "m2_run", compound_public_id="Cpd-A", state_id="EGFR_160-185", box_id="box_1")

    assert result.selection["compound_public_id"] == "Cpd-A"
    assert result.selection["state_id"] == "EGFR_160-185"
    assert result.selection["box_id"] == "box_1"


def test_reference_state_rejected_by_default_and_allowed_diagnostic(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx, state="3GT8_raw", state_role="reference")
    patch_vina_success(monkeypatch, ctx)

    rejected = run_m3_vina_smoke(ctx, "m2_run", mode="run")
    allowed = run_m3_vina_smoke(ctx, "m2_run", mode="run", force=True, allow_reference_smoke=True)

    assert rejected.status == "FAIL"
    assert allowed.status == "WARN"
    assert allowed.m3_t5_allowed is False
    assert allowed.selection["state_role"] == "reference"


def test_missing_required_manifests_fail_run_mode(tmp_path):
    ctx = make_ctx(tmp_path)

    result = run_m3_vina_smoke(ctx, "m2_run", mode="run")

    assert result.status == "FAIL"
    assert any("missing ligand_preparation_manifest" in item for item in result.blockers)
    assert any("missing receptor_preparation_manifest" in item for item in result.blockers)
    assert any("missing docking_box_manifest" in item for item in result.blockers)


def test_m3_t4_not_allowed_blocks_unless_allow_independent(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx, t3_allowed=False)
    patch_vina_success(monkeypatch, ctx)

    blocked = run_m3_vina_smoke(ctx, "m2_run")
    bypass = run_m3_vina_smoke(ctx, "m2_run", force=True, allow_independent=True)

    assert blocked.status == "FAIL"
    assert bypass.status == "WARN"
    assert bypass.m3_t5_allowed is False


def test_m3_t2_handoff_false_blocks_vina_smoke(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx, t2_allowed=False)
    patch_vina_success(monkeypatch, ctx)

    result = run_m3_vina_smoke(ctx, "m2_run")

    assert result.status == "FAIL"
    assert result.vina["invocation_count"] == 0
    assert any("m3_t3_allowed" in blocker for blocker in result.blockers)


def test_ligand_and_receptor_hash_drift_block_selection(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    ligand, _ = write_inputs(ctx)
    ligand.write_text(LIGAND_PDBQT.replace("0.000", "9.000", 1), encoding="utf-8")
    patch_vina_success(monkeypatch, ctx)

    result = run_m3_vina_smoke(ctx, "m2_run")

    assert result.status == "FAIL"
    assert result.vina["invocation_count"] == 0
    assert result.counts["eligible_ligands"] == 0

    ctx2 = make_ctx(tmp_path / "receptor_hash")
    _, receptor = write_inputs(ctx2)
    receptor.write_text(RECEPTOR_PDBQT.replace("1.000", "9.000", 1), encoding="utf-8")
    patch_vina_success(monkeypatch, ctx2)

    result2 = run_m3_vina_smoke(ctx2, "m2_run")

    assert result2.status == "FAIL"
    assert result2.vina["invocation_count"] == 0
    assert result2.counts["eligible_receptor_states"] == 0


def test_bad_box_values_fail(tmp_path, monkeypatch):
    for key, value in [("box_center_x", "nan"), ("box_size_x", "0"), ("non_atp_pass", "false"), ("lower_lateral_pass", "false"), ("dimer_accessibility_pass", "false"), ("traceability_status", "FAIL")]:
        ctx = make_ctx(tmp_path / key)
        write_inputs(ctx, box_updates={key: value})
        patch_vina_success(monkeypatch, ctx)
        result = run_m3_vina_smoke(ctx, "m2_run")
        assert result.status == "FAIL"


def test_paths_outside_run_dir_fail(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx, ligand_outside=True)
    patch_vina_success(monkeypatch, ctx)
    assert run_m3_vina_smoke(ctx, "m2_run").status == "FAIL"

    ctx2 = make_ctx(tmp_path / "receptor")
    write_inputs(ctx2, receptor_outside=True)
    patch_vina_success(monkeypatch, ctx2)
    assert run_m3_vina_smoke(ctx2, "m2_run").status == "FAIL"


def test_vina_missing_timeout_and_nonzero_are_recorded(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path / "missing")
    write_inputs(ctx)
    monkeypatch.setattr("egfr_myo1d.compound.vina_smoke.detect_vina", lambda exe="vina": VinaAvailability(exe, False, None))
    assert run_m3_vina_smoke(ctx, "m2_run").status == "FAIL"

    ctx2 = make_ctx(tmp_path / "timeout")
    write_inputs(ctx2)
    monkeypatch.setattr("egfr_myo1d.compound.vina_smoke.detect_vina", lambda exe="vina": VinaAvailability(exe, True, "test"))
    monkeypatch.setattr("egfr_myo1d.compound.vina_smoke.run_vina_once", lambda argv, stdout, stderr, timeout_seconds, now_iso: VinaRunResult(True, None, True, now_iso(), now_iso(), 1.0))
    timed = run_m3_vina_smoke(ctx2, "m2_run")
    assert timed.status == "FAIL"
    assert timed.vina["timeout"] is True

    ctx3 = make_ctx(tmp_path / "nonzero")
    write_inputs(ctx3)
    monkeypatch.setattr("egfr_myo1d.compound.vina_smoke.detect_vina", lambda exe="vina": VinaAvailability(exe, True, "test"))
    monkeypatch.setattr("egfr_myo1d.compound.vina_smoke.run_vina_once", lambda argv, stdout, stderr, timeout_seconds, now_iso: VinaRunResult(True, 2, False, now_iso(), now_iso(), 0.1))
    assert run_m3_vina_smoke(ctx3, "m2_run").status == "FAIL"


def test_invalid_vina_outputs_block_pass(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    monkeypatch.setattr(
        "egfr_myo1d.compound.vina_smoke.detect_vina",
        lambda exe="vina": VinaAvailability(executable=exe, available=True, version="AutoDock Vina test"),
    )

    def fake_garbage(argv, stdout_file, stderr_file, timeout_seconds, now_iso):
        out = Path(argv[argv.index("--out") + 1])
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("not a pdbqt pose\n", encoding="utf-8")
        stdout_file.write_text("Vina completed\n", encoding="utf-8")
        stderr_file.write_text("", encoding="utf-8")
        return VinaRunResult(True, 0, False, now_iso(), now_iso(), 0.01)

    monkeypatch.setattr("egfr_myo1d.compound.vina_smoke.run_vina_once", fake_garbage)

    result = run_m3_vina_smoke(ctx, "m2_run")

    assert result.status == "FAIL"
    assert result.m3_t5_allowed is False
    assert any("minimal pose" in blocker for blocker in result.blockers)
    assert any("Vina log missing" in blocker for blocker in result.blockers)


def test_no_forbidden_outputs_or_public_coordinate_leaks(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    patch_vina_success(monkeypatch, ctx)

    result = run_m3_vina_smoke(ctx, "m2_run")
    combined = "\n".join(
        path.read_text(encoding="utf-8", errors="ignore")
        for path in [result.manifest_csv, result.qc_csv, result.qc_json, result.report_md, result.phase3_log]
    )

    assert result.status == "PASS"
    assert "SECRET_A" not in combined
    assert "ATOM      1" not in combined
    for forbidden in ["qsub", "docking_inputs", "vina_raw", "tables"]:
        assert not (ctx.run_dir / "phase3_compounds" / forbidden).exists()
