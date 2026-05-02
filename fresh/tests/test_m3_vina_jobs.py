import csv
import json
import sys
from pathlib import Path

from egfr_myo1d.compound.vina_adapter import VinaAvailability
from egfr_myo1d.compound.vina_job_runner import run_vina_chunk
from egfr_myo1d.compound.vina_jobs import run_m3_vina_jobs
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext


LIGAND_PDBQT = "ATOM      1  C   LIG A   1       0.000   0.000   0.000  0.00  0.00     0.000 C\nEND\n"
RECEPTOR_PDBQT = (
    "ATOM      1  CA  GLY A 900       1.000   2.000   3.000  1.00 20.00     0.000 C\n"
    "ATOM      2  CA  GLY B 900       4.000   5.000   6.000  1.00 20.00     0.000 C\nEND\n"
)


def make_ctx(tmp_path, run_id="m3_vina_jobs"):
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


def write_inputs(ctx, states=("EGFR_160-185",), ligands=("Cpd-A",), boxes_per_state=1, include_reference=False, box_updates=None, ligand_outside=False, receptor_outside=False, smoke_allowed=True):
    phase3 = ctx.run_dir / "phase3_compounds"
    manifests = phase3 / "manifests"
    qc = phase3 / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    manifests.mkdir(parents=True, exist_ok=True)
    (qc / "m3_ligand_prep_summary.json").write_text(json.dumps({"overall_status": "PASS", "m3_t3_allowed": True}), encoding="utf-8")
    (qc / "m3_receptor_box_summary.json").write_text(json.dumps({"overall_status": "PASS", "m3_t4_allowed": True}), encoding="utf-8")
    (qc / "vina_smoke_qc.json").write_text(json.dumps({"overall_status": "PASS", "m3_t5_allowed": smoke_allowed}), encoding="utf-8")
    write_csv(manifests / "vina_smoke_manifest.csv", ["smoke_id", "smoke_status"], [{"smoke_id": "smoke_test", "smoke_status": "PASS"}])

    ligand_rows = []
    for cid in ligands:
        ligand_dir = phase3 / "prepared_ligands" / cid
        ligand_dir.mkdir(parents=True, exist_ok=True)
        ligand = (ctx.repo_root / f"outside_{cid}.pdbqt") if ligand_outside else ligand_dir / f"{cid}.pdbqt"
        ligand.write_text(LIGAND_PDBQT, encoding="utf-8")
        ligand_rows.append(
            {
                "compound_public_id": cid,
                "prepared_pdbqt_file": ctx.relative_to_repo(ligand),
                "prepared_pdbqt_sha256": sha(ligand),
                "pdbqt_validation_success": "true",
                "preparation_status": "PASS",
            }
        )
    write_csv(manifests / "ligand_preparation_manifest.csv", ["compound_public_id", "prepared_pdbqt_file", "prepared_pdbqt_sha256", "pdbqt_validation_success", "preparation_status"], ligand_rows)

    receptor_rows = []
    box_rows = []
    all_states = list(states)
    if include_reference:
        all_states.append("3GT8_raw")
    for state in all_states:
        role = "reference" if state == "3GT8_raw" else "primary"
        receptor_dir = phase3 / "receptor_pdbqt" / state
        receptor_dir.mkdir(parents=True, exist_ok=True)
        receptor = (ctx.repo_root / f"outside_{state}.pdbqt") if receptor_outside else receptor_dir / "receptor.pdbqt"
        receptor.write_text(RECEPTOR_PDBQT, encoding="utf-8")
        receptor_rows.append(
            {
                "state_id": state,
                "state_role": role,
                "prepared_receptor_pdbqt_file": ctx.relative_to_repo(receptor),
                "prepared_receptor_pdbqt_sha256": sha(receptor),
                "pdbqt_validation_success": "true",
                "preparation_status": "PASS",
                "allowed_for_vina_smoke": "true",
                "source_receptor_is_old_workflow": "false",
                "source_receptor_is_monomer_only": "false",
            }
        )
        for index in range(1, boxes_per_state + 1):
            family = f"fam_{index}"
            box = {
                "pocket_family_id": family,
                "state_id": state,
                "state_role": role,
                "protomer_id": "A",
                "box_id": f"box_{index}",
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
            box_rows.append(box)
    write_csv(
        manifests / "receptor_preparation_manifest.csv",
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
        receptor_rows,
    )
    write_csv(manifests / "docking_box_manifest.csv", list(box_rows[0].keys()), box_rows)


def patch_vina(monkeypatch):
    monkeypatch.setattr(
        "egfr_myo1d.compound.vina_jobs.detect_vina",
        lambda exe="vina": VinaAvailability(executable=exe, available=True, version="AutoDock Vina test"),
    )


def test_generate_writes_job_and_pbs_manifests_without_invoking_vina_or_qsub(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx, states=("EGFR_160-185", "EGFR_170-200"), ligands=("Cpd-A", "Cpd-B", "Cpd-C"), boxes_per_state=2)
    patch_vina(monkeypatch)

    result = run_m3_vina_jobs(ctx, m2_run_id="m2_run", profile="mini", mode="generate")
    jobs = read_csv(result.job_manifest_csv)
    pbs = read_csv(result.pbs_manifest_csv)

    assert result.status == "PASS"
    assert result.qsub_submission_allowed is True
    assert result.production_submission_allowed is False
    assert len(jobs) == 3 * 4 * 2
    assert {row["pbs_job_name"] for row in pbs} >= {"m3_vina_node04", "m3_vina_node05", "m3_vina_node06"}
    for row in pbs:
        if int(row["assigned_job_count"]) == 0:
            assert row["pbs_generation_status"] == "EMPTY_COMPATIBILITY"
            assert row["qsub_command"] == ""
    assert result.vina["generator_invoked_vina"] is False
    assert result.vina["generator_invoked_qsub"] is False


def test_dry_run_writes_manifest_and_warns_without_vina(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    monkeypatch.setattr("egfr_myo1d.compound.vina_jobs.detect_vina", lambda exe="vina": VinaAvailability(exe, False, None))

    result = run_m3_vina_jobs(ctx, m2_run_id="m2_run", mode="dry-run")

    assert result.status == "WARN"
    assert result.counts["planned_vina_commands"] == 2
    assert result.vina["generator_invoked_vina"] is False
    assert result.job_manifest_csv.is_file()


def test_job_enumeration_is_deterministic_and_job_ids_unique(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx, states=("EGFR_160-185",), ligands=("Cpd-B", "Cpd-A"), boxes_per_state=2)
    patch_vina(monkeypatch)

    first = run_m3_vina_jobs(ctx, m2_run_id="m2_run", mode="generate")
    first_ids = [row["job_id"] for row in read_csv(first.job_manifest_csv)]
    second = run_m3_vina_jobs(ctx, m2_run_id="m2_run", mode="generate", force=True)
    second_ids = [row["job_id"] for row in read_csv(second.job_manifest_csv)]

    assert first_ids == second_ids
    assert len(first_ids) == len(set(first_ids))
    assert first_ids[0].startswith("vina_mini_Cpd-A")


def test_explicit_filters_reduce_job_scope(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx, states=("EGFR_160-185", "EGFR_170-200"), ligands=("Cpd-A", "Cpd-B"), boxes_per_state=2)
    patch_vina(monkeypatch)

    result = run_m3_vina_jobs(ctx, m2_run_id="m2_run", compound_public_id="Cpd-B", state_id="EGFR_170-200", box_id="box_2")
    jobs = read_csv(result.job_manifest_csv)

    assert result.status == "PASS"
    assert len(jobs) == 2
    assert {row["compound_public_id"] for row in jobs} == {"Cpd-B"}
    assert {row["state_id"] for row in jobs} == {"EGFR_170-200"}
    assert {row["box_id"] for row in jobs} == {"box_2"}


def test_reference_state_excluded_by_default_and_allowed_only_with_flags(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx, states=(), ligands=("Cpd-A",), include_reference=True)
    patch_vina(monkeypatch)

    blocked = run_m3_vina_jobs(ctx, m2_run_id="m2_run")
    allowed = run_m3_vina_jobs(ctx, m2_run_id="m2_run", include_3gt8=True, allow_reference_jobs=True, force=True)

    assert blocked.status == "FAIL"
    assert allowed.status == "WARN"
    assert allowed.qsub_submission_allowed is False
    assert allowed.selection["eligible_states"] == ["3GT8_raw"]


def test_missing_inputs_and_m3_t5_gate_fail_generate(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    patch_vina(monkeypatch)
    missing = run_m3_vina_jobs(ctx, m2_run_id="m2_run")
    assert missing.status == "FAIL"
    assert any("ligand_preparation_manifest" in item for item in missing.blockers)

    ctx2 = make_ctx(tmp_path / "gate")
    write_inputs(ctx2, smoke_allowed=False)
    patch_vina(monkeypatch)
    blocked = run_m3_vina_jobs(ctx2, m2_run_id="m2_run")
    bypass = run_m3_vina_jobs(ctx2, m2_run_id="m2_run", force=True, allow_independent=True)
    assert blocked.status == "FAIL"
    assert bypass.status == "WARN"
    assert bypass.qsub_submission_allowed is False


def test_production_profile_requires_mini_completion_or_diagnostic_bypass(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    patch_vina(monkeypatch)

    blocked = run_m3_vina_jobs(ctx, m2_run_id="m2_run", profile="production")
    bypass = run_m3_vina_jobs(ctx, m2_run_id="m2_run", profile="production", force=True, allow_production_plan_without_mini=True)

    assert blocked.status == "FAIL"
    assert bypass.status == "WARN"
    assert bypass.production_submission_allowed is False


def test_bad_box_values_fail_selection(tmp_path, monkeypatch):
    for key, value in [
        ("box_center_x", "not_numeric"),
        ("box_size_x", "0"),
        ("non_atp_pass", "false"),
        ("lower_lateral_pass", "false"),
        ("dimer_accessibility_pass", "false"),
        ("traceability_status", "FAIL"),
    ]:
        ctx = make_ctx(tmp_path / key)
        write_inputs(ctx, box_updates={key: value})
        patch_vina(monkeypatch)
        result = run_m3_vina_jobs(ctx, m2_run_id="m2_run")
        assert result.status == "FAIL"
        assert result.counts["eligible_boxes"] == 0


def test_paths_outside_run_fail(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path / "lig")
    write_inputs(ctx, ligand_outside=True)
    patch_vina(monkeypatch)
    assert run_m3_vina_jobs(ctx, m2_run_id="m2_run").status == "FAIL"

    ctx2 = make_ctx(tmp_path / "rec")
    write_inputs(ctx2, receptor_outside=True)
    patch_vina(monkeypatch)
    assert run_m3_vina_jobs(ctx2, m2_run_id="m2_run").status == "FAIL"


def test_ligand_and_receptor_hash_drift_block_job_planning(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    ligand = ctx.run_dir / "phase3_compounds" / "prepared_ligands" / "Cpd-A" / "Cpd-A.pdbqt"
    ligand.write_text(LIGAND_PDBQT.replace("0.000", "9.000", 1), encoding="utf-8")
    patch_vina(monkeypatch)

    result = run_m3_vina_jobs(ctx, m2_run_id="m2_run")

    assert result.status == "FAIL"
    assert result.counts["eligible_ligands"] == 0

    ctx2 = make_ctx(tmp_path / "receptor_hash")
    write_inputs(ctx2)
    receptor = ctx2.run_dir / "phase3_compounds" / "receptor_pdbqt" / "EGFR_160-185" / "receptor.pdbqt"
    receptor.write_text(RECEPTOR_PDBQT.replace("1.000", "9.000", 1), encoding="utf-8")
    patch_vina(monkeypatch)

    result2 = run_m3_vina_jobs(ctx2, m2_run_id="m2_run")

    assert result2.status == "FAIL"
    assert result2.counts["eligible_receptor_states"] == 0


def test_pbs_scripts_are_concrete_and_hpc_safe(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx, states=("EGFR_160-185", "EGFR_170-200"), boxes_per_state=1)
    patch_vina(monkeypatch)

    result = run_m3_vina_jobs(ctx, m2_run_id="m2_run", profile="mini", mode="generate", ppn=16, cpu_per_vina=1)
    pbs_files = sorted((ctx.run_dir / "phase3_compounds" / "qsub").glob("*.pbs"))
    text = "\n".join(path.read_text(encoding="utf-8") for path in pbs_files)

    assert result.status == "PASS"
    assert "#PBS -o " in text and "$PBS_JOBID" not in text and "$RUN_ID" not in text
    assert 'PYTHONPATH="$REPO_ROOT/fresh/src:${PYTHONPATH:-}"' in text
    for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
        assert f"export {key}=1" in text
    assert "conda activate pyrosetta" in text
    assert max(result.hpc["worker_count_by_chunk"].values()) <= 16


def test_single_node_plan_marks_only_chunk_pbs_submit_ready(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx, states=("EGFR_160-185", "EGFR_170-200"), boxes_per_state=1)
    patch_vina(monkeypatch)

    result = run_m3_vina_jobs(
        ctx,
        m2_run_id="m2_run",
        profile="mini",
        mode="generate",
        nodes=["node04"],
        ppn=32,
        cpu_per_vina=4,
    )
    pbs = read_csv(result.pbs_manifest_csv)
    ready = [row for row in pbs if row["pbs_generation_status"] == "PASS"]
    wrapper = [row for row in pbs if row["pbs_job_name"] == "m3_vina_node04"]

    assert result.status == "PASS"
    assert {row["node"] for row in pbs} == {"node04"}
    assert len(ready) == result.counts["planned_chunks"]
    assert all("chunk" in row["pbs_job_name"] and row["qsub_command"] for row in ready)
    assert wrapper and wrapper[0]["pbs_generation_status"] == "EMPTY_COMPATIBILITY"
    assert wrapper[0]["assigned_job_count"] == "0"
    assert wrapper[0]["qsub_command"] == ""


def test_no_confidential_tokens_or_coordinates_in_public_outputs(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    patch_vina(monkeypatch)

    result = run_m3_vina_jobs(ctx, m2_run_id="m2_run")
    public_text = "\n".join(
        path.read_text(encoding="utf-8", errors="ignore")
        for path in [result.job_manifest_csv, result.pbs_manifest_csv, result.qc_json, result.qc_csv, result.report_md, result.phase3_log]
    )

    assert result.status == "PASS"
    assert "SECRET_A" not in public_text
    assert "SMILES" not in public_text
    assert "ATOM      1" not in public_text
    assert not (ctx.run_dir / "phase3_compounds" / "docking_outputs" / "broad_anchor_scan_optional").exists()
    for forbidden in ["compound_pose_raw.csv", "final_m3_candidate_hypotheses.csv"]:
        assert not (ctx.run_dir / "phase3_compounds" / "tables" / forbidden).exists()


def test_skeleton_gitkeep_outputs_do_not_block_mini_job_plan(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    for rel in [
        "phase3_compounds/docking_inputs/production",
        "phase3_compounds/docking_outputs/focused_pocket_first/production",
        "phase3_compounds/docking_outputs/broad_anchor_scan_optional",
        "phase3_compounds/vina_raw/production",
    ]:
        path = ctx.run_dir / rel
        path.mkdir(parents=True, exist_ok=True)
        (path / ".gitkeep").write_text("", encoding="utf-8")
    patch_vina(monkeypatch)

    result = run_m3_vina_jobs(ctx, m2_run_id="m2_run", profile="mini", mode="generate")

    assert result.status == "PASS"
    assert result.qsub_submission_allowed is True
    assert result.counts["production_outputs_created"] == 0
    assert result.counts["broad_outputs_created"] == 0
    assert not any("production output" in blocker or "broad docking" in blocker for blocker in result.blockers)


def test_previous_smiles_word_blocker_does_not_poison_rerun(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    qc = ctx.run_dir / "phase3_compounds" / "qc"
    reports = ctx.run_dir / "phase3_compounds" / "reports"
    qc.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    (qc / "vina_job_plan_qc.json").write_text(
        json.dumps({"overall_status": "FAIL", "blockers": ["ligand SMILES printed in public outputs/logs/manifests/PBS"]}),
        encoding="utf-8",
    )
    (reports / "m3_task5_vina_job_plan.md").write_text("- ligand SMILES printed in public outputs/logs/manifests/PBS\n", encoding="utf-8")
    patch_vina(monkeypatch)

    result = run_m3_vina_jobs(ctx, m2_run_id="m2_run", profile="mini", mode="generate", force=True)
    qc_summary = json.loads(result.qc_json.read_text(encoding="utf-8"))

    assert result.status == "PASS"
    assert qc_summary["confidentiality"]["smiles_logged"] is False
    assert not any("SMILES" in blocker for blocker in result.blockers)


def test_actual_smiles_value_in_public_output_blocks_job_plan(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    leak = ctx.run_dir / "phase3_compounds" / "reports" / "leaked_structure.txt"
    leak.parent.mkdir(parents=True, exist_ok=True)
    leak.write_text("SMILES: CCO\n", encoding="utf-8")
    patch_vina(monkeypatch)

    result = run_m3_vina_jobs(ctx, m2_run_id="m2_run", profile="mini", mode="generate")
    qc_summary = json.loads(result.qc_json.read_text(encoding="utf-8"))

    assert result.status == "FAIL"
    assert qc_summary["confidentiality"]["smiles_logged"] is True
    assert any("SMILES" in blocker for blocker in result.blockers)


def test_runner_dry_run_executes_only_requested_chunk_contract(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx, states=("EGFR_160-185", "EGFR_170-200"), boxes_per_state=1)
    patch_vina(monkeypatch)
    result = run_m3_vina_jobs(ctx, m2_run_id="m2_run", profile="mini", mode="generate")
    jobs = read_csv(result.job_manifest_csv)
    chunk_id = jobs[0]["chunk_id"]

    runner = run_vina_chunk(ctx, job_manifest=result.job_manifest_csv, chunk_id=chunk_id, max_workers=4, dry_run_runner=True)

    assert runner.status == "PASS"
    assert runner.rows_seen == len([row for row in jobs if row["chunk_id"] == chunk_id])
    assert runner.rows_run == 0


def test_runner_rejects_tampered_manifest_argv_and_disallowed_rows(tmp_path, monkeypatch):
    ctx = make_ctx(tmp_path)
    write_inputs(ctx)
    patch_vina(monkeypatch)
    result = run_m3_vina_jobs(ctx, m2_run_id="m2_run", profile="mini", mode="generate")
    jobs = read_csv(result.job_manifest_csv)
    row = dict(jobs[0])

    row["allowed_for_execution"] = "false"
    tampered = ctx.run_dir / "phase3_compounds" / "manifests" / "tampered_jobs.csv"
    write_csv(tampered, list(row.keys()), [row])
    try:
        run_vina_chunk(ctx, job_manifest=tampered, chunk_id=row["chunk_id"], dry_run_runner=True)
    except ValueError as exc:
        assert "allowed_for_execution" in str(exc)
    else:
        raise AssertionError("runner accepted disallowed row")

    row = dict(jobs[0])
    argv = json.loads(row["vina_argv_json"])
    argv[0] = "python"
    row["vina_argv_json"] = json.dumps(argv)
    tampered2 = ctx.run_dir / "phase3_compounds" / "manifests" / "tampered_jobs_argv.csv"
    write_csv(tampered2, list(row.keys()), [row])
    try:
        run_vina_chunk(ctx, job_manifest=tampered2, chunk_id=row["chunk_id"], dry_run_runner=True)
    except ValueError as exc:
        assert "executable" in str(exc)
    else:
        raise AssertionError("runner accepted tampered argv")


def test_safe_submit_helper_help_is_importable():
    script = Path(__file__).resolve().parents[1] / "scripts" / "submit_m3_compound_docking.py"
    assert script.is_file()
    assert "--i-understand-this-submits-hpc-jobs" in script.read_text(encoding="utf-8")
    assert "default=False" in script.read_text(encoding="utf-8")
