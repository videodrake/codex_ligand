import builtins
import csv
import json
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.validation.ppi_sampling_plan import plan_ppi_sampling


REPO_ROOT = Path(__file__).resolve().parents[2]
TASK4_INPUTS = REPO_ROOT / "fresh" / "tests" / "fixtures" / "task4_inputs"


def unique_run_id(prefix="pytest_task6"):
    return "{}_{}".format(prefix, uuid.uuid4().hex[:12])


def make_tmp_run_context(tmp_path):
    run_id = unique_run_id()
    run_dir = tmp_path / run_id
    return RunContext(
        repo_root=REPO_ROOT,
        fresh_root=REPO_ROOT / "fresh",
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


def copy_task4_inputs(tmp_path):
    root = tmp_path / "task6_inputs"
    shutil.copytree(TASK4_INPUTS, root)
    return root


def load_contract(root):
    return json.loads((root / "ppi_input_contract.json").read_text(encoding="utf-8"))


def write_contract(root, contract):
    (root / "ppi_input_contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def run_plan(tmp_path, root, profile="codex_dev", mode="smoke_input", strict=False):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = plan_ppi_sampling(ctx, mode, profile, root, strict=strict)
    return ctx, report


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def read_jsonl(path):
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def read_csv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def test_cli_help_includes_plan_ppi_sampling():
    env = os.environ.copy()
    env["PYTHONPATH"] = "{}{}{}".format(REPO_ROOT / "fresh" / "src", os.pathsep, env.get("PYTHONPATH", ""))
    proc = subprocess.run(
        [sys.executable, "-m", "egfr_myo1d.cli", "--help"],
        cwd=str(REPO_ROOT),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    assert proc.returncode == 0
    assert "plan-ppi-sampling" in proc.stdout


def test_init_context_then_plan_succeeds_with_fixture_inputs(tmp_path):
    root = copy_task4_inputs(tmp_path)
    ctx, report = run_plan(tmp_path, root)

    assert report["status"] == "PASS_WITH_WARNINGS"
    assert (ctx.run_dir / "prepared" / "ppi" / "ppi_sampling_plan.json").is_file()


def test_outputs_are_written_only_under_run_dir(tmp_path):
    root = copy_task4_inputs(tmp_path)
    ctx, report = run_plan(tmp_path, root)

    assert report["status"] == "PASS_WITH_WARNINGS"
    required = [
        ctx.run_dir / "prepared" / "ppi" / "ppi_sampling_plan.json",
        ctx.run_dir / "prepared" / "ppi" / "ppi_job_specs.jsonl",
        ctx.run_dir / "prepared" / "ppi" / "ppi_job_specs.csv",
        ctx.run_dir / "prepared" / "ppi" / "pose_acceptance_policy.json",
        ctx.manifest_dir / "ppi_sampling_plan_manifest.json",
        ctx.manifest_dir / "ppi_sampling_plan_report.json",
        ctx.qc_dir / "ppi_sampling_plan_audit.csv",
        ctx.qc_dir / "ppi_pose_qc_policy_audit.csv",
        ctx.qc_dir / "ppi_sampling_blockers.csv",
    ]
    for path in required:
        assert path.is_file()
        assert path.resolve().is_relative_to(ctx.run_dir.resolve())


def test_job_specs_contain_primary_and_comparator_roles(tmp_path):
    root = copy_task4_inputs(tmp_path)
    ctx, _report = run_plan(tmp_path, root)
    jobs = read_jsonl(ctx.run_dir / "prepared" / "ppi" / "ppi_job_specs.jsonl")

    primary = [job for job in jobs if job["partner_construct_id"] == "MYO1D_sheet8_9_12_core_955_1001"]
    comparator = [job for job in jobs if job["partner_construct_id"] == "MYO1D_ext_beta_meander_955_1006_tail_masked"]
    assert primary
    assert comparator
    assert all(job["partner_role"] == "production_primary" for job in primary)
    assert all(job["partner_role"] == "noise_monitor" for job in comparator)
    assert not any(job["partner_role"] == "production_primary" for job in comparator)


def test_962_1006_artifact_gets_no_jobs_and_is_quarantined(tmp_path):
    root = copy_task4_inputs(tmp_path)
    ctx, _report = run_plan(tmp_path, root)
    jobs = read_jsonl(ctx.run_dir / "prepared" / "ppi" / "ppi_job_specs.jsonl")
    blockers = read_csv(ctx.qc_dir / "ppi_sampling_blockers.csv")

    assert not any("962_1006" in job["partner_construct_id"] for job in jobs)
    assert any(row["input_id"] == "myo1d_962_1006_terminal_artifact_bad" for row in blockers)
    assert any(row["production_primary_allowed"] == "False" for row in blockers)


def test_active_face_sheet12_and_no_bonus_are_propagated(tmp_path):
    root = copy_task4_inputs(tmp_path)
    ctx, _report = run_plan(tmp_path, root)
    jobs = read_jsonl(ctx.run_dir / "prepared" / "ppi" / "ppi_job_specs.jsonl")
    policy = read_json(ctx.run_dir / "prepared" / "ppi" / "pose_acceptance_policy.json")

    assert jobs[0]["active_face_residues"] == [961, 962, 963, 964, 968, 969, 970, 971, 972]
    assert jobs[0]["sheet12_support_residues"] == [993, 994, 995, 996, 997]
    assert jobs[0]["score_bonus_allowed"] is False
    assert jobs[0]["key_residue_bonus_weight"] == 0.0
    assert policy["score_bonus_allowed"] is False
    assert policy["key_residue_bonus_weight"] == 0.0


def test_manifest_contains_sha256_entries_for_present_inputs(tmp_path):
    root = copy_task4_inputs(tmp_path)
    ctx, _report = run_plan(tmp_path, root)
    manifest = read_json(ctx.manifest_dir / "ppi_sampling_plan_manifest.json")

    assert manifest["source_contract_sha256"]
    assert manifest["job_count_total"] == 12
    for record in manifest["source_file_sha256_entries"].values():
        if record.get("present"):
            assert record["sha256"]


def test_hpc_strict_blocks_v924r_receptor_from_production_plan(tmp_path):
    root = copy_task4_inputs(tmp_path)
    contract = load_contract(root)
    contract["receptor"]["path"] = "egfr_v924r_warn.pdb"
    write_contract(root, contract)
    ctx, report = run_plan(tmp_path, root, profile="hpc_strict")
    jobs = read_jsonl(ctx.run_dir / "prepared" / "ppi" / "ppi_job_specs.jsonl")

    assert report["status"] == "FAIL"
    assert jobs == []
    assert any(item["code"] == "v924r_like_mutation" for item in report["blockers"])


def test_hpc_strict_blocks_ambiguous_single_chain_receptor(tmp_path):
    root = copy_task4_inputs(tmp_path)
    contract = load_contract(root)
    contract["receptor"]["path"] = "egfr_single_chain_X_ambiguous.pdb"
    write_contract(root, contract)
    ctx, report = run_plan(tmp_path, root, profile="hpc_strict")
    jobs = read_jsonl(ctx.run_dir / "prepared" / "ppi" / "ppi_job_specs.jsonl")

    assert report["status"] == "FAIL"
    assert jobs == []
    assert any(item["code"] == "single_chain_dimer_ambiguity" for item in report["blockers"])


def test_invalid_membrane_frame_fails_cleanly(tmp_path):
    root = copy_task4_inputs(tmp_path)
    bad_frame = {
        "receptor_id": "EGFR_synthetic_inactive_dimer_AB",
        "membrane_normal": [0.0, 0.0, 0.0],
        "membrane_plane_point": [0.0, 0.0, 0.0],
        "dimer_axis": [1.0, 0.0, 0.0],
        "protomer_a_centroid": [0.0, 0.0, 0.0],
        "protomer_b_centroid": [20.0, 0.0, 0.0],
        "coordinate_convention": "invalid fixture",
    }
    (root / "membrane_frame_bad.json").write_text(json.dumps(bad_frame), encoding="utf-8")
    contract = load_contract(root)
    contract["membrane_frame"]["path"] = "membrane_frame_bad.json"
    write_contract(root, contract)
    ctx, report = run_plan(tmp_path, root)

    assert report["status"] == "FAIL"
    assert (ctx.manifest_dir / "ppi_sampling_plan_report.json").is_file()
    assert any(item["code"] == "invalid_membrane_frame" for item in report["blockers"])


def test_path_traversal_run_id_rejected():
    env = os.environ.copy()
    env["PYTHONPATH"] = "{}{}{}".format(REPO_ROOT / "fresh" / "src", os.pathsep, env.get("PYTHONPATH", ""))
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "egfr_myo1d.cli",
            "plan-ppi-sampling",
            "--run-id",
            "../bad_run",
            "--mode",
            "smoke_input",
            "--profile",
            "codex_dev",
            "--input-root",
            "fresh/tests/fixtures/task4_inputs",
            "--contract",
            "fresh/tests/fixtures/task4_inputs/ppi_input_contract.json",
        ],
        cwd=str(REPO_ROOT),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    assert proc.returncode != 0
    assert not (REPO_ROOT / "bad_run").exists()


def test_pose_policy_contains_required_future_classes(tmp_path):
    root = copy_task4_inputs(tmp_path)
    ctx, _report = run_plan(tmp_path, root)
    policy = read_json(ctx.run_dir / "prepared" / "ppi" / "pose_acceptance_policy.json")

    required = {
        "accepted_active_8_9_supported_12",
        "accepted_active_8_9_only",
        "review_sheet12_dominant",
        "reject_tail_dominant_artifact",
        "reject_flipped_or_back_face",
        "reject_membrane_proximal",
        "reject_atp_site_overlap",
        "reject_chain_or_residue_reset",
        "reject_quarantined_input",
    }
    assert required <= set(policy["future_pose_classes"])


def test_old_workflow_protected_paths_have_no_diff():
    proc = subprocess.run(
        [
            "git",
            "diff",
            "--name-only",
            "--",
            "run_production.py",
            "main.py",
            "egfr_pipeline/",
            "config/",
            "docs/runbook.md",
            "output/",
            "results_export/",
        ],
        cwd=str(REPO_ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    assert proc.returncode == 0
    assert proc.stdout.strip() == ""


def test_command_does_not_import_or_execute_heavy_tools(tmp_path, monkeypatch):
    root = copy_task4_inputs(tmp_path)
    blocked = {"pyrosetta", "lightdock", "vina", "fpocket", "p2rank", "prank"}
    original_import = builtins.__import__

    def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name.split(".")[0].lower() in blocked:
            raise AssertionError("Task 6 must not import {}".format(name))
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    _ctx, report = run_plan(tmp_path, root)

    assert report["status"] == "PASS_WITH_WARNINGS"
