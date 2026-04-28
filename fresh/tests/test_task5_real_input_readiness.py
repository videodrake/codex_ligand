import json
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.validation.real_inputs import validate_real_inputs


REPO_ROOT = Path(__file__).resolve().parents[2]
TASK4_INPUTS = REPO_ROOT / "fresh" / "tests" / "fixtures" / "task4_inputs"


def unique_run_id(prefix="pytest_task5"):
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
    root = tmp_path / "real_inputs"
    shutil.copytree(TASK4_INPUTS, root)
    return root


def load_contract(root):
    return json.loads((root / "ppi_input_contract.json").read_text(encoding="utf-8"))


def write_contract(root, contract):
    (root / "ppi_input_contract.json").write_text(
        json.dumps(contract, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def run_readiness(tmp_path, root, profile="codex_dev", mode="smoke_input", strict=False):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = validate_real_inputs(ctx, mode, profile, root, strict=strict)
    return ctx, report


def test_cli_help_includes_validate_real_inputs():
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
    assert "validate-real-inputs" in proc.stdout


def test_missing_real_contract_fails_cleanly(tmp_path):
    root = tmp_path / "empty_real_inputs"
    root.mkdir()
    ctx, report = run_readiness(tmp_path, root, profile="hpc_strict")

    assert report["status"] == "FAIL"
    assert (ctx.manifest_dir / "real_input_readiness_report.json").is_file()
    assert (ctx.qc_dir / "real_input_blockers.csv").is_file()
    assert report["qc_summary"]["production_blocker"] == 1


def test_valid_fixture_derived_real_style_contract_warns_in_codex_dev(tmp_path):
    root = copy_task4_inputs(tmp_path)
    ctx, report = run_readiness(tmp_path, root, profile="codex_dev", mode="smoke_input")

    assert report["status"] == "WARN"
    assert report["qc_summary"]["FAIL"] == 0
    assert report["egfr"]["dimer_context_preserved"] is True
    assert (ctx.manifest_dir / "real_input_manifest.json").is_file()


def test_hpc_strict_blocks_single_chain_ambiguous_receptor(tmp_path):
    root = copy_task4_inputs(tmp_path)
    contract = load_contract(root)
    contract["receptor"]["path"] = "egfr_single_chain_X_ambiguous.pdb"
    write_contract(root, contract)
    _ctx, report = run_readiness(tmp_path, root, profile="hpc_strict")

    assert report["status"] == "FAIL"
    assert any(item["code"] == "single_chain_dimer_ambiguity" for item in report["blockers"])
    assert report["qc_summary"]["quarantine"] >= 1


def test_hpc_strict_blocks_v924r_like_receptor(tmp_path):
    root = copy_task4_inputs(tmp_path)
    contract = load_contract(root)
    contract["receptor"]["path"] = "egfr_v924r_warn.pdb"
    write_contract(root, contract)
    _ctx, report = run_readiness(tmp_path, root, profile="hpc_strict")

    assert report["status"] == "FAIL"
    assert any(item["code"] == "v924r_like_mutation" for item in report["blockers"])
    assert "mutations_repaired" in report["egfr"]
    assert report["egfr"]["mutations_repaired"] is False


def test_hpc_strict_blocks_myo1d_962_start_as_primary(tmp_path):
    root = copy_task4_inputs(tmp_path)
    contract = load_contract(root)
    contract["myo1d"]["primary_construct"]["path"] = "myo1d_962_1006_terminal_artifact_bad.pdb"
    contract["myo1d"]["primary_construct"]["biological_residue_range"] = [955, 1001]
    write_contract(root, contract)
    _ctx, report = run_readiness(tmp_path, root, profile="hpc_strict")

    assert report["status"] == "FAIL"
    assert any(item["code"] == "myo1d_terminal_artifact" for item in report["blockers"])


def test_valid_myo1d_955_1001_primary_readiness_is_pass(tmp_path):
    root = copy_task4_inputs(tmp_path)
    _ctx, report = run_readiness(tmp_path, root, profile="codex_dev")

    assert report["myo1d"]["primary"]["status"] == "PASS"
    assert report["myo1d"]["primary"]["residue_range"] == [955, 1001]


def test_955_1006_comparator_remains_noise_monitor_only(tmp_path):
    root = copy_task4_inputs(tmp_path)
    _ctx, report = run_readiness(tmp_path, root, profile="codex_dev")

    assert report["myo1d"]["comparator"]["residue_range"] == [955, 1006]
    assert any(item["code"] == "myo1d_955_1006_comparator_only" for item in report["messages"])


def test_active_face_annotation_only_and_no_score_bonus(tmp_path):
    root = copy_task4_inputs(tmp_path)
    _ctx, report = run_readiness(tmp_path, root, profile="codex_dev")

    assert report["active_face_contract"]["score_bonus_allowed"] is False
    assert report["myo1d"]["score_bonus_allowed"] is False
    assert {961, 962, 963, 964, 968, 969, 970, 971, 972} <= set(report["active_face_contract"]["primary_active_face_residues"])


def test_output_paths_are_constrained_under_run_dir(tmp_path):
    root = copy_task4_inputs(tmp_path)
    ctx, _report = run_readiness(tmp_path, root, profile="codex_dev")
    required = [
        ctx.manifest_dir / "real_input_manifest.json",
        ctx.manifest_dir / "real_input_readiness_report.json",
        ctx.qc_dir / "real_egfr_readiness_audit.csv",
        ctx.qc_dir / "real_myo1d_readiness_audit.csv",
        ctx.qc_dir / "real_membrane_frame_audit.csv",
        ctx.qc_dir / "real_input_blockers.csv",
    ]
    for path in required:
        assert path.is_file()
        assert path.resolve().is_relative_to(ctx.run_dir.resolve())


def test_manifest_contains_sha256_for_all_present_inputs(tmp_path):
    root = copy_task4_inputs(tmp_path)
    ctx, _report = run_readiness(tmp_path, root, profile="codex_dev")
    manifest = json.loads((ctx.manifest_dir / "real_input_manifest.json").read_text(encoding="utf-8"))

    for record in manifest["inputs"].values():
        if record.get("present"):
            assert record["sha256"]


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

