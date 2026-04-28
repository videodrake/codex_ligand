import json
import os
import subprocess
import sys
import uuid
from pathlib import Path

import pytest

from egfr_myo1d.cli import main
from egfr_myo1d.core.logging_utils import append_phase_status, initialize_logs
from egfr_myo1d.core.manifest import initialize_manifests, write_input_manifest
from egfr_myo1d.core.run_context import RunContext, RunContextError
from egfr_myo1d.validation import preflight as preflight_module


REPO_ROOT = Path(__file__).resolve().parents[2]


def unique_run_id(prefix="pytest_task2"):
    return "{}_{}".format(prefix, uuid.uuid4().hex[:12])


def test_run_context_creates_standard_dirs():
    ctx = RunContext.create(unique_run_id(), "smoke_env", repo_root=REPO_ROOT)
    expected_dirs = [
        ctx.manifest_dir,
        ctx.logs_dir,
        ctx.jobs_log_dir,
        ctx.errors_dir,
        ctx.qc_dir,
        ctx.reports_dir,
        ctx.scratch_dir,
        ctx.tmp_dir,
    ]

    assert all(path.is_dir() for path in expected_dirs)
    for path in expected_dirs:
        assert path.resolve().is_relative_to(ctx.run_dir.resolve())


def test_run_id_rejects_path_traversal():
    for bad_run_id in ["../bad_run", "bad/run", "bad\\run", "bad..run"]:
        with pytest.raises(RunContextError):
            RunContext.create(bad_run_id, "smoke_env", repo_root=REPO_ROOT)


def test_run_manifest_written():
    ctx = RunContext.create(unique_run_id(), "smoke_env", repo_root=REPO_ROOT)
    initialize_logs(ctx)
    initialize_manifests(ctx, "smoke_env")

    manifest_path = ctx.manifest_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert manifest["run_id"] == ctx.run_id
    assert manifest["workflow_stage"] == "foundation"
    assert manifest["policy"]["old_workflow_modified"] is False
    assert manifest["policy"]["outputs_restricted_to_run_dir"] is True
    assert manifest["policy"]["docking_enabled"] is False


def test_input_manifest_marks_missing_files_without_failing_smoke_env():
    ctx = RunContext.create(unique_run_id(), "smoke_env", repo_root=REPO_ROOT)
    manifest = write_input_manifest(ctx, "smoke_env")

    assert manifest["mode"] == "smoke_env"
    assert manifest["missing_required_inputs"] == []
    assert {entry["input_id"] for entry in manifest["inputs"]} >= {
        "EGFR_160-185",
        "EGFR_170-200",
        "3GT8_raw",
        "plus10_full_frame",
        "AF-O94832-F1-model_v6",
        "Cpd-A",
        "Cpd-B",
        "Cpd-C",
    }


def test_logging_writes_master_and_phase_status_jsonl():
    ctx = RunContext.create(unique_run_id(), "smoke_env", repo_root=REPO_ROOT)
    initialize_logs(ctx)
    append_phase_status(ctx, "preflight", "WARN", "mock warning", {"x": 1})

    assert (ctx.logs_dir / "master.log").is_file()
    lines = (ctx.logs_dir / "phase_status.jsonl").read_text(encoding="utf-8").splitlines()
    assert lines
    record = json.loads(lines[-1])
    assert record["phase"] == "preflight"
    assert record["status"] == "WARN"
    assert record["run_id"] == ctx.run_id


def test_preflight_smoke_env_writes_environment_report_in_codex_dev(monkeypatch):
    ctx = RunContext.create(unique_run_id(), "smoke_env", repo_root=REPO_ROOT)
    initialize_logs(ctx)

    def fake_import(name, module, required_kind, profile):
        if name in {"numpy", "pandas"}:
            return {
                "category": "python_import",
                "name": name,
                "status": "PASS",
                "present": True,
                "details": {"module": module},
            }
        return {
            "category": "python_import",
            "name": name,
            "status": "WARN",
            "present": False,
            "details": {"module": module, "error": "mock missing"},
        }

    def fake_command(name, commands, required_kind, profile):
        return {
            "category": "external_command",
            "name": name,
            "status": "WARN",
            "present": False,
            "details": {"commands": commands, "error": "mock missing"},
        }

    monkeypatch.setattr(preflight_module, "check_import", fake_import)
    monkeypatch.setattr(preflight_module, "check_command", fake_command)

    report = preflight_module.run_preflight(ctx, "smoke_env", "codex_dev")

    assert report["status"] == "WARN"
    assert report["counts"]["FAIL"] == 0
    assert (ctx.manifest_dir / "environment_report.json").is_file()


def test_status_handles_existing_run(capsys):
    run_id = unique_run_id()
    ctx = RunContext.create(run_id, "smoke_env", repo_root=REPO_ROOT)
    initialize_logs(ctx)
    initialize_manifests(ctx, "smoke_env")
    append_phase_status(ctx, "init_run", "PASS", "ready", {})

    exit_code = main(["status", "--run-id", run_id])
    output = capsys.readouterr().out

    assert exit_code == 0
    assert "run_id={}".format(run_id) in output
    assert "manifest run_manifest.json: present" in output


def test_cli_help_runs():
    env = os.environ.copy()
    env["PYTHONPATH"] = "{}{}{}".format(
        REPO_ROOT / "fresh" / "src",
        os.pathsep,
        env.get("PYTHONPATH", ""),
    )
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
    assert "init-run" in proc.stdout
    assert "preflight" in proc.stdout
