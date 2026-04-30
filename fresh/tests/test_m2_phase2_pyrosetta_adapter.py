import json
import os
import subprocess
import sys
import uuid
from pathlib import Path

import pytest

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.m2.ppi_inputs import generate_m2_1_ppi_inputs
from egfr_myo1d.orchestrator.prepare_inputs import run_prepare_inputs
from egfr_myo1d.ppi.pyrosetta_adapter import generate_pyrosetta_harness
from egfr_myo1d.ppi.run_ppi_job import dry_run_job


REPO_ROOT = Path(__file__).resolve().parents[2]
INTEGRATION_FIXTURE = (
    REPO_ROOT / "fresh" / "tests" / "fixtures" / "m1_phase8_integration"
)


def unique_run_id(prefix="pytest_m2_phase2"):
    return "{0}_{1}".format(prefix, uuid.uuid4().hex[:12])


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


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


@pytest.fixture
def ctx_with_m2_1(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    aggregate = run_prepare_inputs(
        ctx,
        mode="smoke_input",
        profile="codex_dev",
        input_root=INTEGRATION_FIXTURE,
    )
    assert aggregate.status in {"PASS", "PASS_WITH_WARNINGS"}
    m2_1 = generate_m2_1_ppi_inputs(ctx)
    assert m2_1.status in {"PASS", "PASS_WITH_WARNINGS"}
    return ctx


def test_m2_2_harness_generates_smoke_job_manifest(ctx_with_m2_1):
    report = generate_pyrosetta_harness(ctx_with_m2_1, mode="smoke_input")
    assert report.status in {"PASS", "PASS_WITH_WARNINGS"}
    assert len(report.jobs) == 1
    assert report.manifest_path.is_file()
    assert report.job_manifest_csv.is_file()
    assert report.job_manifest_jsonl.is_file()
    assert report.launch_script_path.is_file()

    payload = read_json(report.manifest_path)
    assert payload["execution_allowed"] is False
    assert payload["pyrosetta_imported"] is False
    assert payload["docking_executed"] is False
    assert payload["relaxation_executed"] is False
    assert payload["jobs"][0]["models_per_seed"] == 5
    assert payload["jobs"][0]["seed"] == 0
    stdout_path = payload["jobs"][0]["stdout_path"].replace("\\", "/")
    assert stdout_path.startswith("fresh/runs/") or "logs/jobs" in stdout_path


def test_m2_2_mini_harness_uses_two_states_two_seeds(ctx_with_m2_1):
    report = generate_pyrosetta_harness(ctx_with_m2_1, mode="mini")
    assert report.status in {"PASS", "PASS_WITH_WARNINGS"}
    assert len(report.jobs) == 4
    assert {job.seed for job in report.jobs} == {0, 1}
    assert {job.models_per_seed for job in report.jobs} == {20}
    assert {job.state_id for job in report.jobs} == {"EGFR_160-185", "EGFR_170-200"}


def test_m2_2_dry_run_job_writes_status_without_pyrosetta_execution(ctx_with_m2_1):
    report = generate_pyrosetta_harness(ctx_with_m2_1, mode="smoke_input")
    job_name = report.jobs[0].job_name
    status_path = dry_run_job(ctx_with_m2_1, job_name)
    assert status_path.is_file()
    payload = read_json(status_path)
    assert payload["status"] == "DRY_RUN_PASS"
    assert payload["pyrosetta_imported"] is False
    assert payload["docking_executed"] is False
    assert payload["relaxation_executed"] is False
    status_path.resolve().relative_to(ctx_with_m2_1.run_dir.resolve())


def test_m2_2_missing_m2_1_manifest_fails(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    report = generate_pyrosetta_harness(ctx, mode="smoke_input")
    assert report.status == "FAIL"
    assert report.blockers
    payload = read_json(report.manifest_path)
    assert payload["docking_executed"] is False


def test_m2_2_refuses_real_execution_cli(ctx_with_m2_1):
    report = generate_pyrosetta_harness(ctx_with_m2_1, mode="smoke_input")
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "egfr_myo1d.ppi.run_ppi_job",
            "--run-id",
            ctx_with_m2_1.run_id,
            "--job-name",
            report.jobs[0].job_name,
            "--dry-run",
            "false",
        ],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode != 0
    assert b"refuses real execution" in proc.stderr


def test_m2_2_adapter_source_does_not_import_or_execute_pyrosetta():
    adapter_source = (
        REPO_ROOT / "fresh/src/egfr_myo1d/ppi/pyrosetta_adapter.py"
    ).read_text(encoding="utf-8")
    runner_source = (
        REPO_ROOT / "fresh/src/egfr_myo1d/ppi/run_ppi_job.py"
    ).read_text(encoding="utf-8")
    for source in [adapter_source, runner_source]:
        assert "import pyrosetta" not in source
        assert "os.system" not in source
        assert "qsub " not in source
    # The runner may use subprocess only for Python 3 bootstrap, mirroring cli.py.
    assert "subprocess" not in adapter_source


def _cli_env():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    )
    return env


def test_cli_help_includes_prepare_m2_pyrosetta_harness():
    proc = subprocess.run(
        [sys.executable, "-m", "egfr_myo1d.cli", "--help"],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode == 0
    out = proc.stdout.decode("utf-8", errors="replace") + proc.stderr.decode(
        "utf-8", errors="replace"
    )
    assert "prepare-m2-pyrosetta-harness" in out
