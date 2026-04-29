"""Tests for M1 Phase 6 — PBS generator (handoff §10).

Closes M1 §23 #6: qsub smoke can be generated. Actual `qsub` execution is
HPC-pending (user-side step). These tests verify file generation and
structural correctness.
"""

import json
import os
import re
import subprocess
import sys
import uuid
from pathlib import Path

import pytest

from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext
from egfr_myo1d.hpc.pbs import (
    DEFAULT_CONDA_ENV,
    DEFAULT_CONDA_SH,
    PBSFile,
    THREAD_ENV_KEYS,
    derive_smoke_env_command_lines,
    derive_smoke_input_command_lines,
    generate_pbs,
    load_hpc_gate,
    render_pbs_content,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
GOLDEN_DIR = REPO_ROOT / "fresh" / "tests" / "fixtures" / "m1_phase6_pbs"


def unique_run_id(prefix="pytest_m1_phase6"):
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


# ---------------------------------------------------------------------------
# Pure render tests
# ---------------------------------------------------------------------------


def _render_smoke_env(run_id="RUNID", node="node04", ppn=4, walltime="02:00:00"):
    return render_pbs_content(
        job_name="smoke_env_{0}".format(run_id),
        node=node,
        ppn=ppn,
        walltime=walltime,
        queue="workq",
        repo_root="/test/repo",
        conda_sh=DEFAULT_CONDA_SH,
        conda_env=DEFAULT_CONDA_ENV,
        thread_limits={k: 1 for k in THREAD_ENV_KEYS},
        stdout_path="/test/repo/fresh/runs/{0}/logs/jobs/smoke_env_{0}.stdout".format(run_id),
        stderr_path="/test/repo/fresh/runs/{0}/logs/jobs/smoke_env_{0}.stderr".format(run_id),
        run_id=run_id,
        command_lines=derive_smoke_env_command_lines(run_id),
    )


def test_pbs_uses_concrete_stdout_stderr_absolute_paths():
    content = _render_smoke_env()
    # Absolute paths starting with / (Linux) inside #PBS -o / -e
    assert "#PBS -o /test/repo/fresh/runs/RUNID/logs/jobs/smoke_env_RUNID.stdout" in content
    assert "#PBS -e /test/repo/fresh/runs/RUNID/logs/jobs/smoke_env_RUNID.stderr" in content
    # Must NOT use $PBS_JOBNAME / $PBS_O_WORKDIR variable expansion in -o/-e
    assert "#PBS -o $" not in content
    assert "#PBS -e $" not in content


def test_pbs_exports_pythonpath_to_fresh_src():
    content = _render_smoke_env()
    assert 'export PYTHONPATH="$REPO_ROOT/fresh/src:${PYTHONPATH:-}"' in content


def test_pbs_sets_all_five_thread_limits_to_1():
    content = _render_smoke_env()
    for key in THREAD_ENV_KEYS:
        assert "export {0}=1".format(key) in content


def test_pbs_activates_pyrosetta_conda_env():
    content = _render_smoke_env()
    assert "source {0}".format(DEFAULT_CONDA_SH) in content
    assert "conda activate {0}".format(DEFAULT_CONDA_ENV) in content


def test_pbs_smoke_env_uses_node04_ppn4_walltime_02h():
    content = _render_smoke_env(node="node04", ppn=4, walltime="02:00:00")
    assert "#PBS -l nodes=node04:ppn=4" in content
    assert "#PBS -l walltime=02:00:00" in content


def test_pbs_smoke_env_command_body():
    content = _render_smoke_env(run_id="RUNID")
    assert "python -m egfr_myo1d.cli preflight --run-id RUNID --mode smoke_env --profile hpc_strict" in content
    assert "python -m egfr_myo1d.cli cleanup --run-id RUNID --mode test --dry-run false" in content


def test_pbs_smoke_input_uses_prepare_inputs_command():
    cmds = derive_smoke_input_command_lines("RUNID", input_root="fresh/data/raw")
    body = "\n".join(cmds)
    assert "prepare-inputs" in body
    assert "--run-id RUNID" in body
    assert "--mode smoke_input" in body
    assert "--input-root fresh/data/raw" in body


# ---------------------------------------------------------------------------
# generate_pbs end-to-end (writes file under run_dir)
# ---------------------------------------------------------------------------


def test_generate_pbs_writes_under_run_dir(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    pbs = generate_pbs(ctx, job_name="smoke_env_test", mode="smoke_env", node="node04")
    assert pbs.output_path.is_file()
    pbs.output_path.resolve().relative_to(ctx.run_dir.resolve())
    assert pbs.status == "PASS"
    # File content includes run_id
    content = pbs.output_path.read_text(encoding="utf-8")
    assert ctx.run_id in content


def test_generate_pbs_production_uses_ppn32(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    pbs = generate_pbs(ctx, job_name="prod_test", mode="production", node="node05")
    content = pbs.output_path.read_text(encoding="utf-8")
    assert "#PBS -l nodes=node05:ppn=32" in content
    assert "#PBS -l walltime=999:00:00" in content


def test_generate_pbs_mini_uses_ppn16(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    pbs = generate_pbs(ctx, job_name="mini_test", mode="mini", node="node04")
    content = pbs.output_path.read_text(encoding="utf-8")
    assert "ppn=16" in content
    assert "walltime=12:00:00" in content


def test_generate_pbs_unknown_node_fails(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    with pytest.raises(ValueError):
        generate_pbs(ctx, job_name="bad_node", mode="smoke_env", node="node99")


def test_generate_pbs_unknown_mode_fails(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    with pytest.raises(ValueError):
        generate_pbs(ctx, job_name="bad_mode", mode="nonsense_mode", node="node04")


def test_generate_pbs_invalid_walltime_fails(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    with pytest.raises(ValueError):
        generate_pbs(
            ctx,
            job_name="bad_walltime",
            mode="smoke_env",
            node="node04",
            walltime="2 hours",
        )


def test_generate_pbs_invalid_job_name_fails(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    with pytest.raises(ValueError):
        generate_pbs(ctx, job_name="bad name with spaces", mode="smoke_env", node="node04")


def test_generate_pbs_walltime_override_warns(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    pbs = generate_pbs(
        ctx,
        job_name="walltime_override",
        mode="smoke_env",
        node="node04",
        walltime="01:00:00",
    )
    assert pbs.status == "WARN"
    assert any("walltime_override" in w for w in pbs.warnings)


def test_generate_pbs_appends_job_status(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    generate_pbs(ctx, job_name="status_test", mode="smoke_env", node="node04")
    job_log = ctx.logs_dir / "job_status.jsonl"
    lines = [
        json.loads(line)
        for line in job_log.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    statuses = [r for r in lines if r.get("job_name") == "status_test"]
    assert len(statuses) == 1
    assert statuses[0]["status"] == "GENERATED"
    assert statuses[0]["node"] == "node04"
    assert statuses[0]["ppn"] == 4
    assert statuses[0]["details"]["mode"] == "smoke_env"


def test_generate_pbs_uses_unix_line_endings(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    pbs = generate_pbs(ctx, job_name="lineend_test", mode="smoke_env", node="node04")
    raw = pbs.output_path.read_bytes()
    assert b"\r\n" not in raw, "PBS file must use \\n line endings (HPC is Linux)"
    assert b"\n" in raw


def test_generate_pbs_path_traversal_rejected(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    from egfr_myo1d.core.run_context import RunContextError

    with pytest.raises((RunContextError, ValueError)):
        generate_pbs(
            ctx,
            job_name="traversal_test",
            mode="smoke_env",
            node="node04",
            output_path=Path("/abs/outside/path.pbs"),
        )


# ---------------------------------------------------------------------------
# qsub auto-call audit (must NOT auto-submit)
# ---------------------------------------------------------------------------


def test_no_auto_qsub_call_in_module_or_scripts():
    """The hpc/pbs module must not call qsub; scripts may print 'qsub' as
    instruction text but must not invoke it as a command."""
    pbs_src = (REPO_ROOT / "fresh" / "src" / "egfr_myo1d" / "hpc" / "pbs.py").read_text(
        encoding="utf-8"
    )
    cleaned = re.sub(r'"""[\s\S]*?"""', "", pbs_src)
    cleaned = re.sub(r"#.*", "", cleaned)
    assert "qsub" not in cleaned, "hpc/pbs.py must not invoke qsub"

    # Scripts may mention qsub in heredocs/comments, but no executable line
    # should run `qsub <args>` as a shell-level command. Track heredoc state
    # so heredoc body lines are excluded from the auto-call check.
    for script in ("submit_smoke_env.sh", "submit_smoke_input.sh"):
        text = (REPO_ROOT / "fresh" / "scripts" / script).read_text(encoding="utf-8")
        heredoc_terminator = None
        for line in text.splitlines():
            stripped = line.strip()
            if heredoc_terminator is not None:
                if stripped == heredoc_terminator:
                    heredoc_terminator = None
                continue  # heredoc body lines are literal text, not commands
            heredoc_match = re.search(r"<<-?\s*'?([A-Za-z_][A-Za-z0-9_]*)'?", line)
            if heredoc_match:
                heredoc_terminator = heredoc_match.group(1)
                continue
            if stripped.startswith("#"):
                continue
            if stripped.startswith("qsub "):
                raise AssertionError(
                    "{0} contains a top-level `qsub` invocation: {1}".format(
                        script, line
                    )
                )


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


def _cli_env():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src") + os.pathsep + env.get("PYTHONPATH", "")
    )
    return env


def test_cli_help_includes_prepare_pbs():
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
    assert "prepare-pbs" in out


def test_cli_prepare_pbs_subcommand_help():
    proc = subprocess.run(
        [sys.executable, "-m", "egfr_myo1d.cli", "prepare-pbs", "--help"],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode == 0
    out = proc.stdout.decode("utf-8", errors="replace")
    for flag in ("--run-id", "--job-name", "--mode", "--node", "--ppn", "--walltime"):
        assert flag in out


def test_cli_prepare_pbs_path_traversal_rejected():
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "egfr_myo1d.cli",
            "prepare-pbs",
            "--run-id",
            "../bad_run",
            "--job-name",
            "bad_traversal",
            "--mode",
            "smoke_env",
            "--node",
            "node04",
        ],
        cwd=str(REPO_ROOT),
        env=_cli_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode != 0


# ---------------------------------------------------------------------------
# hpc.yaml integration
# ---------------------------------------------------------------------------


def test_walltime_per_mode_matches_hpc_yaml(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)
    hpc = load_hpc_gate(ctx)
    walltime_map = hpc.get("walltime", {})
    # smoke_env/smoke_input both look up "smoke" in hpc.yaml
    cases = [
        ("smoke_env", "smoke"),
        ("smoke_input", "smoke"),
        ("mini", "mini"),
        ("scaling", "scaling"),
        ("production", "production"),
    ]
    for mode, lookup in cases:
        if lookup not in walltime_map:
            continue
        pbs = generate_pbs(
            ctx, job_name="wt_{0}".format(mode), mode=mode, node="node04"
        )
        content = pbs.output_path.read_text(encoding="utf-8")
        assert "walltime={0}".format(walltime_map[lookup]) in content
