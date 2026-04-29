"""Tests for M1 Phase 1 — cleanup manager (handoff §18).

Closes M1 §23 #8: cleanup is safe and writes cleanup_report.json.
"""

import json
import os
import subprocess
import sys
import uuid
from pathlib import Path

import pytest

from egfr_myo1d.core.cleanup import (
    CleanupReport,
    PRESERVED_DIR_BASENAMES,
    resolve_dry_run_default,
    run_cleanup,
)
from egfr_myo1d.core.logging_utils import initialize_logs
from egfr_myo1d.core.run_context import RunContext, RunContextError


REPO_ROOT = Path(__file__).resolve().parents[2]


def unique_run_id(prefix="pytest_m1_phase1"):
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


def seed_run(ctx):
    """Populate scratch/, tmp/, and a top-level .silent file in run_dir."""
    initialize_logs(ctx)
    # Scratch contents (deletable)
    (ctx.scratch_dir / "pose_001.pdb").write_text("ATOM 1\n")
    (ctx.scratch_dir / "pose_002.pdb").write_text("ATOM 2\n")
    nested = ctx.scratch_dir / "subdir"
    nested.mkdir(parents=True, exist_ok=True)
    (nested / "deeper.tmp").write_text("scratch")
    # Tmp contents (deletable)
    (ctx.tmp_dir / "vina_run.pdbqt.tmp").write_text("VINA")
    (ctx.tmp_dir / "rosetta.silent").write_text("ROSETTA")
    # Top-level silent / tmp at run_dir root (deletable by pattern)
    (ctx.run_dir / "abandoned.silent").write_text("DOCKING_TEMPORARY")
    (ctx.run_dir / "leftover.vina.tmp").write_text("LEFTOVER")
    # Preserved manifest content
    (ctx.manifest_dir / "run_manifest.json").write_text("{}")
    (ctx.qc_dir / "audit.csv").write_text("col1,col2\n")
    (ctx.reports_dir / "summary.md").write_text("# summary\n")


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def assert_under_run_dir(ctx, path):
    try:
        Path(path).resolve().relative_to(ctx.run_dir.resolve())
    except ValueError:
        raise AssertionError("{0} escaped run dir".format(path))


# ---------------------------------------------------------------------------
# Module-level helper tests
# ---------------------------------------------------------------------------


def test_resolve_dry_run_default_test_mode_is_false():
    assert resolve_dry_run_default("test") is False


def test_resolve_dry_run_default_production_mode_is_true():
    assert resolve_dry_run_default("production") is True


# ---------------------------------------------------------------------------
# Core behavior tests
# ---------------------------------------------------------------------------


def test_cleanup_dry_run_makes_no_changes(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    seed_run(ctx)

    report = run_cleanup(ctx, mode="test", dry_run=True, profile="codex_dev")

    assert isinstance(report, CleanupReport)
    assert report.dry_run is True
    assert report.status == "PASS"
    assert report.candidate_count == len(report.deleted_files)
    assert report.deleted_count == 0
    # Files still exist on disk
    assert (ctx.scratch_dir / "pose_001.pdb").exists()
    assert (ctx.tmp_dir / "vina_run.pdbqt.tmp").exists()
    assert (ctx.run_dir / "abandoned.silent").exists()
    # All deleted_files entries report deleted=False
    for record in report.deleted_files:
        assert record["deleted"] is False


def test_cleanup_test_mode_removes_scratch_and_tmp(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    seed_run(ctx)

    report = run_cleanup(ctx, mode="test", dry_run=False, profile="codex_dev")

    assert report.status == "PASS"
    assert report.dry_run is False
    assert report.deleted_count >= 5  # scratch x3 + tmp x2 + top-level x2 = 7 expected
    # Files actually gone
    assert not (ctx.scratch_dir / "pose_001.pdb").exists()
    assert not (ctx.scratch_dir / "subdir" / "deeper.tmp").exists()
    assert not (ctx.tmp_dir / "vina_run.pdbqt.tmp").exists()
    assert not (ctx.tmp_dir / "rosetta.silent").exists()
    assert not (ctx.run_dir / "abandoned.silent").exists()
    assert not (ctx.run_dir / "leftover.vina.tmp").exists()


def test_cleanup_preserves_manifest_logs_qc_reports(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    seed_run(ctx)

    run_cleanup(ctx, mode="test", dry_run=False, profile="codex_dev")

    # Preserved files still exist
    assert (ctx.manifest_dir / "run_manifest.json").exists()
    assert (ctx.qc_dir / "audit.csv").exists()
    assert (ctx.reports_dir / "summary.md").exists()
    # Logs directory and its subtree still exist
    assert (ctx.logs_dir / "master.log").exists()
    assert (ctx.logs_dir / "phase_status.jsonl").exists()
    assert (ctx.logs_dir / "job_status.jsonl").exists()
    assert (ctx.errors_dir / "error_summary.txt").exists()
    assert (ctx.errors_dir / "failed_jobs.csv").exists()
    # Cleanup report itself exists in manifest/
    assert (ctx.manifest_dir / "cleanup_report.json").exists()


def test_cleanup_production_default_is_dry_run(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    seed_run(ctx)

    # Do NOT pass dry_run; rely on production default (True)
    report = run_cleanup(ctx, mode="production", profile="codex_dev")

    assert report.dry_run is True
    assert report.candidate_count == len(report.deleted_files)
    assert report.deleted_count == 0
    # Files NOT actually deleted
    assert (ctx.scratch_dir / "pose_001.pdb").exists()
    assert (ctx.tmp_dir / "vina_run.pdbqt.tmp").exists()


def test_cleanup_production_explicit_dry_run_false_deletes(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    seed_run(ctx)

    report = run_cleanup(ctx, mode="production", dry_run=False, profile="codex_dev")

    assert report.dry_run is False
    assert report.deleted_count > 0
    assert not (ctx.scratch_dir / "pose_001.pdb").exists()


def test_cleanup_writes_cleanup_report_json(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    seed_run(ctx)

    report = run_cleanup(ctx, mode="test", dry_run=True, profile="codex_dev")

    cleanup_report = ctx.manifest_dir / "cleanup_report.json"
    assert cleanup_report.is_file()
    assert_under_run_dir(ctx, cleanup_report)
    payload = read_json(cleanup_report)
    assert payload["run_id"] == ctx.run_id
    assert payload["mode"] == "test"
    assert payload["dry_run"] is True
    assert payload["profile"] == "codex_dev"
    assert payload["status"] == "PASS"
    assert payload["deleted_count"] == report.deleted_count
    assert payload["preserved_count"] == report.preserved_count
    assert "deleted_files" in payload
    assert "preserved_files" in payload
    assert "errors" in payload
    assert "timestamp" in payload


def test_cleanup_appends_phase_status(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    seed_run(ctx)

    run_cleanup(ctx, mode="test", dry_run=False, profile="hpc_strict")

    phase_status_path = ctx.logs_dir / "phase_status.jsonl"
    lines = [
        json.loads(line)
        for line in phase_status_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    cleanup_records = [r for r in lines if r.get("phase") == "cleanup"]
    assert len(cleanup_records) == 1
    record = cleanup_records[0]
    assert record["status"] == "PASS"
    assert record["run_id"] == ctx.run_id
    assert "details" in record
    assert record["details"]["mode"] == "test"
    assert record["details"]["dry_run"] is False
    assert record["details"]["profile"] == "hpc_strict"


def test_cleanup_refuses_outside_run_dir(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)

    # Simulate the candidate enumerator producing an outside path by directly
    # invoking the safety check via _is_under: the public API never actually
    # generates such candidates, but the defensive check inside run_cleanup
    # raises if one slips. We exercise it by patching:
    from egfr_myo1d.core import cleanup as cleanup_module

    outside = tmp_path.parent / "outside_file.tmp"
    outside.write_text("OUTSIDE")

    real_enum = cleanup_module._enumerate_deletion_candidates

    def fake_enum(run_dir):
        return [outside]

    cleanup_module._enumerate_deletion_candidates = fake_enum
    try:
        with pytest.raises(RunContextError):
            run_cleanup(ctx, mode="test", dry_run=False, profile="codex_dev")
    finally:
        cleanup_module._enumerate_deletion_candidates = real_enum
        if outside.exists():
            outside.unlink()


def test_cleanup_invalid_mode_raises(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    initialize_logs(ctx)

    with pytest.raises(ValueError):
        run_cleanup(ctx, mode="bogus", dry_run=False, profile="codex_dev")


def test_cleanup_run_dir_missing_raises(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    # Do NOT call initialize_logs; run_dir does not exist
    with pytest.raises(RunContextError):
        run_cleanup(ctx, mode="test", dry_run=False, profile="codex_dev")


def test_cleanup_records_preserved_files_relative_paths(tmp_path):
    ctx = make_tmp_run_context(tmp_path)
    seed_run(ctx)

    report = run_cleanup(ctx, mode="test", dry_run=True, profile="codex_dev")

    paths = {entry["path"] for entry in report.preserved_files}
    # Forward-slash relative paths under run_dir
    assert any(p.startswith("manifest/") for p in paths)
    assert any(p.startswith("logs/") for p in paths)
    assert any(p.startswith("qc/") for p in paths)
    assert any(p.startswith("reports/") for p in paths)
    # No path escapes run_dir; no absolute paths in the records
    for entry in report.preserved_files:
        assert not entry["path"].startswith("/")
        assert not entry["path"].startswith("\\")
        # Ensure preserved dir basenames cover everything
        first = entry["path"].split("/", 1)[0]
        assert first in PRESERVED_DIR_BASENAMES


# ---------------------------------------------------------------------------
# CLI-level test
# ---------------------------------------------------------------------------


def test_cli_help_includes_cleanup():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src")
        + os.pathsep
        + env.get("PYTHONPATH", "")
    )
    proc = subprocess.run(
        [sys.executable, "-m", "egfr_myo1d.cli", "--help"],
        cwd=str(REPO_ROOT),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode == 0
    out = (proc.stdout or b"").decode("utf-8", errors="replace") + (
        proc.stderr or b""
    ).decode("utf-8", errors="replace")
    assert "cleanup" in out


def test_cli_cleanup_subcommand_help():
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src")
        + os.pathsep
        + env.get("PYTHONPATH", "")
    )
    proc = subprocess.run(
        [sys.executable, "-m", "egfr_myo1d.cli", "cleanup", "--help"],
        cwd=str(REPO_ROOT),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode == 0
    out = (proc.stdout or b"").decode("utf-8", errors="replace")
    assert "--run-id" in out
    assert "--mode" in out
    assert "--dry-run" in out


def test_cli_cleanup_path_traversal_rejected(tmp_path, monkeypatch):
    """Negative-test: ../bad_run must not be accepted as a run_id."""
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT / "fresh" / "src")
        + os.pathsep
        + env.get("PYTHONPATH", "")
    )
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "egfr_myo1d.cli",
            "cleanup",
            "--run-id",
            "../bad_run",
            "--mode",
            "test",
        ],
        cwd=str(REPO_ROOT),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.returncode != 0
