from pathlib import Path

from scripts.precheck_guard import should_enforce_precheck, validate_precheck_status


def test_should_enforce_precheck_for_normal_modes():
    assert should_enforce_precheck("auto", "0") is True
    assert should_enforce_precheck("force", "0") is True
    assert should_enforce_precheck("status", "0") is False
    assert should_enforce_precheck("auto", "1") is False


def test_validate_precheck_status_requires_marker(tmp_path: Path):
    ok, message = validate_precheck_status(tmp_path, mode="auto", skip_precheck_guard="0")

    assert ok is False
    assert "Missing pre-qsub success marker" in message
    assert "run_pre_qsub_checks.pbs" in message


def test_validate_precheck_status_accepts_existing_marker(tmp_path: Path):
    status_dir = tmp_path / "output" / "precheck"
    status_dir.mkdir(parents=True, exist_ok=True)
    (status_dir / "last_pass.json").write_text('{"status":"passed"}', encoding="utf-8")

    ok, message = validate_precheck_status(tmp_path, mode="auto", skip_precheck_guard="0")

    assert ok is True
    assert "Found pre-qsub success marker" in message


def test_validate_precheck_status_rejects_failed_marker(tmp_path: Path):
    status_dir = tmp_path / "output" / "precheck"
    status_dir.mkdir(parents=True, exist_ok=True)
    (status_dir / "last_pass.json").write_text('{"status":"failed"}', encoding="utf-8")

    ok, message = validate_precheck_status(tmp_path, mode="auto", skip_precheck_guard="0")

    assert ok is False
    assert "not marked as passed" in message
