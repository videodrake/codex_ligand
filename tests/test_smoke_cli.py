"""Low-cost CLI smoke tests for pre-qsub validation."""

import os
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _run_main(*args: str) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    env["PYTHONUTF8"] = "1"
    return subprocess.run(
        [sys.executable, "-X", "utf8", "main.py", *args],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=env,
        check=False,
    )


def test_main_help_smoke() -> None:
    result = _run_main("--help")
    assert result.returncode == 0, result.stderr
    assert "vina" in result.stdout
    assert "validate" in result.stdout
    assert "full" in result.stdout


def test_validate_help_smoke() -> None:
    result = _run_main("validate", "--help")
    assert result.returncode == 0, result.stderr
    assert "show this help message" in result.stdout
