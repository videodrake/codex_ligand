"""Small, command-provenance-safe Vina adapter for M3-T4 smoke tests."""

from __future__ import annotations

import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class VinaAvailability:
    executable: str
    available: bool
    version: str | None


@dataclass(frozen=True)
class VinaRunResult:
    invoked: bool
    exit_code: int | None
    timeout: bool
    started_at: str | None
    finished_at: str | None
    duration_seconds: float | None


def detect_vina(vina_executable: str = "vina") -> VinaAvailability:
    """Detect Vina availability using only an allowed version probe."""
    resolved = shutil.which(vina_executable) if not Path(vina_executable).is_file() else vina_executable
    if not resolved:
        return VinaAvailability(executable=vina_executable, available=False, version=None)
    try:
        proc = subprocess.run(
            [resolved, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return VinaAvailability(executable=vina_executable, available=False, version=None)
    text = (proc.stdout or proc.stderr or "").strip()
    version = text.splitlines()[0] if text else None
    return VinaAvailability(executable=resolved, available=proc.returncode == 0 or bool(version), version=version)


def build_vina_argv(
    vina_executable: str,
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: tuple[str, str, str],
    size: tuple[str, str, str],
    output_pdbqt: Path,
    vina_log: Path,
    exhaustiveness: int,
    num_modes: int,
    cpu: int,
    seed: int,
) -> list[str]:
    # This cluster's AutoDock Vina build does not support --log; stdout is captured
    # separately and mirrored to the manifest vina_log_file by smoke/runner code.
    return [
        vina_executable,
        "--receptor",
        str(receptor_pdbqt),
        "--ligand",
        str(ligand_pdbqt),
        "--center_x",
        str(center[0]),
        "--center_y",
        str(center[1]),
        "--center_z",
        str(center[2]),
        "--size_x",
        str(size[0]),
        "--size_y",
        str(size[1]),
        "--size_z",
        str(size[2]),
        "--out",
        str(output_pdbqt),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        str(num_modes),
        "--cpu",
        str(cpu),
        "--seed",
        str(seed),
    ]


def run_vina_once(
    argv: list[str],
    stdout_file: Path,
    stderr_file: Path,
    timeout_seconds: int,
    now_iso,
) -> VinaRunResult:
    """Run exactly one Vina command and capture stdout/stderr separately."""
    stdout_file.parent.mkdir(parents=True, exist_ok=True)
    stderr_file.parent.mkdir(parents=True, exist_ok=True)
    started_at = now_iso()
    started = time.monotonic()
    try:
        proc = subprocess.run(
            argv,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=timeout_seconds,
        )
        stdout_file.write_text(proc.stdout or "", encoding="utf-8")
        stderr_file.write_text(proc.stderr or "", encoding="utf-8")
        finished_at = now_iso()
        return VinaRunResult(
            invoked=True,
            exit_code=proc.returncode,
            timeout=False,
            started_at=started_at,
            finished_at=finished_at,
            duration_seconds=round(time.monotonic() - started, 3),
        )
    except subprocess.TimeoutExpired as exc:
        stdout_file.write_text(exc.stdout or "", encoding="utf-8")
        stderr_file.write_text(exc.stderr or "Vina smoke command timed out\n", encoding="utf-8")
        finished_at = now_iso()
        return VinaRunResult(
            invoked=True,
            exit_code=None,
            timeout=True,
            started_at=started_at,
            finished_at=finished_at,
            duration_seconds=round(time.monotonic() - started, 3),
        )
    except OSError as exc:
        stdout_file.write_text("", encoding="utf-8")
        stderr_file.write_text(f"Vina invocation failed: {exc.__class__.__name__}\n", encoding="utf-8")
        finished_at = now_iso()
        return VinaRunResult(
            invoked=True,
            exit_code=None,
            timeout=False,
            started_at=started_at,
            finished_at=finished_at,
            duration_seconds=round(time.monotonic() - started, 3),
        )
