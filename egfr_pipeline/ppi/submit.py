"""PBS job submission for PPI docking.

Extracted from main.py to enable testing and reuse.
"""
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Tuple


def resolve_pbs_script(
    config_ini: Optional[str],
    repo_root: Path,
    is_test: Optional[bool] = None,
) -> Tuple[Path, str]:
    """Determine PBS script path and log prefix from config.

    Args:
        config_ini: INI config path (used to infer test vs prod).
        repo_root: Repository root for locating PBS scripts.
        is_test: Override test/prod detection. If None, inferred from config_ini.

    Returns:
        (pbs_path, log_prefix)
    """
    if is_test is None:
        is_test = "test" in (config_ini or "")
    pbs_name = "run_ppi_test.pbs" if is_test else "run_ppi_prod.pbs"
    log_prefix = "ppi_test" if is_test else "ppi_prod"
    return repo_root / "config" / pbs_name, log_prefix


def check_qsub_available() -> bool:
    """Check if qsub is available on PATH."""
    return shutil.which("qsub") is not None


def submit_pbs_job(
    pbs_path: Path,
    config_ini: Optional[str] = None,
    run_mode: Optional[str] = None,
    cwd: Optional[Path] = None,
) -> Tuple[bool, str]:
    """Submit a PBS job via qsub.

    Args:
        pbs_path: Path to the PBS script.
        config_ini: Config file to pass via -v CONFIG_FILE=...
        run_mode: Run mode to pass via -v RUN_MODE=...
        cwd: Working directory for qsub.

    Returns:
        (success, message) — message is job_id on success, error on failure.
    """
    cmd = ["qsub"]
    if run_mode:
        cmd += ["-v", f"RUN_MODE={run_mode}"]
    elif config_ini:
        cmd += ["-v", f"CONFIG_FILE={config_ini}"]
    cmd.append(str(pbs_path))

    result = subprocess.run(
        cmd, capture_output=True, text=True,
        cwd=str(cwd) if cwd else None,
    )

    if result.returncode == 0:
        return True, result.stdout.strip()
    else:
        return False, result.stderr.strip()
