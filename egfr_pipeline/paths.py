"""Centralized output path resolution for both Workflow A and B.

Every module that needs an output path should import from here instead
of computing paths locally.  This module is the single source of truth
for the directory layout.

Target layout
=============
output/
├── workflow_a/
│   ├── phase1_vina_docking/{receptor_id}/
│   ├── phase2_ppi_docking/{state}/prod_seed{n}/
│   ├── phase3_ppi_postprocess/
│   ├── phase4_vina_postprocess/
│   ├── phase5_verdict/
│   ├── phase6_report/
│   ├── phase7_validation/
│   └── logs/
├── workflow_b/
│   ├── phase1_ppi_analysis/
│   ├── phase2_pocket_analysis/
│   ├── phase3_focused_docking/
│   └── phase4_scoring/
└── precheck/
"""

from pathlib import Path
from typing import Optional

# Repo root (two levels up from egfr_pipeline/paths.py)
REPO_ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# Core resolution
# ---------------------------------------------------------------------------

def output_root(config: Optional[dict] = None) -> Path:
    """Top-level output directory from config or default."""
    if config:
        return Path(config.get("output_root", "./output"))
    return Path("./output")


def workflow_a_root(config: Optional[dict] = None) -> Path:
    return output_root(config) / "workflow_a"


def workflow_b_root(config: Optional[dict] = None) -> Path:
    return output_root(config) / "workflow_b"


# ---------------------------------------------------------------------------
# Workflow A phase directories
# ---------------------------------------------------------------------------

def wa_phase1_vina_docking(config: Optional[dict] = None) -> Path:
    return workflow_a_root(config) / "phase1_vina_docking"


def wa_phase1_vina_receptor(config: Optional[dict] = None, receptor_id: str = "") -> Path:
    return wa_phase1_vina_docking(config) / receptor_id


def wa_phase2_ppi_docking(config: Optional[dict] = None) -> Path:
    return workflow_a_root(config) / "phase2_ppi_docking"


def wa_phase2_ppi_seed(
    config: Optional[dict] = None,
    state: str = "",
    run_type: str = "prod",
    seed: int = 0,
) -> Path:
    return wa_phase2_ppi_docking(config) / state / f"{run_type}_seed{seed}"


def wa_phase2_runtime_inputs(config: Optional[dict] = None) -> Path:
    return wa_phase2_ppi_docking(config) / "runtime_inputs"


def wa_phase3_ppi_postprocess(config: Optional[dict] = None) -> Path:
    return workflow_a_root(config) / "phase3_ppi_postprocess"


def wa_phase4_vina_postprocess(config: Optional[dict] = None) -> Path:
    return workflow_a_root(config) / "phase4_vina_postprocess"


def wa_phase5_verdict(config: Optional[dict] = None) -> Path:
    return workflow_a_root(config) / "phase5_verdict"


def wa_phase6_report(config: Optional[dict] = None) -> Path:
    return workflow_a_root(config) / "phase6_report"


def wa_phase7_validation(config: Optional[dict] = None) -> Path:
    return workflow_a_root(config) / "phase7_validation"


def wa_logs(config: Optional[dict] = None) -> Path:
    return workflow_a_root(config) / "logs"


# ---------------------------------------------------------------------------
# Workflow B phase directories
# ---------------------------------------------------------------------------

def wb_phase1_ppi_analysis(config: Optional[dict] = None) -> Path:
    return workflow_b_root(config) / "phase1_ppi_analysis"


def wb_phase2_pocket_analysis(config: Optional[dict] = None) -> Path:
    return workflow_b_root(config) / "phase2_pocket_analysis"


def wb_phase3_focused_docking(config: Optional[dict] = None) -> Path:
    return workflow_b_root(config) / "phase3_focused_docking"


def wb_phase4_scoring(config: Optional[dict] = None) -> Path:
    return workflow_b_root(config) / "phase4_scoring"


# ---------------------------------------------------------------------------
# Precheck
# ---------------------------------------------------------------------------

def precheck_dir(config: Optional[dict] = None) -> Path:
    return output_root(config) / "precheck"


def precheck_status_file(config: Optional[dict] = None) -> Path:
    return precheck_dir(config) / "last_pass.json"


# ---------------------------------------------------------------------------
# Backward compatibility — maps old flat project_root to new phase dirs
# ---------------------------------------------------------------------------

# Old layout: output/egfr_myo1d_vina/ (everything flat)
# Old layout: output/phase1_ppi/ (PPI separate)
# These helpers let code migrate incrementally.

def legacy_project_root(config: dict) -> Path:
    """Old-style project root for code not yet migrated."""
    output = Path(config.get("output_root", "./output"))
    project_name = config.get("project_name")
    return output / project_name if project_name else output


def legacy_phase1_ppi_dir() -> Path:
    """Old-style PPI output dir."""
    return REPO_ROOT / "output" / "phase1_ppi"


# ---------------------------------------------------------------------------
# Ensure directories exist
# ---------------------------------------------------------------------------

def ensure_wa_dirs(config: Optional[dict] = None) -> None:
    """Create all Workflow A phase directories."""
    for fn in [
        wa_phase1_vina_docking,
        wa_phase2_ppi_docking,
        wa_phase3_ppi_postprocess,
        wa_phase4_vina_postprocess,
        wa_phase5_verdict,
        wa_phase6_report,
        wa_phase7_validation,
        wa_logs,
    ]:
        fn(config).mkdir(parents=True, exist_ok=True)


def ensure_wb_dirs(config: Optional[dict] = None) -> None:
    """Create all Workflow B phase directories."""
    for fn in [
        wb_phase1_ppi_analysis,
        wb_phase2_pocket_analysis,
        wb_phase3_focused_docking,
        wb_phase4_scoring,
    ]:
        fn(config).mkdir(parents=True, exist_ok=True)
