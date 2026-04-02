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


# ---------------------------------------------------------------------------
# Log index — symlinks to scattered PPI logs for discoverability
# ---------------------------------------------------------------------------

_RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]
_N_SEEDS = 5


def create_log_index(config: Optional[dict] = None) -> Path:
    """Create workflow_a/logs/ with symlinks to PPI per-seed logs + LOG_INDEX.txt."""
    import os
    import shutil

    logs = wa_logs(config)
    logs.mkdir(parents=True, exist_ok=True)

    # PPI seed log symlinks
    ppi_seed_dir = logs / "ppi_seeds"
    if ppi_seed_dir.exists():
        shutil.rmtree(ppi_seed_dir)
    ppi_seed_dir.mkdir(parents=True, exist_ok=True)

    ppi_root = wa_phase2_ppi_docking(config)
    for state in _RECEPTOR_STATES:
        for seed in range(_N_SEEDS):
            seed_dir = ppi_root / state / f"prod_seed{seed}"
            for log_name in ["logs/pipeline.log", "logs/workers.log", "PROGRESS.log"]:
                src = seed_dir / log_name
                link_name = f"{state}_seed{seed}_{log_name.replace('logs/', '').replace('/', '_')}"
                link_path = ppi_seed_dir / link_name
                try:
                    rel = os.path.relpath(src, ppi_seed_dir)
                    if link_path.exists() or link_path.is_symlink():
                        link_path.unlink()
                    os.symlink(rel, link_path)
                except OSError:
                    pass

    # LOG_INDEX.txt
    index_path = logs / "LOG_INDEX.txt"
    lines = [
        "=== EGFR-MYO1D Pipeline Log Index ===",
        "",
        "┌─────────────────────────────────────────────────────────┐",
        "│  문제 발생 시 확인 순서 (Troubleshooting Flowchart)     │",
        "├─────────────────────────────────────────────────────────┤",
        "│                                                         │",
        "│  1. PROGRESS.log  ← 먼저 확인 (단계별 요약 + 에러 수)  │",
        "│     → '[!]' 또는 'ERROR' 검색                          │",
        "│                                                         │",
        "│  2. pipeline.log  ← 필터/클러스터링 상세 (메인 프로세스)│",
        "│     → 'FilterStage', 'FAIL', 'WARNING' 검색            │",
        "│                                                         │",
        "│  3. workers.log   ← 개별 모델 에러 (워커 프로세스)     │",
        "│     → 'ERROR', traceback 검색                          │",
        "│                                                         │",
        "│  요약: PROGRESS → pipeline → workers 순서로 범위 좁히기│",
        "└─────────────────────────────────────────────────────────┘",
        "",
        "Log files:",
        "  PROGRESS.log     실시간 진행 + 에러 요약  (tail -f 용)",
        "  pipeline.log     메인 프로세스 상세 로그   (필터/설정)",
        "  workers.log      워커 프로세스 상세 로그   (모델별 에러)",
        "",
        "PPI docking logs (symlinked):",
        "  logs/ppi_seeds/  Per-seed 로그 심볼릭 링크",
        "",
        "Quick commands:",
        "  # 실시간 모니터링",
        "  tail -f phase2_ppi_docking/{state}/prod_seed{n}/PROGRESS.log",
        "",
        "  # 에러 한 번에 찾기 (모든 시드)",
        "  grep -h 'ERROR\\|CRITICAL\\|\\[!\\]' logs/ppi_seeds/*_PROGRESS.log",
        "",
        "  # 워커 에러 상세 (특정 시드)",
        "  grep 'ERROR' logs/ppi_seeds/{state}_seed{n}_workers.log",
        "",
        "  # 필터 단계 추적",
        "  grep 'FilterStage\\|Result\\|Step' logs/ppi_seeds/*_pipeline.log",
        "",
    ]
    index_path.write_text("\n".join(lines), encoding="utf-8")

    return index_path


def create_results_guide(config: Optional[dict] = None) -> Path:
    """Create workflow_a/RESULTS_GUIDE.txt — file-level navigation for users."""
    wa = _workflow_a_root(config)
    guide_path = wa / "RESULTS_GUIDE.txt"
    lines = [
        "=" * 64,
        "  Workflow A 결과 파일 가이드",
        "=" * 64,
        "",
        "★ 최종 결과 (여기서 시작하세요)",
        "─" * 50,
        "  phase6_report/project_report.txt",
        "    → Vina + PPI 통합 보고서",
        "",
        "  phase5_verdict/valid_sites.csv",
        "    → 후보 결합 포켓 순위 (통합 스코어)",
        "",
        "  phase5_verdict/cross_method_agreement.csv",
        "    → Vina/PPI 교차 검증 결과",
        "",
        "  phase7_validation/validation_summary.txt",
        "    → 출력 품질 검증",
        "",
        "",
        "◎ Phase별 핵심 결과",
        "─" * 50,
        "  [Vina]",
        "  phase4_vina_postprocess/vina_pocket_table.csv",
        "    → 포켓 요약 (위치, 친화도, 안정성)",
        "  phase4_vina_postprocess/vina_drug_pocket_map.csv",
        "    → 리간드→포켓 매핑",
        "",
        "  [PPI]",
        "  phase2_ppi_docking/{state}/prod_seed{n}/RESULTS_GUIDE.txt",
        "    → 시드별 상세 가이드 (시드 폴더 내)",
        "  phase3_ppi_postprocess/ppi_pyrosetta_summary.csv",
        "    → 전 시드 통합 인터페이스 요약",
        "",
        "",
        "◇ 중간 데이터 (일반적으로 불필요)",
        "─" * 50,
        "  phase1_vina_docking/       원본 Vina 포즈 (.pdbqt)",
        "  phase4_vina_postprocess/   Vina 후처리 상세 CSV 10+개",
        "  phase2_ppi_docking/        시드별 PPI 중간 CSV 11+개",
        "",
        "",
        "파일 수 참고:",
        "  Workflow A 전체: ~30+ CSV (대부분 중간 결과)",
        "  사용자가 볼 파일: 위 ★ 4개 + ◎ 3개 = 7개",
        "",
        "─" * 50,
        "  ★=필수  ◎=권장  ◇=선택 (재현성/디버깅용)",
        "",
    ]
    try:
        guide_path.write_text("\n".join(lines), encoding="utf-8")
    except OSError:
        pass
    return guide_path
