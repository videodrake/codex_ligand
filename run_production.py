#!/usr/bin/env python3
"""프로덕션 전체 파이프라인 — Vina + PPI + Verdict + Report.

Usage:
    conda activate pyrosetta
    python run_production.py            # 자동 이어하기 (완료된 Phase 스킵)
    python run_production.py --force    # 전체 재실행 (기존 결과 무시)
    python run_production.py --from 4   # Phase 4부터 실행 (이전 Phase 스킵)

전체 흐름:
  Phase 1: Vina blind docking (3 receptor × 3 ligand, exhaustiveness=128)  ~15분
           결과물: output/{project}/{receptor_id}/{ligand}_blind.pdbqt
  Phase 2: PPI docking — Phase 1 monomer-based (3 states × 5 seeds = 300K models)
           결과물: {docking_dir}/final_result/final_ranking.csv
  Phase 3: PPI postprocess — monomer 타겟은 chain restoration 불필요 (자동 스킵)
  Phase 4: Vina postprocess (parse → contacts → cluster → summarize → compare → bootstrap)
           결과물: output/{project}/vina_pocket_table.csv
  Phase 5: Verdict (3축 통합 scoring)
           결과물: output/{project}/valid_sites.csv
  Phase 6: Report
           결과물: output/{project}/project_report.txt
  Phase 7: Validate (항상 실행)
"""

import argparse
import configparser
import csv
import json
import os
import sys
import time
import traceback
from pathlib import Path
from typing import Callable, List, Optional

from egfr_pipeline.config import load_config as load_pipeline_config
from egfr_pipeline.config import resolve_resource_config
from egfr_pipeline.step_view import (
    build_current_run_manifest,
    record_step1_outputs,
    record_step2_outputs,
    record_step3_outputs,
    record_step4_outputs,
    record_step5_outputs,
    record_step6_outputs,
    record_step7_outputs,
    refresh_root_step_views,
    step_output_view_enabled,
    update_run_overview,
    utc_now_iso,
    write_run_status,
)
from egfr_pipeline.runtime import (
    cap_worker_count,
    lane_is_complete,
    lane_runtime_dir,
    prepare_scratch_workspace,
    reset_lane_completion,
    resolve_runtime_resources,
    sync_tree,
    write_lane_completion,
    write_lane_manifest,
)
REPO_ROOT = Path(__file__).resolve().parent
CONFIG_PATH = REPO_ROOT / "config" / "example-project.yaml"

# PPI 설정 — Phase 1 monomer-based (3 states × 5 seeds)
RECEPTOR_STATES = ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]
PRODUCTION_N_SEEDS = 5  # generate_configs.py와 동일


def _build_ppi_targets() -> List[dict]:
    from egfr_pipeline.paths import wa_phase2_ppi_seed
    targets = []
    for state in RECEPTOR_STATES:
        for seed_idx in range(PRODUCTION_N_SEEDS):
            targets.append({
                "name": f"{state}_seed{seed_idx}",
                "config_ini": f"config/phase1/phase1_prod_{state}_seed{seed_idx}.ini",
                "docking_dir": str(wa_phase2_ppi_seed(state=state, seed=seed_idx)),
                "mapping_csv": "",
                "receptor_id": state,
                "partner_name": "ext_beta_meander",
                "construct_type": "full_kinase_domain",
                "orientation_validation_status": "not_available",
                "seed_index": seed_idx,
                "is_production": True,
            })
    return targets


PPI_TARGETS = _build_ppi_targets()


_FORCE_MODE = False


def _record_step2_outputs_with_targets(config_path, repo_root=None):
    return record_step2_outputs(
        config_path,
        repo_root=repo_root,
        ppi_targets=PPI_TARGETS,
    )


def _record_step7_outputs_with_result(
    config_path,
    repo_root=None,
    validation_result=None,
    error_message=None,
):
    return record_step7_outputs(
        config_path,
        repo_root=repo_root,
        validation_result=validation_result,
        error_message=error_message,
    )


STEP_VIEW_COLLECTORS = {
    1: ("step1_vina_raw", record_step1_outputs),
    2: ("step2_ppi_raw", _record_step2_outputs_with_targets),
    3: ("step3_ppi_postprocess", record_step3_outputs),
    4: ("step4_vina_postprocess", record_step4_outputs),
    5: ("step5_verdict", record_step5_outputs),
    6: ("step6_report", record_step6_outputs),
    7: ("step7_validate", _record_step7_outputs_with_result),
}


def banner(msg: str):
    width = 60
    print()
    print("=" * width)
    print(f"  {msg}")
    print("=" * width)


def run_step(name: str, func, *args, **kwargs):
    banner(name)
    t0 = time.time()
    try:
        result = func(*args, **kwargs)
        elapsed = time.time() - t0
        minutes = int(elapsed // 60)
        seconds = int(elapsed % 60)
        print(f"\n  [OK] {name} 완료 ({minutes}분 {seconds}초)")
        return True, result, None
    except Exception as e:
        elapsed = time.time() - t0
        minutes = int(elapsed // 60)
        seconds = int(elapsed % 60)
        print(f"\n  [FAIL] {name} 실패 ({minutes}분 {seconds}초)")
        print(f"  Error: {e}")
        import traceback
        traceback.print_exc()
        return False, None, e


def _print_completion_summary(
    phase_states: list,
    config: dict,
    hours: int,
    minutes: int,
) -> None:
    """프로덕션 완료 후 결과 요약 블록 출력."""
    from egfr_pipeline.paths import (
        output_root,
        wa_phase4_vina_postprocess, wa_phase5_verdict,
        wa_phase6_report, wa_phase7_validation, wa_logs,
        create_log_index,
    )

    # Log index 생성
    try:
        create_log_index(config)
    except Exception:
        pass

    w = 60
    print()
    print("╔" + "═" * w + "╗")
    print(f"║  {'프로덕션 완료 요약':<{w - 2}}║")
    print("╠" + "═" * w + "╣")

    # Phase별 상태
    print(f"║  {'Phase':<40} {'상태':<8} {'시간':>8} ║")
    print("║  " + "─" * (w - 4) + "  ║")
    for entry in phase_states:
        num = entry.get("phase_number", "?")
        name = entry.get("phase_name", "")
        # 이름 40자 이내로 자르기
        if len(name) > 38:
            name = name[:35] + "..."
        status = entry.get("status", "pending")
        dur = entry.get("duration_seconds", 0)
        dur_min = int(dur // 60) if dur else 0
        dur_sec = int(dur % 60) if dur else 0

        icon = {"completed": "[OK]", "failed": "[FAIL]", "skipped": "[SKIP]"}.get(status, "[--]")
        time_str = f"{dur_min}m{dur_sec}s" if dur else ""
        print(f"║  {icon} {num}. {name:<35} {time_str:>8} ║")

        if status == "failed":
            err = entry.get("last_error", "")
            if err:
                err_line = str(err).split("\n")[0][:50]
                print(f"║      → {err_line:<{w - 8}}║")

    print("╠" + "═" * w + "╣")
    print(f"║  {'총 소요 시간: ' + str(hours) + '시간 ' + str(minutes) + '분':<{w - 2}}║")
    print("╠" + "═" * w + "╣")

    # 핵심 결과 파일
    failed = [e for e in phase_states if e.get("status") == "failed"]

    print(f"║  {'핵심 결과 파일:':<{w - 2}}║")
    result_files = [
        ("보고서", wa_phase6_report(config) / "project_report.txt"),
        ("판정", wa_phase5_verdict(config) / "valid_sites.csv"),
        ("포켓", wa_phase4_vina_postprocess(config) / "vina_pocket_table.csv"),
        ("검증", wa_phase7_validation(config) / "validation_summary.txt"),
    ]
    for label, path in result_files:
        exists = "✓" if path.exists() else "✗"
        print(f"║    {exists} {label}: {path}  ║")

    # F-4.2: Derived step view health
    print("║" + " " * w + "║")
    step_view_dir = output_root(config)
    n_generated = 0
    n_total_steps = 7
    for i in range(1, n_total_steps + 1):
        step_dir = step_view_dir / f"step{i}_*"
        import glob as _glob
        matches = _glob.glob(str(step_view_dir / f"step{i}_*"))
        if matches and Path(matches[0]).is_dir():
            n_generated += 1
    step_icon = "✓" if n_generated == n_total_steps else "✗"
    print(f"║  Derived Step View: {step_icon} {n_generated}/{n_total_steps} generated  ║")
    if n_generated < n_total_steps:
        print(f"║    → python scripts/rebuild_step_views.py  ║")

    print("║" + " " * w + "║")
    print(f"║  {'로그:':<{w - 2}}║")
    print(f"║    {wa_logs(config) / 'LOG_INDEX.txt'}  ║")

    # 다음 행동 안내
    print("╠" + "═" * w + "╣")
    if failed:
        fail_nums = [str(e.get("phase_number", "?")) for e in failed]
        print(f"║  Phase {', '.join(fail_nums)} 실패. 수정 후:")
        first_fail = fail_nums[0]
        print(f"║    python run_production.py --from {first_fail}")
    else:
        print(f"║  {'모든 Phase 성공. 결과 확인 순서:':<{w - 2}}║")
        print(f"║    1. project_report.txt (종합 보고서)")
        print(f"║    2. valid_sites.csv (후보 포켓 순위)")
        print(f"║    3. PyMOL: 1_OVERVIEW_Clusters.pml (시각화)")

    print("╚" + "═" * w + "╝")
    print()


def _load_config():
    config = load_pipeline_config(str(CONFIG_PATH))
    resources, warnings = resolve_resource_config(config)
    config = dict(config)
    config["resources"] = resources
    if warnings:
        config["_resource_warnings"] = warnings
    return config


def _emit_resource_warnings(config: dict) -> None:
    for warning in config.get("_resource_warnings", []):
        print(f"  [WARN] {warning}")


def _select_ppi_targets(
    *,
    state: Optional[str] = None,
    seed: Optional[int] = None,
) -> List[dict]:
    selected: List[dict] = []
    for target in PPI_TARGETS:
        if state and target.get("receptor_id") != state:
            continue
        target_seed = int(target.get("seed_index", 0))
        if seed is not None and target_seed != int(seed):
            continue
        selected.append(target)
    return selected


def _project_root() -> Path:
    """Workflow A root directory."""
    from egfr_pipeline.paths import workflow_a_root
    return workflow_a_root(_load_config())


# F-4.5: Phase completion sanity check — core output file health
PHASE_CORE_OUTPUTS = {
    4: [
        ("Pocket table", "vina_pocket_table.csv"),
        ("Pose table", "vina_pose_table.csv"),
    ],
    5: [
        ("Verdict", "valid_sites.csv"),
        ("Agreement", "cross_method_agreement.csv"),
    ],
    6: [
        ("Report", "project_report.txt"),
    ],
    7: [
        ("Validation", "validation_summary.txt"),
    ],
}


def _sanity_check_phase_outputs(phase_num: int, config: dict) -> None:
    """Print health of core output files after phase completion (F-4.5)."""
    from egfr_pipeline.paths import (
        wa_phase4_vina_postprocess, wa_phase5_verdict,
        wa_phase6_report, wa_phase7_validation,
    )
    phase_dirs = {
        4: wa_phase4_vina_postprocess(config),
        5: wa_phase5_verdict(config),
        6: wa_phase6_report(config),
        7: wa_phase7_validation(config),
    }
    outputs = PHASE_CORE_OUTPUTS.get(phase_num, [])
    base_dir = phase_dirs.get(phase_num)
    if not outputs or not base_dir:
        return

    for label, filename in outputs:
        path = base_dir / filename
        if not path.exists():
            print(f"    [✗] {label}: NOT FOUND")
        elif path.stat().st_size < 100:
            print(f"    [△] {label}: {path.stat().st_size}B (possibly empty)")
        else:
            size = path.stat().st_size
            print(f"    [✓] {label}: {size:,}B")


def _csv_has_rows(path: Path) -> bool:
    """CSV 파일이 존재하고 데이터 행이 1개 이상인지 확인."""
    if not path.exists():
        return False
    try:
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            next(reader, None)  # header
            return next(reader, None) is not None
    except Exception:
        return False


def _ppi_docking_dir(target: dict) -> Optional[Path]:
    """PPI 도킹 결과 디렉토리 경로 반환."""
    direct_dir = target.get("docking_dir")
    if direct_dir:
        path = Path(str(direct_dir))
        return path if path.is_absolute() else REPO_ROOT / path

    config_path = REPO_ROOT / target["config_ini"]
    if not config_path.exists():
        return None

    config = configparser.ConfigParser()
    config.read(config_path, encoding="utf-8")

    input_pdb_rel = config.get("Path", "input_pdb_name", fallback="")
    input_pdb = REPO_ROOT / input_pdb_rel
    if not input_pdb.exists():
        return None

    return input_pdb.parent


SEED_MARKER_FILENAME = "seed_complete.json"


def _seed_marker_path(target: dict) -> Optional[Path]:
    """PPI target의 seed_complete.json 경로."""
    docking_dir = _ppi_docking_dir(target)
    if docking_dir is None:
        return None
    return docking_dir / SEED_MARKER_FILENAME


def _seed_is_complete(target: dict) -> bool:
    """시드 완료 마커 확인. 하위호환: 마커 없으면 final_ranking.csv로 backfill."""
    marker = _seed_marker_path(target)
    if marker is None:
        return False
    # 1) 마커 파일 직접 확인
    if marker.exists():
        try:
            payload = json.loads(marker.read_text(encoding="utf-8"))
            if payload.get("status") == "completed":
                return True
        except (OSError, json.JSONDecodeError):
            pass
    # 2) 하위호환 migration: final_ranking.csv 존재 시 마커 자동 생성
    docking_dir = marker.parent
    ranking = docking_dir / "final_result" / "final_ranking.csv"
    if ranking.exists():
        _write_seed_completion(target, status="completed")
        print(f"  [MIGRATE] {target['name']}: backfilled seed marker")
        return True
    return False


def _write_seed_completion(
    target: dict,
    *,
    status: str,
    duration_seconds: int = 0,
) -> Optional[Path]:
    """seed_complete.json 마커를 출력 디렉토리에 기록."""
    docking_dir = _ppi_docking_dir(target)
    if docking_dir is None:
        return None
    docking_dir.mkdir(parents=True, exist_ok=True)
    marker_path = docking_dir / SEED_MARKER_FILENAME
    payload = {
        "status": status,
        "target_name": target["name"],
        "receptor_id": target.get("receptor_id", ""),
        "seed_index": target.get("seed_index", 0),
        "config_ini": target.get("config_ini", ""),
        "docking_dir": target.get("docking_dir", ""),
        "completed_at": utc_now_iso(),
        "duration_seconds": duration_seconds,
    }
    marker_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return marker_path


def _execution_mode_label(args: argparse.Namespace) -> str:
    if getattr(args, "lane", ""):
        return f"lane:{args.lane}"
    if args.only:
        return f"only:{args.only}"
    if args.from_phase:
        return f"from:{args.from_phase}"
    if args.force:
        return "force"
    return "full"


def _selected_step_numbers(
    args: argparse.Namespace,
    only_phases: set[int],
) -> List[int]:
    if only_phases:
        return sorted(step_num for step_num in only_phases if step_num in STEP_VIEW_COLLECTORS)
    if args.from_phase:
        return sorted(step_num for step_num in STEP_VIEW_COLLECTORS if step_num >= args.from_phase)
    if args.force:
        return sorted(STEP_VIEW_COLLECTORS)
    return []


def _empty_phase_states() -> List[dict]:
    phase_states: List[dict] = []
    for phase_num, name, _ in PHASES:
        phase_states.append(
            {
                "phase_number": phase_num,
                "phase_name": name,
                "status": "pending",
                "started_at": "",
                "completed_at": "",
                "duration_seconds": None,
                "skip_reason": "",
                "last_error": "",
            }
        )
    return phase_states


def _phase_state_entry(phase_states: List[dict], phase_num: int) -> dict:
    for entry in phase_states:
        if int(entry.get("phase_number", -1)) == phase_num:
            return entry
    raise KeyError(f"Unknown phase number: {phase_num}")


def _build_run_status_payload(
    config: dict,
    *,
    execution_mode: str,
    run_id: str,
    started_at: str,
    phase_states: List[dict],
    overall_status: str,
    current_phase_number: Optional[int] = None,
    current_phase_name: str = "",
    completed_at: Optional[str] = None,
    last_error: str = "",
) -> dict:
    project_root = _project_root()
    resolved = sum(
        1
        for entry in phase_states
        if str(entry.get("status", "")) in {"completed", "skipped", "failed"}
    )
    total = len(phase_states)
    return {
        "run_id": run_id,
        "project_name": config.get("project_name", project_root.name),
        "project_root": project_root.as_posix(),
        "source_config": str(CONFIG_PATH.relative_to(REPO_ROOT).as_posix()),
        "execution_mode": execution_mode,
        "overall_status": overall_status,
        "started_at": started_at,
        "completed_at": completed_at or "",
        "updated_at": utc_now_iso(),
        "current_phase_number": current_phase_number,
        "current_phase_name": current_phase_name,
        "last_error": last_error,
        "phase_states": phase_states,
        "summary": {
            "resolved_phases": resolved,
            "total_phases": total,
            "progress_percent": int((resolved / total) * 100) if total else 0,
        },
    }


def _persist_run_status(
    config: dict,
    *,
    execution_mode: str,
    run_id: str,
    started_at: str,
    phase_states: List[dict],
    overall_status: str,
    current_phase_number: Optional[int] = None,
    current_phase_name: str = "",
    completed_at: Optional[str] = None,
    last_error: str = "",
) -> dict:
    payload = _build_run_status_payload(
        config,
        execution_mode=execution_mode,
        run_id=run_id,
        started_at=started_at,
        phase_states=phase_states,
        overall_status=overall_status,
        current_phase_number=current_phase_number,
        current_phase_name=current_phase_name,
        completed_at=completed_at,
        last_error=last_error,
    )
    project_root = _project_root()
    write_run_status(project_root, payload)
    update_run_overview(project_root, run_status=payload)
    return payload


def _pyrosetta_raw_run_paths(config: dict) -> List[str]:
    paths: List[str] = []

    result_dirs = (config.get("ppi") or {}).get("pyrosetta_result_dirs") or {}
    for entries in result_dirs.values():
        if not entries:
            continue
        for entry in entries:
            path = entry.get("path")
            if path and path not in paths:
                paths.append(str(path))

    for target in PPI_TARGETS:
        docking_dir = _ppi_docking_dir(target)
        if docking_dir and str(docking_dir) not in paths:
            paths.append(str(docking_dir))

    return paths


def _refresh_step_view_outputs(
    phase_num: int,
    args: argparse.Namespace,
    config: dict,
    phase_result=None,
    phase_error: Optional[Exception] = None,
    fresh_run: bool = False,
    stale_steps: Optional[List[int]] = None,
) -> None:
    collector_entry = STEP_VIEW_COLLECTORS.get(phase_num)
    if not collector_entry:
        return

    step_label, collector = collector_entry
    try:
        if phase_num == 7:
            step_dir = collector(
                CONFIG_PATH,
                repo_root=REPO_ROOT,
                validation_result=phase_result,
                error_message=str(phase_error) if phase_error else None,
            )
        else:
            if phase_error is not None:
                return
            step_dir = collector(CONFIG_PATH, repo_root=REPO_ROOT)
        manifest_path, index_path = refresh_root_step_views(
            CONFIG_PATH,
            repo_root=REPO_ROOT,
            execution_mode=_execution_mode_label(args),
            fresh_run=fresh_run,
            stale_steps=stale_steps,
            ppi_config_paths=[target["config_ini"] for target in PPI_TARGETS],
            pyrosetta_raw_run_paths=_pyrosetta_raw_run_paths(config),
        )
        print(f"  [STEP] Derived view refreshed: {step_dir}")
        print(f"  [STEP] Run manifest updated: {manifest_path}")
        print(f"  [STEP] Step index updated: {index_path}")
    except Exception as exc:
        print(f"  [WARN] Derived {step_label} refresh failed: {exc}")
        traceback.print_exc()


# ---------------------------------------------------------------------------
# Phase 완료 체크 함수
# ---------------------------------------------------------------------------

def check_phase1() -> List[str]:
    """Phase 1 (Vina): 모든 receptor×ligand .pdbqt 결과 존재 여부.
    결과물: workflow_a/phase1_vina_docking/{receptor_id}/{ligand}_blind.pdbqt
    """
    from egfr_pipeline.paths import wa_phase1_vina_docking
    config = _load_config()
    vina_root = wa_phase1_vina_docking(config)
    mode = config.get("mode", config.get("vina", {}).get("mode", "blind"))
    missing = []

    for rec in config.get("receptors", []):
        for lig in config.get("ligands", []):
            lig_name = Path(lig["pdbqt"]).stem.replace("_ligand", "")
            pdbqt = vina_root / rec["id"] / f"{lig_name}_{mode}.pdbqt"
            if not pdbqt.exists():
                missing.append(f"{rec['id']}/{lig_name}_{mode}.pdbqt")

    return missing


def check_phase2() -> List[str]:
    """Phase 2 (PPI): 각 target의 seed 완료 마커 확인."""
    missing = []
    for target in PPI_TARGETS:
        if not _seed_is_complete(target):
            missing.append(target["name"])
    return missing


def check_phase3() -> List[str]:
    """Phase 3 (PPI postprocess): residue + summary tables exist and are non-empty."""
    from egfr_pipeline.paths import wa_phase3_ppi_postprocess
    ppi_post = wa_phase3_ppi_postprocess(_load_config())
    residue_path = ppi_post / "ppi_pyrosetta_residues.csv"
    summary_path = ppi_post / "ppi_pyrosetta_summary.csv"
    if _csv_has_rows(residue_path) and _csv_has_rows(summary_path):
        return []
    missing = []
    if not _csv_has_rows(residue_path):
        missing.append("ppi_pyrosetta_residues.csv 없음 또는 비어있음")
    if not _csv_has_rows(summary_path):
        missing.append("ppi_pyrosetta_summary.csv 없음 또는 비어있음")
    return missing


def check_phase4() -> List[str]:
    """Phase 4 (Vina Postprocess): 핵심 CSV들 존재+비어있지 않음."""
    from egfr_pipeline.paths import wa_phase4_vina_postprocess
    config = _load_config()
    project_root = wa_phase4_vina_postprocess(config)
    required = [
        "vina_pose_table.csv",
        "vina_pocket_table.csv",
        "vina_drug_pocket_map.csv",
        "vina_pocket_comparison.csv",
        "vina_pocket_bootstrap.csv",
        "vina_postprocess_coverage.csv",
    ]
    missing = []
    for name in required:
        if not _csv_has_rows(project_root / name):
            missing.append(name)
    coverage_path = project_root / "vina_postprocess_coverage.csv"
    if coverage_path.exists():
        with open(coverage_path, newline="", encoding="utf-8") as handle:
            coverage_rows = list(csv.DictReader(handle))
        expected_pairs = len(config.get("receptors", [])) * len(config.get("ligands", []))
        if coverage_rows and len(coverage_rows) != expected_pairs:
            missing.append(f"coverage rows mismatch ({len(coverage_rows)}/{expected_pairs})")
        for row in coverage_rows:
            if str(row.get("status", "")).strip().lower() != "parsed":
                missing.append(
                    f"{row.get('receptor_id', '?')}/{row.get('ligand_id', '?')}={row.get('status', 'unknown')}"
                )
    return missing


def check_phase5() -> List[str]:
    """Phase 5 (Verdict): valid_sites.csv 존재+비어있지 않음."""
    from egfr_pipeline.paths import wa_phase5_verdict
    project_root = wa_phase5_verdict(_load_config())
    missing = []
    for name in ["valid_sites.csv", "cross_method_agreement.csv"]:
        if not _csv_has_rows(project_root / name):
            missing.append(name)
    return missing


def check_phase6() -> List[str]:
    """Phase 6 (Report): project_report.txt 존재 + 내용 검증 (F-4.1)."""
    from egfr_pipeline.paths import wa_phase6_report
    project_root = wa_phase6_report(_load_config())
    missing = []
    report = project_root / "project_report.txt"
    if not report.exists() or report.stat().st_size < 100:
        missing.append("project_report.txt")
    elif report.exists():
        # F-4.1: Detect placeholder reports that look "complete" but are empty shells
        content = report.read_text(encoding="utf-8")
        placeholder_patterns = [
            "No Vina pocket data available",
            "No cross-receptor comparison data available",
            "No PPI auxiliary data available",
            "No verdict data available",
        ]
        if all(pat in content for pat in placeholder_patterns):
            missing.append("project_report.txt (placeholder only)")
    return missing


PHASE_CHECKS = {
    1: ("Vina blind docking 결과 (.pdbqt)", check_phase1),
    2: ("PPI docking 결과 (seed_complete.json)", check_phase2),
    3: ("PPI postprocess (ppi_pyrosetta_residues.csv)", check_phase3),
    4: ("Vina postprocess (pocket_table 등 5개 CSV)", check_phase4),
    5: ("Verdict (valid_sites.csv)", check_phase5),
    6: ("Report (project_report.txt)", check_phase6),
}


def print_status(
    *,
    config: Optional[dict] = None,
    step_view_enabled_flag: bool = True,
):
    """각 Phase별 완료 상태를 출력한다."""
    config = config or _load_config()

    # Step view stale 감지용 manifest (enabled일 때만)
    stale_phases: List[int] = []
    if step_view_enabled_flag:
        manifest = build_current_run_manifest(
            CONFIG_PATH,
            repo_root=REPO_ROOT,
            execution_mode="status",
            ppi_config_paths=[target["config_ini"] for target in PPI_TARGETS],
            pyrosetta_raw_run_paths=_pyrosetta_raw_run_paths(config),
        )

    print("\n  Pipeline 상태:")
    print(f"  {'#':<4} {'상태':<8} {'결과물':<44} {'상세'}")
    print(f"  {'─'*4} {'─'*8} {'─'*44} {'─'*20}")

    for phase_num, (desc, check_fn) in PHASE_CHECKS.items():
        missing = check_fn()
        done = not missing
        if done:
            phase_status = "[DONE]"
            detail = ""
        else:
            phase_status = "[TODO]"
            if phase_num == 2:
                total = len(PPI_TARGETS)
                completed = total - len(missing)
                detail = f"{completed}/{total} seeds 완료"
            else:
                detail = f"누락: {', '.join(missing[:3])}"
                if len(missing) > 3:
                    detail += f" ...+{len(missing)-3}"

        print(f"  {phase_num:<4} {phase_status:<8} {desc:<44} {detail}")

        # Phase 2 시드별 상태 그리드
        if phase_num == 2:
            for state in RECEPTOR_STATES:
                state_targets = [t for t in PPI_TARGETS if t["receptor_id"] == state]
                state_done = sum(1 for t in state_targets if _seed_is_complete(t))
                seeds_str = []
                for t in state_targets:
                    mark = "DONE" if _seed_is_complete(t) else "----"
                    seeds_str.append(f"seed{t['seed_index']}[{mark}]")
                print(f"         {state:<17} ({state_done}/{len(state_targets)})  {' '.join(seeds_str)}")

        # Step view stale 감지
        if step_view_enabled_flag and done and phase_num in STEP_VIEW_COLLECTORS:
            sv = manifest["step_status"].get(f"step{phase_num}", "not_generated")
            if sv in {"not_generated", "incomplete"}:
                stale_phases.append(phase_num)

    print(f"  {7:<4} {'[항상]':<8} {'Validate (매 실행마다 검증)':<44}")

    # Step view stale 경고 (문제가 있을 때만 표시)
    if stale_phases:
        phases_str = ", ".join(str(p) for p in stale_phases)
        print(f"\n  [WARN] Phase {phases_str}의 steps/ 출력이 최신이 아닙니다. 재실행하면 자동 갱신됩니다.")
    print()


# ---------------------------------------------------------------------------
# Phase 실행 함수 (기존과 동일)
# ---------------------------------------------------------------------------

def phase1_vina(
    *,
    allocated_cpus: Optional[int] = None,
    engine: str = "cpu-vina",
):
    """Vina blind docking — receptor별 병렬 실행."""
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from egfr_pipeline.paths import ensure_wa_dirs

    if engine != "cpu-vina":
        raise RuntimeError(f"Unsupported Vina engine for baseline lane: {engine}")

    config = _load_config()
    ensure_wa_dirs(config)
    _emit_resource_warnings(config)
    from egfr_pipeline.vina.vina_executor import ensure_project_config_inputs_ready

    ensure_project_config_inputs_ready(config)
    receptors = config.get("receptors", [])
    ligands = config.get("ligands", [])
    resources_cfg = config.get("resources")
    if not resources_cfg:
        resources_cfg, _warnings = resolve_resource_config(config)
    vina_cfg = config.get("vina", {}) if isinstance(config.get("vina"), dict) else {}
    runtime_allocated = allocated_cpus
    if runtime_allocated is None:
        has_scheduler_env = bool(os.environ.get("PBS_NP")) or any(
            os.environ.get(name)
            for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")
        )
        if not has_scheduler_env:
            runtime_allocated = os.cpu_count() or 1
    runtime = resolve_runtime_resources(allocated_cpus=runtime_allocated)
    cpu_per_job = int(resources_cfg.get("vina", {}).get("cpu_per_job", vina_cfg.get("cpu", 0)))
    configured_parallel = int(
        resources_cfg.get("vina", {}).get("parallel_receptors", len(receptors) or 1)
    )
    max_workers = max(1, min(len(receptors) or 1, configured_parallel))
    if cpu_per_job <= 0 and max_workers > 1:
        print(
            "  [WARN] vina.cpu=0 allows each docking job to use all visible cores; "
            "forcing receptor-level parallelism to 1."
        )
        max_workers = 1

    capped_workers = cap_worker_count(
        requested_workers=max_workers,
        n_jobs=len(receptors) or 1,
        cpu_per_job=cpu_per_job,
        available_cpus=runtime.effective_cpus,
    )
    if capped_workers < max_workers:
        print(
            f"  [WARN] Requested Vina CPU budget exceeds allocated CPUs: "
            f"{max_workers * max(cpu_per_job, 1)} > {runtime.effective_cpus}. "
            f"Capping receptor parallelism {max_workers} -> {capped_workers}."
        )
        max_workers = capped_workers

    print(f"  Receptor {len(receptors)}개 × Ligand {len(ligands)}개 = {len(receptors) * len(ligands)} 도킹 잡")
    print(
        f"  exhaustiveness={vina_cfg.get('exhaustiveness', '?')}, "
        f"n_poses={vina_cfg.get('n_poses', '?')}, "
        f"energy_range={vina_cfg.get('energy_range', '?')}, "
        f"cpu/job={cpu_per_job}, seed={vina_cfg.get('seed', 'auto')}, "
        f"engine={engine}"
    )
    print(
        f"  병렬 실행: {max_workers}개 프로세스 "
        f"(allocated={runtime.allocated_cpus}, effective={runtime.effective_cpus})"
    )

    from egfr_pipeline.vina.vina_executor import dock_one_receptor as _dock_one_receptor

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_dock_one_receptor, rec, ligands, config): rec["id"]
            for rec in receptors
        }
        for future in as_completed(futures):
            rec_id = futures[future]
            try:
                future.result()
                print(f"\n  [OK] {rec_id} 도킹 완료")
            except Exception as e:
                print(f"\n  [FAIL] {rec_id} 도킹 실패: {e}")


def _legacy_phase2_ppi():
    """PyRosetta PPI global blind docking — 순차 실행 (CPU 집약적)."""
    for target in PPI_TARGETS:
        name = target["name"]
        config_ini = target["config_ini"]
        input_pdb = target["input_pdb"]

        pdb_path = Path(input_pdb)
        if not pdb_path.exists():
            alt_path = Path(input_pdb.replace("_wt.pdb", ".pdb"))
            if alt_path.exists():
                pdb_path = alt_path
                print(f"  [INFO] {name}: _wt.pdb 없음, {alt_path.name} 사용")
            else:
                print(f"  [ERROR] {name}: 입력 PDB 없음: {input_pdb}")
                continue

        # 이미 완료된 target은 스킵 (--force 시 무시)
        if not _FORCE_MODE:
            docking_dir = _ppi_docking_dir(target)
            ranking = docking_dir / "final_result" / "final_ranking.csv" if docking_dir else None
            if ranking and ranking.exists():
                print(f"  [SKIP] {name}: 이미 완료 ({ranking})")
                continue

        print(f"\n  --- {name} (20K models) ---")
        print(f"  Config: {config_ini}")
        print(f"  Input:  {pdb_path}")

        from egfr_pipeline.pyrosetta_docking.pipeline_manager import PipelineManager
        PipelineManager(config_ini, str(pdb_path)).execute()


def _legacy_phase3_ppi_postprocess():
    """PPI 결과 chain restoration + residue extraction."""
    config_str = str(CONFIG_PATH)

    for target in PPI_TARGETS:
        name = target["name"]
        mapping_csv = target["mapping_csv"]
        if not mapping_csv:
            print(f"  [SKIP] {name}: monomer 타겟 — chain restoration 불필요")
            continue
        receptor_id = target["receptor_id"]
        partner_name = target["partner_name"]
        construct_type = target.get("construct_type", "full_kinase_domain")
        orientation_validation_status = target.get(
            "orientation_validation_status",
            "not_available",
        )

        input_pdb = Path(target["input_pdb"])
        if not input_pdb.exists():
            input_pdb = Path(target["input_pdb"].replace("_wt.pdb", ".pdb"))
        docking_dir = REPO_ROOT / input_pdb.stem

        if not docking_dir.exists():
            print(f"  [WARN] {name}: 도킹 결과 디렉토리 없음: {docking_dir}")
            continue

        print(f"\n  --- {name} postprocess ---")
        print(f"  Docking dir: {docking_dir}")
        print(f"  Mapping: {mapping_csv}")

        from egfr_pipeline.ppi.postprocess_ppi import postprocess_ppi_results
        postprocess_ppi_results(
            config_path=config_str,
            docking_dir=str(docking_dir),
            mapping_csv=mapping_csv,
            receptor_id=receptor_id,
            partner_name=partner_name,
            construct_type=construct_type,
            orientation_validation_status=orientation_validation_status,
        )


def phase4_vina_postprocess():
    """Vina 후처리 전체 체인."""
    from egfr_pipeline.config import load_config
    from egfr_pipeline.paths import wa_phase4_vina_postprocess, ensure_wa_dirs
    config_str = str(CONFIG_PATH)
    config = load_config(config_str)
    pp = config.get("postprocess") or {}
    ensure_wa_dirs(config)
    project_root = wa_phase4_vina_postprocess(config)

    def _stage_enabled(flag_name: str) -> bool:
        return bool(pp.get(flag_name, True))

    print("\n  --- parse ---")
    if _stage_enabled("parse_results"):
        from egfr_pipeline.vina.parse_poses import build_pose_table_from_config
        out = build_pose_table_from_config(config_str)
        print(f"  → {out}")
        print(f"  → {project_root / 'vina_postprocess_coverage.csv'}")
    else:
        print("  [SKIP] parse_results=false")

    print("\n  --- contacts ---")
    if _stage_enabled("extract_contacts"):
        from egfr_pipeline.vina.pose_contacts import enrich_pose_table_with_contacts
        cutoff = pp.get("contact_cutoff", 4.0)
        out = enrich_pose_table_with_contacts(config_str, cutoff=cutoff)
        print(f"  → {out}")
        print(f"  → {project_root / 'vina_contact_distances.csv'}")
    else:
        print("  [SKIP] extract_contacts=false")

    print("\n  --- cluster ---")
    if _stage_enabled("cluster_pockets"):
        from egfr_pipeline.vina.pocket_cluster import cluster_pose_table
        max_per_pocket = pp.get("max_poses_per_pocket", 0)
        if max_per_pocket > 0:
            print(f"  Pocket cap: max {max_per_pocket} poses per pocket (round-robin by ligand)")
        out = cluster_pose_table(
            config_str,
            cutoff=pp.get("pocket_cutoff", 8.0),
            max_iterations=pp.get("cluster_max_iterations", 10),
            min_pocket_size=pp.get("min_pocket_size", 2),
            merge_by_residue=pp.get("merge_by_residue", True),
            merge_jaccard=pp.get("merge_jaccard", 0.3),
            merge_overlap=pp.get("merge_overlap", 0.5),
            merge_centroid_fallback=pp.get("merge_centroid_fallback", 6.0),
            max_per_pocket=max_per_pocket,
        )
        print(f"  → {out}")
        print(f"  → {project_root / 'vina_clustering_merge_log.csv'}")
        print(f"  → {project_root / 'vina_clustering_parameters.json'}")
        if max_per_pocket > 0:
            print(f"  → {project_root / 'vina_pocket_cap_report.csv'}")
    else:
        print("  [SKIP] cluster_pockets=false")

    print("\n  --- summarize ---")
    if _stage_enabled("summarize_pockets"):
        from egfr_pipeline.vina.pocket_summary import summarize_from_config
        pocket_csv, drug_csv, occupancy_csv = summarize_from_config(config_str)
        print(f"  → {pocket_csv}")
        print(f"  → {drug_csv}")
        print(f"  → {occupancy_csv}")
    else:
        print("  [SKIP] summarize_pockets=false")

    print("\n  --- compare ---")
    if _stage_enabled("compare_pockets"):
        from egfr_pipeline.vina.cross_receptor import compare_from_config
        out = compare_from_config(config_str)
        print(f"  → {out}")
    else:
        print("  [SKIP] compare_pockets=false")

    print("\n  --- bootstrap ---")
    from egfr_pipeline.vina.pocket_stability import bootstrap_from_config
    print("  [INFO] Bootstrap here measures pose-resampling clustering stability, not fresh docking reproducibility.")
    out = bootstrap_from_config(config_str)
    print(f"  → {out}")


def phase5_verdict():
    """3축 증거 통합 판정."""
    from egfr_pipeline.verdict import generate_verdict
    agr_csv, ver_csv = generate_verdict(str(CONFIG_PATH))
    if ver_csv and str(ver_csv) != ".":
        print(f"  → {agr_csv}")
        print(f"  → {ver_csv}")
    else:
        print("  [WARN] pocket table 없음 — verdict 스킵")


def phase6_report():
    """보고서 생성."""
    from egfr_pipeline.report import generate_report
    report_path, csv_path = generate_report(str(CONFIG_PATH))
    print(f"  → {report_path}")
    print(f"  → {csv_path}")


def phase7_validate():
    """출력 검증."""
    from egfr_pipeline.validate import run_validation
    result = run_validation(str(CONFIG_PATH), repo_root=str(REPO_ROOT))
    print(result.summary())
    return result


def _phase3_override_dirs() -> dict:
    overrides = {}
    for target in PPI_TARGETS:
        docking_dir = _ppi_docking_dir(target)
        if not docking_dir or not docking_dir.exists():
            continue
        overrides.setdefault(target["receptor_id"], []).append(
            {
                "path": str(docking_dir),
                "partner": target.get("partner_name", ""),
                "construct_type": target.get("construct_type", "full_kinase_domain"),
                "orientation_validation_status": target.get(
                    "orientation_validation_status",
                    "not_available",
                ),
            }
        )
    return overrides


def phase2_ppi(
    *,
    state: Optional[str] = None,
    seed: Optional[int] = None,
    allocated_cpus: Optional[int] = None,
    scratch_dir: Optional[Path] = None,
):
    """PyRosetta PPI docking using runtime-prepared Phase 1 inputs."""
    from egfr_pipeline.phase1.launch_docking import run_single
    from egfr_pipeline.phase1.prepare_inputs import prepare_phase1_inputs

    prepare_phase1_inputs()
    failures: List[str] = []
    targets = _select_ppi_targets(state=state, seed=seed)
    if not targets:
        selector = f"state={state or '*'}, seed={seed if seed is not None else '*'}"
        raise RuntimeError(f"No PPI targets match selector: {selector}")

    if _FORCE_MODE:
        for target in targets:
            marker = _seed_marker_path(target)
            if marker and marker.exists():
                marker.unlink()
                print(f"  [RESET] {target['name']}: seed marker removed")

    for target in targets:
        name = target["name"]
        config_ini = target["config_ini"]
        config_path = REPO_ROOT / config_ini
        if not config_path.exists():
            print(f"  [ERROR] {name}: config not found: {config_ini}")
            continue

        if not _FORCE_MODE and _seed_is_complete(target):
            print(f"  [SKIP] {name}: seed marker ({_seed_marker_path(target)})")
            continue

        print(f"\n  --- {name} (20K models) ---")
        print(f"  Config: {config_ini}")
        print(f"  Output: {_ppi_docking_dir(target)}")

        run_kwargs = {
            "seed_index": int(target.get("seed_index", 0)),
            "is_production": bool(target.get("is_production", True)),
            "dry_run": False,
        }
        if allocated_cpus is not None:
            run_kwargs["allocated_cpus"] = allocated_cpus
        if scratch_dir is not None:
            run_kwargs["scratch_dir"] = scratch_dir
        t0 = time.time()
        meta = run_single(
            config_path,
            target["receptor_id"],
            **run_kwargs,
        )
        elapsed = int(time.time() - t0)
        if meta.get("status") == "failed":
            failures.append(f"{name}: {meta.get('error', 'unknown error')}")
        else:
            marker = _write_seed_completion(target, status="completed", duration_seconds=elapsed)
            if marker:
                print(f"  [MARKER] {name}: {marker}")

    if failures:
        raise RuntimeError("; ".join(failures))


def phase3_ppi_postprocess():
    """Extract project-level PPI evidence tables from current Phase 2 run dirs."""
    from egfr_pipeline.ppi.pyrosetta_extract import extract_pyrosetta_batch

    override_dirs = _phase3_override_dirs()
    if not override_dirs:
        print("  [WARN] No Phase 2 PyRosetta run directories are available for extraction.")
        return

    residue_csv, summary_csv = extract_pyrosetta_batch(
        str(CONFIG_PATH),
        pyrosetta_result_dirs=override_dirs,
    )
    print(f"  Residues: {residue_csv}")
    print(f"  Summary:  {summary_csv}")
    print(f"  Long-form residues: {residue_csv.parent / 'ppi_pyrosetta_residue_long.csv'}")
    print(f"  Model table:       {residue_csv.parent / 'ppi_pyrosetta_model_table.csv'}")


def _finalize_lane() -> object:
    phase5_verdict()
    phase6_report()
    return phase7_validate()


# ---------------------------------------------------------------------------
# Advanced PPI-First pipeline lane functions (Workflow B)
# ---------------------------------------------------------------------------

def _validate_ko_consistency(patch_csv_path) -> None:
    """Ko et al. active face consistency check (AC-1.3).

    Verifies that Phase 1 hotspot residues include sufficient active face
    (sheet 8/9) residues from Ko et al. experimental data.
    """
    import csv
    import logging
    from pathlib import Path
    from egfr_pipeline.region_definitions import KO_SHEET_8_9, KO_SHEET_10_11, KO_SHEET_12

    log = logging.getLogger(__name__)
    patch_path = Path(patch_csv_path)
    if not patch_path.exists():
        return  # File existence already checked by caller

    with open(patch_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        return

    # Collect hotspot residue numbers (partner-side, i.e. chain B / MYO1D)
    hotspot_resnums = set()
    for row in rows:
        if row.get("is_hotspot_any_state", "").lower() in ("true", "1", "yes"):
            try:
                hotspot_resnums.add(int(row["residue_num"]))
            except (ValueError, KeyError):
                continue

    if not hotspot_resnums:
        return  # No hotspot residues — separate quality guard (Task 3.5)

    active_face_in_hotspot = hotspot_resnums & KO_SHEET_8_9
    sheet_10_11_in_hotspot = hotspot_resnums & KO_SHEET_10_11
    sheet_12_in_hotspot = hotspot_resnums & KO_SHEET_12

    if len(active_face_in_hotspot) < 3:
        raise ValueError(
            f"Ko et al. FAIL: hotspot에 active face (sheet 8/9) 잔기가 "
            f"{len(active_face_in_hotspot)}개뿐입니다 (최소 3개 필요).\n"
            f"Active face 잔기: {sorted(KO_SHEET_8_9)}\n"
            f"Hotspot 잔기: {sorted(hotspot_resnums)}\n"
            f"PPI 예측이 Ko et al. 실험 결과와 괴리됩니다. Workflow B를 중단합니다."
        )

    if sheet_10_11_in_hotspot:
        log.warning(
            "Ko et al. WARNING: hotspot에 sheet 10/11 잔기 감지 (%s). "
            "이 잔기들은 실험적으로 WT-level 기능 유지 — 결합 기여도 낮음.",
            sorted(sheet_10_11_in_hotspot),
        )

    if sheet_12_in_hotspot:
        log.info(
            "Ko et al. INFO: hotspot에 sheet 12 잔기 감지 (%s). "
            "Structural support 영역 — 모니터링 대상.",
            sorted(sheet_12_in_hotspot),
        )


def _validate_adv_handoff(lane: str) -> None:
    """Advanced lane 사전 조건 파일 존재 검증. defense-in-depth."""
    from egfr_pipeline.paths import (
        wa_phase2_ppi_docking, wb_phase1_ppi_analysis,
        wb_phase2_pocket_analysis, wb_phase3_focused_docking,
    )
    config = _load_config()
    ppi_dock_dir = wa_phase2_ppi_docking(config)
    wb_p1 = wb_phase1_ppi_analysis(config)
    wb_p2 = wb_phase2_pocket_analysis(config)
    wb_p3 = wb_phase3_focused_docking(config)

    checks = {
        "adv-phase1": (
            lambda: check_phase2(),
            "PPI 도킹 미완료. Workflow A Phase 2를 먼저 완료하세요.\n"
            "Run: qsub config/run_ppi_state_seed.pbs"
        ),
        "adv-phase2": (
            lambda: [] if (wb_p1 / "phase1_downstream_patch_reference.csv").exists() else ["handoff"],
            f"Phase 1 핸드오프 없음: {wb_p1 / 'phase1_downstream_patch_reference.csv'}\n"
            "Run: qsub config/run_adv_phase1.pbs"
        ),
        "adv-phase3-setup": (
            lambda: [] if (wb_p2 / "phase3_candidate_pocket_reference.csv").exists() else ["handoff"],
            f"Phase 2 핸드오프 없음: {wb_p2 / 'phase3_candidate_pocket_reference.csv'}\n"
            "Run: qsub config/run_adv_phase2.pbs"
        ),
        "adv-phase3-execute": (
            lambda: [] if (wb_p3 / "phase3_docking_job_table.csv").exists() else ["setup"],
            f"Phase 3 job table 없음: {wb_p3 / 'phase3_docking_job_table.csv'}\n"
            "Run: qsub config/run_adv_phase3_setup.pbs"
        ),
        "adv-phase3-post": (
            lambda: [] if _csv_has_rows(wb_p3 / "phase3_round_log.csv") else ["round_log"],
            f"Phase 3 round log 없거나 비어있음.\n"
            "최소 1회 execute round 필요.\n"
            "Run: qsub -v ROUND=0 config/run_adv_phase3_execute.pbs"
        ),
        "adv-phase4": (
            lambda: [] if (wb_p3 / "phase4_docking_evidence_reference.csv").exists() else ["handoff"],
            f"Phase 3 핸드오프 없음: {wb_p3 / 'phase4_docking_evidence_reference.csv'}\n"
            "Run: qsub config/run_adv_phase3_post.pbs"
        ),
    }

    if lane not in checks:
        return

    check_fn, msg = checks[lane]
    missing = check_fn()
    if missing:
        raise FileNotFoundError(msg)

    # Ko et al. consistency check for Phase 1 → Phase 2 handoff (AC-1.3)
    if lane == "adv-phase2":
        _validate_ko_consistency(wb_p1 / "phase1_downstream_patch_reference.csv")

    # ── Data quality guards (AC-3.5) ──────────────────────────────────
    _validate_handoff_data_quality(lane, wb_p1, wb_p2, wb_p3)


def _validate_handoff_data_quality(
    lane: str,
    wb_p1: "Path",
    wb_p2: "Path",
    wb_p3: "Path",
) -> None:
    """Data quality guards for handoff files (AC-3.5).

    Raises FileNotFoundError for FAIL conditions (pipeline must stop).
    Logs warnings for WARNING conditions (pipeline continues).
    """
    import logging
    log = logging.getLogger(__name__)

    # Phase 1 → Phase 2: hotspot 0개 → FAIL
    if lane == "adv-phase2":
        patch_path = wb_p1 / "phase1_downstream_patch_reference.csv"
        if patch_path.exists():
            with open(patch_path, newline="", encoding="utf-8") as f:
                rows = list(csv.DictReader(f))
            hotspot_rows = [
                r for r in rows
                if r.get("is_hotspot_any_state", "").lower() in ("true", "1")
            ]
            if len(hotspot_rows) == 0:
                raise FileNotFoundError(
                    f"Phase 1 데이터 품질 FAIL: hotspot 잔기 0개.\n"
                    f"File: {patch_path}\n"
                    f"Phase 1 결과를 검토하세요 — orientation filter 또는 "
                    f"cluster consensus에 문제가 있을 수 있습니다."
                )

    # Phase 2 → Phase 3: 포켓 0개 → FAIL, 전부 irrelevant → WARNING
    if lane == "adv-phase3-setup":
        pocket_path = wb_p2 / "phase3_candidate_pocket_reference.csv"
        if pocket_path.exists():
            with open(pocket_path, newline="", encoding="utf-8") as f:
                rows = list(csv.DictReader(f))
            if len(rows) == 0:
                raise FileNotFoundError(
                    f"Phase 2 데이터 품질 FAIL: 후보 포켓 0개.\n"
                    f"File: {pocket_path}\n"
                    f"Phase 2 pocket proposal/merge를 검토하세요."
                )
            # Check if all pockets are irrelevant (low_relevance or skip)
            relevant = [
                r for r in rows
                if r.get("docking_priority", "").lower() not in ("skip", "")
                and r.get("relationship_class", "").lower() != "low_relevance_candidate"
            ]
            if len(relevant) == 0:
                log.warning(
                    "Phase 2 WARNING: 모든 후보 포켓이 irrelevant 또는 skip으로 분류됨. "
                    "Phase 3 focused docking 결과가 제한적일 수 있습니다."
                )

    # Phase 3 → Phase 4: 유효 포즈 < 5 → FAIL, 리간드 1종 → WARNING
    if lane == "adv-phase4":
        evidence_path = wb_p3 / "phase4_docking_evidence_reference.csv"
        if evidence_path.exists():
            with open(evidence_path, newline="", encoding="utf-8") as f:
                rows = list(csv.DictReader(f))
            # Count valid poses (rows with non-empty affinity)
            valid_poses = [
                r for r in rows
                if r.get("best_affinity", "").strip() not in ("", "NA", "None")
            ]
            if len(valid_poses) < 5:
                raise FileNotFoundError(
                    f"Phase 3 데이터 품질 FAIL: 유효 포즈 {len(valid_poses)}개 (최소 5개 필요).\n"
                    f"File: {evidence_path}\n"
                    f"Phase 3 focused docking이 충분한 결과를 생성하지 못했습니다."
                )
            # Check ligand diversity
            ligand_ids = {
                r.get("ligand_id", "").strip()
                for r in valid_poses
                if r.get("ligand_id", "").strip()
            }
            if len(ligand_ids) <= 1:
                log.warning(
                    f"Phase 3 WARNING: 리간드 1종({next(iter(ligand_ids), '?')})만 유효 포즈 보유. "
                    f"다양성 부족으로 consensus 해석에 주의가 필요합니다."
                )


def _adv_phase1():
    """Run Phase 1 analysis: TG 1.1.5 (standardize) + TG 1.2 (extract interface)
    + PPI evidence extraction. Assumes PPI docking (Workflow A Phase 2) is complete."""
    _validate_adv_handoff("adv-phase1")
    from egfr_pipeline.phase1.extract_interface import run_extract_interface
    from egfr_pipeline.phase1.standardize_scores import run_standardize_scores

    print("  --- TG 1.1.5: Standardize Scores ---")
    run_standardize_scores()
    print("  --- TG 1.2: Extract Interface ---")
    run_extract_interface()
    # PPI evidence extraction (reuse Workflow A Phase 3 logic)
    print("  --- PPI Evidence Extraction ---")
    phase3_ppi_postprocess()


def _adv_phase2():
    """Run Phase 2 cascade (TG 2.0-2.7). Requires Phase 1 analysis complete."""
    _validate_adv_handoff("adv-phase2")
    from egfr_pipeline.phase2.rerun_cascade import run_phase2_cascade
    run_phase2_cascade(parse_only=True)


def _adv_phase3_setup(allocated_cpus: Optional[int] = None):
    """Run Phase 3 setup (TG 3.0-3.2 + script generation)."""
    _validate_adv_handoff("adv-phase3-setup")
    from egfr_pipeline.phase3.rerun_cascade import run_phase3_cascade
    run_phase3_cascade(mode="setup", allocated_cpus=allocated_cpus)


def _adv_phase3_execute(
    round_id: Optional[int] = None,
    allocated_cpus: Optional[int] = None,
):
    """Run Phase 3 Vina execution for one round."""
    _validate_adv_handoff("adv-phase3-execute")
    from egfr_pipeline.phase3.rerun_cascade import run_phase3_cascade
    run_phase3_cascade(
        mode="execute",
        round_id=round_id,
        allocated_cpus=allocated_cpus,
    )


def _adv_phase3_post():
    """Run Phase 3 post-docking analysis (TG 3.4-3.6)."""
    _validate_adv_handoff("adv-phase3-post")
    from egfr_pipeline.phase3.rerun_cascade import run_phase3_cascade
    run_phase3_cascade(mode="post")


def _adv_phase4():
    """Run Phase 4 cascade (TG 4.0-4.6)."""
    _validate_adv_handoff("adv-phase4")
    from egfr_pipeline.phase4.rerun_cascade import run_phase4_cascade
    run_phase4_cascade()


def _lane_requested_cpus(config: dict, lane: str) -> int:
    resources_cfg = config.get("resources", {})
    vina_resources = resources_cfg.get("vina", {})
    gpu_resources = resources_cfg.get("gpu", {})
    if lane == "vina-cpu":
        return max(
            1,
            int(vina_resources.get("cpu_per_job", 1))
            * int(vina_resources.get("parallel_receptors", 1)),
        )
    if lane == "vina-gpu":
        return max(1, int(gpu_resources.get("cpu_per_job", 1)))
    if lane == "ppi":
        return max(1, int(resources_cfg.get("ppi", {}).get("cpu_per_job", 1)))
    if lane in {"ppi-post", "vina-post"}:
        return 4
    if lane == "finalize":
        return 2
    # Advanced PPI-First pipeline lanes
    if lane == "adv-phase1":
        return 4
    if lane == "adv-phase2":
        return 32
    if lane == "adv-phase3-setup":
        return 2
    if lane == "adv-phase3-execute":
        return max(1, int(vina_resources.get("cpu_per_job", 1)))
    if lane == "adv-phase3-post":
        return 4
    if lane == "adv-phase4":
        return 2
    return 1


def _phase_numbers_for_lane(lane: str) -> List[int]:
    mapping = {
        "vina-cpu": [1],
        "ppi": [2],
        "ppi-post": [3],
        "vina-post": [4],
        "finalize": [5, 6, 7],
        # Advanced lanes don't map to standard phases
        "adv-phase1": [],
        "adv-phase2": [],
        "adv-phase3-setup": [],
        "adv-phase3-execute": [],
        "adv-phase3-post": [],
        "adv-phase4": [],
    }
    return mapping.get(lane, [])


def _refresh_lane_step_views(lane: str, config: dict) -> None:
    dummy_args = argparse.Namespace(
        lane=lane,
        only="",
        from_phase=0,
        force=_FORCE_MODE,
    )
    for phase_num in _phase_numbers_for_lane(lane):
        _refresh_step_view_outputs(
            phase_num,
            dummy_args,
            config,
            phase_result=None,
            phase_error=None,
            fresh_run=not _project_root().exists(),
            stale_steps=[],
        )


def _run_lane(
    lane: str,
    *,
    state: Optional[str] = None,
    seed: Optional[int] = None,
    allocated_cpus: Optional[int] = None,
    engine: str = "cpu-vina",
    gpu_required: bool = False,
    scratch_dir: Optional[Path] = None,
) -> object:
    config = _load_config()
    _emit_resource_warnings(config)
    requested_cpus = _lane_requested_cpus(config, lane)
    runtime = resolve_runtime_resources(
        requested_cpus=requested_cpus,
        allocated_cpus=allocated_cpus,
    )
    project_root = _project_root()
    scratch_workspace = prepare_scratch_workspace(
        project_root,
        lane,
        state=state,
        seed=seed,
        engine=engine,
        scratch_dir=scratch_dir,
    )

    if lane in {"vina-gpu", "phase3-gpu"}:
        gpu_cfg = config.get("resources", {}).get("gpu", {})
        if not bool(gpu_cfg.get("enabled", False)):
            raise RuntimeError(
                f"{lane} lane is disabled. Set resources.gpu.enabled=true to opt in."
            )
        if gpu_required and not runtime.visible_gpu_ids:
            raise RuntimeError("GPU required but CUDA_VISIBLE_DEVICES is empty.")
        raise RuntimeError(
            "GPU lanes are intentionally explicit but still require a site-local engine wrapper."
        )

    if lane == "status":
        print_status(config=config, step_view_enabled_flag=step_output_view_enabled(config))
        return {"status": "completed", "lane": lane}

    if not _FORCE_MODE and lane_is_complete(
        project_root,
        lane,
        state=state,
        seed=seed,
        engine=engine,
    ):
        print(f"  [SKIP] lane {lane} already completed via marker.")
        return {"status": "skipped", "lane": lane}

    reset_lane_completion(
        project_root,
        lane,
        state=state,
        seed=seed,
        engine=engine,
    )
    write_lane_manifest(
        project_root,
        lane,
        resources=runtime,
        state=state,
        seed=seed,
        engine=engine,
        config_path=CONFIG_PATH,
        extra={
            "status": "running",
            "requested_cpus": requested_cpus,
            "gpu_required": gpu_required,
            "scratch_workspace": str(scratch_workspace) if scratch_workspace else "",
        },
    )

    started = time.time()
    try:
        if lane == "vina-cpu":
            result = phase1_vina(
                allocated_cpus=runtime.effective_cpus,
                engine=engine,
            )
        elif lane == "ppi":
            result = phase2_ppi(
                state=state,
                seed=seed,
                allocated_cpus=runtime.effective_cpus,
                scratch_dir=scratch_workspace,
            )
        elif lane == "ppi-post":
            result = phase3_ppi_postprocess()
        elif lane == "vina-post":
            result = phase4_vina_postprocess()
        elif lane == "finalize":
            result = _finalize_lane()
        elif lane == "adv-phase1":
            result = _adv_phase1()
        elif lane == "adv-phase2":
            result = _adv_phase2()
        elif lane == "adv-phase3-setup":
            result = _adv_phase3_setup(runtime.effective_cpus)
        elif lane == "adv-phase3-execute":
            result = _adv_phase3_execute(
                round_id=seed,  # reuse --seed as round selector
                allocated_cpus=runtime.effective_cpus,
            )
        elif lane == "adv-phase3-post":
            result = _adv_phase3_post()
        elif lane == "adv-phase4":
            result = _adv_phase4()
        elif lane == "status":
            print_status(config=config, step_view_enabled_flag=step_output_view_enabled(config))
            result = None
        else:
            raise ValueError(f"Unsupported lane: {lane}")

        if lane != "status":
            _refresh_lane_step_views(lane, config)
        if lane == "ppi":
            try:
                from egfr_pipeline.paths import create_log_index
                create_log_index(config)
            except Exception:
                pass
        write_lane_completion(
            project_root,
            lane,
            status="completed",
            state=state,
            seed=seed,
            engine=engine,
            extra={
                "duration_seconds": max(0, int(time.time() - started)),
                "scratch_workspace": str(scratch_workspace) if scratch_workspace else "",
            },
        )
        return result
    except Exception as exc:
        write_lane_completion(
            project_root,
            lane,
            status="failed",
            state=state,
            seed=seed,
            engine=engine,
            extra={
                "duration_seconds": max(0, int(time.time() - started)),
                "error": str(exc),
                "scratch_workspace": str(scratch_workspace) if scratch_workspace else "",
            },
        )
        raise


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

PHASES = [
    (1, "Phase 1: Vina Blind Docking (production)", phase1_vina),
    (2, f"Phase 2: PPI Global Blind Docking ({len(RECEPTOR_STATES) * PRODUCTION_N_SEEDS * 20}K models, {len(RECEPTOR_STATES)} states x {PRODUCTION_N_SEEDS} seeds)", phase2_ppi),
    (3, "Phase 3: PPI Postprocess (evidence extraction + aggregation)", phase3_ppi_postprocess),
    (4, "Phase 4: Vina Postprocess (전체)", phase4_vina_postprocess),
    (5, "Phase 5: Site Verdict (3축 통합)", phase5_verdict),
    (6, "Phase 6: Report 생성", phase6_report),
    (7, "Phase 7: Validate", phase7_validate),
]


def main():
    parser = argparse.ArgumentParser(description="EGFR-MYO1D Production Pipeline")
    parser.add_argument(
        "--lane",
        choices=[
            "vina-cpu",
            "ppi",
            "ppi-post",
            "vina-post",
            "finalize",
            "status",
            "vina-gpu",
            "phase3-gpu",
            # Advanced PPI-First pipeline lanes (Workflow B)
            "adv-phase1",
            "adv-phase2",
            "adv-phase3-setup",
            "adv-phase3-execute",
            "adv-phase3-post",
            "adv-phase4",
        ],
        default="",
        help="Scheduler lane entry point (recommended PBS interface)",
    )
    parser.add_argument("--state", type=str, default="",
                        help="Single receptor state selector for lane runs")
    parser.add_argument("--seed", type=int, default=None,
                        help="Single seed selector for lane runs")
    parser.add_argument("--allocated-cpus", type=int, default=None,
                        help="Scheduler-allocated CPUs used for safe worker caps")
    parser.add_argument("--engine", choices=["cpu-vina", "gnina"], default="cpu-vina",
                        help="Explicit docking engine for lane runs")
    parser.add_argument("--gpu-required", action="store_true",
                        help="Fail fast when a GPU lane has no visible devices")
    parser.add_argument("--scratch-dir", type=Path, default=None,
                        help="Optional scratch base directory for heavy lanes")
    parser.add_argument("--force", action="store_true",
                        help="전체 재실행 (기존 결과 무시)")
    parser.add_argument("--from", type=int, default=0, dest="from_phase",
                        help="지정 Phase부터 실행 (예: --from 4)")
    parser.add_argument("--only", type=str, default="",
                        help="지정 Phase만 실행 (예: --only 1,4,5,6,7)")
    parser.add_argument("--status", action="store_true",
                        help="각 Phase 완료 상태만 출력")
    parser.add_argument(
        "--disable-step-view",
        action="store_true",
        help="Derived step output view generation을 비활성화",
    )
    args = parser.parse_args()

    global _FORCE_MODE
    _FORCE_MODE = args.force

    if args.lane:
        _run_lane(
            args.lane,
            state=args.state or None,
            seed=args.seed,
            allocated_cpus=args.allocated_cpus,
            engine=args.engine,
            gpu_required=args.gpu_required,
            scratch_dir=args.scratch_dir,
        )
        return

    only_phases = set()
    if args.only:
        only_phases = {int(x.strip()) for x in args.only.split(",")}

    config = _load_config()
    step_view_enabled_flag = step_output_view_enabled(
        config,
        cli_disabled=args.disable_step_view,
    )
    vina_cfg = config.get("vina", {})
    ppi_models = f"{len(RECEPTOR_STATES) * PRODUCTION_N_SEEDS * 20}K"

    banner_width = 54
    vina_ppi_line = (f"Vina (exh={vina_cfg.get('exhaustiveness', '?')}, "
                     f"poses={vina_cfg.get('n_poses', '?')}) + "
                     f"PPI ({ppi_models})")
    print()
    print(f"╔{'═' * banner_width}╗")
    print(f"║  {'EGFR-MYO1D Production Pipeline':<{banner_width - 2}}║")
    print(f"║  {vina_ppi_line:<{banner_width - 2}}║")
    print(f"╚{'═' * banner_width}╝")

    fresh_run = not _project_root().exists()
    execution_mode = _execution_mode_label(args)
    print_status(config=config, step_view_enabled_flag=step_view_enabled_flag)

    if args.status:
        return

    run_started_at = utc_now_iso()
    run_id = f"{config.get('project_name', 'project')}-{run_started_at.replace(':', '').replace('-', '')}"
    phase_states = _empty_phase_states()
    if step_view_enabled_flag:
        _persist_run_status(
            config,
            execution_mode=execution_mode,
            run_id=run_id,
            started_at=run_started_at,
            phase_states=phase_states,
            overall_status="running",
        )

    pending_stale_steps = set(_selected_step_numbers(args, only_phases)) if step_view_enabled_flag else set()
    if step_view_enabled_flag and pending_stale_steps:
        refresh_root_step_views(
            CONFIG_PATH,
            repo_root=REPO_ROOT,
            execution_mode=execution_mode,
            fresh_run=fresh_run,
            stale_steps=sorted(pending_stale_steps),
            ppi_config_paths=[target["config_ini"] for target in PPI_TARGETS],
            pyrosetta_raw_run_paths=_pyrosetta_raw_run_paths(config),
        )

    t_start = time.time()

    for phase_num, name, func in PHASES:
        phase_state = _phase_state_entry(phase_states, phase_num)
        # --only: 지정된 Phase만 실행
        if only_phases and phase_num not in only_phases:
            print(f"\n  [SKIP] {name} (--only {args.only})")
            phase_state["status"] = "skipped"
            phase_state["skip_reason"] = f"Excluded by --only {args.only}"
            phase_state["completed_at"] = utc_now_iso()
            if step_view_enabled_flag:
                _persist_run_status(
                    config,
                    execution_mode=execution_mode,
                    run_id=run_id,
                    started_at=run_started_at,
                    phase_states=phase_states,
                    overall_status="running",
                )
            continue

        # --from: 지정 Phase 이전은 스킵
        if phase_num < args.from_phase:
            print(f"\n  [SKIP] {name} (--from {args.from_phase})")
            phase_state["status"] = "skipped"
            phase_state["skip_reason"] = f"Excluded by --from {args.from_phase}"
            phase_state["completed_at"] = utc_now_iso()
            if step_view_enabled_flag:
                _persist_run_status(
                    config,
                    execution_mode=execution_mode,
                    run_id=run_id,
                    started_at=run_started_at,
                    phase_states=phase_states,
                    overall_status="running",
                )
            continue

        # Phase 7 (Validate)은 항상 실행
        if phase_num < 7 and not args.force and phase_num != args.from_phase:
            if not (only_phases and phase_num in only_phases):
                check_fn = PHASE_CHECKS.get(phase_num, (None, None))[1]
                if check_fn:
                    missing = check_fn()
                    if not missing:
                        phase_state["status"] = "skipped"
                        phase_state["skip_reason"] = "Canonical outputs already exist"
                        phase_state["completed_at"] = utc_now_iso()
                        if step_view_enabled_flag:
                            _persist_run_status(
                                config,
                                execution_mode=execution_mode,
                                run_id=run_id,
                                started_at=run_started_at,
                                phase_states=phase_states,
                                overall_status="running",
                            )
                        print(f"\n  [SKIP] {name} — 결과물 이미 존재")
                        continue

        phase_state["status"] = "running"
        phase_state["started_at"] = utc_now_iso()
        phase_state["completed_at"] = ""
        phase_state["duration_seconds"] = None
        phase_state["skip_reason"] = ""
        phase_state["last_error"] = ""
        if step_view_enabled_flag:
            _persist_run_status(
                config,
                execution_mode=execution_mode,
                run_id=run_id,
                started_at=run_started_at,
                phase_states=phase_states,
                overall_status="running",
                current_phase_number=phase_num,
                current_phase_name=name,
            )

        phase_start_time = time.time()
        phase_success, phase_result, phase_error = run_step(name, func)
        phase_state["completed_at"] = utc_now_iso()
        phase_state["duration_seconds"] = max(0, int(time.time() - phase_start_time))
        phase_state["status"] = "completed" if phase_success else "failed"
        if phase_error is not None:
            phase_state["last_error"] = str(phase_error)
        # F-4.5: Sanity check core outputs after completion
        if phase_success and phase_num in PHASE_CORE_OUTPUTS:
            _sanity_check_phase_outputs(phase_num, config)
        if step_view_enabled_flag and (phase_success or phase_num == 7):
            pending_stale_steps.discard(phase_num)
            _refresh_step_view_outputs(
                phase_num,
                args,
                config,
                phase_result=phase_result,
                phase_error=phase_error,
                fresh_run=fresh_run,
                stale_steps=sorted(pending_stale_steps),
            )
        if step_view_enabled_flag:
            _persist_run_status(
                config,
                execution_mode=execution_mode,
                run_id=run_id,
                started_at=run_started_at,
                phase_states=phase_states,
                overall_status="running",
                last_error=phase_state["last_error"],
            )

    # 요약
    elapsed = time.time() - t_start
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)

    if step_view_enabled_flag:
        for entry in phase_states:
            if entry.get("status") == "pending":
                entry["status"] = "skipped"
                entry["skip_reason"] = "Phase was not executed in this run"
                entry["completed_at"] = utc_now_iso()
        failed_errors = [entry.get("last_error", "") for entry in phase_states if entry.get("status") == "failed"]
        _persist_run_status(
            config,
            execution_mode=execution_mode,
            run_id=run_id,
            started_at=run_started_at,
            phase_states=phase_states,
            overall_status="failed" if failed_errors else "completed",
            completed_at=utc_now_iso(),
            last_error=next((error for error in failed_errors if error), ""),
        )

    # Step-based output organization (F-4.3: collect post-run warnings)
    post_run_warnings: List[str] = []

    try:
        from egfr_pipeline.output_organizer import organize_outputs

        organize_outputs(str(CONFIG_PATH))
    except Exception as exc:
        post_run_warnings.append(f"Output organization failed: {exc}")
        print(f"  [WARN] Output organization failed: {exc}")

    _print_completion_summary(phase_states, config, hours, minutes)

    # F-4.3: Non-zero exit on post-run failures or phase failures
    has_failed_phases = any(
        entry.get("status") == "failed" for entry in phase_states
    )
    if post_run_warnings:
        print("\n  [WARN] 후처리 경고:")
        for w in post_run_warnings:
            print(f"    - {w}")
    if has_failed_phases or post_run_warnings:
        sys.exit(1)


if __name__ == "__main__":
    main()
