#!/usr/bin/env python3
"""프로덕션 전체 파이프라인 — Vina + PPI + Verdict + Report.

Usage:
    conda activate pyrosetta
    python run_production.py            # 자동 이어하기 (완료된 Phase 스킵)
    python run_production.py --force    # 전체 재실행 (기존 결과 무시)
    python run_production.py --from 4   # Phase 4부터 실행 (이전 Phase 스킵)

전체 흐름:
  Phase 1: Vina blind docking (3 receptor × 3 ligand, exhaustiveness=128)  ~15분
           결과물: output/{project}//{receptor_id}/{ligand}_blind.pdbqt
  Phase 2: PPI docking (PyRosetta 20K models × 2 targets)                  ~24-36시간
           결과물: {pdb_stem}/final_result/final_ranking.csv
  Phase 3: PPI postprocess (chain restoration + residue extraction)
           결과물: output/{project}/ppi_pyrosetta_residues.csv
  Phase 4: Vina postprocess (parse → contacts → cluster → summarize → compare → bootstrap)
           결과물: output/{project}/vina_pocket_table.csv
  Phase 5: Verdict (3축 통합 scoring)
           결과물: output/{project}/valid_sites.csv
  Phase 6: Report
           결과물: output/{project}/project_report.txt
  Phase 7: Validate (항상 실행)
"""

import argparse
import csv
import sys
import time
from pathlib import Path
from typing import List, Optional

REPO_ROOT = Path(__file__).resolve().parent
CONFIG_PATH = REPO_ROOT / "config" / "example-project.yaml"

# PPI 설정
PPI_TARGETS = [
    {
        "name": "TH1",
        "config_ini": "config/ppi_prod_TH1.ini",
        "input_pdb": "input/PPI/prepared/EGFR_dimer_TH1_wt.pdb",
        "mapping_csv": "input/PPI/prepared/EGFR_dimer_TH1_mapping.csv",
        "receptor_id": "3GT8_raw",
        "partner_name": "MYO1D_TH1",
    },
    {
        "name": "beta-meander",
        "config_ini": "config/ppi_prod_beta_meander.ini",
        "input_pdb": "input/PPI/prepared/EGFR_dimer_beta_meander_wt.pdb",
        "mapping_csv": "input/PPI/prepared/EGFR_dimer_beta_meander_mapping.csv",
        "receptor_id": "3GT8_raw",
        "partner_name": "MYO1D_beta",
    },
]


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
        return result
    except Exception as e:
        elapsed = time.time() - t0
        minutes = int(elapsed // 60)
        seconds = int(elapsed % 60)
        print(f"\n  [FAIL] {name} 실패 ({minutes}분 {seconds}초)")
        print(f"  Error: {e}")
        import traceback
        traceback.print_exc()
        return None


def _load_config():
    import yaml
    with open(CONFIG_PATH, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _project_root() -> Path:
    config = _load_config()
    output_root = Path(config.get("output_root", "./output"))
    project_name = config.get("project_name", "")
    return output_root / project_name if project_name else output_root


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
    input_pdb = Path(target["input_pdb"])
    if not input_pdb.exists():
        input_pdb = Path(target["input_pdb"].replace("_wt.pdb", ".pdb"))
    return REPO_ROOT / input_pdb.stem


# ---------------------------------------------------------------------------
# Phase 완료 체크 함수
# ---------------------------------------------------------------------------

def check_phase1() -> List[str]:
    """Phase 1 (Vina): 모든 receptor×ligand .pdbqt 결과 존재 여부.
    결과물: output/{project}/{receptor_id}/{ligand}_blind.pdbqt
    """
    config = _load_config()
    project_root = _project_root()
    mode = config.get("mode", config.get("vina", {}).get("mode", "blind"))
    missing = []

    for rec in config.get("receptors", []):
        for lig in config.get("ligands", []):
            lig_name = Path(lig["pdbqt"]).stem.replace("_ligand", "")
            pdbqt = project_root / rec["id"] / f"{lig_name}_{mode}.pdbqt"
            if not pdbqt.exists():
                missing.append(f"{rec['id']}/{lig_name}_{mode}.pdbqt")

    return missing


def check_phase2() -> List[str]:
    """Phase 2 (PPI): 각 target의 final_ranking.csv 존재 여부.
    결과물: {pdb_stem}/final_result/final_ranking.csv
    """
    missing = []
    for target in PPI_TARGETS:
        docking_dir = _ppi_docking_dir(target)
        ranking = docking_dir / "final_result" / "final_ranking.csv" if docking_dir else None
        if not ranking or not ranking.exists():
            missing.append(f"{target['name']}: final_ranking.csv")
    return missing


def check_phase3() -> List[str]:
    """Phase 3 (PPI Postprocess): ppi_pyrosetta_residues.csv 존재+비어있지 않음.
    결과물: output/{project}/ppi_pyrosetta_residues.csv
    """
    path = _project_root() / "ppi_pyrosetta_residues.csv"
    if _csv_has_rows(path):
        return []
    return ["ppi_pyrosetta_residues.csv 없음 또는 비어있음"]


def check_phase4() -> List[str]:
    """Phase 4 (Vina Postprocess): 핵심 CSV들 존재+비어있지 않음.
    결과물: vina_pose_table.csv, vina_pocket_table.csv, vina_drug_pocket_map.csv,
           vina_pocket_comparison.csv, vina_pocket_bootstrap.csv
    """
    project_root = _project_root()
    required = [
        "vina_pose_table.csv",
        "vina_pocket_table.csv",
        "vina_drug_pocket_map.csv",
        "vina_pocket_comparison.csv",
        "vina_pocket_bootstrap.csv",
    ]
    missing = []
    for name in required:
        if not _csv_has_rows(project_root / name):
            missing.append(name)
    return missing


def check_phase5() -> List[str]:
    """Phase 5 (Verdict): valid_sites.csv 존재+비어있지 않음.
    결과물: valid_sites.csv, cross_method_agreement.csv
    """
    project_root = _project_root()
    missing = []
    for name in ["valid_sites.csv", "cross_method_agreement.csv"]:
        if not _csv_has_rows(project_root / name):
            missing.append(name)
    return missing


def check_phase6() -> List[str]:
    """Phase 6 (Report): project_report.txt 존재.
    결과물: project_report.txt, combined_residue_evidence.csv
    """
    project_root = _project_root()
    missing = []
    report = project_root / "project_report.txt"
    if not report.exists() or report.stat().st_size < 100:
        missing.append("project_report.txt")
    return missing


PHASE_CHECKS = {
    1: ("Vina blind docking 결과 (.pdbqt)", check_phase1),
    2: ("PPI docking 결과 (final_ranking.csv)", check_phase2),
    3: ("PPI postprocess (ppi_pyrosetta_residues.csv)", check_phase3),
    4: ("Vina postprocess (pocket_table 등 5개 CSV)", check_phase4),
    5: ("Verdict (valid_sites.csv)", check_phase5),
    6: ("Report (project_report.txt)", check_phase6),
}


def print_status():
    """각 Phase별 완료 상태 출력."""
    print("\n  Phase 상태 점검:")
    print(f"  {'Phase':<8} {'상태':<8} {'결과물':<50} {'상세'}")
    print(f"  {'─'*8} {'─'*8} {'─'*50} {'─'*20}")

    for phase_num, (desc, check_fn) in PHASE_CHECKS.items():
        missing = check_fn()
        if not missing:
            status = "[DONE]"
            detail = ""
        else:
            status = "[TODO]"
            detail = f"누락: {', '.join(missing[:3])}"
            if len(missing) > 3:
                detail += f" ...+{len(missing)-3}"
        print(f"  Phase {phase_num:<3} {status:<8} {desc:<50} {detail}")

    print(f"  Phase 7   [항상]   Validate (매 실행마다 검증)")
    print()


# ---------------------------------------------------------------------------
# Phase 실행 함수 (기존과 동일)
# ---------------------------------------------------------------------------

def phase1_vina():
    """Vina blind docking — receptor별 병렬 실행."""
    from concurrent.futures import ProcessPoolExecutor, as_completed

    config = _load_config()
    receptors = config.get("receptors", [])
    ligands = config.get("ligands", [])

    print(f"  Receptor {len(receptors)}개 × Ligand {len(ligands)}개 = {len(receptors) * len(ligands)} 도킹 잡")
    print(f"  exhaustiveness={config['vina']['exhaustiveness']}, n_poses={config['vina']['n_poses']}")
    print(f"  병렬 실행: {len(receptors)}개 프로세스")

    from run_full_test import _dock_one_receptor

    with ProcessPoolExecutor(max_workers=len(receptors)) as executor:
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


def phase2_ppi():
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

        # 이미 완료된 target은 스킵
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


def phase3_ppi_postprocess():
    """PPI 결과 chain restoration + residue extraction."""
    config_str = str(CONFIG_PATH)

    for target in PPI_TARGETS:
        name = target["name"]
        mapping_csv = target["mapping_csv"]
        receptor_id = target["receptor_id"]
        partner_name = target["partner_name"]

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
        )


def phase4_vina_postprocess():
    """Vina 후처리 전체 체인."""
    from egfr_pipeline.config import load_config
    config_str = str(CONFIG_PATH)
    config = load_config(config_str)

    print("\n  --- parse ---")
    from egfr_pipeline.vina.parse_poses import build_pose_table_from_config
    out = build_pose_table_from_config(config_str)
    print(f"  → {out}")

    print("\n  --- contacts ---")
    from egfr_pipeline.vina.contacts import enrich_pose_table_with_contacts
    cutoff = (config.get("postprocess") or {}).get("contact_cutoff", 4.0)
    out = enrich_pose_table_with_contacts(config_str, cutoff=cutoff)
    print(f"  → {out}")

    print("\n  --- cluster ---")
    from egfr_pipeline.vina.cluster import cluster_pose_table
    pp = config.get("postprocess") or {}
    out = cluster_pose_table(
        config_str,
        cutoff=pp.get("pocket_cutoff", 4.0),
        merge_by_residue=pp.get("merge_by_residue", False),
    )
    print(f"  → {out}")

    print("\n  --- summarize ---")
    from egfr_pipeline.vina.summarize import summarize_from_config
    pocket_csv, drug_csv = summarize_from_config(config_str)
    print(f"  → {pocket_csv}")
    print(f"  → {drug_csv}")

    print("\n  --- compare ---")
    from egfr_pipeline.vina.compare import compare_from_config
    out = compare_from_config(config_str)
    print(f"  → {out}")

    print("\n  --- bootstrap ---")
    from egfr_pipeline.vina.bootstrap import bootstrap_from_config
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

PHASES = [
    (1, "Phase 1: Vina Blind Docking (production)", phase1_vina),
    (2, "Phase 2: PPI Global Blind Docking (50K models)", phase2_ppi),
    (3, "Phase 3: PPI Postprocess (chain restore + extract)", phase3_ppi_postprocess),
    (4, "Phase 4: Vina Postprocess (전체)", phase4_vina_postprocess),
    (5, "Phase 5: Site Verdict (3축 통합)", phase5_verdict),
    (6, "Phase 6: Report 생성", phase6_report),
    (7, "Phase 7: Validate", phase7_validate),
]


def main():
    parser = argparse.ArgumentParser(description="EGFR-MYO1D Production Pipeline")
    parser.add_argument("--force", action="store_true",
                        help="전체 재실행 (기존 결과 무시)")
    parser.add_argument("--from", type=int, default=0, dest="from_phase",
                        help="지정 Phase부터 실행 (예: --from 4)")
    parser.add_argument("--only", type=str, default="",
                        help="지정 Phase만 실행 (예: --only 1,4,5,6,7)")
    parser.add_argument("--status", action="store_true",
                        help="각 Phase 완료 상태만 출력")
    args = parser.parse_args()

    only_phases = set()
    if args.only:
        only_phases = {int(x.strip()) for x in args.only.split(",")}

    config = _load_config()
    vina_cfg = config.get("vina", {})
    ppi_models = "50K"  # from ini files

    print()
    print("╔══════════════════════════════════════════════════════╗")
    print(f"║  EGFR-MYO1D Production Pipeline                     ║")
    print(f"║  Vina (exh={vina_cfg.get('exhaustiveness', '?')}, "
          f"poses={vina_cfg.get('n_poses', '?')}) + "
          f"PPI ({ppi_models}) ║")
    print("╚══════════════════════════════════════════════════════╝")

    print_status()

    if args.status:
        return

    t_start = time.time()

    for phase_num, name, func in PHASES:
        # --only: 지정된 Phase만 실행
        if only_phases and phase_num not in only_phases:
            print(f"\n  [SKIP] {name} (--only {args.only})")
            continue

        # --from: 지정 Phase 이전은 스킵
        if phase_num < args.from_phase:
            print(f"\n  [SKIP] {name} (--from {args.from_phase})")
            continue

        # Phase 7 (Validate)은 항상 실행
        if phase_num < 7 and not args.force and phase_num != args.from_phase:
            if not (only_phases and phase_num in only_phases):
                check_fn = PHASE_CHECKS.get(phase_num, (None, None))[1]
                if check_fn:
                    missing = check_fn()
                    if not missing:
                        print(f"\n  [SKIP] {name} — 결과물 이미 존재")
                        continue

        run_step(name, func)

    # 요약
    elapsed = time.time() - t_start
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)

    banner("프로덕션 완료")
    print(f"  총 소요 시간: {hours}시간 {minutes}분")
    print(f"  출력 디렉토리: {_project_root()}/")
    print(f"  config: {CONFIG_PATH}")
    print()
    print("  확인할 파일:")
    project = _project_root()
    print(f"    {project}/vina_pose_table.csv")
    print(f"    {project}/vina_pocket_table.csv")
    print(f"    {project}/valid_sites.csv")
    print(f"    {project}/project_report.txt")
    print()


if __name__ == "__main__":
    main()
