#!/usr/bin/env python3
"""프로덕션 전체 파이프라인 — Vina + PPI + Verdict + Report.

Usage:
    conda activate pyrosetta
    python run_production.py

전체 흐름:
  Phase 1: Vina blind docking (3 receptor × 3 ligand, exhaustiveness=128)  ~15분
  Phase 2: PPI docking (PyRosetta 20K models × 2 targets)                  ~24-36시간
  Phase 3: PPI postprocess (chain restoration + residue extraction)
  Phase 4: Vina postprocess (parse → contacts → cluster → summarize → compare → bootstrap)
  Phase 5: Verdict (3축 통합 scoring)
  Phase 6: Report + Validate
"""

import sys
import time
from pathlib import Path

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


# ── Phase 1: Vina Docking ──

def phase1_vina():
    """Vina blind docking — receptor별 병렬 실행."""
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import yaml

    with open(CONFIG_PATH, encoding="utf-8") as f:
        config = yaml.safe_load(f)

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


# ── Phase 2: PPI Docking ──

def phase2_ppi():
    """PyRosetta PPI global blind docking — 순차 실행 (CPU 집약적)."""
    for target in PPI_TARGETS:
        name = target["name"]
        config_ini = target["config_ini"]
        input_pdb = target["input_pdb"]

        # _wt.pdb가 없으면 원본 .pdb 사용
        pdb_path = Path(input_pdb)
        if not pdb_path.exists():
            alt_path = Path(input_pdb.replace("_wt.pdb", ".pdb"))
            if alt_path.exists():
                pdb_path = alt_path
                print(f"  [INFO] {name}: _wt.pdb 없음, {alt_path.name} 사용")
            else:
                print(f"  [ERROR] {name}: 입력 PDB 없음: {input_pdb}")
                continue

        print(f"\n  --- {name} (20K models) ---")
        print(f"  Config: {config_ini}")
        print(f"  Input:  {pdb_path}")

        from egfr_pipeline.pyrosetta_docking.pipeline_manager import main as ppi_main
        sys.argv = ["pipeline_manager", config_ini, str(pdb_path)]
        ppi_main()


# ── Phase 3: PPI Postprocess ──

def phase3_ppi_postprocess():
    """PPI 결과 chain restoration + residue extraction."""
    config_str = str(CONFIG_PATH)

    for target in PPI_TARGETS:
        name = target["name"]
        mapping_csv = target["mapping_csv"]
        receptor_id = target["receptor_id"]
        partner_name = target["partner_name"]

        # PPI 도킹 결과 디렉토리 찾기
        # PyRosetta는 input PDB 이름 기반으로 출력 디렉토리 생성
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


# ── Phase 4: Vina Postprocess ──

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


# ── Phase 5: Verdict ──

def phase5_verdict():
    """3축 증거 통합 판정."""
    from egfr_pipeline.verdict import generate_verdict
    agr_csv, ver_csv = generate_verdict(str(CONFIG_PATH))
    if ver_csv and str(ver_csv) != ".":
        print(f"  → {agr_csv}")
        print(f"  → {ver_csv}")
    else:
        print("  [WARN] pocket table 없음 — verdict 스킵")


# ── Phase 6: Report + Validate ──

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


def main():
    print()
    print("╔══════════════════════════════════════════════════════╗")
    print("║  EGFR-MYO1D Production Pipeline                     ║")
    print("║  Vina (128 exh) + PPI (20K) + Verdict + Report      ║")
    print("╚══════════════════════════════════════════════════════╝")

    t_start = time.time()

    # Phase 1: Vina (빠름, ~15분)
    run_step("Phase 1: Vina Blind Docking (production)", phase1_vina)

    # Phase 2: PPI (느림, ~24-36시간)
    run_step("Phase 2: PPI Global Blind Docking (20K models)", phase2_ppi)

    # Phase 3: PPI Postprocess
    run_step("Phase 3: PPI Postprocess (chain restore + extract)", phase3_ppi_postprocess)

    # Phase 4: Vina Postprocess
    run_step("Phase 4: Vina Postprocess (전체)", phase4_vina_postprocess)

    # Phase 5: Verdict
    run_step("Phase 5: Site Verdict (3축 통합)", phase5_verdict)

    # Phase 6: Report
    run_step("Phase 6: Report 생성", phase6_report)

    # Phase 7: Validate
    run_step("Phase 7: Validate", phase7_validate)

    # 요약
    elapsed = time.time() - t_start
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)

    banner("프로덕션 완료")
    print(f"  총 소요 시간: {hours}시간 {minutes}분")
    print(f"  출력 디렉토리: output/egfr_myo1d_vina/")
    print(f"  config: {CONFIG_PATH}")
    print()
    print("  확인할 파일:")
    print("    output/egfr_myo1d_vina/vina_pose_table.csv")
    print("    output/egfr_myo1d_vina/vina_pocket_table.csv")
    print("    output/egfr_myo1d_vina/valid_sites.csv")
    print("    output/egfr_myo1d_vina/project_report.txt")
    print()


if __name__ == "__main__":
    main()
