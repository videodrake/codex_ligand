#!/usr/bin/env python3
"""전체 파이프라인 1회 테스트 — 명령어 한 줄로 실행.

Usage:
    conda activate pyrosetta
    python run_full_test.py

기본 동작:
  1. PDBQT 변환 (receptor PDB → PDBQT, ligand SDF → PDBQT)
  2. Vina blind docking (exhaustiveness=8, n_poses=5 — 빠른 테스트용)
  3. Postprocess 전체 (parse → contacts → cluster → summarize → compare → bootstrap)
  4. Verdict (사이트 판정)
  5. Report (보고서 생성)
  6. Validate (출력 검증)

소요 시간: ~10-30분 (서버 성능에 따라)
"""

import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
CONFIG_TEMPLATE = REPO_ROOT / "config" / "example-project.yaml"

# ── 테스트 설정 (빠른 실행용) ──
TEST_EXHAUSTIVENESS = 8    # 프로덕션: 128, 테스트: 8
TEST_N_POSES = 5           # 프로덕션: 20, 테스트: 5
TEST_BOOTSTRAP_REPS = 10   # 프로덕션: 100, 테스트: 10


def banner(msg: str):
    width = 60
    print()
    print("=" * width)
    print(f"  {msg}")
    print("=" * width)


def run_step(name: str, func, *args, **kwargs):
    """Run a step with timing and error handling."""
    banner(name)
    t0 = time.time()
    try:
        result = func(*args, **kwargs)
        elapsed = time.time() - t0
        print(f"\n  [OK] {name} 완료 ({elapsed:.1f}초)")
        return result
    except Exception as e:
        elapsed = time.time() - t0
        print(f"\n  [FAIL] {name} 실패 ({elapsed:.1f}초)")
        print(f"  Error: {e}")
        return None


def check_prerequisites():
    """필수 도구 확인."""
    banner("Step 0: 환경 점검")

    # Python packages
    required = ["yaml", "numpy", "pandas", "scipy"]
    missing = []
    for pkg in required:
        try:
            __import__(pkg)
        except ImportError:
            missing.append(pkg)
    if missing:
        print(f"  [FAIL] 필수 패키지 없음: {missing}")
        print(f"  → pip install {' '.join(p if p != 'yaml' else 'pyyaml' for p in missing)}")
        return False

    # PDBQT conversion tools
    tools = ["mk_prepare_ligand.py", "prepare_receptor", "obabel"]
    found = []
    for tool in tools:
        if shutil.which(tool):
            found.append(tool)
    if not found:
        print("  [WARN] PDBQT 변환 도구를 찾을 수 없습니다:")
        print("         mk_prepare_ligand.py (meeko), prepare_receptor (ADFR), obabel (OpenBabel)")
        print("         → pip install meeko  또는  conda install -c conda-forge openbabel")
        print("         변환 도구 없이는 Step 1 (Vina 도킹)을 실행할 수 없습니다.")
    else:
        print(f"  변환 도구: {', '.join(found)}")

    # Input files
    receptors = list((REPO_ROOT / "input" / "receptors").glob("*.pdb"))
    ligands = list((REPO_ROOT / "input" / "ligands").glob("*.sdf"))
    print(f"  수용체 PDB: {len(receptors)}개")
    print(f"  리간드 SDF: {len(ligands)}개")

    if not receptors or not ligands:
        print("  [FAIL] input/receptors/*.pdb 또는 input/ligands/*.sdf 파일 없음")
        return False

    print("  [OK] 환경 점검 통과")
    return True


def create_test_config() -> Path:
    """테스트용 config 생성 (낮은 exhaustiveness)."""
    import yaml

    with open(CONFIG_TEMPLATE, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    # 테스트용 설정 오버라이드
    config["project_name"] = "full_test"
    config["output_root"] = "./output"
    config["vina"]["exhaustiveness"] = TEST_EXHAUSTIVENESS
    config["vina"]["n_poses"] = TEST_N_POSES
    # top-level mode도 동기화 (dock.py가 양쪽 다 체크)
    config["mode"] = config["vina"].get("mode", "blind")
    # top-level exhaustiveness/n_poses도 설정 (dock.py fallback용)
    config["exhaustiveness"] = TEST_EXHAUSTIVENESS
    config["n_poses"] = TEST_N_POSES
    config["bootstrap"] = {
        "n_replicates": TEST_BOOTSTRAP_REPS,
        "sample_fraction": 0.8,
        "seed": 42,
    }

    test_config = REPO_ROOT / "config" / "full_test.yaml"
    with open(test_config, "w", encoding="utf-8") as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True)

    print(f"  테스트 config: {test_config}")
    print(f"  exhaustiveness={TEST_EXHAUSTIVENESS}, n_poses={TEST_N_POSES}")
    print(f"  bootstrap replicates={TEST_BOOTSTRAP_REPS}")
    return test_config


def step1_convert_pdbqt():
    """PDB/SDF → PDBQT 변환."""
    from egfr_pipeline.vina.dock import prepare_receptor, prepare_ligand

    receptor_dir = REPO_ROOT / "input" / "receptors"
    ligand_dir = REPO_ROOT / "input" / "ligands"

    # Receptors
    for pdb in sorted(receptor_dir.glob("*.pdb")):
        pdbqt = receptor_dir / f"{pdb.stem}_receptor.pdbqt"
        if pdbqt.exists():
            print(f"  [Skip] {pdbqt.name} (이미 존재)")
        else:
            print(f"  [변환] {pdb.name} → {pdbqt.name}")
            if not prepare_receptor(pdb, pdbqt):
                raise RuntimeError(f"Receptor 변환 실패: {pdb.name}")

    # Ligands
    for sdf in sorted(ligand_dir.glob("*.sdf")):
        pdbqt = ligand_dir / f"{sdf.stem}.pdbqt"
        if pdbqt.exists():
            print(f"  [Skip] {pdbqt.name} (이미 존재)")
        else:
            print(f"  [변환] {sdf.name} → {pdbqt.name}")
            if not prepare_ligand(sdf, pdbqt):
                raise RuntimeError(f"Ligand 변환 실패: {sdf.name}")


def step2_vina_docking(config_path: Path):
    """Vina blind docking."""
    from egfr_pipeline.vina.dock import main as vina_main

    argv_backup = sys.argv[:]
    try:
        sys.argv = ["dock", "--config", str(config_path)]
        vina_main()
    finally:
        sys.argv = argv_backup


def step3_postprocess(config_path: str):
    """전체 후처리 체인."""
    from egfr_pipeline.config import load_config, project_root_from_config

    config = load_config(config_path)
    project_root = project_root_from_config(config)

    # parse
    print("\n  --- parse ---")
    from egfr_pipeline.vina.parse_poses import build_pose_table_from_config
    out = build_pose_table_from_config(config_path)
    print(f"  → {out}")

    # contacts
    print("\n  --- contacts ---")
    from egfr_pipeline.vina.contacts import enrich_pose_table_with_contacts
    cutoff = config.get("postprocess", {}).get("contact_cutoff", 4.0)
    out = enrich_pose_table_with_contacts(config_path, cutoff=cutoff)
    print(f"  → {out}")

    # cluster
    print("\n  --- cluster ---")
    from egfr_pipeline.vina.cluster import cluster_pose_table
    pp = config.get("postprocess", {})
    out = cluster_pose_table(
        config_path,
        cutoff=pp.get("pocket_cutoff", 4.0),
        merge_by_residue=pp.get("merge_by_residue", False),
    )
    print(f"  → {out}")

    # summarize
    print("\n  --- summarize ---")
    from egfr_pipeline.vina.summarize import summarize_from_config
    pocket_csv, drug_csv = summarize_from_config(config_path)
    print(f"  → {pocket_csv}")
    print(f"  → {drug_csv}")

    # compare
    print("\n  --- compare ---")
    from egfr_pipeline.vina.compare import compare_from_config
    out = compare_from_config(config_path)
    print(f"  → {out}")

    # bootstrap
    print("\n  --- bootstrap ---")
    from egfr_pipeline.vina.bootstrap import bootstrap_from_config
    out = bootstrap_from_config(config_path)
    print(f"  → {out}")


def step4_verdict(config_path: str):
    """사이트 판정."""
    from egfr_pipeline.verdict import generate_verdict
    agr_csv, ver_csv = generate_verdict(config_path)
    if ver_csv:
        print(f"  → {agr_csv}")
        print(f"  → {ver_csv}")
    else:
        print("  [WARN] pocket table 없음 — verdict 스킵")


def step5_report(config_path: str):
    """보고서 생성."""
    from egfr_pipeline.report import generate_report
    report_path, csv_path = generate_report(config_path)
    print(f"  → {report_path}")
    print(f"  → {csv_path}")


def step6_validate(config_path: str):
    """출력 검증."""
    from egfr_pipeline.validate import run_validation
    result = run_validation(config_path, repo_root=str(REPO_ROOT))
    result.print_summary()
    return result


def main():
    print()
    print("╔══════════════════════════════════════════════╗")
    print("║  EGFR-MYO1D Full Pipeline Test               ║")
    print("║  exhaustiveness=8, n_poses=5 (빠른 테스트)    ║")
    print("╚══════════════════════════════════════════════╝")

    t_start = time.time()

    # Step 0: 환경 점검
    if not check_prerequisites():
        print("\n환경 점검 실패. 위의 안내에 따라 설치 후 다시 실행하세요.")
        sys.exit(1)

    # 테스트 config 생성
    config_path = run_step("Step 0.5: 테스트 Config 생성", create_test_config)
    if not config_path:
        sys.exit(1)
    config_str = str(config_path)

    # Step 1: PDBQT 변환
    run_step("Step 1: PDBQT 변환 (PDB/SDF → PDBQT)", step1_convert_pdbqt)

    # Step 2: Vina 도킹
    run_step("Step 2: Vina Blind Docking", step2_vina_docking, config_path)

    # Step 3: 후처리
    run_step("Step 3: Postprocess (전체)", step3_postprocess, config_str)

    # Step 4: 판정
    run_step("Step 4: Site Verdict", step4_verdict, config_str)

    # Step 5: 보고서
    run_step("Step 5: Report 생성", step5_report, config_str)

    # Step 6: 검증
    run_step("Step 6: Validate", step6_validate, config_str)

    # 요약
    elapsed = time.time() - t_start
    minutes = int(elapsed // 60)
    seconds = int(elapsed % 60)

    banner("테스트 완료")
    print(f"  총 소요 시간: {minutes}분 {seconds}초")
    print(f"  출력 디렉토리: output/full_test/")
    print(f"  config: {config_path}")
    print()
    print("  확인할 파일:")
    print("    output/full_test/vina_pose_table.csv")
    print("    output/full_test/vina_pocket_table.csv")
    print("    output/full_test/valid_sites.csv")
    print("    output/full_test/project_report.txt")
    print()


if __name__ == "__main__":
    main()
