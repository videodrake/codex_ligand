#!/usr/bin/env python3
"""EGFR-MYO1D Docking Pipeline — Unified Interactive CLI.

Usage:
  python main.py              # Interactive menu
  python main.py --help       # Show all CLI options
  python main.py vina         # Run Vina docking directly
  python main.py postprocess  # Run Vina postprocessing
  python main.py report       # Generate report
  python main.py validate     # Validate outputs
  python main.py full         # Run full pipeline (vina → postprocess → report → validate)
"""

import argparse
import os
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = REPO_ROOT / "config" / "example-project.yaml"

MENU = """
╔══════════════════════════════════════════════════════════╗
║     EGFR-MYO1D Docking Pipeline  (통합 실행 메뉴)        ║
╠══════════════════════════════════════════════════════════╣
║                                                          ║
║  [1] Vina Docking            (AutoDock Vina 실행)        ║
║  [2] Vina Postprocess        (결과 파싱/클러스터링)      ║
║  [3] PPI Docking             (PyRosetta PPI 도킹)        ║
║      3a. PDB 준비 (dimer + MYO1D 합치기)                 ║
║      3b. 도킹 실행                                       ║
║      3c. 결과 원복 (chain 번호 정상화)                   ║
║  [4] MD Analysis             (GROMACS 분석)              ║
║  [5] Generate Report         (종합 보고서 생성)          ║
║  [6] Validate Outputs        (출력 검증)                 ║
║  [7] Full Pipeline           (1→2→5→6 자동 실행)         ║
║                                                          ║
║  [q] Quit                                                ║
║                                                          ║
╚══════════════════════════════════════════════════════════╝
"""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def ask_choice(prompt: str, choices: list, allow_multiple: bool = False) -> str:
    """Display numbered choices and return selected value(s)."""
    for i, c in enumerate(choices, 1):
        print(f"  [{i}] {c}")
    while True:
        raw = input(f"\n{prompt}: ").strip()
        if not raw:
            continue
        if allow_multiple:
            parts = [p.strip() for p in raw.replace(",", " ").split()]
            selected = []
            for p in parts:
                if p.isdigit() and 1 <= int(p) <= len(choices):
                    selected.append(choices[int(p) - 1])
            if selected:
                return selected
        else:
            if raw.isdigit() and 1 <= int(raw) <= len(choices):
                return choices[int(raw) - 1]
        print("  잘못된 입력입니다. 다시 선택해주세요.")


def ask_yes_no(prompt: str, default: bool = True) -> bool:
    suffix = " [Y/n]" if default else " [y/N]"
    raw = input(f"{prompt}{suffix}: ").strip().lower()
    if not raw:
        return default
    return raw in ("y", "yes", "ㅇ", "ㅇㅇ")


def ask_input(prompt: str, default: str = "") -> str:
    suffix = f" [{default}]" if default else ""
    raw = input(f"{prompt}{suffix}: ").strip()
    return raw if raw else default


def find_config() -> str:
    """Interactively find or select a config file."""
    configs = sorted(REPO_ROOT.glob("config/*.yaml")) + sorted(REPO_ROOT.glob("config/*.yml"))
    configs += sorted(REPO_ROOT.glob("smoke_test/config*.yaml"))

    if not configs:
        path = ask_input("Config 파일 경로를 입력하세요")
        return path

    print("\n사용 가능한 config 파일:")
    for i, c in enumerate(configs, 1):
        print(f"  [{i}] {c.relative_to(REPO_ROOT)}")
    print(f"  [{len(configs) + 1}] 직접 입력")

    while True:
        raw = input("\n선택: ").strip()
        if raw.isdigit():
            idx = int(raw)
            if 1 <= idx <= len(configs):
                return str(configs[idx - 1])
            if idx == len(configs) + 1:
                return ask_input("Config 파일 경로")
        print("  잘못된 입력입니다.")


# ---------------------------------------------------------------------------
# Task runners
# ---------------------------------------------------------------------------

def run_vina(config_path: str = None, **kwargs):
    """Run Vina docking."""
    print("\n" + "=" * 50)
    print("  [1] AutoDock Vina Docking")
    print("=" * 50)

    from egfr_pipeline.vina.dock import main as vina_main

    # Build sys.argv for the existing dock.py CLI
    argv_backup = sys.argv[:]
    try:
        if config_path:
            sys.argv = ["dock", "--config", config_path]
        else:
            # Let dock.py run its own interactive mode
            sys.argv = ["dock"]
        vina_main()
    finally:
        sys.argv = argv_backup


def run_postprocess(config_path: str = None):
    """Run Vina postprocessing chain."""
    print("\n" + "=" * 50)
    print("  [2] Vina Postprocess")
    print("=" * 50)

    if not config_path:
        config_path = find_config()

    from egfr_pipeline.config import load_config, project_root_from_config

    config = load_config(config_path)
    project_root = project_root_from_config(config)

    steps = [
        ("결과 파싱 (parse_poses)", "parse"),
        ("접촉 잔기 추출 (contacts)", "contacts"),
        ("포켓 클러스터링 (cluster)", "cluster"),
        ("포켓 요약 (summarize)", "summarize"),
        ("교차 비교 (compare)", "compare"),
        ("PPI 잔기 추출 (ppi)", "ppi"),
    ]

    print("\n후처리 단계:")
    for i, (name, _) in enumerate(steps, 1):
        print(f"  [{i}] {name}")
    print(f"  [a] 전체 실행")

    sel = input("\n실행할 단계 선택 (번호/a): ").strip().lower()
    if sel == "a":
        selected = [s[1] for s in steps]
    else:
        indices = [int(x) for x in sel.replace(",", " ").split() if x.isdigit()]
        selected = [steps[i - 1][1] for i in indices if 1 <= i <= len(steps)]

    if not selected:
        print("선택된 단계가 없습니다.")
        return

    receptors = config.get("receptors", [])
    receptor_ids = [r["id"] for r in receptors]

    for step in selected:
        print(f"\n--- {step} ---")
        try:
            _run_postprocess_step(step, config_path, config, project_root, receptor_ids)
        except Exception as e:
            print(f"  [ERROR] {step}: {e}")
            if not ask_yes_no("계속 진행하시겠습니까?"):
                return


def _run_postprocess_step(step, config_path, config, project_root, receptor_ids):
    """Execute a single postprocess step."""
    if step == "parse":
        from egfr_pipeline.vina.parse_poses import build_pose_table_from_config
        print("  Parsing Vina results...")
        out = build_pose_table_from_config(config_path)
        print(f"  → {out}")

    elif step == "contacts":
        from egfr_pipeline.vina.contacts import enrich_pose_table_with_contacts
        cutoff = config.get("postprocess", {}).get("contact_cutoff", 4.0)
        print(f"  Extracting contacts (cutoff={cutoff}Å)...")
        out = enrich_pose_table_with_contacts(config_path, cutoff=cutoff)
        print(f"  → {out}")

    elif step == "cluster":
        from egfr_pipeline.vina.cluster import cluster_pose_table
        cutoff = config.get("postprocess", {}).get("pocket_cutoff", 4.0)
        print(f"  Clustering pockets (cutoff={cutoff}Å)...")
        out = cluster_pose_table(config_path, cutoff=cutoff)
        print(f"  → {out}")

    elif step == "summarize":
        from egfr_pipeline.vina.summarize import summarize_from_config
        print("  Summarizing pockets...")
        pocket_csv, drug_csv = summarize_from_config(config_path)
        print(f"  → {pocket_csv}")
        print(f"  → {drug_csv}")

    elif step == "compare":
        from egfr_pipeline.vina.compare import compare_from_config
        print("  Comparing pockets across receptors...")
        out = compare_from_config(config_path)
        print(f"  → {out}")

    elif step == "ppi":
        from egfr_pipeline.ppi.pyrosetta_extract import extract_pyrosetta_batch
        print("  Extracting PPI residues...")
        extract_pyrosetta_batch(config_path)


def run_pyrosetta(config_ini: str = None):
    """Run PyRosetta PPI docking (sub-menu)."""
    print("\n" + "=" * 50)
    print("  [3] PPI Docking (PyRosetta)")
    print("=" * 50)

    print("\n  [a] PDB 준비 — dimer + MYO1D 도메인 합치기")
    print("  [b] 도킹 실행 — PyRosetta PPI global docking")
    print("  [c] 결과 원복 — chain 번호 정상화 (도킹 후)")
    print("  [d] PBS 스크립트 확인")

    sel = input("\n선택 [a/b/c/d]: ").strip().lower()

    if sel == "a":
        _ppi_prepare()
    elif sel == "b":
        _ppi_run_docking(config_ini)
    elif sel == "c":
        _ppi_restore()
    elif sel == "d":
        _ppi_show_pbs()
    else:
        print("잘못된 입력입니다.")


def _ppi_prepare():
    """Prepare dimer + partner PDB for PPI docking."""
    print("\n--- PDB 준비 (dimer + partner 합치기) ---")

    # Dimer PDB
    dimer_default = "smoke_test/input/original/3gt8_dimer_anp.pdb"
    dimer_path = ask_input("EGFR dimer PDB 경로", dimer_default)

    # Partner selection
    partner_dir = REPO_ROOT / "input" / "PPI"
    partners = sorted(partner_dir.glob("*.pdb"))
    if partners:
        print("\nMYO1D 도메인 파일:")
        for i, p in enumerate(partners, 1):
            print(f"  [{i}] {p.name}")
        raw = input("\n선택: ").strip()
        if raw.isdigit() and 1 <= int(raw) <= len(partners):
            partner_path = str(partners[int(raw) - 1])
        else:
            partner_path = ask_input("Partner PDB 경로")
    else:
        partner_path = ask_input("Partner PDB 경로")

    partner_name = Path(partner_path).stem.replace(" ", "_")
    output_path = f"input/PPI/prepared/EGFR_dimer_{partner_name}.pdb"

    from egfr_pipeline.ppi.prepare_dimer_pdb import prepare_dimer_partner
    prepare_dimer_partner(
        Path(dimer_path), Path(partner_path), Path(output_path),
        partner_name=partner_name,
    )


def _ppi_run_docking(config_ini: str = None):
    """Run PyRosetta PPI docking."""
    print("\n--- PPI 도킹 실행 ---")

    if not config_ini:
        inis = sorted(REPO_ROOT.glob("config/ppi_*.ini"))
        if not inis:
            inis = sorted(REPO_ROOT.glob("config/*.ini"))
        if inis:
            print("\nPPI config 파일:")
            for i, p in enumerate(inis, 1):
                print(f"  [{i}] {p.relative_to(REPO_ROOT)}")
            print(f"  [{len(inis) + 1}] 두 도메인 모두 (순차)")
            print(f"  [{len(inis) + 2}] 직접 입력")
            raw = input("\n선택: ").strip()
            if raw.isdigit() and 1 <= int(raw) <= len(inis):
                config_ini = str(inis[int(raw) - 1])
            elif raw.isdigit() and int(raw) == len(inis) + 1:
                config_ini = "__both__"
            else:
                config_ini = ask_input("Config .ini 경로")
        else:
            config_ini = ask_input("Config .ini 경로")

    print("\n실행 방식:")
    print("  [1] PBS 제출 (qsub) — HPC 서버 권장")
    print("  [2] 직접 실행 (foreground) — 테스트용")

    run_mode = input("\n선택 [1/2]: ").strip()

    if run_mode == "1":
        _ppi_submit_pbs(config_ini)
    elif run_mode == "2":
        if config_ini == "__both__":
            print("\n  ⚠ '두 도메인 모두'는 PBS 제출만 지원합니다.")
            _ppi_submit_pbs(config_ini)
            return
        print(f"\n  Config: {config_ini}")
        if ask_yes_no("실행하시겠습니까? (PyRosetta 필요)", default=False):
            from egfr_pipeline.pyrosetta_docking.pipeline_manager import main as pm_main
            sys.argv = ["pipeline_manager", config_ini]
            try:
                pm_main()
            finally:
                sys.argv = sys.argv[:1]


def _ppi_submit_pbs(config_ini: str = None):
    """Submit PPI docking job via qsub."""
    import shutil
    import subprocess as sp

    # Determine test vs prod PBS
    is_test = "test" in (config_ini or "")
    is_both = config_ini == "__both__"

    if is_both:
        print("\n  실행 규모:")
        print("  [1] 테스트 (1K 모델, ~3-4시간)")
        print("  [2] 프로덕션 (20K 모델, ~24-36시간)")
        scale = input("\n  선택 [1/2]: ").strip()
        is_test = scale != "2"

    pbs_path = REPO_ROOT / "config" / ("run_ppi_test.pbs" if is_test else "run_ppi_prod.pbs")
    log_prefix = "ppi_test" if is_test else "ppi_prod"

    if not shutil.which("qsub"):
        print("\n  ⚠ qsub를 찾을 수 없습니다. PBS가 설치된 서버에서 실행하세요.")
        print(f"  수동 실행: qsub {pbs_path}")
        return

    if is_both:
        cmd = ["qsub", "-v", "RUN_MODE=both", str(pbs_path)]
    else:
        cmd = ["qsub", "-v", f"CONFIG_FILE={config_ini}", str(pbs_path)]

    print(f"\n  PBS: {pbs_path.name}")
    print(f"  실행: {' '.join(cmd)}")

    if ask_yes_no("제출하시겠습니까?"):
        result = sp.run(cmd, capture_output=True, text=True, cwd=str(REPO_ROOT))
        if result.returncode == 0:
            job_id = result.stdout.strip()
            print(f"\n  제출 완료: {job_id}")
            print(f"  로그: {log_prefix}.o / {log_prefix}.e")
            print(f"  상태 확인: qstat {job_id}")
        else:
            print(f"\n  제출 실패:")
            print(f"    {result.stderr.strip()}")


def _ppi_restore():
    """Restore chain numbering after docking."""
    print("\n--- 결과 원복 (chain 번호 정상화) ---")

    print("\n원복 대상:")
    print("  [1] PDB 파일 (도킹 결과 구조)")
    print("  [2] CSV 파일 (final_ranking.csv 등)")
    sel = input("\n선택 [1/2]: ").strip()

    # Find mapping files
    mapping_dir = REPO_ROOT / "input" / "PPI" / "prepared"
    mappings = sorted(mapping_dir.glob("*_mapping.csv"))
    if mappings:
        print("\nMapping 파일:")
        for i, p in enumerate(mappings, 1):
            print(f"  [{i}] {p.name}")
        raw = input("\n선택: ").strip()
        if raw.isdigit() and 1 <= int(raw) <= len(mappings):
            mapping_path = str(mappings[int(raw) - 1])
        else:
            mapping_path = ask_input("Mapping CSV 경로")
    else:
        mapping_path = ask_input("Mapping CSV 경로")

    input_path = ask_input("입력 파일 경로")
    output_path = ask_input("출력 파일 경로")

    from egfr_pipeline.ppi.prepare_dimer_pdb import restore_chains, restore_csv

    if sel == "1":
        restore_chains(Path(input_path), Path(mapping_path), Path(output_path))
    elif sel == "2":
        restore_csv(Path(input_path), Path(mapping_path), Path(output_path))


def _ppi_show_pbs():
    """Show PBS script usage."""
    print("\n--- PBS 스크립트 사용법 ---")
    print()
    print("  [테스트] 1K 모델, ~3-4시간")
    print("  qsub config/run_ppi_test.pbs                                         # beta-meander")
    print("  qsub -v CONFIG_FILE=config/ppi_test_TH1.ini config/run_ppi_test.pbs  # TH1")
    print("  qsub -v RUN_MODE=both config/run_ppi_test.pbs                        # 둘 다")
    print()
    print("  [프로덕션] 20K 모델, ~24-36시간")
    print("  qsub config/run_ppi_prod.pbs                                         # beta-meander")
    print("  qsub -v CONFIG_FILE=config/ppi_prod_TH1.ini config/run_ppi_prod.pbs  # TH1")
    print("  qsub -v RUN_MODE=both config/run_ppi_prod.pbs                        # 둘 다")
    print()
    print("  # 로그 확인")
    print("  tail -f ppi_test.o   # 테스트")
    print("  tail -f ppi_prod.o   # 프로덕션")


def run_md():
    """Run MD (GROMACS) analysis."""
    print("\n" + "=" * 50)
    print("  [4] MD GROMACS Analysis")
    print("=" * 50)

    print("\nMD 분석 모듈:")
    print("  [1] GROMACS trajectory analysis")
    print("  [2] Ligand contact analysis")

    sel = input("\n선택: ").strip()

    if sel == "1":
        print("\n  GROMACS 분석을 위해 다음 파일이 필요합니다:")
        print("  - TPR 파일, XTC 파일, GRO 파일")
        print("\n  직접 실행:")
        print("  python -m egfr_pipeline.md.gromacs_analysis <args>")
        print("\n  서버 환경에서 GROMACS 및 MDAnalysis가 설치된 상태로 실행하세요.")

    elif sel == "2":
        print("\n  리간드 접촉 분석을 위해 다음 파일이 필요합니다:")
        print("  - TPR 파일, XTC 파일")
        print("\n  직접 실행:")
        print("  python -m egfr_pipeline.md.ligand_contacts <args>")
        print("\n  서버 환경에서 MDAnalysis가 설치된 상태로 실행하세요.")


def run_report(config_path: str = None):
    """Generate the combined report."""
    print("\n" + "=" * 50)
    print("  [5] Generate Report")
    print("=" * 50)

    if not config_path:
        config_path = find_config()

    from egfr_pipeline.report import generate_report

    output_dir = ask_input("출력 디렉토리 (기본: project root)", "")
    report_path, csv_path = generate_report(
        config_path,
        output_dir=output_dir if output_dir else None,
    )
    print(f"\n  Report: {report_path}")
    print(f"  CSV:    {csv_path}")


def run_validate(config_path: str = None):
    """Run output validation."""
    print("\n" + "=" * 50)
    print("  [6] Validate Outputs")
    print("=" * 50)

    if not config_path:
        config_path = find_config()

    from egfr_pipeline.validate import run_validation

    repo_root = ask_input("Repository root", str(REPO_ROOT))
    result = run_validation(config_path, repo_root=repo_root)
    result.print_summary()
    return result


def run_full(config_path: str = None):
    """Run full pipeline: Vina → Postprocess → Report → Validate."""
    print("\n" + "=" * 50)
    print("  [7] Full Pipeline (Vina → Postprocess → Report → Validate)")
    print("=" * 50)

    if not config_path:
        config_path = find_config()

    print("\n실행 순서:")
    print("  1. Vina Docking")
    print("  2. Postprocess (전체)")
    print("  3. Report 생성")
    print("  4. Output 검증")

    if not ask_yes_no("\n진행하시겠습니까?"):
        return

    # Step 1: Vina
    print("\n" + "━" * 40)
    print("Step 1/4: Vina Docking")
    print("━" * 40)
    run_vina(config_path)

    # Step 2: Postprocess all
    print("\n" + "━" * 40)
    print("Step 2/4: Postprocess")
    print("━" * 40)
    from egfr_pipeline.config import load_config, project_root_from_config
    config = load_config(config_path)
    project_root = project_root_from_config(config)
    receptor_ids = [r["id"] for r in config.get("receptors", [])]
    for step in ["parse", "contacts", "cluster", "summarize", "compare", "ppi"]:
        try:
            _run_postprocess_step(step, config_path, config, project_root, receptor_ids)
        except Exception as e:
            print(f"  [WARN] {step}: {e}")

    # Step 3: Report
    print("\n" + "━" * 40)
    print("Step 3/4: Generate Report")
    print("━" * 40)
    from egfr_pipeline.report import generate_report
    try:
        report_path, csv_path = generate_report(config_path)
        print(f"  Report: {report_path}")
    except Exception as e:
        print(f"  [WARN] Report generation: {e}")

    # Step 4: Validate
    print("\n" + "━" * 40)
    print("Step 4/4: Validate")
    print("━" * 40)
    from egfr_pipeline.validate import run_validation
    result = run_validation(config_path, repo_root=str(REPO_ROOT))
    result.print_summary()


# ---------------------------------------------------------------------------
# CLI argument parsing
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="EGFR-MYO1D Docking Pipeline — Unified CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", "-c", help="Project config file (.yaml)")

    sub = parser.add_subparsers(dest="command")

    sub.add_parser("vina", help="Run AutoDock Vina docking")
    sub.add_parser("postprocess", help="Run Vina postprocessing chain")
    sub.add_parser("pyrosetta", help="Run PyRosetta PPI docking")
    sub.add_parser("md", help="Run MD GROMACS analysis")
    sub.add_parser("report", help="Generate combined report")
    sub.add_parser("validate", help="Validate pipeline outputs")
    sub.add_parser("full", help="Run full pipeline (vina→postprocess→report→validate)")

    return parser


# ---------------------------------------------------------------------------
# Interactive menu loop
# ---------------------------------------------------------------------------

def interactive_menu():
    """Main interactive menu."""
    print(MENU)

    while True:
        sel = input("선택 [1-7, 3a/3b/3c, q]: ").strip().lower()

        if sel == "q":
            print("종료합니다.")
            break
        elif sel == "1":
            run_vina()
        elif sel == "2":
            run_postprocess()
        elif sel == "3":
            run_pyrosetta()
        elif sel == "3a":
            _ppi_prepare()
        elif sel == "3b":
            _ppi_run_docking()
        elif sel == "3c":
            _ppi_restore()
        elif sel == "4":
            run_md()
        elif sel == "5":
            run_report()
        elif sel == "6":
            run_validate()
        elif sel == "7":
            run_full()
        else:
            print("잘못된 입력입니다. 1-7 또는 q를 선택하세요.\n")
            continue

        # After task, show menu again
        print(MENU)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    # If no args, interactive mode
    if len(sys.argv) == 1:
        interactive_menu()
        return

    parser = build_parser()
    args = parser.parse_args()
    config = args.config

    dispatch = {
        "vina": lambda: run_vina(config),
        "postprocess": lambda: run_postprocess(config),
        "pyrosetta": lambda: run_pyrosetta(),
        "md": lambda: run_md(),
        "report": lambda: run_report(config),
        "validate": lambda: run_validate(config),
        "full": lambda: run_full(config),
    }

    if args.command in dispatch:
        dispatch[args.command]()
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
