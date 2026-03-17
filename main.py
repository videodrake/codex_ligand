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
║  [7] Site Verdict            (유효 사이트 자동 판정)     ║
║  [8] PPI Postprocess         (PPI 후처리 자동화)         ║
║  [9] Full Pipeline           (1→2→7→5→6 자동 실행)       ║
║ [10] LightDock Validation    (Secondary PPI 검증)         ║
║ [11] Phase 2 Cascade         (Pocket 분석 재실행)         ║
║ [12] Organize Outputs        (결과물 Step별 정리)         ║
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
        ("부트스트랩 안정성 분석 (bootstrap)", "bootstrap"),
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
        pp = config.get("postprocess", {})
        cutoff = pp.get("pocket_cutoff", 8.0)
        max_iter = pp.get("cluster_max_iterations", 10)
        min_size = pp.get("min_pocket_size", 2)
        merge = pp.get("merge_by_residue", True)
        merge_j = pp.get("merge_jaccard", 0.3)
        merge_oc = pp.get("merge_overlap", 0.5)
        merge_cf = pp.get("merge_centroid_fallback", 6.0)
        print(f"  Clustering pockets (cutoff={cutoff}Å, merge={merge}, min_size={min_size})...")
        out = cluster_pose_table(
            config_path, cutoff=cutoff,
            max_iterations=max_iter,
            min_pocket_size=min_size,
            merge_by_residue=merge,
            merge_jaccard=merge_j,
            merge_overlap=merge_oc,
            merge_centroid_fallback=merge_cf,
        )
        print(f"  → {out}")

    elif step == "summarize":
        from egfr_pipeline.vina.summarize import summarize_from_config
        print("  Summarizing pockets...")
        pocket_csv, drug_csv, occupancy_csv = summarize_from_config(config_path)
        print(f"  → {pocket_csv}")
        print(f"  → {drug_csv}")
        print(f"  → {occupancy_csv}")

    elif step == "compare":
        from egfr_pipeline.vina.compare import compare_from_config
        print("  Comparing pockets across receptors...")
        out = compare_from_config(config_path)
        print(f"  → {out}")

    elif step == "ppi":
        from egfr_pipeline.ppi.pyrosetta_extract import extract_pyrosetta_batch
        print("  Extracting PPI residues...")
        extract_pyrosetta_batch(config_path)

    elif step == "bootstrap":
        from egfr_pipeline.vina.bootstrap import bootstrap_from_config
        bs = config.get("bootstrap", {})
        n = bs.get("n_replicates", 100)
        print(f"  Bootstrap stability analysis ({n} replicates)...")
        out = bootstrap_from_config(config_path)
        print(f"  → {out}")


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
        inis = sorted(REPO_ROOT.glob("config/ppi_*.ini")) + sorted(REPO_ROOT.glob("config/phase1/*.ini"))
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
    from egfr_pipeline.ppi.submit import (
        resolve_pbs_script, check_qsub_available, submit_pbs_job,
    )

    is_both = config_ini == "__both__"
    is_test = "test" in (config_ini or "")

    if is_both:
        print("\n  실행 규모:")
        print("  [1] 테스트 (1K 모델, ~3-4시간)")
        print("  [2] 프로덕션 (20K 모델, ~24-36시간)")
        scale = input("\n  선택 [1/2]: ").strip()
        is_test = scale != "2"

    pbs_path, log_prefix = resolve_pbs_script(config_ini, REPO_ROOT, is_test=is_test)

    if not check_qsub_available():
        print("\n  ⚠ qsub를 찾을 수 없습니다. PBS가 설치된 서버에서 실행하세요.")
        print(f"  수동 실행: qsub {pbs_path}")
        return

    run_mode = "both" if is_both else None
    cfg_arg = None if is_both else config_ini

    print(f"\n  PBS: {pbs_path.name}")
    if ask_yes_no("제출하시겠습니까?"):
        success, message = submit_pbs_job(
            pbs_path, config_ini=cfg_arg, run_mode=run_mode, cwd=REPO_ROOT,
        )
        if success:
            print(f"\n  제출 완료: {message}")
            print(f"  로그: {log_prefix}.o / {log_prefix}.e")
            print(f"  상태 확인: qstat {message}")
        else:
            print(f"\n  제출 실패:")
            print(f"    {message}")


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
    print("  [Phase 1 LightDock 검증]")
    print("  qsub config/run_lightdock.pbs                          # 전체 state")
    print("  qsub -v STATE=3GT8_raw config/run_lightdock.pbs        # 단일 state")
    print()
    print("  [Phase 1 LightDock 테스트]")
    print("  qsub config/run_lightdock_test.pbs                     # 3GT8_raw")
    print()
    print("  [프로덕션 전체 파이프라인]")
    print("  PRECHECK_JOB=$(qsub config/run_pre_qsub_checks.pbs)")
    print("  qsub -W depend=afterok:${PRECHECK_JOB} config/run_production.pbs")
    print()
    print("  # Phase 1 PyRosetta 설정: config/phase1/*.ini")
    print("  # 로그 확인")
    print("  tail -f lightdock.o")


def run_md(config_path: str = None):
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
    print(result.summary())
    return result


def run_verdict(config_path: str = None):
    """Run site verdict (automated pocket validation)."""
    print("\n" + "=" * 50)
    print("  [7] Site Verdict")
    print("=" * 50)

    if not config_path:
        config_path = find_config()

    from egfr_pipeline.verdict import generate_verdict

    output_dir = ask_input("출력 디렉토리 (기본: project root)", "")
    agreement_csv, verdict_csv = generate_verdict(
        config_path,
        output_dir=output_dir if output_dir else None,
    )

    if agreement_csv and verdict_csv:
        print(f"\n  Agreement: {agreement_csv}")
        print(f"  Verdict:   {verdict_csv}")
    else:
        print("\n  판정 실행 실패 — pocket table이 없습니다.")


def run_ppi_postprocess(config_path: str = None):
    """Run PPI post-processing automation."""
    print("\n" + "=" * 50)
    print("  [8] PPI Postprocess")
    print("=" * 50)

    if not config_path:
        config_path = find_config()

    from egfr_pipeline.ppi.postprocess_ppi import postprocess_ppi_results

    docking_dir = ask_input("도킹 결과 디렉토리 경로")
    if not docking_dir:
        print("  경로가 입력되지 않았습니다.")
        return

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

    receptor_id = ask_input("Receptor ID (예: 3GT8_raw)")
    partner_name = ask_input("Partner 이름 (예: beta_meander)", "")

    postprocess_ppi_results(
        config_path=config_path,
        docking_dir=docking_dir,
        mapping_csv=mapping_path,
        receptor_id=receptor_id,
        partner_name=partner_name,
    )


def run_lightdock():
    """Run LightDock secondary validation (sub-menu)."""
    print("\n" + "=" * 50)
    print("  [10] LightDock Secondary Validation")
    print("=" * 50)

    print("\n  [a] Setup — 실행 스크립트 생성 (전체 state)")
    print("  [b] Run — LightDock 실행 (LightDock 설치 필요)")
    print("  [c] Extract — 인터페이스 추출 + orientation filter")
    print("  [d] Convergence — PyRosetta 교차 검증")
    print("  [e] Full — a + c + d (setup + extract + convergence)")
    print("  [f] 가용성 확인 (LightDock 설치 상태)")
    print("  [g] PBS 스크립트 안내")

    sel = input("\n선택 [a-g]: ").strip().lower()

    from egfr_pipeline.phase1 import lightdock_validation as ldv

    if sel == "a":
        for state in ldv.RECEPTOR_STATES:
            print(f"\n{state}:")
            ldv.generate_lightdock_setup(state, ldv.PHASE1_OUTPUT_DIR)
    elif sel == "b":
        state = input("  State [all/3GT8_raw/EGFR_160-185/EGFR_170-200]: ").strip()
        states = ldv.RECEPTOR_STATES if state in ("", "all") else [state]
        for s in states:
            ldv.run_lightdock_scripts(s, ldv.PHASE1_OUTPUT_DIR)
    elif sel == "c":
        for state in ldv.RECEPTOR_STATES:
            print(f"\n{state}:")
            ldv.extract_lightdock_interfaces(state, ldv.PHASE1_OUTPUT_DIR)
    elif sel == "d":
        for state in ldv.RECEPTOR_STATES:
            print(f"\n{state}:")
            ldv.compute_cross_method_convergence(state, ldv.PHASE1_OUTPUT_DIR)
    elif sel == "e":
        sys.argv = ["lightdock_validation", "--all"]
        ldv.main()
    elif sel == "f":
        available, msg = ldv.check_lightdock_available()
        print(f"\n  {msg}")
    elif sel == "g":
        print("\n--- PBS 스크립트 사용법 ---")
        print()
        print("  [프로덕션] 400 swarms, ~12-24시간/state")
        print("  qsub config/run_lightdock.pbs                          # 전체 state")
        print("  qsub -v STATE=3GT8_raw config/run_lightdock.pbs        # 단일 state")
        print()
        print("  [테스트] 50 swarms, ~1-2시간")
        print("  qsub config/run_lightdock_test.pbs                     # 3GT8_raw")
        print("  qsub -v STATE=EGFR_160-185 config/run_lightdock_test.pbs")
        print()
        print("  # 로그 확인")
        print("  tail -f lightdock.o")
    else:
        print("잘못된 입력입니다.")


def _run_phase2_cmd(args):
    """Dispatch for `python main.py phase2` CLI subcommand."""
    from egfr_pipeline.phase2.rerun_cascade import run_phase2_cascade
    run_phase2_cascade(
        parse_only=args.parse_only,
        from_tg=args.from_tg,
        ftmap_dir=args.ftmap_dir,
    )


def run_phase2_cascade_menu():
    """Run Phase 2 cascade (TG 2.0 -> 2.7)."""
    print("\n" + "=" * 50)
    print("  [11] Phase 2 Cascade (Pocket Analysis)")
    print("=" * 50)

    print("\n  [a] 전체 실행 (setup + parse + 분석)")
    print("  [b] Parse-only (fpocket 완료 후 분석만)")
    print("  [c] 특정 TG부터 재실행")
    print("  [d] PBS 스크립트 안내")

    sel = input("\n선택 [a-d]: ").strip().lower()

    if sel == "d":
        print("\n--- PBS 스크립트 사용법 ---")
        print()
        print("  qsub config/run_phase2_cascade.pbs                       # 전체")
        print("  qsub -v SKIP_TOOLS=1 config/run_phase2_cascade.pbs       # fpocket 건너뛰기")
        print("  qsub -v FROM_TG=2.5 config/run_phase2_cascade.pbs        # 특정 TG부터")
        return

    from egfr_pipeline.phase2.rerun_cascade import run_phase2_cascade

    parse_only = sel in ("b", "c")
    from_tg = "2.0"
    ftmap_dir = None

    if sel == "c":
        from_tg = ask_input("시작 TG (예: 2.3)", "2.0")

    ftmap_raw = ask_input("FTMap 디렉토리 (없으면 Enter)", "")
    if ftmap_raw:
        ftmap_dir = Path(ftmap_raw)

    run_phase2_cascade(
        parse_only=parse_only,
        from_tg=from_tg,
        ftmap_dir=ftmap_dir,
    )


def run_organize(config_path: str = None):
    """Organize outputs into step-based directories for easy browsing."""
    print("\n" + "=" * 50)
    print("  Output Organization (Step-based)")
    print("=" * 50)

    if not config_path:
        config_path = find_config()

    from egfr_pipeline.organize import organize_outputs

    organize_outputs(config_path, repo_root=REPO_ROOT)


def run_full(config_path: str = None):
    """Run full pipeline: Vina → Postprocess → Verdict → Report → Validate."""
    print("\n" + "=" * 50)
    print("  [9] Full Pipeline (Vina → Postprocess → Verdict → Report → Validate)")
    print("=" * 50)

    if not config_path:
        config_path = find_config()

    print("\n실행 순서:")
    print("  1. Vina Docking")
    print("  2. Postprocess (전체)")
    print("  3. Site Verdict (유효 사이트 판정)")
    print("  4. Report 생성")
    print("  5. Output 검증")

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

    # Step 3: Verdict
    print("\n" + "━" * 40)
    print("Step 3/5: Site Verdict")
    print("━" * 40)
    try:
        from egfr_pipeline.verdict import generate_verdict
        agr_csv, ver_csv = generate_verdict(config_path)
        if ver_csv:
            print(f"  Verdict: {ver_csv}")
    except Exception as e:
        print(f"  [WARN] Verdict: {e}")

    # Step 4: Report
    print("\n" + "━" * 40)
    print("Step 4/5: Generate Report")
    print("━" * 40)
    from egfr_pipeline.report import generate_report
    try:
        report_path, csv_path = generate_report(config_path)
        print(f"  Report: {report_path}")
    except Exception as e:
        print(f"  [WARN] Report generation: {e}")

    # Step 5: Validate
    print("\n" + "━" * 40)
    print("Step 5/5: Validate")
    print("━" * 40)
    from egfr_pipeline.validate import run_validation
    result = run_validation(config_path, repo_root=str(REPO_ROOT))
    print(result.summary())


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
    sub.add_parser("verdict", help="Run site verdict (automated pocket validation)")
    sub.add_parser("ppi-postprocess", help="Run PPI post-processing automation")
    sub.add_parser("lightdock", help="LightDock secondary validation")
    sub.add_parser("full", help="Run full pipeline (vina→postprocess→verdict→report→validate)")
    sub.add_parser("organize", help="Organize outputs into step-based directories for easy browsing")
    p2 = sub.add_parser("phase2", help="Run Phase 2 cascade (pocket analysis TG 2.0→2.7)")
    p2.add_argument("--parse-only", action="store_true",
                     help="Skip fpocket/P2Rank setup (already executed)")
    p2.add_argument("--from-tg", default="2.0",
                     help="Start from this TG (e.g. '2.5')")
    p2.add_argument("--ftmap-dir", type=Path, default=None,
                     help="FTMap output directory")

    return parser


# ---------------------------------------------------------------------------
# Interactive menu loop
# ---------------------------------------------------------------------------

def interactive_menu():
    """Main interactive menu."""
    print(MENU)

    while True:
        sel = input("선택 [1-12, 3a/3b/3c, q]: ").strip().lower()

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
            run_verdict()
        elif sel == "8":
            run_ppi_postprocess()
        elif sel == "9":
            run_full()
        elif sel == "10":
            run_lightdock()
        elif sel == "11":
            run_phase2_cascade_menu()
        elif sel == "12":
            run_organize()
        else:
            print("잘못된 입력입니다. 1-12 또는 q를 선택하세요.\n")
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
        "pyrosetta": lambda: run_pyrosetta(config),
        "md": lambda: run_md(config),
        "report": lambda: run_report(config),
        "validate": lambda: run_validate(config),
        "verdict": lambda: run_verdict(config),
        "ppi-postprocess": lambda: run_ppi_postprocess(config),
        "lightdock": lambda: run_lightdock(),
        "full": lambda: run_full(config),
        "organize": lambda: run_organize(config),
        "phase2": lambda: _run_phase2_cmd(args),
    }

    if args.command in dispatch:
        dispatch[args.command]()
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
