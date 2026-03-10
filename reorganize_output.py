#!/usr/bin/env python3
"""프로덕션 완료 후 출력 폴더 재구성 스크립트.

기존 플랫 구조를 계층적 구조로 정리한다.
원본 파일을 이동(mv)하며, 심볼릭 링크로 하위 호환을 유지한다.

Usage:
    python reorganize_output.py --dry-run            # 미리보기 (필수 먼저 실행)
    python reorganize_output.py                      # 실제 실행
    python reorganize_output.py --project-dir output/my_project
    python reorganize_output.py --no-symlinks        # 심볼릭 링크 없이

정리 전:
    output/{project}/
    ├── 3GT8_raw/173940_blind.pdbqt       ← Vina + 분석 CSV 혼재
    ├── vina_pose_table.csv
    ├── valid_sites.csv
    └── ...
    EGFR_dimer_TH1_wt/                    ← PPI 결과가 repo 루트에
    EGFR_dimer_beta_meander_wt/

정리 후:
    output/{project}/
    ├── results/                          ← ★ 최종 결과
    ├── vina/                             ← Vina 세부
    │   └── docking/{receptor_id}/
    ├── ppi/                              ← PPI 세부
    │   ├── TH1/
    │   └── beta_meander/
    └── validation/
"""

import argparse
import os
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

# ── PPI 타겟 정의 (run_production.py와 동기화) ──────────────────────
PPI_TARGETS = [
    {"name": "TH1", "pdb_stem": "EGFR_dimer_TH1_wt"},
    {"name": "beta_meander", "pdb_stem": "EGFR_dimer_beta_meander_wt"},
]

# ── 파일 분류 ────────────────────────────────────────────────────────

RESULTS_FILES = [
    "valid_sites.csv",
    "project_report.txt",
    "combined_residue_evidence.csv",
    "cross_method_agreement.csv",
    "vina_consensus_sites.csv",        # 조건부 (합의 부위가 있을 때만 생성)
]

VINA_CSV_FILES = [
    "vina_pose_table.csv",
    "vina_pocket_table.csv",
    "vina_drug_pocket_map.csv",
    "vina_pocket_comparison.csv",
    "vina_pocket_bootstrap.csv",
    "vina_pocket_cutoff_sweep.csv",    # 선택적 (sweep 실행 시만)
]

RECEPTOR_IDS = ["3GT8_raw", "3GT8_cl38_48", "3GT8_cl85_100"]

PPI_CSV_FILES = [
    "ppi_pyrosetta_residues.csv",
    "ppi_pyrosetta_summary.csv",
    "ppi_afm_residues.csv",            # 선택적 (AFM 데이터 있을 때만)
    "ppi_afm_summary.csv",             # 선택적
]


# ── 안전한 유틸리티 함수 ─────────────────────────────────────────────

def safe_rel(path: Path, base: Path) -> str:
    """base 기준 상대 경로. base 하위가 아니면 절대 경로 반환."""
    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)


class Reorganizer:
    def __init__(self, project_dir: Path, dry_run: bool = False,
                 create_symlinks: bool = True):
        self.pd = project_dir.resolve()
        self.dry_run = dry_run
        self.create_symlinks = create_symlinks
        self.n_moved = 0
        self.n_skipped = 0
        self.errors = []

    # ── 로깅 ─────────────────────────────────────────────────────

    def _log(self, action: str, src: str, dst: str = ""):
        tag = "[DRY-RUN] " if self.dry_run else ""
        if dst:
            print(f"  {tag}{action}: {src} → {dst}")
        else:
            print(f"  {tag}{action}: {src}")

    def _warn(self, msg: str):
        print(f"  [WARN] {msg}")

    def _error(self, msg: str):
        print(f"  [ERROR] {msg}")
        self.errors.append(msg)

    # ── 안전한 파일 연산 ─────────────────────────────────────────

    def _mkdir(self, path: Path):
        """디렉토리 생성 (이미 존재하면 무시)."""
        if path.exists():
            return
        self._log("mkdir", safe_rel(path, self.pd))
        if not self.dry_run:
            path.mkdir(parents=True, exist_ok=True)

    def _safe_move_file(self, src: Path, dst: Path):
        """단일 파일 이동. src 없으면 skip, dst 있으면 error."""
        # src가 심볼릭 링크이면 건너뜀 (이미 정리된 상태)
        if src.is_symlink():
            self.n_skipped += 1
            return False
        if not src.exists():
            self.n_skipped += 1
            return False
        if not src.is_file():
            self._error(f"파일이 아님 (디렉토리?): {safe_rel(src, self.pd)}")
            return False
        if dst.exists() and not dst.is_symlink():
            self._error(f"대상 이미 존재: {safe_rel(dst, self.pd)}")
            return False

        self._mkdir(dst.parent)
        self._log("move", safe_rel(src, self.pd), safe_rel(dst, self.pd))
        if not self.dry_run:
            # dst가 심볼릭 링크이면 제거 후 이동
            if dst.is_symlink():
                dst.unlink()
            shutil.move(str(src), str(dst))
        self.n_moved += 1
        return True

    def _safe_move_dir(self, src: Path, dst: Path, src_label: str = ""):
        """디렉토리 이동. src 없으면 skip, dst 있으면 error."""
        label = src_label or safe_rel(src, self.pd)

        if src.is_symlink():
            self.n_skipped += 1
            return False
        if not src.exists():
            self.n_skipped += 1
            return False
        if not src.is_dir():
            self._error(f"디렉토리가 아님: {label}")
            return False
        if dst.exists() and not dst.is_symlink():
            self._error(f"대상 디렉토리 이미 존재: {safe_rel(dst, self.pd)}")
            return False

        self._mkdir(dst.parent)
        self._log("move_dir", label, safe_rel(dst, self.pd))
        if not self.dry_run:
            if dst.is_symlink():
                dst.unlink()
            shutil.move(str(src), str(dst))
        self.n_moved += 1
        return True

    def _make_symlink(self, link_path: Path, target: Path, link_label: str = ""):
        """심볼릭 링크 생성. 이미 존재하면 건너뜀."""
        if not self.create_symlinks:
            return
        label = link_label or safe_rel(link_path, self.pd)

        if link_path.exists() or link_path.is_symlink():
            return  # 이미 존재 → 조용히 건너뜀

        try:
            rel_target = os.path.relpath(target, link_path.parent)
        except ValueError:
            rel_target = str(target)

        self._log("symlink", label, rel_target)
        if not self.dry_run:
            try:
                link_path.symlink_to(rel_target)
            except OSError as e:
                self._warn(f"심볼릭 링크 실패: {e}")

    # ── 사전 검증 ────────────────────────────────────────────────

    def _preflight_check(self) -> bool:
        """실행 전 상태 검증. 문제 있으면 False 반환."""
        pd = self.pd

        if not pd.exists():
            self._error(f"프로젝트 디렉토리 없음: {pd}")
            return False

        # 최소 1개 이상의 대상 파일이 존재하는지
        has_any = False
        for fname in RESULTS_FILES + VINA_CSV_FILES + PPI_CSV_FILES:
            p = pd / fname
            if p.exists() and not p.is_symlink():
                has_any = True
                break
        for rec_id in RECEPTOR_IDS:
            p = pd / rec_id
            if p.exists() and p.is_dir() and not p.is_symlink():
                has_any = True
                break
        for target in PPI_TARGETS:
            p = REPO_ROOT / target["pdb_stem"]
            if p.exists() and p.is_dir() and not p.is_symlink():
                has_any = True
                break

        if not has_any:
            # results/ 등이 이미 존재하면 이미 정리된 것
            if (pd / "results").is_dir() or (pd / "vina").is_dir():
                print("[SKIP] 이미 정리된 구조이거나 이동할 파일이 없습니다.")
                return False
            self._warn("이동할 파일이 없습니다. 프로덕션이 완료되었는지 확인하세요.")
            return False

        return True

    # ── 메인 실행 ────────────────────────────────────────────────

    def run(self) -> bool:
        pd = self.pd

        print(f"\n프로젝트 디렉토리: {pd}")
        print(f"모드: {'미리보기 (dry-run)' if self.dry_run else '실행'}")
        print(f"심볼릭 링크: {'생성' if self.create_symlinks else '생성 안 함'}")
        print()

        if not self._preflight_check():
            return len(self.errors) == 0

        # ── Step 1: 디렉토리 골격 생성 ────────────────────────────
        print("Step 1. 디렉토리 구조 생성")
        for subdir in ["results", "vina", "vina/docking", "ppi", "validation"]:
            self._mkdir(pd / subdir)

        # ── Step 2: 핵심 결과 → results/ ──────────────────────────
        print("\nStep 2. 핵심 결과 → results/")
        for fname in RESULTS_FILES:
            src = pd / fname
            dst = pd / "results" / fname
            if self._safe_move_file(src, dst):
                self._make_symlink(src, dst)

        # ── Step 3: Vina CSV → vina/ ─────────────────────────────
        print("\nStep 3. Vina 분석 CSV → vina/")
        for fname in VINA_CSV_FILES:
            src = pd / fname
            dst = pd / "vina" / fname
            if self._safe_move_file(src, dst):
                self._make_symlink(src, dst)

        # ── Step 4: 수용체 PDBQT 디렉토리 → vina/docking/ ────────
        print("\nStep 4. Vina 도킹 결과 → vina/docking/")
        for rec_id in RECEPTOR_IDS:
            src = pd / rec_id
            dst = pd / "vina" / "docking" / rec_id
            if self._safe_move_dir(src, dst):
                self._make_symlink(src, dst)

        # ── Step 5: PPI CSV → ppi/ ───────────────────────────────
        print("\nStep 5. PPI 잔기 CSV → ppi/")
        for fname in PPI_CSV_FILES:
            src = pd / fname
            dst = pd / "ppi" / fname
            if self._safe_move_file(src, dst):
                self._make_symlink(src, dst)

        # ── Step 6: PPI 도킹 결과 → ppi/{name}/ ─────────────────
        print("\nStep 6. PPI 도킹 결과 → ppi/")
        for target in PPI_TARGETS:
            pdb_stem = target["pdb_stem"]
            name = target["name"]
            dst = pd / "ppi" / name

            # 경로 후보: repo 루트 또는 project_dir 안
            src_repo = REPO_ROOT / pdb_stem
            src_proj = pd / pdb_stem

            # 실제 존재하는 쪽을 선택
            if src_repo.exists() and src_repo.is_dir() and not src_repo.is_symlink():
                src = src_repo
                src_label = pdb_stem + "/ (repo root)"
            elif src_proj.exists() and src_proj.is_dir() and not src_proj.is_symlink():
                src = src_proj
                src_label = pdb_stem + "/ (project dir)"
            else:
                self._warn(f"PPI 결과 없음: {pdb_stem} (아직 Phase 2 미완료?)")
                self.n_skipped += 1
                continue

            if self._safe_move_dir(src, dst, src_label=src_label):
                # repo 루트에 있었으면 거기에 심볼릭 링크
                if src == src_repo:
                    self._make_symlink(
                        src_repo, dst,
                        link_label=f"{pdb_stem}/ (repo root symlink)",
                    )
                # project_dir에 있었으면 거기에 심볼릭 링크
                if src == src_proj:
                    self._make_symlink(src_proj, dst)

        # ── Step 7: relaxed_cache → ppi/relaxed_cache ─────────────
        print("\nStep 7. relaxed_cache 이동")
        cache_src = REPO_ROOT / "relaxed_cache"
        cache_dst = pd / "ppi" / "relaxed_cache"
        if cache_src.exists() and cache_src.is_dir() and not cache_src.is_symlink():
            if self._safe_move_dir(cache_src, cache_dst,
                                   src_label="relaxed_cache/ (repo root)"):
                self._make_symlink(
                    cache_src, cache_dst,
                    link_label="relaxed_cache/ (repo root symlink)",
                )
        else:
            self.n_skipped += 1

        # ── Step 8: README 생성 ──────────────────────────────────
        print("\nStep 8. README 생성")
        readme_path = pd / "results" / "README.txt"
        if not readme_path.exists():
            self._log("create", "results/README.txt")
            if not self.dry_run:
                readme_path.write_text(
                    "EGFR-MYO1D 파이프라인 결과\n"
                    "========================\n\n"
                    "핵심 결과 파일:\n"
                    "  valid_sites.csv               포켓별 증거 판정 (STRONG/MODERATE/WEAK)\n"
                    "  project_report.txt             통합 분석 리포트\n"
                    "  cross_method_agreement.csv     Vina-PPI 교차 검증\n"
                    "  combined_residue_evidence.csv  잔기별 다중 증거 통합\n\n"
                    "세부 데이터:\n"
                    "  ../vina/                       Vina 도킹 + 클러스터링 상세\n"
                    "  ../ppi/                        PyRosetta PPI 도킹 상세\n"
                    "  ../validation/                 출력 검증 결과\n",
                    encoding="utf-8",
                )

        # ── 요약 ─────────────────────────────────────────────────
        print(f"\n{'='*50}")
        print(f"  이동: {self.n_moved}건")
        print(f"  건너뜀: {self.n_skipped}건 (없거나 이미 처리)")
        print(f"  오류: {len(self.errors)}건")
        if self.errors:
            print("\n  오류 목록:")
            for e in self.errors:
                print(f"    - {e}")
        print(f"{'='*50}")

        if self.dry_run:
            print("\n※ --dry-run 모드입니다. 실제 파일은 이동되지 않았습니다.")
            print("  실행하려면: python reorganize_output.py")

        return len(self.errors) == 0


def main():
    parser = argparse.ArgumentParser(
        description="프로덕션 출력 폴더 재구성",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "사용법:\n"
            "  python reorganize_output.py --dry-run   # 먼저 미리보기\n"
            "  python reorganize_output.py             # 실제 실행\n"
        ),
    )
    parser.add_argument(
        "--project-dir", default=None,
        help="프로젝트 출력 디렉토리 (기본: config에서 자동 감지)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="미리보기만 (실제 파일 이동 안 함)",
    )
    parser.add_argument(
        "--no-symlinks", action="store_true",
        help="호환용 심볼릭 링크 생성 안 함",
    )
    args = parser.parse_args()

    # 프로젝트 디렉토리 결정
    if args.project_dir:
        project_dir = Path(args.project_dir)
    else:
        config_path = REPO_ROOT / "config" / "example-project.yaml"
        if config_path.exists():
            import yaml
            with open(config_path, encoding="utf-8") as f:
                config = yaml.safe_load(f)
            output_root = config.get("output_root", "./output")
            project_name = config.get("project_name", "")
            if not project_name:
                print("[ERROR] config에 project_name이 없습니다.")
                sys.exit(1)
            project_dir = REPO_ROOT / output_root / project_name
        else:
            print("[ERROR] config 파일 없음. --project-dir를 지정하세요.")
            sys.exit(1)

    if not project_dir.is_absolute():
        project_dir = REPO_ROOT / project_dir

    reorg = Reorganizer(
        project_dir,
        dry_run=args.dry_run,
        create_symlinks=not args.no_symlinks,
    )

    success = reorg.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
