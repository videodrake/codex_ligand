#!/usr/bin/env python3
"""CSV 스키마 변경 감지 및 영향 범위 표시.

Claude Code hook: 커밋 전 실행.
변경된 CSV 출력 컬럼이 있으면 해당 CSV를 읽는 하위 모듈 목록을 표시.
"""
import re
import subprocess
from pathlib import Path

# 핸드오프 CSV → 소비 모듈 매핑
HANDOFF_MAP = {
    "phase1_downstream_patch_reference.csv": [
        "egfr_pipeline/phase2/patch_ingestion.py",
    ],
    "phase3_candidate_pocket_reference.csv": [
        "egfr_pipeline/phase3/pocket_reference_ingestion.py",
    ],
    "phase4_docking_evidence_reference.csv": [
        "egfr_pipeline/phase4/",  # 전체 Phase 4
    ],
}


def get_changed_files():
    result = subprocess.run(
        ["git", "diff", "--cached", "--name-only"],
        capture_output=True,
        text=True,
    )
    return result.stdout.strip().split("\n")


def check_csv_columns(filepath):
    """변경된 파일에서 CSV 컬럼 정의 변경을 감지."""
    result = subprocess.run(
        ["git", "diff", "--cached", "-U0", filepath],
        capture_output=True,
        text=True,
    )
    patterns = [r"columns\s*=", r"\.to_csv", r"\.rename\(", r"\.drop\(.*columns"]
    for pattern in patterns:
        if re.search(pattern, result.stdout):
            return True
    return False


def main():
    changed = get_changed_files()
    warnings = []

    for f in changed:
        if not f.endswith(".py"):
            continue
        if check_csv_columns(f):
            for csv_name, consumers in HANDOFF_MAP.items():
                try:
                    content = Path(f).read_text(errors="ignore")
                except FileNotFoundError:
                    continue
                if csv_name in content:
                    warnings.append(
                        f"⚠️  {f}에서 {csv_name}의 컬럼이 변경된 것 같습니다.\n"
                        f"   소비 모듈: {', '.join(consumers)}\n"
                        f"   이 모듈들도 수정했는지 확인하세요."
                    )

    if warnings:
        print("\n".join(warnings))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
