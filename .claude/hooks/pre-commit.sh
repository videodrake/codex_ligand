#!/bin/bash
# .claude/hooks/pre-commit.sh
# Claude Code hook: 커밋 전 실행

set -e

# 1. smoke test
pytest tests/ -m smoke --tb=short -q 2>&1 | tail -5

# 2. validate.py 빠른 검증 (파일 존재 + 스키마 체크만)
python -c "
from egfr_pipeline.validate import run_quick_checks
result = run_quick_checks()
if result > 0:
    print(f'WARNING: validate returned code {result}')
" 2>/dev/null || true

# 3. paths.py 변경 감지
if git diff --cached --name-only | grep -q "paths.py"; then
    echo "⚠️  paths.py가 변경되었습니다. 전체 smoke test를 실행합니다."
    pytest tests/ -m smoke --tb=short
fi
