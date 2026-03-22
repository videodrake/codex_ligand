# 3 File System — Vibe Coding Agent Generator

기획서(plan.md) → 3파일(PRD + tasks + tech-stack) → 에이전트 시스템 → 개발 시작. 이 레포는 변환 템플릿과 프로젝트 작업 공간을 관리한다.

## 핵심 구조

- `templates/stage1.md` — 기획서 → 3파일 변환 지시서 (시스템 핵심. 삭제 금지.)
- `templates/stage2.md` — 3파일 → 에이전트 시스템 변환 지시서 (시스템 핵심. 삭제 금지.)
- `projects/{name}/` — 각 프로젝트 작업 공간 (삭제 가능)

## 새 프로젝트 생성 방법

1. `projects/{name}/docs/` 폴더 생성
2. 사용자의 기획 문서를 `projects/{name}/docs/`에 저장 (단일 파일이든 여러 파일이든 상관없음)
3. `templates/stage1.md`와 `templates/stage2.md`를 `projects/{name}/`에 복사
4. 프로젝트 폴더로 이동하여 Stage 1 → Stage 2 순서로 실행

## 파이프라인

```
기획서 (plan.md)
    │  stage1.md 실행
    ▼
PRD + tasks.md + tech-stack.md    ← 3 File System
    │  stage2.md 실행
    ▼
CLAUDE.md + 에이전트 4개 + 스킬 + CONTEXT.md
    │  /recover → /execute → /review → /test
    ▼
완성된 프로젝트
```

## 실험적 근거

EGFR-MYO1D 파이프라인에서 확인된 실험적 사실:

1. **ATP 결합 유지 + 활성 소실**: 약물 처리 시 ATP 결합은 유지되면서 kinase 활성만 소실됨. 따라서 ATP binding site 포켓은 파이프라인에서 false positive로 자동 배제 (`is_atp_site = True`, `exclusion_reason = "ATP_site_experimental"`).
2. **Ko et al. alanine substitution**: MYO1D beta-meander sheet 8/9 잔기가 EGFR과의 직접 결합면(active face). 이 잔기가 PPI hotspot에 3개 미만이면 Workflow B를 중단.
3. **리간드 다양성**: 3종 리간드(173940, 97806, VAX-C12_0) 쌍별 Tanimoto < 0.4 → 화학적으로 다양, cross-chemical consensus 유효.

## 규칙

- templates/ 폴더는 절대 삭제하지 않는다. 프로젝트 삭제 시 projects/{name}/만 삭제.
- 각 Stage는 코드를 작성하지 않고 문서만 생성한다.
- Stage 1 완료 후 사용자 검토를 받고 Stage 2를 진행한다.
