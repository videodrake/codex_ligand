# EGFR-MYO1D Pipeline — 하네스 엔지니어링 적용 설계서

> **프로젝트**: EGFR-MYO1D PPI Docking Pipeline
> **하네스 목적**: AI 코딩 에이전트(Claude Code)가 이 파이프라인을 안정적으로 개발·유지·확장할 수 있도록 제어 환경을 구축
> **작성일**: 2026-04-06

---

## 1. 현황 분석

### 1.1 이미 있는 것 (기존 하네스 요소)

| 요소 | 현재 상태 | 평가 |
|------|-----------|------|
| CLAUDE.md | 존재 — 실험적 근거 3개, 참조 문서 맵, 규칙 포함. 단, 3 File System 관련 내용은 폐기 대상 | ⚠️ 과학 컨텍스트 양호, 에이전트 행동 제약 보강 필요 |
| 커맨드 시스템 | /recover, /execute, /review, /test 4개 | ❌ 3 File System(prd.md, tasks.md) 의존 — 전부 폐기 |
| 문서 | 14개 이상 — PIPELINE_ARCHITECTURE_REPORT, 설계의도, runbook 등 | ✅ 풍부. 에이전트 문서 라우팅만 추가하면 됨 |
| 테스트 | tests/ + validate.py (4그룹 8함수, 종료코드 0/1/2) | ✅ 검증 프레임워크 존재 |
| 코드 내 안전장치 | `_validate_adv_handoff()` (Phase 간 핸드오프 사전 검증), Vina 가용성 가드, Phase 3 모드별 사전조건 | ✅ defense-in-depth 이미 구현됨 |
| 경로 관리 | `egfr_pipeline/paths.py` 중앙 관리 | ✅ 존재 — 하네스에서 보호 필요 |
| 가드레일 | 없음 | ❌ |
| 스킬 | 없음 | ❌ |

### 1.2 폐기 대상 (하네스 설계에서 제외)

- **3 File System** — templates/stage1.md, stage2.md, projects/ 폴더 전체
- **3 File System 산출물** — docs/prd.md, docs/tasks.md (3 File System이 생성하던 파일)
- **커맨드 4개** — .claude/commands/recover.md, execute.md, review.md, test.md (prd.md, tasks.md 의존)
- **nightly_review** — scripts/nightly_review.py, docs/nightly_review_automation.md, docs/nightly_incremental_improvement_automation.md
- **CLAUDE_org.md** — 이전 버전 CLAUDE.md 백업
- **docs/CONTEXT.md** — 기존 내용 폐기, 하네스 세션 메모리 용도로 재작성하여 유지

### 1.3 프로젝트 특수성

- **HPC 환경**: 네트워크 차단, PBS/qsub 기반 — 에이전트가 직접 실행 불가, 스크립트 생성만 가능
- **이중 워크플로우**: Workflow A(blind) + B(PPI-first) — **같은 "Phase 1"이 다른 것을 가리킴** (아래 상세)
- **과학적 제약**: 실험적 사실이 코드 로직에 직접 반영 — 위반 시 과학적 오류
- **다단계 의존**: Workflow B에서 Phase 간 순차 의존 + 핸드오프 CSV 계약
- **역사적 버그**: DockingSlideIntoContact 누락 시 dG=0.0 등 — 스킬에 반드시 포함
- **경로 중앙 관리**: paths.py 수정 = 전체 파이프라인 영향

### 1.4 Phase 번호 혼동 위험 (핵심)

이 프로젝트의 가장 큰 혼동 요소:

```
Workflow A                          Workflow B
─────────                          ─────────
Phase 1 = Vina Blind               Phase 1 = PPI 분석 (TG 1.0~1.6)
Phase 2 = PPI Global Blind         Phase 2 = Pocket Analysis (TG 2.0~2.7)
Phase 3 = PPI Postprocess          Phase 3 = Focused Vina (TG 3.0~3.6)
Phase 4 = Vina Postprocess         Phase 4 = Perturbation Scoring (TG 4.0~4.6)
Phase 5 = Verdict
Phase 6 = Report
Phase 7 = Validate
```

**코드 디렉토리 매핑도 직관적이지 않음**:
- `egfr_pipeline/phase1/` = Workflow **B**의 Phase 1 (PPI 분석)
- `egfr_pipeline/vina/` = Workflow **A**의 Phase 1 + Phase 4
- `egfr_pipeline/ppi/` = Workflow **A**의 Phase 2~3

→ 에이전트가 "Phase 1 수정해줘"라고 들으면 **반드시 어떤 워크플로우인지 먼저 확인**해야 함

---

## 2. 하네스 아키텍처 설계

```
┌──────────────────────────────────────────────────────────┐
│                   CLAUDE.md (Entry Point)                  │
│  ┌────────────┐  ┌────────────┐  ┌──────────────┐        │
│  │ 프로젝트    │  │ 절대 규칙   │  │ 문서 라우팅   │        │
│  │ 한 줄 요약  │  │ 7개        │  │ (스킬 안내)   │        │
│  └────────────┘  └────────────┘  └──────────────┘        │
└──────────────────────┬───────────────────────────────────┘
                       │
       ┌───────────────┼───────────────┐
       ▼               ▼               ▼
┌──────────────┐ ┌──────────────┐ ┌──────────────┐
│   Skills/    │ │   Agents/    │ │   Hooks/     │
│  도메인 지식  │ │  전문 에이전트 │ │  자동 검증   │
├──────────────┤ ├──────────────┤ ├──────────────┤
│ vina-docking │ │ pipeline-dev │ │ pre-commit   │
│ ppi-analysis │ │ reviewer     │ │ csv-schema   │
│ hpc-ops      │ │ science-qa   │ │              │
│ phase-deps   │ │              │ │              │
│ scoring      │ │              │ │              │
│ testing      │ │              │ │              │
│ bug-history  │ │              │ │              │
└──────────────┘ └──────────────┘ └──────────────┘
       │               │               │
       └───────────────┼───────────────┘
                       ▼
             ┌──────────────────┐
             │  docs/CONTEXT.md │
             │  (세션 간 메모리) │
             └──────────────────┘
```

---

## 3. CLAUDE.md 재설계

CLAUDE.md는 백과사전이 아니라 짧은 지도 역할을 한다. 실제 source of truth는 docs/와 스킬 파일에 둔다.

### 3.1 구조 (60줄 이내)

```markdown
# EGFR-MYO1D Docking Pipeline

EGFR-MYO1D PPI 교란 약물 포켓 탐색 파이프라인.
Vina(소분자) + PyRosetta(PPI) 이중 증거 통합, 3개 EGFR 구조 상태 교차 비교.

## 절대 규칙

[아래 3.2 참조]

## 워크플로우 구분

[Phase 번호 매핑표 — 1.4 내용 축약]

## Definition of Done

코드 변경을 완료로 간주하려면:
- pytest -m smoke 전체 통과
- 변경된 Phase의 validate.py 검증 통과
- CSV 스키마 변경 시 하위 ingestion 코드 동시 수정 확인
- paths.py 변경 시 전체 smoke test 통과

## 스킬 (.claude/skills/)

작업 시 관련 스킬을 먼저 읽는다:
- vina-docking/ — Vina 관련 작업 시
- ppi-analysis/ — PPI/PyRosetta 작업 시
- hpc-operations/ — PBS/qsub/서버 관련 시
- phase-dependencies/ — Phase 간 전환 시
- scoring-system/ — Verdict/스코어링 변경 시
- testing/ — 테스트 작성/실행 시
- bug-history/ — PyRosetta/Vina 코드 수정 전 반드시 확인

## 독립 모듈

- egfr_pipeline/md/ — GROMACS MD 분석 (선택적, MDAnalysis 필요). 핵심 파이프라인과 독립.

## 참조 문서 (docs/ = source of truth)

- PIPELINE_ARCHITECTURE_REPORT.md — 전체 아키텍처, 모듈별 입출력
- docs/runbook.md — 실행 가이드
- docs/methodology_limitations.md — 방법론 한계 5개 섹션
- docs/workflow_comparison_guide.md — Workflow A↔B 불일치 해석
- docs/phase4_A3_axis_specification.md — Phase 4 A3 축 계산 로직
- docs/data_inventory.md — 입출력 데이터 인벤토리
- 투두리스트.md — 미구현 항목 추적
```

### 3.2 절대 규칙 (7개)

```markdown
## 절대 규칙

1. **HPC 전용 실행**: 도킹/연산 코드를 직접 실행하지 않는다.
   PBS 스크립트를 생성하고 qsub 명령을 안내한다.

2. **ATP 포켓 배제**: is_atp_site=True인 포켓을 STRONG 판정하면 안 된다.
   ATP 결합 유지 + 활성 소실이 실험적 사실이다.

3. **워크플로우 확인 의무**: "Phase N 수정" 요청 시 반드시 Workflow A/B 중
   어느 쪽인지 확인한다. 같은 번호가 다른 모듈을 가리킨다.
   확인 없이 작업을 시작하지 않는다.

4. **디렉토리-워크플로우 매핑**: 코드 수정 시 이 매핑을 따른다:
   - egfr_pipeline/vina/ → Workflow A (Phase 1 + 4)
   - egfr_pipeline/ppi/ → Workflow A (Phase 2 + 3)
   - egfr_pipeline/phase1/ → Workflow B Phase 1
   - egfr_pipeline/phase2/ → Workflow B Phase 2
   - egfr_pipeline/phase3/ → Workflow B Phase 3
   - egfr_pipeline/phase4/ → Workflow B Phase 4
   - verdict.py, report.py, validate.py → Workflow A Phase 5~7

5. **CSV 스키마 보존**: 기존 CSV 출력의 컬럼명·타입을 변경할 때는
   하위 Phase의 ingestion 코드와 validate.py를 함께 수정한다.

6. **paths.py 보호**: egfr_pipeline/paths.py를 수정하면
   전체 Phase의 경로 해석이 바뀐다.
   수정 후 반드시 모든 Phase의 smoke test를 실행한다.

7. **잔기 번호 체계**: PDB 잔기 번호(author numbering)를 사용한다.
   내부 인덱스(0-based)와 혼동하면 전체 분석이 틀어진다.
```

---

## 4. 스킬 시스템 설계 (.claude/skills/)

### 4.1 skill-vina-docking/SKILL.md

```
name: vina-docking
description: Vina 도킹 관련 작업 시 로딩. 트리거 — Vina, 도킹, 리간드, PDBQT, exhaustiveness, affinity 언급 시. 비트리거 — PPI/PyRosetta 작업, 스코어링 로직만 변경할 때는 이 스킬이 아님.

## 핵심 정보
- vina_executor.py: prepare_receptor → prepare_ligand → run_vina
- exhaustiveness=384 (기본값 8의 48배, blind docking 70Å+ box 대응)
- 3종 리간드: 173940, 97806, VAX-C12_0 (쌍별 Tanimoto < 0.4)
- 에너지 단위: 항상 kcal/mol

## Workflow별 역할
- Workflow A: Phase 1 (blind, 전체 표면) + Phase 4 (postprocess)
- Workflow B: Phase 3 (focused, 포켓별 집중, budget-aware)

## PDBQT 변환 fallback 순서
Meeko → ADFR → MGLTools → OpenBabel (코드: vina_executor.py)

## 주의사항
- Workflow A blind box는 EGFR 전체를 감싸는 70Å+
- Workflow B focused box는 Phase 2에서 정의된 포켓별 좌표 사용
- 두 모드의 결과를 직접 비교하면 안 됨 (탐색 범위가 다름)

## 이 스킬을 쓰지 말아야 할 때
- Verdict/스코어링 점수 체계만 바꿀 때 → scoring-system 스킬
- Phase 간 CSV 핸드오프 문제 → phase-dependencies 스킬
- PBS 스크립트만 수정할 때 → hpc-operations 스킬
```

### 4.2 skill-ppi-analysis/SKILL.md

```
name: ppi-analysis
description: PPI/PyRosetta 관련 작업 시 로딩. 트리거 — PPI, PyRosetta, interface, patch, orientation filter, MYO1D 언급 시. 비트리거 — Vina 소분자 도킹만 다룰 때, MD 분석 시에는 이 스킬이 아님.

## 핵심 정보
- PyRosetta PPI 도킹: 3 states × 5 seeds = 15 세트, 각 20K 모델 = 300K 총
- orientation_filter: sheet 8/9의 active-face normal vs receptor 방향 dot product
  - 양수 = pass (active face가 receptor 향함)
  - 음수 = fail (뒤집힘)
  - consensus 계산에는 pass 모델만 사용

## 실험 데이터 매핑
- Ko et al. alanine substitution: sheet 8/9 잔기 = active face
- PPI hotspot에 이 잔기가 3개 미만이면 Workflow B 중단 조건

## LightDock 검증
- PyRosetta와 독립적인 교차 검증 수단
- LightDock 결과가 PyRosetta와 일치하면 신뢰도 상승
- 단, LightDock 자체의 한계 있음 (PIPELINE_ARCHITECTURE_REPORT.md TG 1.4 참조)

## 위험한 코드 변경
- DockingSlideIntoContact 누락 → 모든 dG가 0.0 (V1.0 역사적 버그)
- FoldTree: 역직렬화 후 반드시 setup_foldtree 재설정 필요
- excluded_residues_A: 막면/다이머 인터페이스 금지 구역, hard filter
- key_residues_B: 실험 데이터 기반, soft bonus (adjusted_dG)
- enable_early_rejection: DockMCMProtocol 전 금지구역 접촉 검사 (연산 절약)
```

### 4.3 skill-hpc-operations/SKILL.md

```
name: hpc-operations
description: PBS/qsub/서버 관련 작업 시 로딩. 트리거 — qsub, PBS, 서버, HPC, 프로덕션 실행 언급 시. 비트리거 — 로컬에서 도는 단위 테스트, 코드 로직만 바꿀 때는 이 스킬이 아님.

## 대원칙
모든 도킹/연산은 반드시 qsub로 HPC 서버에서 실행한다.
에이전트는 PBS 스크립트를 생성/수정만 하고, 직접 실행하지 않는다.

## 실행 방법
- Workflow A: qsub config/run_production.pbs (또는 lane별 병렬)
- Workflow B: qsub config/run_advanced_pipeline.pbs
- 개별 lane: run_production.py --lane {lane_name}

## Lane 목록 (14개)
vina-cpu, ppi, ppi-post, vina-post, finalize, status,
vina-gpu, phase3-gpu,
adv-phase1, adv-phase2, adv-phase3-setup, adv-phase3-execute,
adv-phase3-post, adv-phase4

## 모드 구분
- dry-run: 검증만 수행 (실제 도킹 안 함)
- setup: 서버 실행용 스크립트만 생성
- execute: 실제 도킹 실행

## Phase 3 cascade (Workflow B)
rerun_cascade.py가 setup → execute → post 3모드를 순차 관리
모드별 사전조건: setup → job table, post → round log
```

### 4.4 skill-phase-dependencies/SKILL.md

```
name: phase-dependencies
description: Phase 간 전환, 핸드오프, 의존 관계 작업 시 로딩. 트리거 — Phase 간 전환, 핸드오프, 다음 Phase, 입력 참조 언급 시. 비트리거 — 단일 Phase 내부 로직만 바꿀 때, 스코어링 점수 체계만 변경할 때.

## Workflow A 의존 그래프 (독립 → 병합)
Phase 1(Vina Blind) ──────┐
                           ├→ Phase 5(Verdict) → Phase 6(Report) → Phase 7(Validate)
Phase 2(PPI Blind) ───────┘
  └→ Phase 3(PPI Post) ───┘
Phase 4(Vina Post) ────────┘

핵심: Phase 1과 Phase 2는 서로의 결과를 사용하지 않는다.
Phase 5(Verdict)에서 처음으로 병합.

## Workflow B 순차 의존 그래프
Phase 1 → Phase 2 → Phase 3 → Phase 4
전제 조건: Workflow A의 PPI 도킹(Phase 2)이 완료되어야 시작 가능

## 핸드오프 CSV 파일 (Workflow B)
| 구간 | 파일 | 생성 TG |
|------|------|---------|
| Phase 1 → 2 | phase1_downstream_patch_reference.csv | TG 1.6 |
| Phase 2 → 3 | phase3_candidate_pocket_reference.csv | TG 2.6 |
| Phase 3 → 4 | phase4_docking_evidence_reference.csv | TG 3.6 |

## 기존 안전장치 (코드 내장)
- _validate_adv_handoff(): 각 Phase 시작 전 핸드오프 파일 존재 사전 검증
- Vina 가용성 가드: silent all-skip 방지
- Phase 3 cascade 모드별 사전조건 검증

→ 이 안전장치들을 건드리지 않는다. 추가 검증이 필요하면 기존 패턴을 따른다.
```

### 4.5 skill-scoring-system/SKILL.md

```
name: scoring-system
description: Verdict/스코어링 관련 작업 시 로딩. 트리거 — Verdict, 점수, 축, STRONG/MODERATE/WEAK, 스코어링 언급 시. 비트리거 — Vina 도킹 파라미터만 바꿀 때, PPI 도킹 자체를 수정할 때.

## Workflow A: 3축 체계 (verdict.py)
- 축 1 Vina Quality (50점): affinity + convergence + consensus + stability + diversity
- 축 2 PPI Spatial (20점): spatial + overlap + reproducibility
- 축 3 Cross-Receptor (30점): 다중 구조 상태 일관성
- PPI 없을 시: 60 + 0 + 40 = 100으로 적응적 재배분
- STRONG ≥ 55 (최소 2축에서 의미 있는 점수 필요)

## Workflow B: 4축 체계 (phase4/)
- A1: PPI interface 관계 (orthosteric/rim/allosteric/irrelevant)
- A2: Druggability
- A3: Perturbation relevance (PIPELINE_ARCHITECTURE_REPORT.md Phase 4 섹션 참조)
- A4: State robustness
- A1+A3 합산 가중치 60% → affinity만 좋고 MYO1D 무관한 포켓은 상위 불가

## 판정 원칙
"증거 분류이지 타당성 판정이 아니다."
STRONG도 PyMOL 시각 검증 필수. WEAK도 cryptic pocket 가능성 있음.
```

### 4.6 skill-testing/SKILL.md

```
name: testing
description: 테스트 작성/실행/검증 관련 작업 시 로딩. 트리거 — 테스트, 검증, validate, 회귀, smoke test 언급 시. 비트리거 — 기능 구현 자체만 할 때 (구현 후 테스트는 별도 단계로).

## validate.py (Workflow A 출력 검증)
4개 그룹, 8개 함수:
- 그룹 1: 파일 존재 + ID 일관성 + 추적성 + coverage
- 그룹 2: CSV 스키마 회귀 검사
- 그룹 3: 잔기 번호 일관성 + 알려진 변이 확인
- 그룹 4: 핸드오프 준비 확인
종료 코드: 0(통과) / 1(경고) / 2(실패)

## pytest 마커
- smoke: 빠른 기본 검증 (도킹 없이 실행 가능)
- full: PyRosetta/Vina 필요한 통합 테스트

## mock 전략
- PyRosetta 없는 환경: mock pose 객체 사용
- Vina 없는 환경: 사전 생성된 결과 파일로 후처리 테스트

## paths.py 수정 후
반드시 pytest tests/ -m smoke --tb=short 전체 실행
```

### 4.7 skill-bug-history/SKILL.md

```
name: bug-history
description: PyRosetta/Vina 코드 수정 전 반드시 확인. 트리거 — PyRosetta 코드 수정, 도킹 엔진 변경, 스코어 계산 변경 시. 비트리거 — 문서만 수정할 때, PBS 스크립트만 바꿀 때, 테스트만 추가할 때.

## 역사적 버그 — 반복하면 안 되는 실수들

### BUG-001: DockingSlideIntoContact 누락 (V1.0, 심각도: Critical)
- 증상: 모든 dG 값이 0.0
- 원인: DockMCMProtocol 호출 전 DockingSlideIntoContact 누락
- 교훈: PyRosetta 도킹 프로토콜 수정 시 반드시 dG 값이 비-제로인지 확인

### BUG-002: FoldTree 미재설정 (심각도: High)
- 증상: 역직렬화된 pose에서 도킹 결과가 비정상
- 원인: setup_foldtree(pose, "A_B", movable_jumps) 호출 누락
- 교훈: pose 로드/역직렬화 후 반드시 FoldTree 재설정

### BUG-003: stdout/stderr 리다이렉트 깨짐 (심각도: Medium)
- 증상: PyRosetta 배너가 출력에 섞임
- 원인: pyrosetta_init.py의 리다이렉트 구조 변경
- 교훈: pyrosetta_init.py 수정 시 배너 억제 동작 확인

### BUG-004: Beta-meander 방향 미검증 (심각도: High)
- 증상: back face가 receptor를 향한 모델이 consensus에 포함
- 원인: orientation_filter 미적용 또는 통과 기준 오류
- 교훈: PPI 결과 처리 시 항상 orientation pass/fail 상태 확인

## 코드 수정 전 체크리스트
□ 이 파일의 관련 버그를 확인했는가?
□ 수정 후 dG 값이 합리적인가? (0.0이 아닌가?)
□ FoldTree가 올바르게 설정되어 있는가?
□ orientation_filter가 정상 동작하는가?
```

---

## 5. 에이전트 정의 (.claude/agents/)

이 프로젝트는 "새 앱을 처음부터 빌드"하는 게 아니라 기존 과학 파이프라인의 유지보수·확장이다. Planner→Generator→Evaluator 생성 루프보다는, 역할별 전문 에이전트가 더 적합하다.

### 5.0 권한 경계 원칙

에이전트가 자율적으로 처리할 수 있는 것과, 반드시 사람의 승인이 필요한 것을 구분한다.

**에이전트 자율 가능:**
- 단일 Phase 내부의 버그 수정
- 테스트 추가/수정
- 문서 업데이트
- PBS 스크립트 생성 (실행은 사람)
- 코드 스타일/린트 수정

**반드시 사람 승인 필요:**
- Phase 간 핸드오프 CSV 스키마 변경
- 스코어링 축/가중치 변경 (과학적 판단)
- 새 Phase 추가 또는 워크플로우 구조 변경
- paths.py 수정
- 실험적 사실(ATP 배제, Ko et al.) 관련 로직 변경
- 프로덕션 output/ 결과에 영향을 주는 변경

### 5.1 pipeline-dev.md (기본 개발 에이전트)

```
역할: 파이프라인 모듈 개발 및 수정

접근 가능: 모든 소스 코드, 테스트, 설정

작업 전 필수:
- "어떤 워크플로우의 몇 번째 Phase인가?" 확인 (규칙 3)
- 해당 디렉토리-워크플로우 매핑 확인 (규칙 4)
- 관련 스킬 로딩 (bug-history는 PyRosetta/Vina 수정 시 필수)

제약:
- Phase 간 CSV 스키마 변경 시 하위 ingestion 코드 + validate.py 동시 수정
- paths.py 수정 후 전체 smoke test 필수
- 새 Phase 추가 시 data_inventory.md 업데이트

금지:
- qsub 직접 실행 (스크립트 생성만)
- output/ 디렉토리의 프로덕션 결과 삭제
- input/ 디렉토리의 구조 파일(PDB/SDF) 수정
- _validate_adv_handoff() 등 기존 안전장치 제거
```

### 5.2 reviewer.md (코드 리뷰 에이전트)

```
역할: 구현된 코드의 과학적·구조적 정합성 검토

5개 관점으로 검토:
1. 과학적 정합성: 실험적 사실 3개(ATP, Ko et al., 리간드 다양성)와 모순 없는가?
2. Phase 의존성: 핸드오프 CSV 계약을 깨지 않는가? _validate_adv_handoff()가 여전히 통과하는가?
3. 잔기 번호: PDB author numbering이 일관되는가?
4. 에너지 단위: kcal/mol 일관성
5. 버그 히스토리: skill-bug-history의 패턴을 반복하지 않는가?

출력: PASS / NEEDS_FIX (수정 사항 + 관련 스킬 참조)
```

### 5.3 science-qa.md (과학 검증 에이전트)

```
역할: 파라미터·알고리즘 변경의 과학적 타당성 질문

트리거: 파라미터 변경, 새 필터 추가, 스코어링 변경, 임계값 조정

질문 체크리스트:
□ blind docking의 편향 없는 탐색 원칙을 위반하는가?
□ cross-chemical consensus (3종 리간드)에 영향을 주는가?
□ ATP binding site 배제 로직이 여전히 작동하는가?
□ Ko et al. 실험 데이터와 모순되는가?
□ 3개 receptor state 비교에 영향을 주는가?
□ Workflow A의 "독립 증거원 → 사후 병합" 원칙을 깨는가?
□ Workflow B의 순차 의존 구조를 깨는가?
```

---

## 6. 검증 계층 설계

이 프로젝트에 맞는 검증은 3계층으로 충분하다. 브라우저 E2E나 LLM-as-judge는 HPC 과학 파이프라인에는 해당하지 않는다.

```
계층 1: 기계적 검증 (매 커밋)
  - pytest -m smoke (빠른 기본 검증, 도킹 없이)
  - csv-schema-guard.py (핸드오프 CSV 컬럼 변경 감지)
  - paths.py 변경 감지 → 전체 smoke test

계층 2: 계약 검증 (Phase 수정 시)
  - validate.py 4그룹 8함수 (파일 존재 + 스키마 + 잔기 일관성 + 핸드오프)
  - _validate_adv_handoff() (Phase 간 핸드오프 파일 사전 검증, 코드 내장)
  - 해당 Phase의 단위 테스트

계층 3: 과학적 검증 (파라미터/알고리즘 변경 시)
  - reviewer 에이전트의 5개 관점 검토
  - science-qa 에이전트의 질문 체크리스트
  - PyMOL 시각 검증 (사람이 수행)
```

---

## 7. 훅 시스템 설계 (.claude/hooks/)

훅은 2개로 구성한다. output 보호는 에이전트 금지 규칙으로, Phase 간 검증은 코드 내 `_validate_adv_handoff()`로 이미 담당하므로 별도 훅이 불필요하다.

### 7.1 pre-commit (실행 가능한 bash 스크립트)

```bash
#!/bin/bash
# .claude/hooks/pre-commit.sh
# Claude Code hook: PreToolUse(Bash) 또는 커밋 전 실행

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
```

### 7.2 csv-schema-guard (실행 가능한 Python 스크립트)

```python
#!/usr/bin/env python3
"""CSV 스키마 변경 감지 및 영향 범위 표시.

Claude Code hook: 커밋 전 실행.
변경된 CSV 출력 컬럼이 있으면 해당 CSV를 읽는 하위 모듈 목록을 표시.
"""
import subprocess
import re
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
        capture_output=True, text=True
    )
    return result.stdout.strip().split("\n")

def check_csv_columns(filepath):
    """변경된 파일에서 CSV 컬럼 정의 변경을 감지."""
    result = subprocess.run(
        ["git", "diff", "--cached", "-U0", filepath],
        capture_output=True, text=True
    )
    # pd.DataFrame 컬럼 정의, CSV 헤더 쓰기 등의 변경 감지
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
            # 이 파일이 핸드오프 CSV를 생성하는지 확인
            for csv_name, consumers in HANDOFF_MAP.items():
                content = Path(f).read_text(errors="ignore")
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
```

---

## 8. 피드백 루프 설계

### 8.1 실패 학습 메커니즘

```
에이전트 실패 발생
    ↓
실패 유형 분류:
  - 과학적 오류 → skill-bug-history에 추가
  - Phase 혼동 → CLAUDE.md 매핑표 보강
  - CSV 스키마 깨짐 → csv-schema-guard 매핑 추가
  - 일반 코드 오류 → docs/CONTEXT.md에 기록
    ↓
다음 세션 시작 시 에이전트가 CLAUDE.md + docs/CONTEXT.md를 읽으며 자동 반영
```

### 8.2 Drift 정리

코드가 시간이 지나며 나빠지는 것을 감지하기 위해, 주기적으로(또는 큰 변경 후) 에이전트에게 다음을 요청한다:

- 미사용 import / dead code 탐지
- validate.py 실행 → 경고/실패 항목 보고
- 투두리스트.md의 미구현 placeholder가 여전히 유효한지 확인
- docs/와 실제 코드의 불일치 탐지 (문서 drift)

### 8.3 docs/CONTEXT.md 구조

```markdown
# 프로젝트 컨텍스트

## 현재 작업 상태
- 워크플로우: A / B
- 현재 Phase: N
- 마지막 완료 태스크: ...

## 최근 결정 사항
- [날짜] 결정 내용 + 이유

## 발견된 이슈 (미해결)
- [ ] 이슈 내용

## 실패 패턴 (반복 방지)
- [날짜] 실패 유형 + 원인 + 해결 방법
```

---

## 9. 최종 디렉토리 구조

```
.claude/
├── agents/
│   ├── pipeline-dev.md
│   ├── reviewer.md
│   └── science-qa.md
├── skills/
│   ├── vina-docking/SKILL.md
│   ├── ppi-analysis/SKILL.md
│   ├── hpc-operations/SKILL.md
│   ├── phase-dependencies/SKILL.md
│   ├── scoring-system/SKILL.md
│   ├── testing/SKILL.md
│   └── bug-history/SKILL.md
└── hooks/
    ├── pre-commit.sh
    └── csv-schema-guard.py

CLAUDE.md                            ← 재작성 (절대 규칙 7개 + 스킬 안내 + md/ 모듈 안내)
docs/CONTEXT.md                      ← 세션 간 메모리 (기존 내용 폐기 후 재작성)
```

---

## 10. 구현 우선순위

### Phase 0: 정리 (Claude Code에게 시킬 작업)
0. **폐기 파일 삭제 + 참조 정리** — 아래 섹션 11의 목록대로 실행

### Phase 1: 즉시 (가장 높은 ROI)
1. **CLAUDE.md 재작성** — 절대 규칙 7개, 워크플로우 매핑표, 스킬 안내, md/ 모듈 안내
2. **skill-phase-dependencies** — Phase 혼동이 가장 위험한 실패
3. **skill-bug-history** — 역사적 버그 반복 방지

### Phase 2: 안전망
4. **skill-ppi-analysis** (DockingSlideIntoContact 등 핵심 위험 포함)
5. **skill-vina-docking**
6. **csv-schema-guard.py** 훅

### Phase 3: 완성
7. **나머지 스킬** (hpc-ops, scoring, testing)
8. **에이전트 정의** (pipeline-dev, reviewer, science-qa)
9. **pre-commit.sh** 훅
10. **docs/CONTEXT.md** 재작성

---

## 11. Claude Code 정리 작업 목록

하네스 구축 전에 Claude Code에게 시킬 레거시 정리 작업.
"이 설계서의 섹션 11을 실행해줘"로 지시 가능.

### 11.1 삭제할 파일

```
# 3 File System
templates/stage1.md
templates/stage2.md
templates/                          ← 폴더 자체 (비면)
projects/                           ← 폴더 전체 (하위 포함)

# 3 File System 산출물
docs/prd.md
docs/tasks.md

# 커맨드 (3 File System 의존)
.claude/commands/recover.md
.claude/commands/execute.md
.claude/commands/review.md
.claude/commands/test.md
.claude/commands/                    ← 폴더 자체 (비면)

# nightly_review
scripts/nightly_review.py
docs/nightly_review_automation.md
docs/nightly_incremental_improvement_automation.md

# 이전 버전 CLAUDE.md
CLAUDE_org.md
```

### 11.2 참조 정리 (삭제 후 깨지는 참조 수정)

- **README.md** 문서 안내 테이블에서 nightly_review 2개 문서 행 삭제
- **CLAUDE.md** 전체 재작성 (3 File System 내용 제거 → 본 설계서 섹션 3 기반)
- **테스트**: docs/prd.md, docs/tasks.md 존재를 확인하는 테스트가 있으면 삭제 또는 수정
  - `tests/test_e2e_group7.py`의 `TestDocumentExistence`에서 "docs/prd.md", "docs/tasks.md" 항목 제거

### 11.3 Claude Code가 확인해야 할 것

- `docs/archive/` 폴더 내용 스캔 — 현재 코드와 맞지 않는 문서가 있으면 보고
- `input/PPI/prepared/` — 다른 코드에서 아직 참조하는지 grep으로 확인, 미참조면 삭제 가능 보고
- `paths.py`의 DEPRECATED 함수 — 다른 모듈에서 호출하는지 grep으로 확인, 미호출이면 삭제 가능 보고
- `egfr_pipeline/md/` 모듈 — CLAUDE.md에 "선택적, MDAnalysis 필요, 핵심 파이프라인과 독립" 한 줄 추가
