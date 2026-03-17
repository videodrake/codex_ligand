# EGFR-MYO1D Docking Pipeline

EGFR 키나아제 도메인과 MYO1D beta-meander 간 단백질-단백질 상호작용(PPI)을 교란할 수 있는 **약물 결합 포켓**을 탐색하는 통합 파이프라인.

- AutoDock Vina (소분자 blind docking) + PyRosetta (PPI global blind docking) 이중 증거 통합
- 3가지 EGFR 구조 상태 (3GT8_raw, EGFR_160-185, EGFR_170-200) 교차 비교
- 모든 도킹/연산은 HPC 서버에서 `qsub`로 실행

## 빠른 시작

```bash
# Workflow A: Standard Production (전체 자동)
qsub config/run_pre_qsub_checks.pbs
qsub config/run_production.pbs

# Workflow B: Advanced PPI-First (PPI 결과 기반 정밀 탐색)
qsub config/run_advanced_pipeline.pbs
```

상세 실행 방법: [docs/runbook.md](docs/runbook.md)

## 문서 안내

| 문서 | 내용 | 대상 |
|------|------|------|
| [docs/PROJECT_USAGE_OVERVIEW.md](docs/PROJECT_USAGE_OVERVIEW.md) | 프로젝트 사용 개요 — 각 도구가 뭘 하는지, 어떤 정보를 얻는지 | 처음 읽는 사람 |
| [PIPELINE_ARCHITECTURE_REPORT.md](PIPELINE_ARCHITECTURE_REPORT.md) | 전체 아키텍처 상세 보고서 — 모든 모듈, 알고리즘, 입출력 | 구조 파악 |
| [docs/runbook.md](docs/runbook.md) | 실행 가이드 — qsub 명령, 순서, 결과 확인 | 실행할 때 |
| [docs/environment_setup.md](docs/environment_setup.md) | 환경 설정 — conda, PyRosetta, 서버 설정 | 최초 설정 |
| [config/README.md](config/README.md) | Config 파일 의미 — YAML, INI, PBS 설명 | 설정 변경 |
| [docs/manual_vina.md](docs/manual_vina.md) | AutoDock Vina 상세 매뉴얼 | Vina 참고 |
| [docs/manual_pyrosetta.md](docs/manual_pyrosetta.md) | PyRosetta PPI 도킹 상세 매뉴얼 | PPI 참고 |
| [docs/phase1_notes.md](docs/phase1_notes.md) | Phase 1 참고 노트 — 실행, 샘플링, 필터, LightDock | Phase 1 상세 |
| [docs/data_inventory.md](docs/data_inventory.md) | 입출력 데이터 인벤토리 | 데이터 추적 |
| [docs/pre_qsub_test_line.md](docs/pre_qsub_test_line.md) | 사전 제출 테스트 절차 | 테스트 |
| [docs/test_suite_triage.md](docs/test_suite_triage.md) | 테스트 분류 가이드 | 테스트 |
| [docs/nightly_review_automation.md](docs/nightly_review_automation.md) | 자동 리뷰 스크립트 | 자동화 |
| [CLAUDE.md](CLAUDE.md) | Claude Code AI 컨텍스트 | AI 세션 |

## 두 가지 워크플로우

### Workflow A: Standard Production

Vina blind + PPI blind → 독립 결과 통합. `run_production.py` 자동화.

```
Vina (약물 포켓) ──┐
                    ├→ Verdict (겹치면 = 후보) → Report → Validate
PPI (MYO1D 결합) ──┘
```

### Workflow B: Advanced PPI-First

PPI 결과 기반으로 포켓을 좁혀가며 정밀 탐색. 순차 의존.

```
PPI 도킹 → 포켓 분석(fpocket/P2Rank) → Focused Vina → 4축 스코어링
```

상세 설명: [docs/PROJECT_USAGE_OVERVIEW.md](docs/PROJECT_USAGE_OVERVIEW.md)

## 레포 구조

```
main.py                    # 대화형 CLI
run_production.py          # 프로덕션 자동화 (lane 기반)
egfr_pipeline/             # 코어 구현
  pyrosetta_docking/       #   PPI 도킹 엔진
  vina/                    #   Vina 도킹 + 분석
  phase1/ ~ phase4/        #   Phase별 모듈
  ppi/                     #   PPI 준비/후처리
config/                    # YAML, INI, PBS 스크립트
docs/                      # 문서
input/                     # 수용체/리간드 입력
output/                    # 결과 출력
tests/                     # 테스트
scripts/                   # 유틸리티 스크립트
```

## 결과 확인

프로덕션 완료 후:

1. `output/{project}/step_index.md` — 진입점
2. `output/{project}/project_report.txt` — 종합 보고서
3. `output/{project}/valid_sites.csv` — 후보 포켓 판정
4. `output/{project}/vina_pocket_table.csv` — 포켓별 상세
5. PyMOL: `1_OVERVIEW_Clusters.pml` — 시각화
