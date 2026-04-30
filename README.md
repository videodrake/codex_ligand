# EGFR-MYO1D Docking Pipeline

> AI/agent routing note: the active implementation target is the fresh workflow
> under `fresh/`. Start with `AI_START_HERE.md`, then
> `fresh/docs/AGENT_DOC_ROUTING.md`. The older links below are legacy context and
> must not be used as controlling scientific, workflow, or environment specs.

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

### 구조 이해
| 문서 | 내용 |
|------|------|
| [docs/architecture.md](docs/architecture.md) | 전체 아키텍처 (모듈/알고리즘/입출력) |
| [docs/PROJECT_USAGE_OVERVIEW.md](docs/PROJECT_USAGE_OVERVIEW.md) | 프로젝트 개요 |
| [docs/data_inventory.md](docs/data_inventory.md) | 입출력 인벤토리 |
| [docs/harness_design.md](docs/harness_design.md) | Claude Code 하네스 엔지니어링 설계 |
| [docs/harness_execution.md](docs/harness_execution.md) | 하네스 실행 가이드 |

### 실행
| 문서 | 내용 |
|------|------|
| [docs/runbook.md](docs/runbook.md) | 실행 가이드 |
| [docs/environment_setup.md](docs/environment_setup.md) | 환경 설정 |
| [docs/manual_vina.md](docs/manual_vina.md) | Vina 매뉴얼 |
| [docs/manual_pyrosetta.md](docs/manual_pyrosetta.md) | PyRosetta 매뉴얼 |

### 테스트/검증
| 문서 | 내용 |
|------|------|
| [docs/test_suite_triage.md](docs/test_suite_triage.md) | 테스트 분류 |
| [docs/pre_qsub_test_line.md](docs/pre_qsub_test_line.md) | 사전 제출 테스트 |

### 설정
| 문서 | 내용 |
|------|------|
| [config/README.md](config/README.md) | Config 파일 의미 |

### AI 에이전트
| 문서 | 내용 |
|------|------|
| [CLAUDE.md](CLAUDE.md) | Claude Code 컨텍스트 |

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

1. `output/workflow_a/phase6_report/project_report.txt` — 종합 보고서
2. `output/workflow_a/phase5_verdict/valid_sites.csv` — 후보 포켓 판정
3. `output/workflow_a/phase4_vina_postprocess/vina_pocket_table.csv` — 포켓별 상세
4. `output/workflow_a/phase2_ppi_docking/{state}/` — PPI 도킹 결과
5. PyMOL: `1_OVERVIEW_Clusters.pml` — 시각화
