# 실행 가이드 (Runbook)

Last updated: 2026-03-17

이 문서는 파이프라인 실행에 필요한 모든 절차와 명령어를 통합한 운영 가이드입니다.
모든 무거운 연산(도킹, PPI)은 HPC 서버에서 qsub로 제출합니다. 로컬에서 도킹을 실행하지 마십시오.

참고 문서:
- [architecture.md](architecture.md): 데이터 흐름 및 핸드오프 구조
- [output_artifact_map.md](output_artifact_map.md): 산출물 의미 및 우선순위
- [data_inventory.md](data_inventory.md): 물리적 입출력 위치

---

## 환경 설정

### 서버 접속 및 환경 활성화

```bash
cd ~/codex_ligand
conda activate pyrosetta
```

- 프로덕션과 precheck 모두 동일한 `pyrosetta` conda 환경을 사용합니다.
- 별도의 테스트 전용 환경을 가정하지 마십시오.

### 설정 파일 확인

기본 프로젝트 설정 파일:

```bash
config/example-project.yaml
```

CLI 실행 시 항상 `-c` 옵션을 명시합니다. `-c/--config`는 `main.py`의 최상위 옵션이며 서브커맨드 **앞에** 위치해야 합니다:

```bash
python main.py -c config/example-project.yaml <command>
```

### 사전 확인 사항

제출 전 반드시 확인:

1. 활성 config가 의도한 프로젝트 설정인지 확인
2. 3개 receptor state (`3GT8_raw`, `EGFR_160-185`, `EGFR_170-200`)가 등록되어 있는지 확인
3. 필요한 리간드와 준비된 입력 파일이 해당 lane에 사용 가능한지 확인
4. worker 수가 routine safe bound(16)를 초과하지 않는지 확인

---

## Workflow A: Standard Pipeline (run_production.py)

Vina 중심 리간드 증거 흐름 + PPI 인터페이스 증거를 통합 실행합니다.

### Phase 구성

| Phase | 내용 | 주요 산출물 |
|-------|------|-------------|
| 1 | Vina blind docking (3 receptor x 3 ligand) | `output/workflow_a/phase1_vina_docking/{receptor_id}/` |
| 2 | PPI docking (3 states x 5 seeds = 300K models) | `output/workflow_a/phase2_ppi_docking/{state}/prod_seed{n}/final_ranking.csv` |
| 3 | PPI postprocess | `output/workflow_a/phase3_ppi_postprocess/` |
| 4 | Vina postprocess | `output/workflow_a/phase4_vina_postprocess/vina_pocket_table.csv` |
| 5 | Verdict (3축 통합 scoring) | `output/workflow_a/phase5_verdict/valid_sites.csv` |
| 6 | Report | `output/workflow_a/phase6_report/project_report.txt` |
| 7 | Validate | `output/workflow_a/phase7_validation/` |

### Step 1: Pre-qsub Validation

항상 무거운 작업 전에 precheck를 먼저 실행합니다.

```bash
qsub config/run_pre_qsub_checks.pbs
```

필수 체크포인트: `output/precheck/last_pass.json`가 존재하고 pass 상태여야 합니다.
통과하지 못하면 **중단**하고 설정/입력/환경 문제를 해결합니다.

### Step 2: Production 제출

#### 안전한 체인 제출 (권장)

precheck 성공 후 자동으로 production 실행:

```bash
PRECHECK_JOB=$(qsub config/run_pre_qsub_checks.pbs)
qsub -W depend=afterok:${PRECHECK_JOB} config/run_production.pbs
```

#### 직접 제출

```bash
qsub config/run_production.pbs                                 # 자동 이어하기
qsub -v MODE=force config/run_production.pbs                   # 전체 재실행
qsub -v MODE=from,FROM=4 config/run_production.pbs             # Phase 4부터
qsub -v MODE=status config/run_production.pbs                  # 상태 확인
qsub -v MODE=vina-only config/run_production.pbs               # Vina lane만 (Phase 1,4,5,6,7)
qsub -v MODE=ppi-only config/run_production.pbs                # PPI lane만 (Phase 2,3)
qsub -v MODE=post-only config/run_production.pbs               # Phase 4부터 (후처리)
qsub -v SKIP_PRECHECK_GUARD=1 config/run_production.pbs        # precheck 가드 우회 (의도적으로만)
```

#### Lane별 PBS 스크립트 (개별 제출)

`run_production.pbs` 대신 lane별로 개별 제출할 수 있습니다:

```bash
qsub config/run_vina_cpu.pbs              # Vina CPU docking
qsub config/run_ppi_state_seed.pbs        # PPI docking (state/seed별)
qsub config/run_ppi_postprocess.pbs       # PPI 후처리
qsub config/run_vina_postprocess.pbs      # Vina 후처리
qsub config/run_finalize.pbs              # Verdict + Report + Validate
```

#### LightDock 2차 검증 (Phase 1)

PyRosetta가 Phase 1 primary evidence이며, LightDock은 secondary validation입니다:

```bash
qsub config/run_lightdock.pbs             # 전체 state
qsub config/run_lightdock_test.pbs        # 테스트
```

### Step 3: run_production.py 직접 실행 (서버 셸에서)

PBS 없이 직접 실행하는 경우 (디버깅/개발용):

```bash
python run_production.py                  # 자동 이어하기 (완료된 Phase 스킵)
python run_production.py --force          # 전체 재실행
python run_production.py --from 4         # Phase 4부터 실행
python run_production.py --status         # 상태 확인만 (실행 안 함)
python run_production.py --only 2,3       # 특정 Phase만 실행
```

### Step-Folder 모드 의미

- `--status`: 읽기 전용. step 폴더를 재생성하지 않음.
- `--from N`: Phase N 이상 재실행. 이전 step 폴더는 유지.
- `--only N[,M]`: 지정된 Phase만 재실행.
- `--force`: 기존 출력이 있어도 재실행 후 step view 갱신.
- 클린 재실행: `python scripts/reset_production_outputs.py --execute` 후 재제출.

---

## Workflow B: Advanced PPI-First Pipeline

PPI 증거를 먼저 확보한 후 포켓 제안 -> 다양성 도킹 -> 교란 스코어링까지 수행하는 고급 워크플로우입니다.

**전제 조건**: Workflow A의 Phase 2 (PPI docking) 완료 (`seed_complete.json` 존재)

### 자동 체인 제출

```bash
qsub config/run_advanced_pipeline.pbs                                          # 전체 (Phase 1~4)
qsub -v ADV_FROM=2 config/run_advanced_pipeline.pbs                            # Phase 2부터
qsub -v ADV_ROUNDS=5 config/run_advanced_pipeline.pbs                          # Phase 3를 5 rounds
qsub -v ADV_FROM=3,ADV_ROUNDS=2 config/run_advanced_pipeline.pbs              # Phase 3부터, 2 rounds
```

### Phase 구성

| Phase | PBS 스크립트 | 내용 |
|-------|-------------|------|
| 1 | `run_adv_phase1.pbs` | PPI 분석 |
| 2 | `run_adv_phase2.pbs` | 포켓 분석 |
| 3-setup | `run_adv_phase3_setup.pbs` | 다양성 도킹 준비 |
| 3-execute | `run_adv_phase3_execute.pbs` (x N rounds) | 다양성 도킹 실행 |
| 3-post | `run_adv_phase3_post.pbs` | 다양성 도킹 후처리 |
| 4 | `run_adv_phase4.pbs` | 교란 관련성 스코어링 |

### 개별 Phase 수동 제출

의존성 체인을 수동으로 구성할 때:

```bash
ADV1=$(qsub config/run_adv_phase1.pbs)
ADV2=$(qsub -W depend=afterok:$ADV1 config/run_adv_phase2.pbs)
ADV3S=$(qsub -W depend=afterok:$ADV2 config/run_adv_phase3_setup.pbs)
ADV3E0=$(qsub -v ROUND=0 -W depend=afterok:$ADV3S config/run_adv_phase3_execute.pbs)
ADV3E1=$(qsub -v ROUND=1 -W depend=afterok:$ADV3E0 config/run_adv_phase3_execute.pbs)
ADV3E2=$(qsub -v ROUND=2 -W depend=afterok:$ADV3E1 config/run_adv_phase3_execute.pbs)
ADV3P=$(qsub -W depend=afterok:$ADV3E2 config/run_adv_phase3_post.pbs)
qsub -W depend=afterok:$ADV3P config/run_adv_phase4.pbs
```

---

## 개별 명령 참고

### main.py CLI 명령어

| Command | 용도 | 주요 출력 위치 |
|---------|------|---------------|
| `vina` | Vina blind docking | `output/workflow_a/phase1_vina_docking/` |
| `postprocess` | 포즈 파싱, 포켓 요약, receptor state 비교 | `output/workflow_a/phase4_vina_postprocess/` |
| `pyrosetta` | PyRosetta Phase 1 실행 | `output/workflow_a/phase2_ppi_docking/` |
| `ppi-postprocess` | PPI 후처리 출력 재생성/추출 | `output/workflow_a/phase3_ppi_postprocess/` |
| `verdict` | 사이트 판정 | `output/workflow_a/phase5_verdict/valid_sites.csv` |
| `report` | 텍스트 및 통합 증거 보고서 | `output/workflow_a/phase6_report/project_report.txt` |
| `validate` | 현재 출력 상태 검증 | `output/workflow_a/phase7_validation/` |
| `full` | 기본 통합 CLI 경로 (Vina 중심 baseline) | `output/workflow_a/` |
| `organize` | Step별 출력 정리 | `output/workflow_a/` (phase별 디렉토리) |

참고: `full`은 Vina 중심 routine baseline에 해당하며, 과학적 Phase 1->4 전체 경로를 의미하지 않습니다.

### 수동 명령 시퀀스

**Routine baseline (Vina 중심):**

```bash
python main.py -c config/example-project.yaml vina
python main.py -c config/example-project.yaml postprocess
python main.py -c config/example-project.yaml verdict
python main.py -c config/example-project.yaml report
python main.py -c config/example-project.yaml validate
```

**Phase 1 집중 (receptor-side evidence):**

```bash
python main.py -c config/example-project.yaml pyrosetta
python main.py -c config/example-project.yaml ppi-postprocess
```

**단일 명령 전체 실행:**

```bash
python main.py -c config/example-project.yaml full
```

### 인터랙티브 모드

```bash
python main.py
```

사용 가능한 명령 목록을 확인할 수 있습니다.

### 기타 PBS 스크립트

```bash
qsub config/run_vina_gpu.pbs              # Vina GPU docking
qsub config/run_phase2_cascade.pbs        # Phase 2 cascade
qsub config/run_phase3_gpu_round.pbs      # Phase 3 GPU round
qsub config/run_phase3_cpu_round.pbs      # Phase 3 CPU round
qsub config/run_full_test.pbs             # 전체 테스트
qsub config/run_production_fresh.pbs      # 클린 프로덕션
```

---

## 로그 확인 및 모니터링

### 로그 파일 체계

파이프라인은 3-계층 로거(`pipeline`, `pipeline.worker`, `pipeline.vina`)를 사용합니다.
설정 모듈: `egfr_pipeline/pyrosetta_docking/logging_config.py`

| 로그 파일 | 위치 | 레벨 | 용도 |
|-----------|------|------|------|
| `PROGRESS.log` | `phase2_ppi_docking/{state}/prod_seed{n}/` | INFO | 실시간 진행 상황 (`tail -f`) |
| `pipeline.log` | `phase2_ppi_docking/{state}/prod_seed{n}/logs/` | DEBUG | 상세 파이프라인 로그 |
| `workers.log` | `phase2_ppi_docking/{state}/prod_seed{n}/logs/` | DEBUG | 멀티프로세스 워커 로그 |
| `LOG_INDEX.txt` | `workflow_a/logs/` | — | 로그 인덱스 + 읽기 가이드 |
| PBS stdout | `production.o{jobid}` (PBS 작업 디렉토리) | — | 전체 잡 출력 + 완료 요약 |

### 실시간 모니터링

```bash
# 특정 시드 실시간 추적
tail -f output/workflow_a/phase2_ppi_docking/3GT8_raw/prod_seed0/PROGRESS.log

# 모든 시드 진행 상황 한 눈에
for f in output/workflow_a/phase2_ppi_docking/*/prod_seed*/PROGRESS.log; do
  echo "=== $(basename $(dirname $(dirname $f)))/$(basename $(dirname $f)) ==="
  tail -1 "$f"
done
```

### 에러 진단

```bash
# 로그 인덱스 확인 (심볼릭 링크 기반)
cat output/workflow_a/logs/LOG_INDEX.txt

# 모든 시드에서 에러 검색
grep 'ERROR\|CRITICAL' output/workflow_a/logs/ppi_seeds/*_workers.log

# 필터 단계 추적
grep 'FilterStage' output/workflow_a/logs/ppi_seeds/*_pipeline.log

# Vina 모듈 로그 (프로덕션 모드에서 pipeline.vina 로거로 기록)
grep 'pipeline.vina' output/workflow_a/logs/*.log 2>/dev/null
```

### PBS 잡 출력 확인

```bash
# 완료 요약 (run_production.py가 종료 시 출력)
tail -50 production.o*

# 잡 상태 확인
qstat -u $USER
```

---

## 결과 확인

### 해석 시작점

프로덕션 완료 후 `output/workflow_a/phase6_report/project_report.txt`부터 시작합니다.

권장 읽기 순서:

1. `output/workflow_a/phase6_report/project_report.txt`
2. `output/workflow_a/phase5_verdict/valid_sites.csv`
3. `output/workflow_a/phase4_vina_postprocess/vina_pocket_table.csv`
4. `output/workflow_a/phase3_ppi_postprocess/ppi_pyrosetta_residues.csv`
5. `output/workflow_a/phase2_ppi_docking/{state}/prod_seed{n}/final_ranking.csv`

Canonical 출력은 `output/workflow_a/` 아래 phase별 디렉토리에 있으며 source of truth입니다.

### 출력 위치 매핑

| 명령 계열 | 확인 위치 |
|-----------|----------|
| `vina` | `output/workflow_a/phase1_vina_docking/` |
| `postprocess` | `output/workflow_a/phase4_vina_postprocess/` |
| `verdict`, `report`, `validate`, `full` | `output/workflow_a/phase5_verdict/`, `phase6_report/`, `phase7_validation/` |
| `pyrosetta`, `ppi-postprocess` | `output/workflow_a/phase2_ppi_docking/`, `phase3_ppi_postprocess/` |
| pre-qsub PBS | `output/precheck/` |
| advanced phases | `output/workflow_b/phase2_pocket_analysis/`, `phase3_focused_docking/`, `phase4_scoring/` (해당 lane이 범위 내일 때만) |

### 참조 출력 레이아웃

```text
output/
├── workflow_a/
│   ├── phase1_vina_docking/{receptor_id}/
│   ├── phase2_ppi_docking/{state}/prod_seed{n}/
│   ├── phase3_ppi_postprocess/
│   │   ├── ppi_pyrosetta_residues.csv
│   │   └── ppi_pyrosetta_summary.csv
│   ├── phase4_vina_postprocess/
│   │   ├── vina_pose_table.csv
│   │   ├── vina_pocket_table.csv
│   │   └── vina_drug_pocket_map.csv
│   ├── phase5_verdict/
│   │   ├── valid_sites.csv
│   │   └── cross_method_agreement.csv
│   ├── phase6_report/
│   │   ├── project_report.txt
│   │   └── combined_residue_evidence.csv
│   ├── phase7_validation/
│   └── logs/
├── workflow_b/
│   ├── phase1_ppi_analysis/
│   ├── phase2_pocket_analysis/
│   ├── phase3_focused_docking/
│   └── phase4_scoring/
└── precheck/
```

### 필수 체크포인트

| 단계 | 필수 체크포인트 | 중단 조건 |
|------|---------------|----------|
| Precheck | `output/precheck/last_pass.json` 존재 및 pass | 미존재 또는 실패 |
| Routine Vina postprocess | 핵심 Vina 테이블이 생성됨 | pose/pocket 요약 테이블 미존재 |
| Routine integration | `valid_sites.csv`, `cross_method_agreement.csv`, `project_report.txt` 존재 | 최종 판정 파일 미존재 또는 오래됨 |
| Phase 1 PyRosetta | Phase 1 residue, patch, review 출력 존재 | 구조화된 Phase 1 export 없음 |
| Phase 1 LightDock | LightDock 요청 시 cross-method 수렴 출력 존재 | LightDock 미완료인데 증거로 인용 |

### 해석 규칙

- Receptor state 분리를 모든 리뷰 단계에서 유지
- PyRosetta = Phase 1 primary evidence
- LightDock = independent secondary validation (standalone primary가 아님)
- AFM = 비활성 (명시적 재활성화 전까지)
- `verdict`, `report`, `validate` = routine baseline의 최종 해석 레이어
- Advanced Phase 4 perturbation 출력을 기본 최종 레이어로 승격하지 말 것 (명시적 Phase 4 작업이 아닌 한)

---

## 흔한 실수

- `-c config/example-project.yaml`을 빠뜨리고 의도하지 않은 config로 실행
- `full`이 과학적 4-Phase 전체 계획을 실행한다고 오해
- AFM이 비활성인데 AFM 의존 작업 실행
- 머신 코어 수를 보고 routine safe worker bound(16)를 초과
- `output/workflow_a/`의 pointer stub 파일을 실제 payload로 착각
- 과학적 Phase 번호와 production stage 번호를 혼동
- 3개 receptor state를 요약/핸드오프 해석에서 너무 일찍 병합

## 중단 및 에스컬레이션 조건

다음 중 하나라도 해당되면 중단하고 문제를 해결합니다:

- pre-qsub validation 실패
- receptor state 등록이 불완전하거나 모호
- 현재 미설정된 AFM 입력에 의존하는 실행
- PyRosetta 지원 없이 LightDock만으로 primary evidence 인용
- 예상 결과 파일이 없는데 하위 해석이 이미 시작됨
- 출력 파일이 현재 payload가 아닌 오래된 pointer로 보임
