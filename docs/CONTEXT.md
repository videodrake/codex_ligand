# 프로젝트 컨텍스트

## 현재 작업 상태
- 워크플로우: Workflow A (Phase 2 PPI seed 5~9 HPC 실행 중)
- 현재 작업: HPC 도킹 결과 대기 중
- 다음 작업: **HPC 결과 수집 → Phase 3 postprocess → Phase 4 vina postprocess → Phase 5 verdict**

## 다음 세션에서 해야 할 일 (에이전트 필독)

### 1단계: HPC 결과 확인 — 사용자에게 물어볼 것

사용자에게 아래 명령어 결과를 요청한다:

```bash
# HPC에서 실행 (codex_ligand2 디렉토리)
# PPI seed 완료 상태 확인
find output/workflow_a/phase2_ppi_docking -name "seed_complete.json" | sort

# Vina 결과 유무 확인
ls output/workflow_a/phase1_vina_docking/3GT8_raw/
```

**완료 기대치:**
- PPI: 3 상태 × 10 seeds = 30개 seed_complete.json (기존 seed 0~4 = 16개 완료, seed 5~9 = 15개 실행 중. 3GT8_raw seed5는 _capped 폐기 후 재실행)
- Vina: Phase 1 결과 유무 확인 필요 (아직 미확인)

### 2단계: 결과 가져오기

PPI 결과가 완료되면, 이 환경으로 CSV만 가져온다:

```bash
# HPC에서 실행
cd /work4/hwang/onepack/my_second_project/codex_ligand2
tar czf /tmp/ppi_results.tar.gz \
  $(find output/workflow_a/phase2_ppi_docking -name "*.csv" -o -name "*.json")
```

파일 크기 수 MB 이내. PDB 파일(수 GB)은 불필요.

### 3단계: 결과 수집 후 진행할 파이프라인

| 순서 | Phase | 명령어 | 선행 조건 |
|------|-------|--------|----------|
| 1 | Phase 3: PPI Postprocess | `python run_production.py --only 3` | PPI seed 전부 완료 |
| 2 | Phase 4: Vina Postprocess | `python run_production.py --only 4` | Phase 1 Vina 완료 |
| 3 | Phase 5: Verdict | `python run_production.py --only 5` | Phase 3 + 4 완료 |
| 4 | Phase 6: Report | `python run_production.py --only 6` | Phase 5 완료 |
| 5 | Phase 7: Validate | `python run_production.py --only 7` | Phase 6 완료 |

**Vina가 아직 안 돌았다면**: PBS 스크립트 생성 필요 (`config/run_vina_cpu.pbs`)

### 4단계: 결과 분석 후 처리할 투두 항목

결과가 나온 후 진행 가능한 항목 (투두리스트.md 참조):
- Orientation filter AMBIGUOUS_BAND 검증 (dot product 분포 분석)
- Sheet 12 인터페이스 출현 빈도 재확인
- `pending_*` 리간드 지지 수준 재점수 (Vina 결과 필요)
- Ko et al. sheet 8/9 잔기 3개 이상 검증
- Cross-method Jaccard 리뷰
- Centroid 거리 편향 임계값 리뷰

## 작업 로그
- [2026-04-06] Phase 0-1: 레거시 파일 삭제 완료 — 삭제 3개, 미존재 12개
- [2026-04-06] Phase 0-2: 레거시 참조 정리 — README.md 2개 행 삭제, docs/README.md 2개 행 삭제, tests/test_nightly_review.py 삭제, test_e2e_group7.py 미존재 확인
- [2026-04-06] Phase 0-3: docs/ 스캔 완료 — 활성 7개, 아카이브 35개, paths.py DEPRECATED 없음, input/PPI/prepared/ 활성 참조 4개소 확인
- [2026-04-06] Phase 0-4: CLAUDE.md 재작성 — 설계서 섹션 3 기반 68줄, 3 File System 내용 전부 제거
- [2026-04-06] Phase 1-1: CLAUDE.md 재작성 완료
- [2026-04-06] Phase 1-2: CONTEXT.md 재작성 완료
- [2026-04-06] Phase 2-1: 스킬 생성 — phase-dependencies, bug-history (2/7)
- [2026-04-06] Phase 2-2: 스킬 생성 — ppi-analysis, vina-docking, hpc-operations, scoring-system, testing (7/7 완료)
- [2026-04-06] Phase 3-1: 에이전트 3개 생성 — pipeline-dev, reviewer, science-qa (권한 경계 원칙 포함)
- [2026-04-06] Phase 3-2: 훅 2개 생성 — pre-commit.sh, csv-schema-guard.py
- [2026-04-06] Phase 4: 하네스 구축 완료 — README.md 역할별 재구성, 참조 무결성 검증 통과, 스킬 7 + 에이전트 3 + 훅 2 확인

## 최근 결정 사항
- [2026-04-06] 하네스 엔지니어링 적용 시작. 설계서: harness_engineering_design.md
- [2026-04-06] PRODUCTION_N_SEEDS 5→10 확장. generate_configs.py와 run_production.py 양쪽 동기화 필수
- [2026-04-06] HPC 환경을 codex_ligand2로 이전. output/workflow_a와 input/은 codex_ligand(원본)에서 심볼릭 링크로 연결

## 발견된 이슈 (미해결)
- (없음)

## 해결된 이슈
- [2026-04-06] input/PPI/prepared/ 활성 참조 4개소 정리 완료 — register_pilot_data() 제거(dead code), main.py/postprocess_ppi.py/test 경로 업데이트, data_inventory.md 반영

## HPC 배포 주의사항
- codex_ligand2 (최신 repo)에서 실행 시 input/과 output/workflow_a를 codex_ligand(원본)에서 심볼릭 링크해야 함. `mkdir -p output && ln -s ...`
- `_capped` PDB는 일반 `ext_beta_meander` PDB와 **다른 구조** (Chain B 원자 수 상이: 5362 vs 5653). 혼용하면 결과 비일관. 3GT8_raw seed5 `_capped` 결과 폐기 후 일반 PDB로 재실행
- PRODUCTION_N_SEEDS 변경 시 generate_configs.py와 run_production.py **양쪽 모두** 수정해야 함. 한쪽만 바꾸면 config INI 부재로 RuntimeError 발생
- seed별 PBS 제출: `qsub -v STATE=$STATE,SEED=$SEED config/run_ppi_state_seed.pbs` (node 지정: `-l nodes=node06:ppn=16`)

## 실패 패턴 (반복 방지)
- [2026-04-06] CSV 컬럼 추가 시 하위 전파 판단 누락 + DoD 미확인 커밋
- [2026-04-06] BUG-001 관련 변경 시 버그 번호 명시적 인용 누락 + 디렉토리 매핑 누락(pyrosetta_docking/)
- [2026-04-06] 스코어링 가중치 변경을 사람 승인 없이 코드 수정 + 커밋. 원인: 권한 경계 원칙이 scoring-system 스킬에 명시되지 않았음. 해결: pipeline-dev, scoring-system, CLAUDE.md 3곳에 승인 필요 명시
- [2026-04-06] DoD smoke test를 '영향 없음' 판단으로 생략. 원인: DoD 문구가 조건부 생략을 허용하는 것처럼 읽힘. 해결: 무조건 실행으로 문구 강화
- [2026-04-06] PRODUCTION_N_SEEDS=5 상태에서 seed 5~9 실행 시도 → config INI 부재로 "completed without valid outputs" RuntimeError. 원인: generate_configs.py와 run_production.py의 N_SEEDS 불일치. 해결: 양쪽 10으로 통일 + seed 5~9 INI 생성
- [2026-04-06] codex_ligand2 클론 후 input/ 미링크 상태에서 PBS 제출 → FileNotFoundError (TH1 domain.pdb). 원인: git repo에 대용량 PDB 미포함. 해결: 원본 repo에서 심볼릭 링크
- [2026-04-06] 3GT8_raw seed5에 _capped PDB(5653 atoms) 사용 → 일반 PDB(5362 atoms)와 비일관. 원인: 이전 코드에서 다른 입력 PDB 이름 사용. 해결: seed5 결과 폐기 후 일반 PDB로 재실행
