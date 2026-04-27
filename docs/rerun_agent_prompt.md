# EGFR-MYO1D PPI 도킹 파이프라인 — 재실행 설계 에이전트 프롬프트

## 역할

너는 EGFR-MYO1D protein–protein interaction 연구의 **computational experiment designer**야. 이 프로젝트의 전체 도킹 파이프라인을 처음부터 다시 실행할 수 있도록 워크플로우를 완벽하게 다듬어야 해.

## 프로젝트 배경

EGFR kinase domain (residues 699-1007)과 MYO1D TH1 domain (residues 955-1006, AlphaFold predicted beta-meander)의 PPI 결합 부위를 규명하고, 그 근처에서 druggable pocket을 찾는 연구야.

### 핵심 실험 도구
- **PyRosetta**: 전역 PPI 도킹 (RosettaDock, REF2015 score function)
- **LightDock**: 독립 교차검증 (DFIRE2 scoring)
- **fpocket**: 표면 포켓 탐지
- **AutoDock Vina**: 소분자 도킹 (blind + focused)
- **GROMACS**: MD 시뮬레이션 (CHARMM36, 200ns) — 이미 완료됨, 재실행 불필요

### 3개 EGFR 구조 상태
1. `3GT8_raw` — X-ray crystal structure (PDB: 3GT8, heavy atoms only)
2. `EGFR_160-185` — MD cluster representative (38-48ns window)
3. `EGFR_170-200` — MD cluster representative (85-100ns window)

### MYO1D 파트너
- AlphaFold model, TH1 domain
- Extended beta-meander: residues 955-1006 (52 residues)
- Active face: sheets 8 (961-964) + 9 (968-972) — Ko et al. alanine scanning으로 실험적 검증됨
- Sheet 12 (993-997): structural support (직접 결합 아님)

### 리간드 3종
- 173940, 97806, VAX-C12_0
- Pairwise Tanimoto < 0.4 (구조 다양성)
- 교수님 제공, 임의 교체 금지

## 현재 문제점 (재실행 이유)

### 실험 설계 결함
기존 실행에서 seed 0-4와 seed 5-9가 **두 변수가 동시에 달랐음**:

| Seed | Receptor 구조 | Partner 범위 | construct_type |
|------|--------------|-------------|----------------|
| 0–4 | EGFR **dimer** (2량체→1체인 병합, chain B +1000 offset) | 960–1006 | dimer_offset |
| 5–9 | EGFR **monomer** (단량체) | **955**–1006 (extended) | full_kinase_domain |

**문제**: receptor construct(dimer vs monomer)와 partner range(960 vs 955)가 동시에 바뀌어서, 결과 차이의 원인을 분리할 수 없음. 이건 실험 설계의 기본 원칙 위반이야.

### 해결 방향
**모든 seed에서 partner를 동일하게 extended beta-meander (955-1006)로 통일**하고, dimer/monomer만 변수로 남겨야 함.

## 사용 가능한 컴퓨팅 자원

- **노드 3개**, 각 32 cores
- PBS/qsub 기반 job submission
- 노드당 2 jobs 동시 실행 가능 (job당 16 cores)
- 동시 6 jobs → 30 jobs 기준 5 라운드

## 코드베이스 구조

```
codex_ligand/
├── egfr_pipeline/
│   ├── pyrosetta_docking/     # PyRosetta 도킹 엔진 (공유)
│   │   ├── pipeline_manager.py  # 메인 도킹 파이프라인
│   │   ├── movers.py            # RosettaDock movers
│   │   ├── run_metadata.py      # construct_type 추론
│   │   └── pyrosetta_init.py
│   ├── ppi/
│   │   ├── prepare_dimer_pdb.py # EGFR dimer → 1체인 병합 (CHAIN_B_OFFSET=1000)
│   │   ├── postprocess_ppi.py   # PPI 결과 후처리
│   │   └── pyrosetta_extract.py # 인터페이스 잔기 추출
│   ├── phase1/                  # WF-B Phase 1: PPI 분석
│   │   ├── orientation_filter.py  # PCA + dot product 필터
│   │   ├── cluster_consensus.py   # hotspot 동정 (threshold ≥0.50)
│   │   └── launch_docking.py
│   ├── phase2/                  # WF-B Phase 2: Pocket 분석
│   │   └── patch_relationship.py  # pocket-PPI 관계 분류
│   ├── phase3/                  # WF-B Phase 3: Focused Vina
│   ├── phase4/                  # WF-B Phase 4: Perturbation scoring
│   ├── vina/                    # WF-A: Blind Vina
│   ├── verdict.py               # WF-A Phase 5
│   ├── report.py                # WF-A Phase 6
│   ├── validate.py              # WF-A Phase 7
│   └── paths.py                 # 전체 경로 관리 (수정 주의)
├── config/phase1/               # 도킹 INI configs
├── input/
│   ├── receptors/               # 3GT8_raw.pdb, EGFR_160-185.pdb, EGFR_170-200.pdb
│   └── PPI/                     # MYO1D PDB, docking pair metadata
├── scripts/                     # PBS scripts, extraction scripts
└── docs/                        # 문서, 논문 초안
```

### 핵심 파일들

1. **`egfr_pipeline/ppi/prepare_dimer_pdb.py`** — dimer PDB 준비
   - `prepare_dimer_partner()`: EGFR dimer chain A+B → chain A 하나로 병합 (+1000 offset on chain B) → MYO1D를 chain B로 추가
   - `restore_dimer_result()`: 도킹 결과에서 원래 2체인 복원
   - `CHAIN_B_OFFSET = 1000`

2. **`config/phase1/phase1_prod_{state}_seed{n}.ini`** — 도킹 config
   - `[Docking] total_global_models = 20000`
   - `[System] n_cpus = 16`
   - `[Constraints] excluded_residues_a = 709-720,724-731,736-739,747,783-785,799-805,871-873,917-921`

3. **`egfr_pipeline/phase1/orientation_filter.py`** — 후처리 필터
   - `ACTIVE_FACE_RESIDUES = sheets 8+9 (961-972)`
   - `AMBIGUOUS_BAND = 0.10`

4. **`input/PPI/phase1/docking_pair_metadata.csv`** — 도킹 쌍 정의

## 해야 할 작업 목록

### Phase 0: 입력 구조 준비
1. **Monomer input PDB 검증** — 기존 `input/receptors/*.pdb` 3개가 정상인지 확인
2. **Extended MYO1D (955-1006) PDB 확보** — `input/PPI/` 에서 정확한 파일 찾기/검증
3. **Dimer input PDB 생성** — `prepare_dimer_pdb.py`로 각 state의 dimer + extended partner PDB 생성
4. **Monomer input PDB 생성** — monomer + extended partner PDB 생성 (chain A: EGFR monomer, chain B: MYO1D)

### Phase 1: Config 생성
5. **INI config 30개 생성** — 3 states × 10 seeds
   - Seed 0–4: dimer construct → dimer input PDB
   - Seed 5–9: monomer construct → monomer input PDB
   - **Partner는 모두 동일 (955-1006)**
   - Random seed: 재현 가능한 deterministic 값
   - `n_cpus = 16`

### Phase 2: PBS 스크립트
6. **PBS 스크립트 생성** — 3 노드에 최적 배분
   - 노드당 2 jobs 동시 (16 cores × 2 = 32)
   - 5 라운드로 30 jobs 완료
   - 의존성 관리: 라운드 간 순서 보장 또는 수동 실행 가이드
   - 에러 시 개별 재실행 가능하도록 모듈화

### Phase 3: 후처리 파이프라인 검증
7. **orientation_filter.py** — AMBIGUOUS_BAND=0.10 검증 계획 포함
8. **cluster_consensus.py** — hotspot threshold 0.50 유지
9. **postprocess_ppi.py** — construct_type별 분리 분석 + 통합 분석 모두 가능하도록
10. **LightDock 교차검증** — 동일 입력 구조로 재실행 스크립트

### Phase 4: Downstream 워크플로우 (PPI 완료 후)
11. **fpocket** — 새 결과 기반 포켓 탐지
12. **Vina blind (WF-A)** — 재실행 또는 기존 결과 재활용 판단
13. **Vina focused (WF-B)** — 새 PPI patch 기반
14. **Perturbation scoring** — 4축 스코어 재계산

### Phase 5: 검증 + 문서화
15. **validate.py** 실행
16. **construct_type별 비교 분석** — dimer vs monomer 결과 차이 정량화
17. **논문 Methods 섹션 업데이트** — 새 실험 설계 반영

## 절대 규칙

1. **ATP 포켓 배제**: is_atp_site=True → STRONG 판정 금지
2. **paths.py 수정 금지**: 수정하면 전체 Phase 경로가 바뀜
3. **스코어링 가중치 변경 승인**: verdict.py, score_framework.py 가중치는 사람 승인 필요
4. **잔기 번호**: PDB author numbering 사용 (0-based 금지)
5. **리간드 교체 금지**: 3종 리간드는 교수님 지정
6. **CSV 스키마 보존**: 기존 CSV 컬럼명 변경 시 하위 Phase ingestion 코드 동시 수정
7. **output/ 경로 접근 불가**: 이 환경에서는 HPC output에 접근 불가. 스크립트 생성 + 명령어 안내만 가능

## 기대 결과물

1. **실행 계획서** — 전체 워크플로우 단계별 문서
2. **입력 PDB 준비 스크립트** — dimer/monomer × 3 states
3. **Config INI 30개** — 일관된 설정
4. **PBS 스크립트** — 3 노드 최적 배분, 라운드별 실행
5. **후처리 실행 가이드** — PPI 결과 나온 후 순서대로 실행할 명령어 목록
6. **검증 체크리스트** — 각 단계에서 확인해야 할 항목

## HPC 환경 정보

- 경로: `/work4/hwang/onepack/my_second_project/codex_ligand`  (또는 codex_ligand2)
- PBS/qsub, conda activate pyrosetta
- 노드: 3개 × 32 cores
- Python 3.9, PyRosetta 4, AutoDock Vina 1.2, fpocket 4.0, LightDock 0.9

## 참고 문서 위치

- `docs/phase1_notes.md` — Phase 1 설계 근거, sampling rationale, orientation filter
- `docs/CONTEXT.md` — 버그 이력, 결정 로그
- `docs/runbook.md` — 실행 가이드
- `docs/methodology_limitations.md` — 방법론 한계점
- `CLAUDE.md` — 프로젝트 규칙, 워크플로우 구분
- `.claude/skills/` — 도메인별 스킬 파일
