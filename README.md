# EGFR–MYO1D Pipeline

EGFR kinase domain(C-lobe) 표적 탐색을 위한 통합 계산 파이프라인.
세 가지 수용체 상태(3GT8_raw, 3GT8_cl38_48, 3GT8_cl85_100)에 대해 리간드 도킹, PPI 도킹, 포켓 클러스터링, 잔기 비교를 수행한다.

---

## 빠른 시작

```bash
# 인터랙티브 메뉴
python main.py

# CLI 서브커맨드
python main.py vina -c config/example-project.yaml
python main.py postprocess -c config/example-project.yaml
python main.py report -c config/example-project.yaml
python main.py validate -c config/example-project.yaml
python main.py full -c config/example-project.yaml    # 전체 파이프라인
```

---

## 디렉토리 구조

```
codex_ligand/
├── main.py                           # 통합 CLI 진입점 (인터랙티브 메뉴 + 서브커맨드)
├── CLAUDE.md                         # Claude Code 프로젝트 지침
├── README.md
│
├── egfr_pipeline/                    # 핵심 패키지
│   ├── config.py                     # 공유 설정 로딩 (YAML/JSON)
│   ├── residue_utils.py              # 잔기 정규화 (HSD→HIS 등)
│   ├── report.py                     # 종합 보고서 생성
│   ├── validate.py                   # 출력 검증 (스키마, ID, 잔기 일관성)
│   │
│   ├── vina/                         # AutoDock Vina 리간드 도킹
│   │   ├── dock.py                   # 도킹 실행 (blind/focused)
│   │   ├── parse_poses.py            # 결과 파싱 → vina_pose_table.csv
│   │   ├── contacts.py               # 수용체 접촉 잔기 추출
│   │   ├── cluster.py                # 포켓 클러스터링 (centroid greedy)
│   │   ├── summarize.py              # 포켓/리간드 요약 테이블
│   │   └── compare.py                # 교차 수용체 포켓 비교
│   │
│   ├── ppi/                          # PPI 잔기 표준화
│   │   ├── pyrosetta_extract.py      # PyRosetta 결과에서 인터페이스 잔기 추출
│   │   └── afm_extract.py            # AlphaFold-Multimer 결과 처리
│   │
│   ├── pyrosetta_docking/            # PyRosetta PPI 글로벌 도킹 (v2.0)
│   │   ├── pipeline_manager.py       # 7단계 파이프라인 오케스트레이터
│   │   ├── docking.py                # Relax, Global Docking, Refinement 워커
│   │   ├── analysis.py               # Scoring, RMSD, Interface 분석
│   │   └── common.py                 # PyRosetta 초기화, Pose↔String 변환
│   │
│   └── md/                           # GROMACS MD 분석
│       ├── gromacs_analysis.py        # 궤적 분석
│       └── ligand_contacts.py         # 리간드 접촉 분석
│
├── config/                           # 프로젝트 설정 파일
│   └── example-project.yaml          # Vina 프로젝트 config 예시
│
├── input/                            # 실제 입력 데이터
│   ├── receptors/                    # 수용체 PDB (3개)
│   │   ├── 3GT8_raw.pdb              # Crystal (chain A, 699-1007)
│   │   ├── 3GT8_cl38_48.pdb          # MD cluster 38-48ns (chain X, 634-1014)
│   │   └── 3GT8_cl85_100.pdb         # MD cluster 85-100ns (chain X, 634-1014)
│   └── ligands/                      # 리간드 SDF (3개)
│       ├── 173940_ligand.sdf
│       ├── 97806_ligand.sdf
│       └── VAX-C12_0_ligand.sdf
│
├── output/                           # 파이프라인 출력 (런타임 생성)
├── docs/                             # 기획 문서, 매뉴얼, 가이드
├── legacy/                           # 리팩터링 이전 원본 스크립트 보관
├── smoke_test/                       # 기능 검증용 테스트 데이터
└── tests/                            # 테스트 코드
```

---

## 파이프라인 구성

### 1. AutoDock Vina 리간드 도킹 (`egfr_pipeline/vina/`)

수용체-리간드 도킹 및 후처리 체인:

```
Vina Docking (dock.py)
  → 결과 파싱 (parse_poses.py) → vina_pose_table.csv
  → 접촉 잔기 추출 (contacts.py)
  → 포켓 클러스터링 (cluster.py)
  → 포켓 요약 (summarize.py) → vina_pocket_table.csv, vina_drug_pocket_map.csv
  → 교차 수용체 비교 (compare.py) → vina_pocket_comparison.csv
```

- **Blind / Focused** 도킹 모드 지원
- 수용체 순차 / 리간드 병렬 (max_workers 설정)
- 포켓 클러스터링: centroid 기반 greedy assignment
- 교차 비교: Jaccard, overlap coefficient, centroid 거리, same_patch_candidate 판정

### 2. PyRosetta PPI 글로벌 도킹 (`egfr_pipeline/pyrosetta_docking/`)

2-chain PDB에 대한 7단계 PPI 도킹:

```
Relax → Global Docking → Scoring & Filtering (v2.0/v1.0 자동 분기)
  → L_RMSD Clustering → Refinement → Final Scoring → Visualization
```

- v2.0: 2-Pass 설계, Mini Refinement, Graduated Fallback (Level 0~3)
- v1.0: 레거시 호환 (config에 `[FilterStage1]` 없으면 자동)
- 실행: `python main.py pyrosetta` 또는 직접 `python -m egfr_pipeline.pyrosetta_docking.pipeline_manager config.ini input.pdb`

### 3. PPI 잔기 표준화 (`egfr_pipeline/ppi/`)

PyRosetta / AlphaFold-Multimer 결과에서 인터페이스 잔기를 추출하고 정규화:

- CHARMM 잔기명 변환 (HSD/HSE/HSP→HIS, CYX→CYS)
- Chain ID 차이 정규화 (A/B vs X)
- 출력: `ppi_pyrosetta_residues.csv`, `ppi_pyrosetta_summary.csv`, `ppi_afm_residues.csv`

### 4. 보고서 & 검증

- **Report** (`egfr_pipeline/report.py`): Vina 포켓 + PPI 잔기 증거를 결합한 종합 보고서
  - `project_report.txt` + `combined_residue_evidence.csv`
- **Validate** (`egfr_pipeline/validate.py`): 출력 무결성 검증
  - 파일 존재, CSV 스키마, ID 일관성, 잔기 번호, handoff 준비 상태
  - Exit code: 0=PASS, 1=WARN, 2=FAIL

---

## 설정

프로젝트 설정은 `config/example-project.yaml` 참조:

```yaml
project_name: egfr_myo1d_vina
output_root: ./output
max_workers: 16

receptors:
  - id: 3GT8_raw
    pdb: input/receptors/3GT8_raw.pdb
    pdbqt: input/3GT8_raw_receptor.pdbqt
  # ...

ligands:
  - id: ligand_a
    pdbqt: input/ligand_a_ligand.pdbqt
  # ...

vina:
  mode: blind
  exhaustiveness: 128
  n_poses: 20

postprocess:
  contact_cutoff: 4.0
  pocket_cutoff: 4.0
```

PyRosetta 도킹은 별도 `.ini` 설정 파일 사용 (legacy/ 참조).

---

## 출력 데이터 계약

### Vina 출력

| 파일 | 내용 |
|------|------|
| `vina_pose_table.csv` | 포즈별 affinity, centroid, contact_residues, pocket_id |
| `vina_pocket_table.csv` | 포켓별 통계 (best/mean affinity, union residues, top residues) |
| `vina_drug_pocket_map.csv` | 리간드→포켓 매핑 (dominant pocket, multimodal binding 판정) |
| `vina_pocket_comparison.csv` | 교차 수용체 비교 (Jaccard, overlap, shared residues, same_patch_candidate) |

### PPI 출력

| 파일 | 내용 |
|------|------|
| `ppi_pyrosetta_residues.csv` | 인터페이스 잔기 (정규화 ID, frequency, occupancy, ΔE) |
| `ppi_pyrosetta_summary.csv` | 수용체별 요약 (n_models, n_clusters, top_residues) |
| `ppi_afm_residues.csv` | AlphaFold-Multimer 잔기 (min_ca_distance) |

### 통합 출력

| 파일 | 내용 |
|------|------|
| `project_report.txt` | 종합 텍스트 보고서 |
| `combined_residue_evidence.csv` | Vina + PPI 증거 통합 (evidence_sources 필드) |

---

## 수용체 정보

| ID | 출처 | Chain | 잔기 범위 | 설명 |
|----|------|-------|-----------|------|
| 3GT8_raw | Crystal | A | 699–1007 | EGFR kinase domain 원본 |
| 3GT8_cl38_48 | MD cluster | X | 634–1014 | 38–48ns 대표 구조 |
| 3GT8_cl85_100 | MD cluster | X | 634–1014 | 85–100ns 대표 구조 |

> 잔기 번호 699–1007 구간은 세 수용체 모두 겹치며, chain ID만 다름 (A vs X).
> CHARMM 잔기명(HSD/HSE/HSP)은 파이프라인 내에서 자동 정규화됨.

---

## 환경 요구사항

```bash
# 최소 (Vina 후처리)
pip install pyyaml numpy pandas

# Vina 도킹
pip install vina rdkit matplotlib

# PyRosetta PPI 도킹
# PyRosetta 라이선스 필요 (conda install pyrosetta)

# MD 분석
pip install MDAnalysis
```

- **실행 환경**: Linux (Ubuntu HPC 권장), 16+ CPU cores
- **Python**: 3.9+

---

## 문서

| 문서 | 위치 |
|------|------|
| PyRosetta 매뉴얼 | `docs/manual_pyrosetta.md` |
| Vina 매뉴얼 | `docs/manual_vina.md` |
| PyRosetta 필터링 가이드 | `docs/pyrosetta_docking_fitering_guide.md` |
| PRD | `docs/prd_egfr_myo_1_d_pipeline.md` |
| Task breakdown | `docs/tasks_egfr_myo_1_d_pipeline.md` |
| Runbook | `docs/runbook.md` |
| 프로젝트 컨텍스트 | `docs/project_context.md` |
| Codex 인수인계 | `docs/codex_handoff_egfr_myo_1_d_pipeline_spec.md` |
| 세션 연속 노트 | `docs/CODEX_CONTINUATION_2026-03-09.md` |

---

## 스코어링 메트릭 참고

### Vina

| Metric | 의미 | 좋은 값 |
|--------|------|---------|
| affinity (kcal/mol) | 결합 친화도 | 낮을수록 강한 결합 |

### PyRosetta PPI

| Metric | 의미 | 좋은 값 |
|--------|------|---------|
| dG_separated (REU) | 인터페이스 결합 에너지 | < -10 |
| dSASA (Å²) | 매몰 표면적 | > 800 |
| sc_value | Shape Complementarity | > 0.65 |
| packstat | 원자 패킹 밀도 | > 0.65 |
| dG_density | 에너지 밀도 (dG/dSASA×100) | < -1.5 |
| delta_unsatHbonds | 미충족 수소결합 | < 5 |
| nres_int | 인터페이스 잔기 수 | > 15 |
| hbonds_int | 인터페이스 수소결합 | ≥ 1 |
