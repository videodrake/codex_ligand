# 실행 매뉴얼 — EGFR-MYO1D Pipeline

> 최종 업데이트: 2026-03-10

---

## 0. 환경 준비

conda 환경 `pyrosetta`를 사용한다:

```bash
cd ~/codex_ligand
conda activate pyrosetta
```

### 필수 패키지 설치 (처음 1회)

```bash
pip install pyyaml numpy pandas scipy matplotlib pytest
```

### 선택 패키지 (필요 시)

```bash
# Vina 도킹 (SDF→PDBQT 변환 포함)
pip install vina rdkit

# MD 궤적 분석
pip install MDAnalysis
```

> PyRosetta는 `pyrosetta` 환경에 이미 설치되어 있다.
> AutoDock Vina 바이너리는 별도 설치 필요 (https://vina.scripps.edu/).

---

## 1. Config 파일 작성

`config/example-project.yaml`을 복사해서 실제 경로로 수정:

```yaml
project_name: egfr_myo1d_vina
output_root: ./output
max_workers: 16

receptors:
  - id: 3GT8_raw
    pdb: input/receptors/3GT8_raw.pdb          # CA 좌표 (verdict centroid용)
    pdbqt: input/receptors/3GT8_raw_receptor.pdbqt  # Vina 도킹용
    chain: A
  - id: 3GT8_cl38_48
    pdb: input/receptors/3GT8_cl38_48.pdb
    pdbqt: input/receptors/3GT8_cl38_48_receptor.pdbqt
    chain: A
  - id: 3GT8_cl85_100
    pdb: input/receptors/3GT8_cl85_100.pdb
    pdbqt: input/receptors/3GT8_cl85_100_receptor.pdbqt
    chain: A

ligands:
  - id: 173940
    pdbqt: input/ligands/173940_ligand.pdbqt
  - id: 97806
    pdbqt: input/ligands/97806_ligand.pdbqt
  - id: VAX-C12_0
    pdbqt: input/ligands/VAX-C12_0_ligand.pdbqt

vina:
  mode: blind                # blind=전체 표면, focused=특정 좌표
  exhaustiveness: 128
  n_poses: 20

postprocess:
  contact_cutoff: 4.0       # 접촉 잔기 거리 컷오프
  pocket_cutoff: 4.0        # 포켓 클러스터링 centroid 컷오프
  merge_by_residue: false
  keep_chain: false

bootstrap:
  n_replicates: 100
  sample_fraction: 0.8
  seed: 42

experimental:                # 실험 데이터 있으면 주석 해제
  # known_binding_residues: [744, 752, 831]
  # known_non_binding_residues: [899, 950]
  # source: "Kim et al. 2023"

ppi:
  pyrosetta_result_dirs:     # PyRosetta 도킹 결과 경로
    # 3GT8_raw: /path/to/results/3GT8_raw
  afm_models:                # AlphaFold-Multimer 모델 PDB
    # 3GT8_raw: /path/to/afm_model_3GT8_raw.pdb
  afm_settings:
    receptor_chain: A
    partner_chain: B
    contact_cutoff: 8.0
```

---

## 2. 실행 방법 (2가지)

### 방법 A: 인터랙티브 메뉴

```bash
python main.py
```

메뉴가 뜨면 번호 선택:

| 번호 | 기능 | 설명 |
|------|------|------|
| 1 | Vina Docking | AutoDock Vina 실행 (Vina 바이너리 필요) |
| 2 | Vina Postprocess | 결과 파싱/클러스터링 (단계별 선택 가능) |
| 3 | PPI Docking | PyRosetta 도킹 (HPC 서버에서) |
| 7 | Site Verdict | 증거 기반 사이트 판정 |
| 5 | Generate Report | 종합 보고서 |
| 6 | Validate Outputs | 출력 검증 |
| 9 | Full Pipeline | 1→2→7→5→6 자동 실행 |

### 방법 B: CLI 서브커맨드 (비대화식)

```bash
python main.py vina -c config/my-project.yaml
python main.py postprocess -c config/my-project.yaml
python main.py verdict -c config/my-project.yaml
python main.py report -c config/my-project.yaml
python main.py validate -c config/my-project.yaml
python main.py full -c config/my-project.yaml    # 전체 자동 실행
```

---

## 3. 실행 순서 (상세)

### Step 1: Vina 리간드 도킹

```bash
python main.py vina -c config/my-project.yaml
```

- 각 수용체 × 리간드 조합에 대해 blind docking 실행
- **필요**: Vina 바이너리 + `.pdbqt` 파일
- **출력**: `output/{project_name}/{receptor_id}/{ligand_id}/` 하위에 `.pdbqt` 포즈 파일

### Step 2: 후처리 (Postprocess)

```bash
python main.py postprocess -c config/my-project.yaml
```

인터랙티브에서 `a` (전체 실행) 선택하면 순서대로:

| 단계 | 하는 일 | 출력 |
|------|---------|------|
| parse | Vina 결과 → CSV | `vina_pose_table.csv` |
| contacts | 수용체-리간드 접촉 잔기 추출 | pose_table에 contact_residues 추가 |
| cluster | 포즈를 포켓으로 클러스터링 | pose_table에 pocket_id 추가 |
| summarize | 포켓별 통계 | `vina_pocket_table.csv`, `vina_drug_pocket_map.csv` |
| compare | 교차 수용체 포켓 비교 | `vina_pocket_comparison.csv` |
| ppi | PyRosetta 결과에서 잔기 추출 | `ppi_pyrosetta_residues.csv` |
| bootstrap | 포켓 안정성 분석 | `vina_pocket_bootstrap.csv` |

### Step 3: PPI 도킹 (별도 — HPC 서버)

**3a. Dimer PDB 준비** (로컬):

```bash
python main.py   # 메뉴 [3] → [a]
# EGFR dimer PDB + MYO1D 도메인 PDB 합치기
```

**3b. 도킹 실행** (HPC):

```bash
# 테스트 (1K 모델, 3-4시간)
qsub config/run_ppi_test.pbs

# 프로덕션 (20K 모델, 24-36시간)
qsub config/run_ppi_prod.pbs
```

**3c. 결과 원복** (도킹 후):

```bash
python main.py   # 메뉴 [8] PPI Postprocess
# chain 번호 정상화 + 인터페이스 잔기 추출
```

### Step 4: AlphaFold-Multimer (별도 — 로컬)

AFM 예측을 외부에서 실행한 후, 결과 PDB 경로를 config에 기재:

```yaml
ppi:
  afm_models:
    3GT8_raw: /path/to/afm_model.pdb
```

AFM 잔기 추출은 verdict 실행 시 자동으로 로딩된다.
직접 추출도 가능:

```bash
python -c "
from egfr_pipeline.ppi.afm_extract import extract_afm_batch
extract_afm_batch('config/my-project.yaml')
"
```

### Step 5: 사이트 판정 (Verdict)

```bash
python main.py verdict -c config/my-project.yaml
```

- Vina 포켓 + PPI(PyRosetta+AFM) + 교차 수용체 증거를 결합
- 3축 점수 합산 (100점 만점) → STRONG / MODERATE / WEAK 분류
- **출력**: `cross_method_agreement.csv`, `valid_sites.csv`

### Step 6: 보고서 생성

```bash
python main.py report -c config/my-project.yaml
```

- **출력**: `project_report.txt` (텍스트), `combined_residue_evidence.csv`

### Step 7: 검증

```bash
python main.py validate -c config/my-project.yaml
```

- CSV 스키마, ID 일관성, 잔기 번호 범위 확인
- Exit code: 0=PASS, 1=WARN, 2=FAIL

---

## 4. 한 줄 요약 (전체 자동)

Vina 도킹이 이미 완료되어 결과 파일이 있다면:

```bash
python main.py full -c config/my-project.yaml
```

→ postprocess → verdict → report → validate 자동 실행

---

## 5. 테스트

```bash
# 전체 테스트 (54개, ~1초)
pytest tests/ -v

# 특정 테스트만
pytest tests/ -k afm          # AFM 통합 테스트
pytest tests/ -k bootstrap     # Bootstrap 테스트
pytest tests/ -k e2e           # End-to-end 테스트
pytest tests/ -k schema        # 스키마 일치
```

---

## 6. 출력 파일 위치

모든 출력은 `{output_root}/{project_name}/` 아래:

```
output/egfr_myo1d_vina/
├── vina_pose_table.csv              ← 포즈별 상세 (Step 2)
├── vina_pocket_table.csv            ← 포켓별 통계 (Step 2)
├── vina_drug_pocket_map.csv         ← 리간드→포켓 매핑 (Step 2)
├── vina_pocket_comparison.csv       ← 교차 수용체 비교 (Step 2)
├── vina_pocket_bootstrap.csv        ← 포켓 안정성 (Step 2)
├── ppi_pyrosetta_residues.csv       ← PyRosetta PPI 잔기 (Step 3)
├── ppi_afm_residues.csv             ← AFM PPI 잔기 (Step 4)
├── cross_method_agreement.csv       ← Vina↔PPI 공간 비교 (Step 5)
├── valid_sites.csv                  ← 사이트 판정 결과 (Step 5)
├── vina_consensus_sites.csv         ← 교차 수용체 합의 사이트 (Step 5)
├── project_report.txt               ← 종합 보고서 (Step 6)
└── combined_residue_evidence.csv    ← Vina+PPI+AFM 통합 증거 (Step 6)
```

---

## 7. 핵심 포인트

- **Vina 도킹 자체**만 외부 바이너리(AutoDock Vina)가 필요. 나머지는 전부 Python만으로 실행
- **PPI 도킹**은 HPC에서 PyRosetta로 별도 실행 → 결과 경로를 config에 기재
- **AFM**도 별도 실행 → 모델 PDB 경로를 config에 기재
- 이미 Vina 결과가 있으면 `postprocess`부터 시작 가능 (Step 1 스킵)
- PPI/AFM 데이터 없어도 Vina + Cross-receptor만으로 verdict 가능 (adaptive scoring)

---

## 8. Verdict 점수 체계

3축 합산 100점:

| 축 | PPI 있을 때 | PPI 없을 때 | 기준 |
|----|------------|------------|------|
| Vina Quality | 50점 | 60점 | affinity, pose 수렴, 리간드 합의 |
| PPI Proximity | 20점 | 0점 | pocket↔PPI 3D 거리 |
| Cross-receptor | 30점 | 40점 | 다른 수용체에서 같은 포켓 발견 |

판정:
- **STRONG** (≥55): 우선 검토 대상
- **MODERATE** (≥30): 참고 대상
- **WEAK** (<30): 낮은 증거

> 이것은 유효성 판정이 아니라 **증거 강도 분류**이다. STRONG 포켓도 연구자의 시각적 확인이 필요하다.

---

## 9. 수용체 정보

| ID | 출처 | Chain | 잔기 범위 | 설명 |
|----|------|-------|-----------|------|
| 3GT8_raw | Crystal | A | 699–1007 | EGFR kinase domain 원본 |
| 3GT8_cl38_48 | MD cluster | X | 634–1014 | 38–48ns 대표 구조 |
| 3GT8_cl85_100 | MD cluster | X | 634–1014 | 85–100ns 대표 구조 |

- 잔기 번호 699–1007 구간은 세 수용체 모두 겹침
- CHARMM 잔기명(HSD/HSE/HSP)은 파이프라인 내에서 자동 정규화

---

## 10. 문제 해결

| 증상 | 원인 | 해결 |
|------|------|------|
| `RuntimeError: pyyaml is required` | pyyaml 미설치 | `pip install pyyaml` |
| `vina_pose_table.csv` 없음 | Vina 도킹 미실행 또는 parse 미실행 | Step 1, 2 순서 확인 |
| PPI 잔기 번호 >1700 | chain B offset 미원복 | 메뉴 [8] PPI Postprocess 실행 |
| verdict에서 PPI 점수 0 | PPI 데이터 경로 미설정 | config의 `ppi` 섹션 확인 |
| bootstrap 결과 없음 | postprocess에서 bootstrap 단계 스킵 | postprocess → bootstrap 실행 |
| validate FAIL | CSV 스키마 변경 | `tests/test_pipeline.py` SchemaConsistency 테스트 확인 |
