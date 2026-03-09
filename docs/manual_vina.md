# AutoDock Vina Docking Pipeline - Manual

## 1. 개요

이 파이프라인은 AutoDock Vina Python API를 사용하여 단백질-리간드 도킹을 수행합니다.
하나의 통합 스크립트(`run_docking.py`)로 **Blind Docking**과 **Focused Docking**을 모두 지원합니다.

**두 가지 실행 방식:**
- `python run_docking.py` — 인터랙티브 모드 (단계별 메뉴 선택)
- `python run_docking.py --mode blind ...` — CLI 모드 (명령줄 인자)

---

## 2. 폴더 구조

```
vina_docking/
  run_docking.py              # 통합 도킹 스크립트
  MANUAL.md                   # 이 매뉴얼
  input/                      # 입력 파일 (receptor + ligand)
    3gt8_monomer_anp.pdb      # receptor PDB
    3gt8_monomer_anp_receptor.pdbqt  # 준비된 receptor PDBQT
    97806_ligand.pdbqt         # ligand PDBQT
    97806_ligand.sdf           # ligand 원본 SDF
    ...
  output/                     # 출력 결과
    YYYY-MM-DD_<label>/       # 실행별 폴더 (날짜_receptor이름)
      <ligand>_blind.pdbqt    # 도킹 결과 (모드명 포함)
      <ligand>_blind_pockets.pdbqt       # 완성 포켓 best pose
      <ligand>_blind_all_pockets.pdbqt   # 완성 포켓 전체 포즈
      <ligand>_blind_partial_pockets.pdbqt  # 미완성 포켓 (참고)
      config.yaml             # 실행 설정 (재현용)
    selected_poses/           # 수동 선별한 최종 pose
```

### 자동 파일 인식

**어떤 이름이든 상관없이** `input/`에 넣으면 자동으로 분류됩니다:

| 분류 기준 | 판정 | 예시 |
|-----------|------|------|
| 파일명에 `_receptor` 포함 | Receptor | `3gt8_receptor.pdbqt` |
| 파일명에 `_ligand` 포함 | Ligand | `97806_ligand.sdf` |
| `.sdf`, `.mol2` 확장자 | Ligand (소분자) | `drug.sdf`, `compound.mol2` |
| `.smi` 확장자 | SMILES | `compounds.smi` |
| `.pdb`/`.cif` + 원자 500개 초과 | Receptor (단백질) | `protein.pdb` |
| `.pdb` + 원자 500개 이하 | Ligand (소분자) | `small_mol.pdb` |

**네이밍 규칙 없이도 작동합니다.** 그냥 파일을 넣으면 됩니다.

기존 네이밍 규칙(`_receptor`, `_ligand`)도 여전히 우선 인식됩니다.

---

## 3. 필수 환경

### Python 패키지
```bash
pip install numpy vina pyyaml
```

| 패키지 | 필수 여부 | 용도 |
|--------|----------|------|
| `numpy` | 필수 | 좌표 계산 |
| `vina` | 필수 | AutoDock Vina Python API |
| `pyyaml` | 권장 | config 저장/로드 (미설치 시 JSON fallback) |

### Receptor 변환 도구 (--prepare-receptor 사용 시)
아래 중 하나 이상 설치 필요 (순서대로 시도):
1. **ADFR Suite** — `prepare_receptor` 명령
2. **MGLTools** — `prepare_receptor4.py` 스크립트
3. **OpenBabel** — `obabel` 명령

---

## 4. 사용법

### 4.0 인터랙티브 모드 (권장)

인자 없이 실행하면 단계별 메뉴가 표시됩니다:

```bash
python run_docking.py
```

```
============================================================
  AutoDock Vina Docking Pipeline
  인터랙티브 모드
============================================================

--- 모드 선택 ---
  [1] Blind Docking    (단백질 전체 탐색)
  [2] Focused Docking  (특정 포켓 집중)
  [3] Config 재실행    (이전 설정 로드)
  [q] 종료

선택: _
```

이후 Receptor, Ligand, 영역, 파라미터를 순서대로 선택하고 최종 확인 후 실행됩니다.
- Receptor/Ligand가 1개뿐이면 자동 선택
- Ligand는 복수 선택 가능 (예: `1,2` 또는 `all`)
- Config 재실행: output/ 내 이전 config.yaml을 선택하여 동일 조건 재실행

### 4.1 CLI 모드 — 기본 사용 (자동 감지)

`input/` 폴더에 파일을 넣으면 자동으로 감지합니다.

```bash
# Blind docking — 단백질 전체 표면 탐색
python run_docking.py --mode blind

# Focused docking — 등록된 포켓 프리셋 사용
python run_docking.py --mode focused --region clobe
```

자동 감지 규칙:
- `*_receptor.pdbqt`가 1개면 자동 선택, 복수면 에러 (--receptor로 지정 필요)
- `*_ligand.pdbqt`는 전부 선택
- receptor PDB는 PDBQT 이름에서 유추

### 4.2 특정 파일 지정

```bash
# receptor 수동 지정
python run_docking.py --mode blind \
  --receptor input/3gt8_monomer_anp_receptor.pdbqt \
  --receptor-pdb input/3gt8_monomer_anp.pdb

# 특정 리간드만 도킹
python run_docking.py --mode blind --ligands input/97806_ligand.pdbqt

# 복수 리간드 지정
python run_docking.py --mode blind \
  --ligands input/97806_ligand.pdbqt input/173940_ligand.pdbqt
```

### 4.3 Focused Docking — 좌표 직접 지정

```bash
python run_docking.py --mode focused \
  --center -38.2 0.2 -73.8 \
  --box-size 30 30 30
```

### 4.4 파라미터 조정

```bash
# exhaustiveness 변경 (높을수록 정확, 느림)
python run_docking.py --mode blind --exhaustiveness 256

# pose 수 변경
python run_docking.py --mode blind --n-poses 50

# blind mode box 패딩/최소 크기 변경
python run_docking.py --mode blind --padding 10.0 --min-box 80.0
```

### 4.5 Config 파일로 재현

매 실행마다 `config.yaml`이 자동 저장됩니다. 이를 재사용하면 동일 조건으로 재실행 가능:

```bash
python run_docking.py --config output/2026-02-26_3gt8_monomer_anp/config.yaml
```

config.yaml 예시:
```yaml
mode: focused
receptor: input/3gt8_monomer_anp_receptor.pdbqt
receptor_pdb: input/3gt8_monomer_anp.pdb
ligands:
  - input/97806_ligand.pdbqt
  - input/173940_ligand.pdbqt
region: clobe
center: [-38.2, 0.2, -73.8]
box_size: [30, 30, 30]
exhaustiveness: 128
n_poses: 20
```

### 4.6 Receptor PDB → PDBQT 변환

새 receptor PDB를 준비할 때:

```bash
python run_docking.py --mode blind \
  --receptor-pdb input/3gt8_dimer_anp.pdb \
  --prepare-receptor
```

변환된 PDBQT는 `input/` 폴더에 저장됩니다.

### 4.7 출력 디렉토리 직접 지정

```bash
python run_docking.py --mode blind --output-dir output/my_custom_run

# 또는 라벨만 지정 (날짜는 자동)
python run_docking.py --mode blind --label test_run
# → output/2026-03-05_test_run/
```

---

## 5. 전체 인자 목록

| 인자 | 기본값 | 설명 |
|------|--------|------|
| `--mode` | (필수) | `blind` 또는 `focused` |
| `--receptor` | auto | receptor PDBQT 경로 |
| `--receptor-pdb` | auto | receptor PDB 경로 (blind 모드 box 계산용) |
| `--ligands` | auto | ligand PDBQT 파일 (복수 가능) |
| `--smiles` | - | SMILES 문자열 직접 입력 (자동 3D→PDBQT 변환) |
| `--center X Y Z` | - | focused 모드 box 중심 좌표 |
| `--box-size X Y Z` | 30 30 30 | box 크기 (Angstrom) |
| `--region` | - | 이름 기반 포켓 프리셋 |
| `--exhaustiveness` | 128 | 탐색 강도 (높을수록 정확, 느림) |
| `--n-poses` | 20 | 생성할 pose 수 |
| `--padding` | 5.0 | blind 모드 box 여유 (Angstrom) |
| `--min-box` | 70.0 | blind 모드 최소 box 크기 (Angstrom) |
| `--output-dir` | auto | 출력 디렉토리 직접 지정 |
| `--label` | auto | 출력 디렉토리 라벨 |
| `--config` | - | YAML/JSON config 파일 경로 |
| `--prepare-receptor` | False | PDB → PDBQT 변환 수행 |
| `--n-pockets N` | - | 다중 포켓 탐색: N개 포켓 발견 |
| `--max-per-pocket M` | 3 | 포켓당 최대 포즈 수 |
| `--cluster-radius R` | 5.0 | 같은 포켓 판정 반경 (Å) |
| `--exclude-zone X Y Z R` | - | 해당 좌표 반경 R Å 이내 포즈 제외 (복수 가능) |

### 4.8 다중 포켓 탐색 (--n-pockets)

블라인드 도킹에서 여러 결합 부위를 자동으로 발견합니다.

```bash
# 5개 포켓 탐색 (기본: 포켓당 3개 포즈)
python run_docking.py --mode blind --n-pockets 5

# 포켓당 5개 포즈, 판정 반경 8Å
python run_docking.py --mode blind --n-pockets 5 --max-per-pocket 5 --cluster-radius 8
```

**작동 원리:**
포즈를 에너지 순서대로 하나씩 처리하면서 포켓에 배정합니다:
1. 포즈 수 자동 증가 (충분한 다양성 확보)
2. 각 포즈의 centroid → 가장 가까운 열린 포켓에 배정
3. 포켓이 `max-per-pocket`개 차면 **닫힘** → 다음 포즈는 다른 포켓으로
4. 새 위치의 포즈는 **새 포켓** 생성
5. N개 포켓이 모두 채워지면 종료

**예시 흐름 (5개 포켓, 포켓당 3개):**
```
pose 1 → Pocket #1 (1/3)
pose 2 → Pocket #1 (2/3)
pose 3 → Pocket #2 (1/3)
pose 4 → Pocket #1 (3/3) → Pocket #1 완성!
pose 5 → Pocket #2 (2/3)
pose 6 → Pocket #3 (1/3)
pose 7 → Pocket #3 (2/3)
pose 8 → Pocket #2 (3/3) → Pocket #2 완성!
pose 9 → [Pocket #1 닫힘, skip] → Pocket #4 (1/3)
...
→ 5개 포켓 모두 완성!
```

**파라미터:**
| 인자 | 기본값 | 설명 |
|------|--------|------|
| `--n-pockets N` | - | 탐색할 포켓 수 |
| `--max-per-pocket M` | 3 | 포켓당 포즈 수 (가득 차면 닫힘) |
| `--cluster-radius R` | 5.0 | 같은 포켓 판정 반경 (Å) |

**출력 파일:**
- `<ligand>_blind.pdbqt` — 전체 포즈 (기존)
- `<ligand>_blind_pockets.pdbqt` — 각 포켓의 best pose만 (포켓 수만큼)
- `<ligand>_blind_all_pockets.pdbqt` — 모든 포켓의 전체 포즈

### 4.9 Exclusion Zone (수동 좌표 제외)

좌표를 직접 지정하여 특정 영역의 포즈를 제외할 수도 있습니다:

```bash
python run_docking.py --mode blind --exclude-zone -38.2 0.2 -73.8 10
```

인터랙티브 모드에서는 "포즈 필터링" 메뉴에서 자동/수동을 선택할 수 있습니다.

---

## 6. 통합 전처리 시스템

스크립트 실행 시 `input/` 폴더의 **모든** 분자 파일을 자동으로 스캔 → 분류 → 변환합니다.
**이미 PDBQT가 존재하면 자동으로 스킵합니다.**

실행 시 전처리 요약이 자동 출력됩니다:
```
══════════════════════════════════════════════════
  전처리 완료 — input/ 현황
══════════════════════════════════════════════════
  Receptor: 1개
    > 3gt8_monomer_anp_receptor.pdbqt
  Ligand:   2개
    > 97806_ligand.pdbqt
    > drug_ligand.pdbqt
  변환: receptor 0개, ligand 1개
══════════════════════════════════════════════════
```

### Receptor (PDB/MOL2/CIF → PDBQT)
```
input/3gt8_monomer_anp.pdb  →  input/3gt8_monomer_anp_receptor.pdbqt
input/protein.mol2          →  input/protein_receptor.pdbqt
```
변환 도구 시도 순서: ADFR → MGLTools → OpenBabel

### Ligand (SDF/MOL2 → PDBQT)
```
input/97806_ligand.sdf      →  input/97806_ligand.pdbqt
input/173940_ligand.mol2    →  input/173940_ligand.pdbqt
```
변환 도구 시도 순서: Meeko → ADFR → MGLTools → OpenBabel

### SMILES 입력 (.smi 파일)

SMILES 문자열이 담긴 `.smi` 파일을 `input/`에 넣으면 자동으로 3D 구조 생성 → PDBQT 변환됩니다.

**SMI 파일 포맷** (`input/compounds.smi`):
```
# 주석 줄 (무시됨)
CC(=O)Oc1ccccc1C(=O)O  aspirin
c1ccc2[nH]c(-c3ccccn3)nc2c1  benzimidazole
CCO
```
- 각 줄: `SMILES [이름]` (탭 또는 공백 구분)
- 이름을 생략하면 `{파일명}_0`, `{파일명}_1`, ... 자동 부여
- `#`으로 시작하는 줄과 빈 줄은 무시

**변환 흐름:**
```
compounds.smi 의 "CC(=O)Oc1ccccc1C(=O)O  aspirin"
    → input/aspirin_ligand.sdf    (3D 좌표 생성)
    → input/aspirin_ligand.pdbqt  (도킹용 변환)
```

**CLI에서 직접 SMILES 입력:**
```bash
# 단일 SMILES
python run_docking.py --mode blind --smiles "CC(=O)Oc1ccccc1C(=O)O"

# 복수 SMILES
python run_docking.py --mode blind --smiles "CCO" "CC(=O)Oc1ccccc1C(=O)O"
```

**3D 구조 생성 도구** (하나 이상 필요):
| 도구 | 설치 | 비고 |
|------|------|------|
| RDKit | `conda install -c conda-forge rdkit` | ETKDGv3 + MMFF 최적화 (권장) |
| OpenBabel | `conda install -c conda-forge openbabel` | `--gen3d --minimize` |

### 사용 흐름
1. `input/`에 receptor PDB와 ligand 파일을 넣는다 (SDF/MOL2/SMI 중 아무거나)
2. `python run_docking.py` 실행
3. 자동으로 SMILES→SDF→PDBQT 변환 → 도킹 실행

변환 도구 설치:
```bash
pip install meeko                       # 리간드 SDF→PDBQT 변환 (권장)
conda install -c conda-forge openbabel  # 범용 변환 + SMILES 3D 생성
conda install -c conda-forge rdkit      # SMILES 3D 생성 (선택)
```

---

## 7. Blind vs Focused Docking

### 6.1 알고리즘

두 모드 모두 동일한 AutoDock Vina 엔진을 사용합니다.
차이점은 **탐색 공간(search box)의 크기와 위치**뿐입니다.

**Vina 탐색 알고리즘: Iterated Local Search (ILS)**
1. 리간드의 위치/방향/torsion을 랜덤하게 변형 (mutation)
2. 변형 후 BFGS quasi-Newton 방법으로 local optimization
3. Metropolis criterion으로 수락/거부 결정
4. exhaustiveness 횟수만큼 독립적으로 병렬 실행
5. 모든 결과에서 에너지 순 정렬 + RMSD 클러스터링 → 상위 pose 출력

**Vina Scoring Function 구성 요소:**
- van der Waals 상호작용
- 수소결합
- 소수성 효과
- 회전 가능 결합 수에 따른 엔트로피 페널티

### 6.2 비교

| | Blind | Focused |
|---|---|---|
| Box 크기 | 70+ Angstrom (단백질 전체) | 20~40 Angstrom (포켓) |
| Box 중심 | PDB 좌표에서 자동 계산 | 사용자 지정 or 프리셋 |
| 탐색 범위 | 전체 표면 | 특정 결합 부위 |
| 장점 | 편향 없는 탐색, 새 포켓 발견 | 높은 sampling density, 정확도 |
| 단점 | 동일 exhaustiveness에서 정확도 낮음 | 사전지식 필요 |
| 용도 | 결합 부위 탐색, 초기 스크리닝 | 특정 포켓 결합 최적화 |

### 6.3 권장 워크플로우

```
1. Blind docking     → 전체 표면에서 결합 부위 후보 발견
2. 결과 분석          → PyMOL/Chimera에서 상위 pose 확인
3. Focused docking   → 유의미한 포켓에 대해 집중 도킹
4. 결과 비교          → affinity, pose 안정성 평가
```

---

## 7. Focused Docking 사전지식 입력 방법

Focused docking의 핵심은 **어디에 box를 놓을 것인가**입니다.
box 중심 좌표를 결정하는 방법:

### 7.1 Region 프리셋 (스크립트 내장)

`run_docking.py` 내의 `REGION_PRESETS` 딕셔너리에 등록된 좌표:

```python
REGION_PRESETS = {
    "clobe": {
        "center": [-38.2, 0.2, -73.8],
        "box_size": [30, 30, 30],
        "description": "C-lobe pocket (catalytic loop 808-810 + A-loop PRO848 + AP-2 helix 987-991)"
    },
}
```

새 포켓을 발견하면 여기에 추가:
```python
    "dimer_interface": {
        "center": [x, y, z],
        "box_size": [25, 25, 25],
        "description": "Dimer interface pocket"
    },
```

사용: `python run_docking.py --mode focused --region clobe`

### 7.2 구조 기반 수동 지정

PyMOL이나 UCSF Chimera에서 관심 잔기의 좌표를 확인하고 직접 입력:

**PyMOL에서 좌표 확인:**
```
# 관심 잔기 선택
select pocket, resi 808-810+848+987-991
# 중심 좌표 출력
centerofmass pocket
```

**Chimera에서 좌표 확인:**
```
select :808-810,848,987-991
measure center sel
```

사용: `python run_docking.py --mode focused --center -38.2 0.2 -73.8 --box-size 30 30 30`

### 7.3 Blind Docking 결과에서 추출

1. Blind docking 실행
2. 결과 PDBQT를 PyMOL에서 로드
3. 유의미한 pose의 리간드 중심 좌표 확인
4. 그 좌표를 focused docking의 center로 사용

### 7.4 포켓 예측 도구 활용

외부 도구로 결합 포켓을 자동 예측한 후 좌표를 입력:

| 도구 | 방식 | 설치 |
|------|------|------|
| **fpocket** | Geometry 기반 (Voronoi tessellation) | `conda install -c conda-forge fpocket` |
| **P2Rank** | 머신러닝 기반 | https://github.com/rdk/p2rank |
| **DoGSiteScorer** | Grid 기반 | https://proteins.plus/ (웹) |

fpocket 예시:
```bash
fpocket -f input/3gt8_monomer_anp.pdb
# 결과: 3gt8_monomer_anp_out/ 폴더에 포켓 정보 출력
# 각 포켓의 중심 좌표를 --center로 사용
```

---

## 9. ANP (ATP 포켓 차단)

receptor PDB에 ANP(adenylyl-imidodiphosphate)가 포함되어 있으면,
ATP 결합 포켓이 ANP로 점유되어 리간드가 해당 포켓에 결합하지 못합니다.

- 파일명에 "anp"가 포함된 receptor PDB를 사용하면 자동으로 ANP 원자 수를 확인
- ANP가 없으면 경고 출력 (ATP 포켓이 열려있어 리간드가 잘못 결합할 수 있음)

---

## 10. 출력 결과 해석

### 9.1 결과 테이블

```
============================================================
  97806 Blind Docking Results
  Center: [-28.5, 5.3, -62.1]
  Box: [70.0, 70.0, 85.2]
============================================================
 Mode   Affinity  RMSD_lb  RMSD_ub
    1      -8.52    0.000    0.000
    2      -8.31    2.145    4.892
    3      -7.98    1.853    3.201
    ...
```

| 컬럼 | 의미 |
|------|------|
| Mode | pose 번호 (에너지 순 정렬) |
| Affinity | 결합 에너지 (kcal/mol, 음수일수록 강한 결합) |
| RMSD_lb | 최적 pose 대비 RMSD 하한 (lower bound) |
| RMSD_ub | 최적 pose 대비 RMSD 상한 (upper bound) |

### 9.2 해석 가이드

| Affinity 범위 | 해석 |
|---------------|------|
| -10 이하 | 매우 강한 결합 (드묾) |
| -8 ~ -10 | 강한 결합 |
| -6 ~ -8 | 중간 결합 |
| -4 ~ -6 | 약한 결합 |
| -4 이상 | 유의미하지 않음 |

### 9.3 결과 파일

- `<ligand>_<mode>.pdbqt` — 모든 pose가 포함된 PDBQT 파일
- `config.yaml` — 실행 파라미터 (재현용)

### PyMOL 시각화 스크립트 자동 생성

도킹 완료 후 `output/.../view_results.pml` 자동 생성.
PyMOL에서 실행하면 receptor 색칠 + 도킹 결과를 한눈에 볼 수 있음:

```bash
# PyMOL 실행 (스크립트 로드)
pymol output/2026-03-06_3gt8/view_results.pml

# 또는 PyMOL 내에서
@output/2026-03-06_3gt8/view_results.pml
```

자동 포함되는 내용:
- Receptor PDB 로드 (cartoon, gray 기본)
- 사용자 색상 스크립트 적용 (있으면)
- 도킹 결과 PDBQT 로드 (리간드별 다른 색상)
- 포켓/필터 결과 자동 포함 (있으면)

### 사용자 색상 스크립트 (.pml)

`input/`에 `.pml` 파일을 넣으면 receptor 색칠에 자동 적용됨.
`receptor`라는 객체 이름을 기준으로 작성:

```pymol
# 예시: EGFR 색상 분석 (input/3gt8.pml)
color gray, receptor
color red, receptor and resi 705-800        # N-lobe
color orangered, receptor and resi 728-733  # P-loop
color limegreen, receptor and resi 856-979  # C-lobe
color yellow, receptor and resi 997-1002    # AP2 helix
```

매칭 규칙: `3gt8.pml` → `3gt8_receptor.pdbqt` (파일명 stem 매칭).
`.pml`이 1개면 모든 receptor에 적용.

---

## 11. 새 리간드/receptor 추가하기

### 새 리간드 추가

1. SDF 또는 MOL2 파일을 PDBQT로 변환:
   ```bash
   obabel new_compound.sdf -O input/new_compound_ligand.pdbqt --gen3d
   ```
   또는 ADFR Suite:
   ```bash
   prepare_ligand -l new_compound.mol2 -o input/new_compound_ligand.pdbqt
   ```

2. `input/` 폴더에 `<ID>_ligand.pdbqt` 네이밍으로 저장

3. 도킹 실행 (자동 감지됨):
   ```bash
   python run_docking.py --mode blind
   ```

### 새 receptor 추가

1. PDB 파일을 `input/`에 저장
2. PDBQT 변환:
   ```bash
   python run_docking.py --mode blind \
     --receptor-pdb input/new_receptor.pdb \
     --prepare-receptor
   ```
3. 변환된 `new_receptor_receptor.pdbqt`가 `input/`에 생성됨

---

## 12. 문제 해결

| 증상 | 원인 | 해결 |
|------|------|------|
| "복수 receptor 감지됨" | input/에 *_receptor.pdbqt가 2개 이상 | `--receptor`로 하나 지정 |
| "Blind 모드에는 receptor PDB가 필요합니다" | PDB 파일을 자동 유추할 수 없음 | `--receptor-pdb`로 지정 |
| "Focused 모드에는 --center 지정 필수" | center 좌표 누락 | `--center X Y Z` 또는 `--region` 사용 |
| "pyyaml 미설치" | PyYAML 없음 | `pip install pyyaml` 또는 JSON fallback 사용 |
| "receptor PDBQT 변환 실패" | 변환 도구 미설치 | ADFR, MGLTools, 또는 OpenBabel 설치 |
| ANP WARNING | receptor에 ANP 미포함 | ANP 포함 PDB 사용 또는 의도된 경우 무시 |
| "SMILES→SDF 변환 실패" | RDKit/OpenBabel 미설치 | `conda install -c conda-forge openbabel` 또는 `rdkit` 설치 |
| "Invalid SMILES" | SMILES 문법 오류 | SMILES 문자열 확인 (https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html) |
