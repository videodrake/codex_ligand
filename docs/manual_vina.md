> Status note (2026-03-12): Read `docs/current_pipeline_status.md` first.
> Vina remains active, but current Phase 1 secondary validation uses LightDock, not AFM.

# AutoDock Vina Docking Pipeline - Manual

> **Deprecation Notice**: 이전 버전의 `run_docking.py` 단독 실행은 더 이상 지원되지 않습니다.
> 모든 기능은 통합 CLI인 `main.py`를 통해 접근합니다.
> 이 저장소에서는 통합 CLI(`main.py`)를 기준으로 사용하세요.
> 구버전 단독 스크립트는 이 저장소에 없으며, 필요하면 `../autodot_vina/`를 참고하세요.

## 1. 개요

이 파이프라인은 AutoDock Vina Python API를 사용하여 단백질-리간드 도킹을 수행합니다.
통합 CLI(`main.py`)를 통해 **Vina 도킹, 후처리, 보고서 생성, 검증, 사이트 판정**을 모두 관리합니다.

**두 가지 실행 방식:**
- `python main.py` -- 인터랙티브 모드 (메뉴 [1]-[9] 선택)
- `python main.py -c config/example-project.yaml vina` -- CLI 모드 (서브커맨드)

**사용 가능한 서브커맨드:**

| 서브커맨드 | 설명 | 메뉴 번호 |
|------------|------|-----------|
| `vina` | AutoDock Vina 도킹 실행 | [1] |
| `postprocess` | Vina 결과 파싱/클러스터링 | [2] |
| `pyrosetta` | PyRosetta PPI 도킹 | [3] |
| `md` | GROMACS MD 분석 | [4] |
| `report` | 종합 보고서 생성 | [5] |
| `validate` | 출력 검증 | [6] |
| `verdict` | 유효 사이트 자동 판정 | [7] |
| `ppi-postprocess` | PPI 후처리 자동화 | [8] |
| `full` | 전체 파이프라인 (1->2->7->5->6) | [9] |

---

## 2. 폴더 구조

```
codex_ligand/
  main.py                          # 통합 CLI 진입점
  config/
    example-project.yaml           # 프로젝트 설정 (YAML)
    phase1/                        # Phase 1 PyRosetta 설정 (18 INI files)
    run_lightdock.pbs              # Phase 1 LightDock PBS
    run_lightdock_test.pbs         # Phase 1 LightDock 테스트 PBS
  egfr_pipeline/
    config.py                      # Config 로드 유틸리티
    report.py                      # 종합 보고서 생성
    validate.py                    # 출력 검증
    verdict.py                     # 사이트 판정 (STRONG/MODERATE/WEAK)
    vina/
      vina_executor.py              # Vina 도킹 코어 (구 run_docking.py)
      parse_poses.py               # 결과 파싱 (구 parse_vina_results.py)
      pose_contacts.py             # 접촉 잔기 추출 (구 extract_contacts.py)
      pocket_cluster.py            # 포켓 클러스터링 (구 cluster_pockets.py)
      pocket_summary.py            # 포켓 요약 (구 summarize_pockets.py)
      cross_receptor.py            # 교차 비교 (구 compare_pockets.py)
    ppi/
      prepare_dimer_pdb.py         # Dimer PDB 준비 + chain 원복
      pyrosetta_extract.py         # PPI 잔기 추출 (구 extract_ppi_residues.py)
      postprocess_ppi.py           # PPI 후처리 자동화
      afm_extract.py               # AlphaFold-Multimer 파서 (stub)
    pyrosetta_docking/
      pipeline_manager.py          # PyRosetta PPI 도킹 오케스트레이터
      movers.py                    # PyRosetta 워커 (Relax, Docking, Refinement)
      scoring.py                   # 스코어링, RMSD, Interface 분석
      pyrosetta_init.py            # PyRosetta 유틸리티
  # 구버전 단독 Vina 스크립트는 이 워크스페이스에 포함되지 않음
  # 필요 시 ../autodot_vina/run_docking.py 참고

  input/                           # 입력 파일
  output/                          # 출력 결과
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
| `pyyaml` | 필수 | 프로젝트 config 로드 |
| `pandas` | 필수 | 후처리 CSV 작업 |
| `scipy` | 필수 | 클러스터링 (scipy.cluster.hierarchy) |

### Receptor 변환 도구 (--prepare-receptor 사용 시)
아래 중 하나 이상 설치 필요 (순서대로 시도):
1. **ADFR Suite** -- `prepare_receptor` 명령
2. **MGLTools** -- `prepare_receptor4.py` 스크립트
3. **OpenBabel** -- `obabel` 명령

---

## 4. 사용법

### 4.0 인터랙티브 모드 (권장)

인자 없이 실행하면 통합 메뉴가 표시됩니다:

```bash
python main.py
```

```
╔══════════════════════════════════════════════════════════╗
║     EGFR-MYO1D Docking Pipeline  (통합 실행 메뉴)        ║
╠══════════════════════════════════════════════════════════╣
║                                                          ║
║  [1] Vina Docking            (AutoDock Vina 실행)        ║
║  [2] Vina Postprocess        (결과 파싱/클러스터링)      ║
║  [3] PPI Docking             (PyRosetta PPI 도킹)        ║
║      3a. PDB 준비 (dimer + MYO1D 합치기)                 ║
║      3b. 도킹 실행                                       ║
║      3c. 결과 원복 (chain 번호 정상화)                   ║
║  [4] MD Analysis             (GROMACS 분석)              ║
║  [5] Generate Report         (종합 보고서 생성)          ║
║  [6] Validate Outputs        (출력 검증)                 ║
║  [7] Site Verdict            (유효 사이트 자동 판정)     ║
║  [8] PPI Postprocess         (PPI 후처리 자동화)         ║
║  [9] Full Pipeline           (1→2→7→5→6 자동 실행)       ║
║                                                          ║
║  [q] Quit                                                ║
║                                                          ║
╚══════════════════════════════════════════════════════════╝
```

Vina 도킹을 실행하려면 `1`을 선택합니다. 이후 `vina_executor.py`의 인터랙티브 모드가 시작되어
Receptor, Ligand, 영역, 파라미터를 순서대로 선택하고 최종 확인 후 실행됩니다.

- Receptor/Ligand가 1개뿐이면 자동 선택
- Ligand는 복수 선택 가능 (예: `1,2` 또는 `all`)

### 4.1 CLI 모드 -- 서브커맨드 실행

YAML 프로젝트 설정 파일(`config/example-project.yaml`)을 사용합니다.

```bash
# Vina 도킹
python main.py -c config/example-project.yaml vina

# 후처리 (파싱 + 클러스터링 + 비교)
python main.py -c config/example-project.yaml postprocess

# 사이트 판정
python main.py -c config/example-project.yaml verdict

# 보고서 생성
python main.py -c config/example-project.yaml report

# 출력 검증
python main.py -c config/example-project.yaml validate

# 전체 파이프라인 (vina → postprocess → verdict → report → validate)
python main.py -c config/example-project.yaml full
```

### 4.2 Config 파일

프로젝트 설정은 `config/` 디렉토리의 YAML 파일로 관리합니다.

```bash
config/example-project.yaml   # 메인 프로젝트 설정
```

Config 파일에는 receptor 목록, ligand 경로, 도킹 파라미터, 후처리 설정 등이 포함됩니다.
인터랙티브 모드에서는 사용 가능한 config 파일 목록이 자동으로 표시됩니다.

### 4.3 Focused Docking -- 좌표 직접 지정

`vina_executor.py`의 인터랙티브 모드(메뉴 [1])를 통해 focused docking 실행 시
center 좌표와 box size를 직접 입력할 수 있습니다.

### 4.4 파라미터 조정

파라미터는 config YAML 파일에서 설정하거나, 인터랙티브 모드에서 단계별로 지정합니다.

주요 파라미터:

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `exhaustiveness` | 128 | 탐색 강도 (높을수록 정확, 느림) |
| `n_poses` | 20 | 생성할 pose 수 |
| `padding` | 5.0 | blind 모드 box 여유 (Angstrom) |
| `min_box` | 70.0 | blind 모드 최소 box 크기 (Angstrom) |
| `contact_cutoff` | 4.0 | 접촉 잔기 판정 거리 (Angstrom) |
| `pocket_cutoff` | 4.0 | 포켓 클러스터링 거리 (Angstrom) |

### 4.5 Receptor PDB -> PDBQT 변환

인터랙티브 모드(메뉴 [1])에서 새 receptor PDB를 넣으면 자동 변환을 제안합니다.

### 4.6 출력 디렉토리

출력은 config에 정의된 프로젝트 루트 아래에 생성됩니다.
각 실행의 설정은 output 디렉토리 내에 `config.yaml`로 자동 저장됩니다.

---

## 5. 후처리 체인 (Postprocess)

Vina 도킹 완료 후, 결과를 구조화된 CSV로 변환하는 6단계 후처리 체인이 있습니다.

```
  parse_poses → contacts → cluster → summarize → compare → ppi
```

### 5.1 각 단계 설명

| 순서 | 모듈 | 입력 | 출력 | 설명 |
|------|------|------|------|------|
| 1 | `parse_poses.py` | Vina PDBQT 결과 | `vina_pose_table.csv` | 포즈별 affinity, RMSD, 좌표 파싱 |
| 2 | `pose_contacts.py` | pose table + receptor PDB | pose table (enriched) | 각 포즈의 접촉 잔기 추출 |
| 3 | `pocket_cluster.py` | enriched pose table | pose table (clustered) | 공간 기반 포켓 클러스터링 |
| 4 | `pocket_summary.py` | clustered pose table | `vina_pocket_table.csv`, `vina_drug_pocket_map.csv` | 포켓별 요약 통계 |
| 5 | `cross_receptor.py` | pocket tables (다중 receptor) | `vina_pocket_comparison.csv` | 교차 receptor 포켓 비교 |
| 6 | PPI 잔기 추출 | PPI 도킹 결과 | PPI 잔기 CSV | PyRosetta 결과에서 인터페이스 잔기 표준화 |

### 5.2 실행 방법

```bash
# CLI 모드 — 전체 후처리
python main.py -c config/example-project.yaml postprocess

# 인터랙티브 모드 — 메뉴 [2] 선택 후 단계별 또는 전체(a) 실행
python main.py
# → [2] 선택 → config 파일 선택 → [a] 전체 실행
```

### 5.3 사이트 판정 (Verdict)

후처리 완료 후, `verdict` 모듈이 Vina 포켓과 PPI 결과를 종합하여 각 사이트의 증거 수준을 자동 판정합니다.

- **출력**: `cross_method_agreement.csv`, `valid_sites.csv`
- **판정 등급**: STRONG (>=55점) / MODERATE (>=30점) / WEAK (<30점)
- **적응적 점수 배분**: PPI 데이터가 있으면 Vina(50)+PPI(20)+Cross(30)=100, 없으면 Vina(60)+Cross(40)=100

```bash
python main.py -c config/example-project.yaml verdict
```

---

## 6. Blind vs Focused Docking

### 6.1 알고리즘

두 모드 모두 동일한 AutoDock Vina 엔진을 사용합니다.
차이점은 **탐색 공간(search box)의 크기와 위치**뿐입니다.

**Vina 탐색 알고리즘: Iterated Local Search (ILS)**
1. 리간드의 위치/방향/torsion을 랜덤하게 변형 (mutation)
2. 변형 후 BFGS quasi-Newton 방법으로 local optimization
3. Metropolis criterion으로 수락/거부 결정
4. exhaustiveness 횟수만큼 독립적으로 병렬 실행
5. 모든 결과에서 에너지 순 정렬 + RMSD 클러스터링 -> 상위 pose 출력

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

`egfr_pipeline/vina/vina_executor.py` 내의 `REGION_PRESETS` 딕셔너리에 등록된 좌표:

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

인터랙티브 모드(메뉴 [1])에서 focused docking 선택 시 프리셋 목록이 자동 표시됩니다.

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

인터랙티브 모드에서 center 좌표와 box size를 입력합니다.

### 7.3 Blind Docking 결과에서 추출

1. Blind docking 실행 (메뉴 [1] 또는 `python main.py vina`)
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
# 각 포켓의 중심 좌표를 focused docking에서 사용
```

---

## 8. 다중 포켓 탐색 (--n-pockets)

블라인드 도킹에서 여러 결합 부위를 자동으로 발견합니다.
인터랙티브 모드(메뉴 [1])에서 "포켓 탐색" 옵션을 선택하거나, vina_executor.py 내부 설정으로 지정합니다.

**작동 원리:**
포즈를 에너지 순서대로 하나씩 처리하면서 포켓에 배정합니다:
1. 포즈 수 자동 증가 (충분한 다양성 확보)
2. 각 포즈의 centroid -> 가장 가까운 열린 포켓에 배정
3. 포켓이 `max-per-pocket`개 차면 **닫힘** -> 다음 포즈는 다른 포켓으로
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

→ 5개 포켓 모두 완성!
```

**파라미터:**
| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `n_pockets` | - | 탐색할 포켓 수 |
| `max_per_pocket` | 3 | 포켓당 포즈 수 (가득 차면 닫힘) |
| `cluster_radius` | 5.0 | 같은 포켓 판정 반경 (Angstrom) |

**출력 파일:**
- `<ligand>_blind.pdbqt` -- 전체 포즈 (기존)
- `<ligand>_blind_pockets.pdbqt` -- 각 포켓의 best pose만 (포켓 수만큼)
- `<ligand>_blind_all_pockets.pdbqt` -- 모든 포켓의 전체 포즈

---

## 9. ANP (ATP 포켓 차단)

receptor PDB에 ANP(adenylyl-imidodiphosphate)가 포함되어 있으면,
ATP 결합 포켓이 ANP로 점유되어 리간드가 해당 포켓에 결합하지 못합니다.

- 파일명에 "anp"가 포함된 receptor PDB를 사용하면 자동으로 ANP 원자 수를 확인
- ANP가 없으면 경고 출력 (ATP 포켓이 열려있어 리간드가 잘못 결합할 수 있음)

---

## 10. 통합 전처리 시스템

도킹 실행 시 `input/` 폴더의 **모든** 분자 파일을 자동으로 스캔 -> 분류 -> 변환합니다.
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

### Receptor (PDB/MOL2/CIF -> PDBQT)
```
input/3gt8_monomer_anp.pdb  →  input/3gt8_monomer_anp_receptor.pdbqt
input/protein.mol2          →  input/protein_receptor.pdbqt
```
변환 도구 시도 순서: ADFR -> MGLTools -> OpenBabel

### Ligand (SDF/MOL2 -> PDBQT)
```
input/97806_ligand.sdf      →  input/97806_ligand.pdbqt
input/173940_ligand.mol2    →  input/173940_ligand.pdbqt
```
변환 도구 시도 순서: Meeko -> ADFR -> MGLTools -> OpenBabel

### SMILES 입력 (.smi 파일)

SMILES 문자열이 담긴 `.smi` 파일을 `input/`에 넣으면 자동으로 3D 구조 생성 -> PDBQT 변환됩니다.

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

인터랙티브 모드(메뉴 [1])에서 SMILES 입력 옵션을 사용합니다.

**3D 구조 생성 도구** (하나 이상 필요):
| 도구 | 설치 | 비고 |
|------|------|------|
| RDKit | `conda install -c conda-forge rdkit` | ETKDGv3 + MMFF 최적화 (권장) |
| OpenBabel | `conda install -c conda-forge openbabel` | `--gen3d --minimize` |

### 사용 흐름
1. `input/`에 receptor PDB와 ligand 파일을 넣는다 (SDF/MOL2/SMI 중 아무거나)
2. `python main.py` 실행, 메뉴 [1] 선택
3. 자동으로 SMILES->SDF->PDBQT 변환 -> 도킹 실행

변환 도구 설치:
```bash
pip install meeko                       # 리간드 SDF→PDBQT 변환 (권장)
conda install -c conda-forge openbabel  # 범용 변환 + SMILES 3D 생성
conda install -c conda-forge rdkit      # SMILES 3D 생성 (선택)
```

---

## 11. 출력 결과 해석

### 11.1 결과 테이블

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

```

| 컬럼 | 의미 |
|------|------|
| Mode | pose 번호 (에너지 순 정렬) |
| Affinity | 결합 에너지 (kcal/mol, 음수일수록 강한 결합) |
| RMSD_lb | 최적 pose 대비 RMSD 하한 (lower bound) |
| RMSD_ub | 최적 pose 대비 RMSD 상한 (upper bound) |

### 11.2 해석 가이드

| Affinity 범위 | 해석 |
|---------------|------|
| -10 이하 | 매우 강한 결합 (드묾) |
| -8 ~ -10 | 강한 결합 |
| -6 ~ -8 | 중간 결합 |
| -4 ~ -6 | 약한 결합 |
| -4 이상 | 유의미하지 않음 |

### 11.3 결과 파일

- `<ligand>_<mode>.pdbqt` -- 모든 pose가 포함된 PDBQT 파일
- `config.yaml` -- 실행 파라미터 (재현용)

### 11.4 후처리 출력 파일

| 파일 | 생성 단계 | 설명 |
|------|----------|------|
| `vina_pose_table.csv` | parse_poses | 전체 포즈 테이블 (affinity, 좌표, 접촉잔기, 포켓 ID) |
| `vina_pocket_table.csv` | summarize | 포켓별 요약 (평균 affinity, 잔기 목록) |
| `vina_drug_pocket_map.csv` | summarize | 리간드-포켓 매핑 |
| `vina_pocket_comparison.csv` | compare | 교차 receptor 포켓 비교 |
| `cross_method_agreement.csv` | verdict | Vina-PPI 교차 방법 일치도 |
| `valid_sites.csv` | verdict | 유효 사이트 판정 결과 (STRONG/MODERATE/WEAK) |
| `project_report.txt` | report | 종합 보고서 |

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

매칭 규칙: `3gt8.pml` -> `3gt8_receptor.pdbqt` (파일명 stem 매칭).
`.pml`이 1개면 모든 receptor에 적용.

---

## 12. 새 리간드/receptor 추가하기

### 새 리간드 추가

1. SDF 또는 MOL2 파일을 `input/`에 넣으면 자동 변환됩니다.
   또는 수동 변환:
   ```bash
   obabel new_compound.sdf -O input/new_compound_ligand.pdbqt --gen3d
   ```

2. `input/` 폴더에 `<ID>_ligand.pdbqt` 네이밍으로 저장

3. 도킹 실행 (자동 감지됨):
   ```bash
   python main.py -c config/example-project.yaml vina
   # 또는 인터랙티브: python main.py → [1]
   ```

### 새 receptor 추가

1. PDB 파일을 `input/`에 저장
2. 인터랙티브 모드(메뉴 [1])에서 receptor를 선택하면 자동 변환을 제안합니다
3. 변환된 `new_receptor_receptor.pdbqt`가 `input/`에 생성됨

---

## 13. 전체 파이프라인 워크플로우

권장되는 전체 실행 흐름:

```
python main.py -c config/example-project.yaml full
```

이 명령은 다음을 순차 실행합니다:

```
Step 1: Vina Docking          → 모든 receptor x ligand 조합 도킹
Step 2: Postprocess           → parse → contacts → cluster → summarize → compare → ppi
Step 3: Site Verdict          → Vina + PPI 교차 판정 (STRONG/MODERATE/WEAK)
Step 4: Generate Report       → 종합 보고서 (project_report.txt)
Step 5: Validate              → 출력 파일 무결성 검증
```

인터랙티브 모드에서는 메뉴 [9]를 선택하면 동일하게 실행됩니다.

---

## 14. 문제 해결

| 증상 | 원인 | 해결 |
|------|------|------|
| "복수 receptor 감지됨" | input/에 *_receptor.pdbqt가 2개 이상 | config에서 receptor를 지정 |
| "Blind 모드에는 receptor PDB가 필요합니다" | PDB 파일을 자동 유추할 수 없음 | config에서 receptor_pdb 경로 지정 |
| "Focused 모드에는 center 지정 필수" | center 좌표 누락 | 인터랙티브 모드에서 좌표 입력 또는 region 프리셋 사용 |
| "pyyaml 미설치" | PyYAML 없음 | `pip install pyyaml` |
| "receptor PDBQT 변환 실패" | 변환 도구 미설치 | ADFR, MGLTools, 또는 OpenBabel 설치 |
| ANP WARNING | receptor에 ANP 미포함 | ANP 포함 PDB 사용 또는 의도된 경우 무시 |
| "SMILES->SDF 변환 실패" | RDKit/OpenBabel 미설치 | `conda install -c conda-forge openbabel` 또는 `rdkit` 설치 |
| config 파일을 못 찾음 | YAML 경로 오류 | `config/example-project.yaml` 경로 확인 |
| 서브커맨드 인식 안됨 | 오타 | `python main.py --help`로 사용 가능한 서브커맨드 확인 |

---

## 15. 구버전(legacy) 명령 대응표

이전 `run_docking.py` 기반 명령은 다음과 같이 변경되었습니다:

| 구버전 | 신버전 |
|--------|--------|
| `python run_docking.py` | `python main.py` (메뉴 [1]) |
| `python run_docking.py --mode blind` | `python main.py -c config/example-project.yaml vina` |
| `python run_docking.py --config output/.../config.yaml` | `python main.py -c config/example-project.yaml vina` |
| `python parse_vina_results.py ...` | `python main.py -c config/example-project.yaml postprocess` (단계 1) |
| `python extract_contacts.py ...` | `python main.py -c config/example-project.yaml postprocess` (단계 2) |
| `python cluster_pockets.py ...` | `python main.py -c config/example-project.yaml postprocess` (단계 3) |
| `python summarize_pockets.py ...` | `python main.py -c config/example-project.yaml postprocess` (단계 4) |
| `python compare_pockets.py ...` | `python main.py -c config/example-project.yaml postprocess` (단계 5) |
| `python generate_report.py ...` | `python main.py -c config/example-project.yaml report` |
| `python validate_outputs.py ...` | `python main.py -c config/example-project.yaml validate` |

이 저장소에는 `legacy/` 디렉토리가 없습니다.
구버전 단독 실행 흐름이 필요하면 `../autodot_vina/`를 참고하세요.
