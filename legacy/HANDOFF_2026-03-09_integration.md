# 인수인계: 데이터 통합 프로세스 구현 계획

> 작성일: 2026-03-09
> 이전 세션에서 완료된 작업과 다음 단계를 상세히 기술한다.
> 새 세션에서 이 파일과 README.md, CLAUDE.md를 먼저 읽을 것.

---

## 1. 현재 상태 요약

### 1.1 완료된 작업

#### 리팩터링 (이번 세션)
- 21개 루트 Python 파일 → `egfr_pipeline/` 패키지로 통합
- `main.py` 통합 CLI 생성 (인터랙티브 메뉴 + 서브커맨드)
- 14개 문서 → `docs/`로 정리, 원본 28개 → `legacy/`로 보관
- 6개 입력 파일 → `input/receptors/` + `input/ligands/`로 정리
- README.md 전면 재작성

#### PPI 도킹 세팅 (이번 세션)
- `egfr_pipeline/ppi/prepare_dimer_pdb.py` 생성
  - EGFR dimer (chain A+B) → chain A로 병합 (B에 +1000 offset)
  - MYO1D 도메인 → chain B로 추가
  - 도킹 후 chain 원복 기능 (restore_chains, restore_csv)
- 준비된 PDB 2개:
  - `input/PPI/prepared/EGFR_dimer_beta_meander.pdb` (616 + 47 res)
  - `input/PPI/prepared/EGFR_dimer_TH1.pdb` (616 + 206 res)
- 막 인접 제외 잔기: 사용자 제공 72개 잔기
- config 분리: test (1K) / prod (20K)
- PBS 스크립트: `run_ppi_test.pbs` / `run_ppi_prod.pbs`
  - `#PBS -q workq`, `nodes=node05:ppn=16`, `walltime=999:00:00`
  - `source ~/.bashrc && conda activate pyrosetta` 필수

#### 현재 테스트 진행 중
- 서버 `/work4/hwang/onepack/codex_ligand/`에서 PPI 테스트 실행 중
- `qsub config/run_ppi_test.pbs` (beta-meander, 1K 모델)
- 결과는 `EGFR_dimer_beta_meander/` 디렉토리에 생성될 예정

### 1.2 Git 상태

```
최신 커밋: ca4d9ad (Reduce PBS ppn from 32 to 16)
브랜치: main
원격: https://github.com/videodrake/codex_ligand
서버 경로: /work4/hwang/onepack/codex_ligand/
```

### 1.3 파이프라인 자동화 수준 (현재)

| 단계 | 자동화 | 비고 |
|------|--------|------|
| 포즈 찾기 (Vina) | 95% | 완전 자동 |
| 사이트 찾기 (클러스터링) | 90% | 완전 자동 |
| 교차 분석 (Vina↔PPI) | 40% | 정량적 교차 검증 없음 |
| 유효 사이트 판정 | 0% | 자동 판정 모듈 없음 |
| PPI 후처리 통합 | 50% | chain 원복, 추출 트리거 수동 |

---

## 2. 구현 계획: 데이터 통합 프로세스

### 2.1 전체 아키텍처 (목표)

```
Phase 1: 도킹 실행
  ├── Vina: 3 receptors × 3 ligands → 포즈
  ├── PyRosetta PPI: EGFR dimer × beta-meander → PPI 포즈
  └── PyRosetta PPI: EGFR dimer × TH1 → PPI 포즈

Phase 2: 사이트 발견
  ├── Vina 후처리: parse → contacts → cluster → summarize
  ├── PPI 후처리: chain 원복(자동) → 잔기 추출(자동)
  └── 교차 수용체 비교: Vina 포켓 간 비교

Phase 3: 교차 검증 ← 새로 구현
  ├── Vina 포켓 잔기 ↔ PPI 인터페이스 잔기 매칭
  ├── 교차 방법 합의도 점수 계산
  └── 교차 수용체 일관성 평가

Phase 4: 자동 판정 ← 새로 구현
  ├── 포켓별 신뢰도 점수 (다중 증거 기반)
  ├── VALID / UNCERTAIN / INVALID 판정
  └── valid_sites.csv 출력
```

### 2.2 새로 만들 모듈

#### 모듈 1: `egfr_pipeline/verdict.py` (핵심)

**목적**: 모든 증거를 종합하여 각 포켓의 유효성을 자동 판정

**입력 파일:**
- `vina_pocket_table.csv` — Vina 포켓 (affinity, pose 수, 잔기)
- `vina_pocket_comparison.csv` — 교차 수용체 비교
- `vina_drug_pocket_map.csv` — 리간드→포켓 매핑
- `ppi_pyrosetta_residues.csv` — PPI 인터페이스 잔기 (있으면)
- `ppi_pyrosetta_summary.csv` — PPI 요약 (있으면)

**출력 파일:**
- `cross_method_agreement.csv` — 포켓별 Vina↔PPI 잔기 일치도
- `valid_sites.csv` — 최종 판정 테이블

**판정 로직 (초안):**

```python
def score_pocket(pocket_row, ppi_residues, comparison_rows):
    score = 0.0
    reasons = []

    # 1. Vina 포켓 품질 (40점)
    #    - best_affinity < -7.0 kcal/mol → +15
    #    - n_pose >= 5 → +10
    #    - n_ligand >= 2 → +15 (다수 리간드가 같은 포켓)

    # 2. PPI 잔기 일치도 (30점) — PPI 데이터 있을 때만
    #    - Vina 포켓 union_contact_residues와 PPI 잔기의 Jaccard
    #    - jaccard >= 0.3 → +15
    #    - overlap_coeff >= 0.5 → +15

    # 3. 교차 수용체 일관성 (30점)
    #    - same_patch_candidate가 2개 이상 receptor에서 발견 → +15
    #    - 3개 receptor 모두에서 발견 → +30

    # 판정:
    #   score >= 60 → VALID
    #   score >= 30 → UNCERTAIN
    #   score < 30  → INVALID

    return score, verdict, reasons
```

**`cross_method_agreement.csv` 스키마:**

```
receptor_id, pocket_id, n_vina_residues, n_ppi_residues,
n_shared_residues, jaccard, overlap_coeff,
shared_residue_list, ppi_mean_occupancy_of_shared,
vina_best_affinity, ppi_best_dg, agreement_level
```

**`valid_sites.csv` 스키마:**

```
receptor_id, pocket_id, verdict, confidence_score,
vina_quality_score, ppi_agreement_score, cross_receptor_score,
best_affinity, n_pose, n_ligand, n_shared_with_ppi,
cross_receptor_matches, reasons
```

**함수 구조:**

```python
# egfr_pipeline/verdict.py

def load_all_evidence(config_path: str) -> dict:
    """모든 증거 파일 로드."""

def compute_cross_method_agreement(
    pocket_table: List[dict],
    ppi_residues: List[dict],
) -> List[dict]:
    """Vina 포켓 잔기와 PPI 인터페이스 잔기의 일치도 계산."""

def compute_cross_receptor_consistency(
    comparison_rows: List[dict],
) -> Dict[str, List[str]]:
    """같은 포켓이 몇 개 receptor에서 나타나는지."""

def score_pocket(
    pocket: dict,
    ppi_agreement: Optional[dict],
    cross_receptor: List[str],
) -> Tuple[float, str, List[str]]:
    """포켓 점수 계산 → (score, verdict, reasons)."""

def generate_verdict(
    config_path: str,
    output_dir: Optional[str] = None,
) -> Tuple[Path, Path]:
    """메인 함수: 전체 판정 실행."""
    # 1. load_all_evidence()
    # 2. compute_cross_method_agreement()
    # 3. compute_cross_receptor_consistency()
    # 4. 각 포켓에 score_pocket() 적용
    # 5. cross_method_agreement.csv 저장
    # 6. valid_sites.csv 저장
    # 7. 콘솔에 요약 출력
```

#### 모듈 2: PPI 후처리 자동화

**현재 문제:** PyRosetta 도킹 완료 후 다음 단계가 전부 수동
- chain 원복 (restore_chains, restore_csv)
- PPI 잔기 추출 (extract_pyrosetta_batch)
- 리포트 재생성

**해결: `egfr_pipeline/ppi/postprocess_ppi.py`**

```python
def postprocess_ppi_results(
    config_path: str,
    docking_output_dir: str,    # e.g., "EGFR_dimer_beta_meander/"
    mapping_csv: str,           # e.g., "input/PPI/prepared/*_mapping.csv"
    partner_name: str,          # e.g., "beta_meander"
):
    """도킹 완료 후 자동 후처리 파이프라인."""

    # Step 1: chain 원복 — 모든 final_result PDB + final_ranking.csv
    restored_dir = restore_all_results(docking_output_dir, mapping_csv)

    # Step 2: 원복된 결과를 ppi.pyrosetta_result_dirs 형식으로 정리
    organize_for_extraction(restored_dir, config_path, partner_name)

    # Step 3: PPI 잔기 추출 자동 실행
    extract_pyrosetta_batch(config_path)

    # Step 4: 리포트 재생성
    generate_report(config_path)

    print(f"PPI 후처리 완료: {restored_dir}")
```

**`restore_all_results()` 상세:**

```python
def restore_all_results(docking_dir: Path, mapping_csv: Path) -> Path:
    """도킹 출력 디렉토리 전체를 원복."""
    restored_dir = docking_dir / "restored"

    # 1. final_result/*.pdb → restored/final_result/*.pdb
    for pdb in (docking_dir / "final_result").glob("Rank*.pdb"):
        restore_chains(pdb, mapping_csv, restored_dir / "final_result" / pdb.name)

    # 2. final_ranking.csv → restored/final_ranking.csv
    restore_csv(
        docking_dir / "final_ranking.csv", mapping_csv,
        restored_dir / "final_ranking.csv",
        columns=["binding_residues_A"]
    )

    # 3. cluster_results/cluster_summary.csv → restored/
    if (docking_dir / "cluster_results" / "cluster_summary.csv").exists():
        restore_csv(...)

    return restored_dir
```

### 2.3 main.py 업데이트 계획

**현재 Full Pipeline (옵션 7):**
```
1. Vina Docking
2. Vina Postprocess (전체)
3. Report 생성
4. Validate
```

**변경 후:**
```
1. Vina Docking
2. Vina Postprocess (전체)
3. PPI 후처리 (chain 원복 + 잔기 추출) — 결과 있으면 자동
4. 교차 검증 + 판정 (verdict.py) ← 신규
5. Report 생성 (verdict 결과 포함)
6. Validate
```

**메뉴 업데이트:**
```
║  [1] Vina Docking            (AutoDock Vina 실행)        ║
║  [2] Vina Postprocess        (결과 파싱/클러스터링)      ║
║  [3] PPI Docking             (PyRosetta PPI 도킹)        ║
║  [4] MD Analysis             (GROMACS 분석)              ║
║  [5] Generate Report         (종합 보고서 생성)          ║
║  [6] Validate Outputs        (출력 검증)                 ║
║  [7] Site Verdict            (유효 사이트 자동 판정) ← 신규║
║  [8] Full Pipeline           (1→2→7→5→6 자동 실행)       ║
```

### 2.4 report.py 업데이트 계획

**현재:** 4개 섹션 (Vina 포켓, 교차 비교, 리간드 매핑, PPI 보조)
**추가:** Section 5 — Site Verdict 요약

```
Section 5: Automated Site Verdict
=================================
Pocket P001 (3GT8_raw):  VALID    (score=72/100)
  - Vina: affinity=-8.2, 12 poses, 3 ligands
  - PPI: 8/15 interface residues overlap (jaccard=0.42)
  - Cross-receptor: same patch in 3GT8_cl38_48 (P003), 3GT8_cl85_100 (P002)

Pocket P002 (3GT8_raw):  UNCERTAIN (score=45/100)
  - Vina: affinity=-6.1, 3 poses, 1 ligand
  - PPI: no data available
  - Cross-receptor: no match found

Pocket P005 (3GT8_cl38_48): INVALID (score=15/100)
  - Vina: affinity=-4.3, 1 pose, 1 ligand
  - PPI: 0 interface residues overlap
  - Cross-receptor: no match found
```

---

## 3. 구현 순서 (우선순위)

### Step 1: 테스트 결과 확인 (서버에서)
```bash
# 테스트 완료 확인
ls EGFR_dimer_beta_meander/final_ranking.csv
cat EGFR_dimer_beta_meander/final_result/docking_validation_report.txt
```

### Step 2: PPI 후처리 자동화 (`postprocess_ppi.py`)
- chain 원복 자동화
- 잔기 추출 자동 트리거
- 파일 위치: `egfr_pipeline/ppi/postprocess_ppi.py`

### Step 3: verdict.py 구현
- cross_method_agreement.csv 생성
- valid_sites.csv 생성
- 파일 위치: `egfr_pipeline/verdict.py`

### Step 4: main.py + report.py 업데이트
- 메뉴에 [7] Site Verdict 추가
- Full Pipeline에 verdict 포함
- report에 Section 5 추가

### Step 5: TH1 도메인 테스트
```bash
qsub -v CONFIG_FILE=config/ppi_test_TH1.ini config/run_ppi_test.pbs
```

### Step 6: 프로덕션 실행
```bash
qsub -v RUN_MODE=both config/run_ppi_prod.pbs
```

---

## 4. 기술 참조

### 4.1 핵심 파일 경로

| 파일 | 역할 |
|------|------|
| `main.py` | 통합 CLI 진입점 |
| `egfr_pipeline/config.py` | 공유 설정 로딩 |
| `egfr_pipeline/residue_utils.py` | 잔기 정규화 |
| `egfr_pipeline/vina/compare.py` | 교차 수용체 비교 |
| `egfr_pipeline/ppi/prepare_dimer_pdb.py` | dimer 병합/원복 |
| `egfr_pipeline/ppi/pyrosetta_extract.py` | PPI 잔기 추출 |
| `egfr_pipeline/report.py` | 보고서 생성 |
| `egfr_pipeline/validate.py` | 출력 검증 |
| `egfr_pipeline/verdict.py` | **신규** — 자동 판정 |
| `egfr_pipeline/ppi/postprocess_ppi.py` | **신규** — PPI 후처리 자동화 |

### 4.2 데이터 계약 (CSV 스키마)

**기존 (변경 없음):**
- `vina_pose_table.csv`: receptor_id, ligand_id, pose_rank, affinity, centroid_xyz, pocket_id, contact_residues
- `vina_pocket_table.csv`: receptor_id, pocket_id, n_pose, n_ligand, best/mean_affinity, union_contact_residues, top_residues
- `vina_drug_pocket_map.csv`: receptor_id, ligand_id, dominant_pocket_id, best_affinity
- `vina_pocket_comparison.csv`: receptor_a/b, pocket_a/b, centroid_dist, residue_jaccard, same_patch_candidate
- `ppi_pyrosetta_residues.csv`: receptor_id, residue_id, frequency, occupancy
- `combined_residue_evidence.csv`: residue_id, evidence_sources (vina+ppi)

**신규:**
- `cross_method_agreement.csv`: receptor_id, pocket_id, n_vina_residues, n_ppi_residues, n_shared, jaccard, overlap_coeff, agreement_level
- `valid_sites.csv`: receptor_id, pocket_id, verdict, confidence_score, reasons

### 4.3 PPI dimer 관련 참조

**Dimer 병합 규칙:**
- EGFR chain A: 잔기 699-1007 (원래 번호 유지)
- EGFR chain B: 잔기 701-1007 → +1000 offset → 1701-2007
- MYO1D partner: chain B (원래 번호 유지)
- ANP 리간드 (chain C, D): 제거

**막 인접 제외 잔기 (사용자 제공):**
```
Chain A 원본: 709-720,724-731,736-739,747,783-785,799-805,871-873,917-921
Chain B(+1000): 1713-1720,1726-1729,1799-1804,1868-1874,1917-1920
```

**Mapping 파일 위치:**
- `input/PPI/prepared/EGFR_dimer_beta_meander_mapping.csv`
- `input/PPI/prepared/EGFR_dimer_TH1_mapping.csv`

**원복 규칙:**
- PDB: A:699-1007 → A, A:1701-2007 → B, B → C (MYO1D partner)
- CSV: 잔기번호 1701+ → B:701+ 형태로 변환

### 4.4 PyRosetta 도킹 출력 구조 (참고)

```
EGFR_dimer_beta_meander/
├── filter_passed/           # 필터 통과 구조
├── cluster_results/
│   ├── C01_M01_S-15.23.pdb  # 클러스터 대표
│   └── cluster_summary.csv
├── final_result/
│   ├── Rank01_*.pdb          # 최종 랭킹 구조
│   ├── Rank01_*_Energies.csv
│   ├── view_results.pml
│   ├── 2_DETAIL_C01.pml
│   └── docking_validation_report.txt
├── final_ranking.csv         # 종합 랭킹
├── energy_funnel.png
└── 1_OVERVIEW_Clusters.pml
```

### 4.5 서버 환경 정보

```
서버: node04/node05
경로: /work4/hwang/onepack/codex_ligand/
Python: conda환경 "pyrosetta"
PBS: workq 큐, ppn=16
GitHub: https://github.com/videodrake/codex_ligand
git pull 시 gnome-ssh-askpass 경고 발생 → 무시 가능, username/password 수동 입력
로컬 수정 있으면 git stash → pull → stash drop
```

---

## 5. 새 세션 시작 프롬프트

```
이 파일들을 순서대로 읽어:
1. README.md
2. CLAUDE.md
3. docs/HANDOFF_2026-03-09_integration.md

현재 상태: PPI 테스트 도킹이 서버에서 실행 중 또는 완료됨.
다음 작업: 테스트 결과 확인 후 verdict.py (자동 사이트 판정 모듈) 구현.
서버 경로: /work4/hwang/onepack/codex_ligand/
```

---

## 6. 주의사항

1. **PyRosetta 코드 수정 금지**: `egfr_pipeline/pyrosetta_docking/` 내부 파일의 try-except, hasattr 체크, pool.imap 패턴 절대 변경하지 말 것
2. **잔기 정규화 필수**: Vina와 PPI 잔기를 비교할 때 반드시 `residue_utils.normalize_residue_id()` 사용 (HSD→HIS 등)
3. **chain 원복 전 분석 금지**: PPI 도킹 결과의 잔기 1701+는 원래 chain B의 701+임. 원복 없이 Vina 잔기와 비교하면 매칭 안 됨
4. **PPI 데이터 없어도 verdict 작동**: PPI 결과가 없으면 Vina-only 점수로 판정 (PPI 점수 = 0, 나머지로 판단)
5. **config에 verdict 설정 추가 필요**: 판정 임계값을 config.yaml에 넣어야 사용자가 조정 가능
