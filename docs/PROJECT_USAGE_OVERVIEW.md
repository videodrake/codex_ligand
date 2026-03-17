# EGFR-MYO1D 프로젝트 사용 개요

## 한 줄 요약

EGFR과 MYO1D의 단백질-단백질 상호작용(PPI)을 **약물로 교란할 수 있는 포켓**을 찾는 파이프라인.

---

## 사용하는 도구와 각 도구에서 얻는 정보

| 도구 | 하는 일 | 얻는 정보 |
|------|---------|----------|
| **AutoDock Vina** | 소분자를 EGFR 표면에 도킹 | "EGFR 표면 어디에 약물이 잘 붙는가" (포켓 위치 + 결합 에너지) |
| **PyRosetta** | MYO1D 단백질을 EGFR에 도킹 | "EGFR 표면 어디에 MYO1D가 붙는가" (PPI 인터페이스 잔기) |
| **fpocket / P2Rank** | EGFR 표면의 포켓 형태 분석 | "구조적으로 약물이 들어갈 만한 움푹한 곳이 어디인가" |
| **LightDock** | PyRosetta와 다른 방식으로 단백질 도킹 | "PyRosetta 결과가 다른 방법으로도 재현되는가" (교차 검증) |

---

## 두 가지 워크플로우

### Workflow A: Standard Production — "일단 다 돌려보고 비교"

```
Vina (약물 포켓 탐색)  ──┐
                          ├──→  겹치는 곳 = 약물로 PPI 교란 가능한 후보
PPI (MYO1D 결합 부위)  ──┘
```

**목적**: 약물이 잘 붙는 곳(Vina)과 MYO1D가 붙는 곳(PPI)이 **공간적으로 겹치는지** 확인.
겹치면 → 그 포켓에 약물을 넣으면 MYO1D가 못 붙을 수 있다.

**실행**:
```bash
qsub config/run_production.pbs
```

**결과물**:
- `valid_sites.csv` — 각 포켓의 증거 강도 (STRONG / MODERATE / WEAK)
- `project_report.txt` — 종합 보고서
- `vina_pocket_table.csv` — 포켓별 약물 결합 에너지, 위치, 수렴도
- `ppi_pyrosetta_summary.csv` — PPI 인터페이스 잔기 요약

### Workflow B: Advanced PPI-First — "PPI 결과를 기반으로 좁혀가며 탐색"

```
PPI 도킹 → PPI 인터페이스 정의
              ↓
         포켓 분석 (fpocket/P2Rank으로 약물 가능 포켓 탐지)
              ↓
         Focused Vina (발견된 포켓에만 집중 도킹)
              ↓
         통합 스코어링 (4축: PPI 신뢰도 + 약물성 + 교란 가능성 + 상태 강건성)
```

**목적**: PPI 인터페이스 근처 포켓만 골라서 **정밀 탐색** → 최종 후보 순위.

**실행**:
```bash
# Workflow A의 PPI 도킹 완료 후
qsub config/run_advanced_pipeline.pbs
```

**결과물**:
- `phase4_final_review_table.csv` — 포켓 × 리간드 최종 순위
- 기계적 분류: orthosteric(직접 차단) / rim(가장자리 교란) / allosteric(원격 조절)

---

## 왜 3개 EGFR 구조를 쓰는가?

| 구조 | 설명 |
|------|------|
| 3GT8_raw | X-ray 결정 구조 (정적) |
| EGFR_160-185 | MD 시뮬레이션 클러스터 1 (동적 상태) |
| EGFR_170-200 | MD 시뮬레이션 클러스터 2 (동적 상태) |

단백질은 움직이므로 한 구조만으로는 부족하다. **3개 상태 모두에서 나타나는 포켓/인터페이스가 가장 신뢰할 수 있는 후보**이다.

---

## 실행 환경

- **모든 도킹/연산**: HPC 서버 (`qsub`로 제출)
- **코드 수정/테스트**: 이 Codespace (개발 환경)
- 서버에서는 `conda activate pyrosetta` 환경 사용

---

## 결과 확인 순서

1. **Validation Report** (`docking_validation_report.txt`) → 품질 체크 PASS/FAIL 확인
2. **valid_sites.csv** 또는 **phase4_final_review_table.csv** → 후보 포켓 순위 확인
3. **PyMOL로 시각화** → 후보 포켓이 생물학적으로 타당한 위치인지 확인
   ```bash
   scp -r user@server:/path/to/output/ ~/Desktop/
   pymol 1_OVERVIEW_Clusters.pml
   ```

---

## 요약 다이어그램

```
┌─────────────────────────────────────────────────────┐
│                    입력                               │
│  EGFR PDB (3개) + MYO1D + 소분자 리간드 라이브러리    │
└───────┬─────────────────────────┬───────────────────┘
        │                         │
   ┌────▼─────┐             ┌─────▼──────┐
   │   Vina   │             │  PyRosetta │
   │ 약물이   │             │ MYO1D가    │
   │ 어디에   │             │ 어디에     │
   │ 붙는가?  │             │ 붙는가?    │
   └────┬─────┘             └─────┬──────┘
        │                         │
        └────────┬────────────────┘
                 │
        ┌────────▼────────┐
        │   둘이 겹치면    │
        │   = 약물 후보    │
        └────────┬────────┘
                 │
        ┌────────▼─────────────┐
        │  fpocket / P2Rank    │
        │  그 자리가 정말      │
        │  약물이 들어갈 만한  │
        │  구조인가?           │
        └────────┬─────────────┘
                 │
        ┌────────▼────────┐
        │  Focused Vina   │
        │  그 포켓에 약물  │
        │  정밀 도킹       │
        └────────┬────────┘
                 │
        ┌────────▼────────┐
        │  최종 순위       │
        │  + 기계적 분류   │
        └─────────────────┘
```
