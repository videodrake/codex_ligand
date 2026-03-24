# EGFR-MYO1D Pipeline: Human Decision Points

파이프라인은 증거를 정리하고, 사람이 판단한다.
이 문서는 자동화로 대체할 수 없는 인간 판단 지점을 단계별로 정리한다.

---

## Phase 1 — Vina Blind Docking

| ID | 검토 항목 | 왜 사람이 봐야 하는가 | 확인 방법 | 판단 기준 |
|----|----------|---------------------|----------|----------|
| H-1.1 | 리간드 화학 다양성 | 3종 리간드가 실제로 다른 화학 공간을 대표하는지는 Tanimoto 수치만으로 판단 불가 | RDKit fingerprint + MW/LogP/TPSA 비교 | Tanimoto < 0.5이면 양호. 0.7 이상이면 consensus 해석 주의 |
| H-1.2 | Search box 적절성 | Box가 전체 kinase domain을 포함하는지 | PyMOL에서 box 시각화 | 모든 C-lobe 표면이 box 안에 있어야 함 |

---

## Phase 2 — PPI Global Blind Docking

| ID | 검토 항목 | 왜 사람이 봐야 하는가 | 확인 방법 | 판단 기준 |
|----|----------|---------------------|----------|----------|
| H-2.1 | Beta-meander 배향 | Orientation filter가 active face를 올바르게 식별하는지 | PyMOL에서 top-ranked 모델 3-5개의 sheet 8/9 side chain 방향 확인 | Side chain이 receptor를 향하면 pass |
| H-2.2 | 비물리적 도킹 배향 | Membrane face나 N-lobe에만 붙은 모델이 결과에 포함되는지 | `1_OVERVIEW_Clusters.pml`에서 클러스터별 위치 확인 | C-lobe surface 클러스터만 의미 있음 |
| H-2.3 | Dimer vs Monomer 입력 확인 | 입력 구조가 실제로 의도한 oligomeric state인지 | ChimeraX에서 입력 PDB 열어서 chain 수 확인 | chain X에 dimer가 숨어있을 수 있음 |
| H-2.4 | 3GT8 V924R mutation 영향 | 결정화 mutation이 PPI interface 에너지에 영향을 주는지 | 3GT8 단독 포켓이 다른 state에서도 나타나는지 비교 | 3GT8에서만 나타나면 mutation artifact 가능성 고려 |

---

## Phase 3 — PPI Postprocess

| ID | 검토 항목 | 왜 사람이 봐야 하는가 | 확인 방법 | 판단 기준 |
|----|----------|---------------------|----------|----------|
| H-3.1 | PPI patch 범위 합리성 | 계산만으로는 spurious contact를 걸러낼 수 없음 | PyMOL에서 top occupancy 잔기 위치 확인 | C-lobe beta-sheet 영역(940-1006)에 집중되어야 함 |
| H-3.2 | Cross-seed 수렴 | 특정 seed에서만 나타나는 patch가 있는지 | `ppi_pyrosetta_residues.csv`의 `frac_runs_supporting` 확인 | 1.0 = 5/5 seed (강건). 0.2 = 1/5 (seed-specific) |
| H-3.3 | Occupancy 분포 해석 | Broad interface가 생물학적 실재인지 sampling 부족인지 | 전체 occupancy 분포 확인 | flat partner(beta-meander)에서는 broad interface가 정상일 수 있음 |

---

## Phase 4 — Vina Postprocess

| ID | 검토 항목 | 왜 사람이 봐야 하는가 | 확인 방법 | 판단 기준 |
|----|----------|---------------------|----------|----------|
| H-4.1 | ATP site 겹침 판정 | 포켓이 ATP pocket 자체인지 인접 영역인지 자동 구분 불가 | PyMOL에서 P010 contact residue를 ATP site 37개 잔기와 시각 비교 | 포켓 중심이 ATP site이면 false positive |
| H-4.2 | 클러스터 병합 적절성 | 잔기 공유로 병합된 포켓이 실제로 같은 cavity인지 | PyMOL에서 병합 포켓의 pose들을 동시에 표시 | 2개 별도 cavity가 병합됐으면 over-merge |
| H-4.3 | Ligand-pocket 배정 | Dominant pocket 배정이 실제 선호와 일치하는지 | `vina_drug_pocket_map.csv`의 dominant_ligand_fraction 확인 | Fraction < 0.3이면 "dominant" 표기가 오해 유발 |

---

## Phase 5 — Site Verdict

| ID | 검토 항목 | 왜 사람이 봐야 하는가 | 확인 방법 | 판단 기준 |
|----|----------|---------------------|----------|----------|
| H-5.1 | STRONG 포켓 시각적 검증 | STRONG은 증거 수렴이지 생물학적 타당성이 아님. 모든 STRONG은 반드시 눈으로 확인 | PyMOL에서 각 STRONG 포켓의 pose + PPI residue 동시 표시 | 접근 가능한 표면 cavity에 있어야 함. 내부 void는 artifact |
| H-5.2 | PPI proximity 해석 | 15A이 소분자로 실제 교란 가능한 거리인지 | PyMOL에서 포켓~PPI surface 사이 구조 확인 | 연결 cavity가 있으면 교란 가능. protein core가 가로막으면 allosteric만 가능 |
| H-5.3 | Cross-state 불일치 해석 | 3GT8 단독 STRONG이 mutation artifact인지 | 3GT8 V924R(PDB 948) 위치와 해당 포켓의 공간 관계 확인 | Mutation에서 10A 이내면 artifact 가능성 높음 |
| H-5.4 | Verdict 가중치 민감도 | 50:20:30 배분에서 순위가 안정적인지 | Plan D-1 민감도 시뮬레이션 결과 검토 | Top 3가 6가지 조합에서 모두 유지되면 robust |

---

## Phase 6-7 — Report & Validation

| ID | 검토 항목 | 왜 사람이 봐야 하는가 | 확인 방법 | 판단 기준 |
|----|----------|---------------------|----------|----------|
| H-6.1 | Known mutation warning 해석 | 경고가 결과에 얼마나 영향을 주는지 | validation_summary.txt 검토 + 해당 잔기가 STRONG 포켓에 포함되는지 확인 | 핵심 잔기에 해당하면 주의 |
| H-6.2 | Optional file 부재 판단 | Optional 파일 누락이 의도인지 실수인지 | Validation warning의 optional 누락 목록 확인 | AFM 미실행이 의도적이면 OK |

---

## Workflow B — Advanced Pipeline

| ID | 검토 항목 | 왜 사람이 봐야 하는가 | 확인 방법 | 판단 기준 |
|----|----------|---------------------|----------|----------|
| H-B.1 | Phase 1->2 핸드오프 품질 | Patch가 pocket proposal에 적절한 입력인지 | PyMOL에서 patch residue와 fpocket 포켓 위치 비교 | Patch 근처에 포켓이 제안되어야 함 |
| H-B.2 | Phase 3 search budget 분배 | Budget이 특정 포켓에 과도 집중되지 않았는지 | phase3_budget_tracking.csv 검토 | 한 포켓이 70% 이상 소모하면 diversity 문제 |
| H-B.3 | Workflow A<->B 결과 비교 | 두 워크플로우가 같은 사이트를 지목하는지 | 비교 모듈 결과 + PyMOL 시각 비교 | 불일치 시 각 워크플로우의 bias 분석 필요 |

---

## 전체 프로젝트 수준

| ID | 검토 항목 | 왜 사람이 봐야 하는가 | 확인 방법 | 판단 기준 |
|----|----------|---------------------|----------|----------|
| H-P.1 | Rigid-body 한계 판단 | CS001이 induced-fit 포켓인지 | MD simulation 또는 문헌 기반 | 현재 파이프라인으로는 답 불가. MD 추가 여부 결정 |
| H-P.2 | 실험 검증 우선순위 | STRONG 5개 중 어떤 것을 먼저 테스트할지 | 전체 결과 + 실험 feasibility 종합 | CS001이 최우선이나 실험 조건에 따라 변경 가능 |
| H-P.3 | 최종 결론 범위 | 무엇을 주장할 수 있는지 | 결과의 한계를 정직하게 평가 | "계산적 증거가 수렴한다"까지만 가능. "교란한다"는 실험 후에만 |
