# Workflow A↔B 비교 모듈 설계 (AC-5.1)

## 1. 목적

Workflow A(Verdict 3축)와 Workflow B(Phase 4 Perturbation 4축)의 포켓 결과를 체계적으로 비교하여 Consensus/A-only/B-only/Conflict로 분류한다.

## 2. 입력 데이터

### Workflow A (Verdict)
| 파일 | 역할 | 핵심 컬럼 |
|------|------|-----------|
| `valid_sites.csv` | 포켓 판정 | pocket_id, receptor_id, verdict, confidence_score, vina_quality_score, ppi_proximity_score, is_atp_site, allosteric_candidate |
| `vina_pocket_table.csv` | 포켓 centroid + 잔기 | centroid_x/y/z, union_contact_residues |

### Workflow B (Phase 4)
| 파일 | 역할 | 핵심 컬럼 |
|------|------|-----------|
| `phase4_final_review_table.csv` | 포켓 순위 | pocket_id, receptor_id, mechanistic_class, perturbation_score, rank |
| `candidate_pockets.csv` (Phase 2) | 포켓 centroid + 잔기 | centroid_x/y/z, residue_ids |

## 3. 매칭 방법

**선택: 옵션 1 (centroid + Jaccard)**

```
match = centroid_distance < 8.0 Å  AND  residue_jaccard ≥ 0.3
```

- centroid_distance: 유클리드 거리 (vina_pocket_table centroid vs candidate_pockets centroid)
- residue_jaccard: |A ∩ B| / |A ∪ B| (union_contact_residues vs residue_ids)

### 매칭 절차
1. Workflow A 포켓 중 `is_atp_site=True`인 포켓 제외
2. Workflow B 포켓은 모두 포함 (irrelevant 포함 — A-only allosteric 식별용)
3. A×B 전체 쌍에 대해 centroid 거리 + Jaccard 계산
4. 가장 가까운 매칭 우선 (greedy, 거리 기준)
5. 매칭 실패(EC-5.1: 8Å 이내 없음) → "미매칭" 분류

## 4. 분류 기준

| 분류 | A 조건 | B 조건 | 해석 |
|------|--------|--------|------|
| **Consensus** | STRONG 또는 MODERATE | rank ≤ 상위 50% 또는 perturbation_score > 0.5 | 최강 계산적 증거, 실험 최우선 |
| **A-only** | STRONG 또는 MODERATE | 미매칭 또는 irrelevant (low_relevance/uncertain) | Druggable하나 PPI 무관 → allosteric 검토 |
| **B-only** | WEAK 또는 미매칭 | rank ≤ 상위 50% | PPI 근처 얕은 포켓 → blind docking 편향 가능 |
| **Conflict** | STRONG + B=irrelevant, 또는 A=WEAK + B=상위 | 반대 판정 | 메트릭 특이성 → PyMOL 수동 확인 |

### "상위" 정의
- Workflow B에서 `perturbation_score > 0.5` 또는 `rank ≤ ceil(n_pockets / 2)`

### 플래그
- A-only 포켓: `allosteric_candidate=True`이면 "allosteric 후보?" 플래그
- B-only 포켓: "blind docking 편향?" 플래그

## 5. 제외 조건

- `is_atp_site=True` 포켓: Workflow A에서 제외 (실험적 근거, ATP binding 유지)
- Workflow B irrelevant 포켓: 제외하지 않음 (A-only 식별에 필요)

## 6. 출력

| 파일 | 내용 |
|------|------|
| `workflow_comparison.csv` | pocket_id, receptor_id, workflow_a_verdict, workflow_b_class, comparison_category, allosteric_flag, bias_flag, centroid_dist, residue_jaccard |
| `workflow_comparison_report.md` | 분류별 요약, Consensus 목록, Conflict 포켓 상세 |

## 7. 엣지 케이스

| EC | 조건 | 처리 |
|----|------|------|
| EC-5.1 | 8Å 이내 매칭 없음 | "미매칭"으로 분류, 각 워크플로우에서 독립 해석 |
| EC-5.2 | allosteric 후보 0개 | 정상, 보고서에 "allosteric 후보 미식별" 기록 |
| 동일 포켓 다중 매칭 | A 1개가 B 2개에 8Å 이내 | 가장 가까운 1개만 매칭 (greedy) |
| 출력 CSV 미존재 | 서버 미실행 | 합성 데이터로 동작, 빈 결과 반환 |
