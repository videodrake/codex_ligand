# Skill: 워크플로우 비교 및 문서 정비

## 적용 범위
- Group: 5
- Features: F-5

## 도메인 지식

### 두 워크플로우의 출력 구조
| Workflow A | Workflow B |
|------------|------------|
| `valid_sites.csv` — Verdict 3축 판정 | `phase4_final_review_table.csv` — 4축 최종 순위 |
| 포켓 ID, centroid, 잔기 목록, 총점, 판정 | 포켓 ID, 관계 분류, A1-A4 점수, 최종 순위 |
| STRONG/MODERATE/WEAK | orthosteric/rim/allosteric/irrelevant |

### 포켓 매칭 방법
- **옵션 1** (권장): centroid 거리 < 8Å AND 잔기 Jaccard ≥ 0.3
  - cross_receptor.py의 same-patch 로직 재활용 가능
  - 높은 정확도, 구현 복잡도 중간
- **옵션 2**: centroid 거리만 (< 8Å)
  - 단순하지만 다른 포켓이 근접할 때 오매칭 가능

### 4가지 분류와 해석
| 분류 | 조건 | 해석 | 후속 |
|------|------|------|------|
| Consensus | A=STRONG/MOD + B=상위 | 최강 계산적 증거 | 실험 최우선 후보 |
| A-only | A=STRONG/MOD + B=무관 | Druggable하나 PPI 무관 | Allosteric 후보 검토 |
| B-only | B=상위 + A=WEAK/무관 | 얕은 PPI 근처 포켓 | Blind docking 편향 때문 |
| Conflict | 반대 판정 | 메트릭 특이성 | PyMOL 수동 확인 필수 |

### Allosteric 후보 분류
- 조건: Vina 축 ≥ 35점(만점 50의 70%) AND PPI 축 ≤ 5점(만점 20의 25%)
- 의미: 약물이 잘 결합하지만 PPI 패치와 공간적으로 떨어져 있음
- 가능한 메커니즘: 약물 결합 → 구조 변화 유도 → 간접적 PPI 교란
- 현재 파이프라인으로는 이 메커니즘을 직접 검증할 수 없음 → MD 시뮬레이션 필요

### 방법론적 한계 5개 섹션
1. **Rigid-body docking**: induced fit 미반영, 3 MD 클러스터로 부분 보완
2. **LightDock 독립성**: "method-diverse"이지 "method-independent"가 아님 (같은 입력 구조)
3. **입력 구조 공유**: Vina/PyRosetta/LightDock 모두 동일 3개 PDB → 공통 맹점 가능
4. **Solvent 효과**: implicit solvation만, water-mediated H-bond 무시
5. **Vina scoring 편향**: 소수성 과대평가, 극성 PPI interface 근처 과소평가

### 면책 조항 삽입 방법
- CSV 파일: 첫 행에 `# DISCLAIMER: ...` 주석
- TXT 파일: 최상단에 면책 문단
- downstream 파서가 `#` 주석 행을 무시하는지 확인 필요 → 아니면 별도 metadata 파일

## 주의 사항
- 비교 모듈은 Verdict 가중치 확정(Group 4) 후 실행해야 최종 결과 기반 비교가 됨
- ATP site 포켓(is_atp_site=True)은 비교에서 제외 (실험적 근거)
- Workflow B irrelevant 포켓도 비교에 포함해야 A-only allosteric 후보를 식별할 수 있음
- paths.py deprecation은 주석 추가만 — 함수 자체를 삭제하지 않음
