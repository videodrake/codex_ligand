name: scoring-system
description: Verdict/스코어링 관련 작업 시 로딩. 트리거 — Verdict, 점수, 축, STRONG/MODERATE/WEAK, 스코어링 언급 시. 비트리거 — Vina 도킹 파라미터만 바꿀 때, PPI 도킹 자체를 수정할 때.

## Workflow A: 3축 체계 (verdict.py)
- 축 1 Vina Quality (30점): affinity + convergence + consensus + stability + diversity
- 축 2 PPI Spatial (40점): spatial + overlap + reproducibility
- 축 3 Cross-Receptor (30점): 다중 구조 상태 일관성
- PPI 없을 시: 60 + 0 + 40 = 100으로 적응적 재배분
- STRONG ≥ 55 (최소 2축에서 의미 있는 점수 필요)

## Workflow B: 4축 체계 (phase4/)
- A1: PPI interface 관계 (orthosteric/rim/allosteric/irrelevant)
- A2: Druggability
- A3: Perturbation relevance
- A4: State robustness
- A1+A3 합산 가중치 60% → affinity만 좋고 MYO1D 무관한 포켓은 상위 불가

## 판정 원칙
"증거 분류이지 타당성 판정이 아니다."
STRONG도 PyMOL 시각 검증 필수. WEAK도 cryptic pocket 가능성 있음.

## 상세 참조
- PIPELINE_ARCHITECTURE_REPORT.md — Phase 4 Scoring + Verdict 섹션
- 설계의도.md — 3축/4축 설계 의도
