# Workflow A↔B 불일치 해석 가이드

> workflow_comparison.csv의 결과를 해석할 때 참조한다.
> 불일치 시나리오별 원인과 후속 조치를 기술한다.

## 배경

- **Workflow A** (Verdict): Vina blind docking → 3축 스코어링 (Vina 50 + PPI 20 + Cross 30)
- **Workflow B** (Phase 4): PPI-first → Pocket proposal → Focused docking → 4축 Perturbation scoring

두 워크플로우는 서로 다른 가설 공간을 탐색한다:
- A는 "어디에 약물이 결합할 수 있는가?" (druggability-first)
- B는 "어디에 결합하면 MYO1D를 교란할 수 있는가?" (perturbation-first)

## 불일치 시나리오

### 시나리오 1: A=STRONG, B=irrelevant

**상황:** Workflow A에서 높은 점수를 받았으나, Workflow B에서 PPI와 무관하다고 분류된 포켓.

**가능한 원인:**
- 포켓이 PPI 인터페이스에서 떨어져 있지만 약물 결합이 우수함
- Vina blind docking이 소수성 포켓을 과대평가 (방법론적 한계 #5)
- 포켓이 EGFR 기능(kinase 활성)에는 중요하나 MYO1D 결합과 무관

**후속 조치:**
1. `allosteric_candidate` 플래그 확인 — True이면 allosteric 메커니즘 가능성 검토
2. PyMOL에서 포켓 위치와 PPI 패치 거리를 시각적으로 확인
3. 구조적으로 allosteric 경로가 있는지 MD 시뮬레이션으로 검증 검토

### 시나리오 2: B=상위, A=WEAK

**상황:** Workflow B에서 높은 perturbation score를 받았으나, Workflow A에서 WEAK 판정.

**가능한 원인:**
- PPI 패치 근처 얕은 포켓 — PPI 교란 잠재력은 있으나 약물 결합이 약함
- Vina blind docking에서 탐색되지 않은 포켓 (blind docking 편향)
- 포켓 크기가 작거나 접근성이 낮아 Vina 포즈 수가 부족

**후속 조치:**
1. `bias_flag=True` 여부 확인 — B-only 포켓은 blind docking 편향 가능성
2. 해당 포켓의 druggability tier 확인 (Phase 2)
3. 약물 결합 개선 가능성 검토 (fragment-based 접근 등)

### 시나리오 3: 둘 다 무관심 (A=WEAK, B=하위)

**상황:** 두 워크플로우 모두에서 낮은 점수.

**가능한 원인:**
- 실제로 약물 결합도 PPI 교란도 어려운 사이트
- 두 방법의 공통 맹점에 의해 과소평가된 사이트 (입력 구조 편향)

**후속 조치:**
1. 일반적으로 추가 조사 불필요
2. 단, 생물학적으로 관심이 있는 잔기(실험적으로 알려진 기능성 잔기)가 포함되어 있다면 공통 맹점 가능성 검토

## 일반 원칙

| 일치도 | 신뢰 수준 | 후속 |
|--------|-----------|------|
| **Consensus** (A+B 모두 높음) | 최고 | 실험 최우선 후보 |
| **A-only** (A 높음, B 낮음) | 중간 | Allosteric 여부 검토 |
| **B-only** (B 높음, A 낮음) | 중간-낮음 | Blind docking 편향 감안 |
| **Conflict** (반대 판정) | 최저 | PyMOL 수동 검토 필수 |
| **양쪽 낮음** | 해당 없음 | 일반적으로 제외 |
