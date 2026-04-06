# Methodology Limitations

EGFR-MYO1D PPI 교란 약물 포켓 탐색 파이프라인의 방법론적 한계를 정리한다.
결과 해석 시 반드시 이 문서를 참조한다.

## 1. 구조적 한계

### 1.1 수용체 상태 다양성
- 3개 수용체 상태(3GT8_raw, EGFR_160-185, EGFR_170-200) 중 **실험 결정 구조는 3GT8_raw 1개**뿐이다.
- 나머지 2개는 MD 클러스터 대표 구조로, 독립적 실험 증거가 아니다.
- 3개 상태만으로는 교차 상태 견고성(cross-state robustness) 평가의 일반화가 제한된다.

### 1.2 파일럿 데이터 아티팩트
- 레거시 파일럿 데이터(pilot_data_reference.csv)는 N-lobe 부재, VAL962 N-terminal artifact, orientation filter 미적용 조건에서 생성되었다.
- 현재 방법론과 직접 비교가 불가능하다.

### 1.3 백본 유연성 미반영
- PyRosetta 글로벌 도킹은 fixed-backbone을 가정한다.
- MYO1D beta-meander의 동적 유연성이 모델링되지 않는다.
- 도킹 궤적 수를 늘려도 backbone mismatch 한계는 해소되지 않는다 (`docs/phase1_notes.md:145-150`).

## 2. 도킹 방법론 한계

### 2.1 Vina
- **Blind search box (70 A+)**: 극히 넓은 탐색 공간. exhaustiveness=384(기본값의 48배)로도 약한 결합 포즈를 놓칠 수 있다.
- **친화도 임계값 보정 범위**: C-lobe 표면 포켓(-5 ~ -8 kcal/mol) 기준으로 보정되었다. ATP 포켓(-9 ~ -12 kcal/mol) 수준의 깊은 포켓은 대상이 아니다 (`verdict.py:144-150`).
- **Workflow A/B 결과 직접 비교 불가**: Blind box(A)와 focused box(B)는 탐색 범위가 다르므로 점수를 직접 비교할 수 없다.

### 2.2 PyRosetta PPI
- **DockingSlideIntoContact 누락 시 dG=0.0**: V1.0 역사적 버그. 빠뜨리면 전체 에너지 스코어가 무효화된다.
- **FoldTree 복잡성**: 역직렬화 후 FoldTree 재설정 오류가 하위 스코어링 전체에 전파될 수 있다.
- **Hard mask 적용**: dimer interface 및 ATP site 배제가 절대 필터(soft penalty 아님)로 동작한다.

### 2.3 LightDock (보조 검증)
- **2차 증거 전용**: PyRosetta를 대체하지 않는다. Phase 1 수신 측 패치 리뷰 지원만 담당한다.
- **스코어링 함수 비호환**: PyRosetta REU(dG_separated)와 LightDock DFIRE2는 단위가 다르므로 점수 직접 비교가 불가능하다.
- **Chain 재할당 위험**: LightDock 출력 PDB의 잔기 재할당으로 orientation filter가 실패할 수 있다 → `insufficient_data` fallback.
- **LightDock 단독 잔기 배제**: PyRosetta 뒷받침 없는 LightDock-only 잔기는 Phase 2 패치 입력으로 사용할 수 없다.

## 3. 실험 근거 제약

### 3.1 Ko et al. alanine substitution
- **Sheet 8/9만 실험 검증됨**: Sheet 10-12의 기여는 추론이다.
- **PPI hotspot 임계값**: Ko et al. sheet 8/9 잔기 3개 미만이면 Workflow B 중단 — 이진 판정이며 그래디언트가 없다.

### 3.2 리간드 다양성
- **3종 고정**: 173940, 97806, VAX-C12_0 (쌍별 Tanimoto < 0.4).
- **SAR 구축 불가**: 3개 분자만으로는 구조-활성 관계 도출이나 scaffold 간 경향성 검증이 불가능하다.
- 구조적 편향 방지 목적으로 임의 교체가 금지되어 있다.

### 3.3 ATP 결합 유지
- ATP 포켓은 MYO1D 결합 후에도 유지된다는 실험적 사실에 기반하여, ATP site를 교란 포켓으로 판정하지 않는다.
- 이 가정이 틀릴 경우 잠재적 교란 포켓 하나를 놓치게 된다.

## 4. 스코어링 체계 주의사항

### 4.1 Centroid 거리 편향
- Vina centroid는 포켓 공동(cavity) 내부(표면 아래 ~3-5 A)에 위치한다.
- 실제 포켓 입구-PPI 표면 거리보다 체계적으로 과대 추정된다.
- 임계값은 보수적으로 설정되어 있다 (`verdict.py:1337-1341`).

### 4.2 상태 특이적 포켓 처리
- 한 상태에서만 검출된 포켓은 0.4점(0.0이 아님)으로 처리된다.
- 수용체 상태가 1개뿐인 경우, 진정한 상태 특이성과 과소 샘플링을 구분할 수 없다 (`phase4/state_interpretation.py:59-68`).

### 4.3 PPI 데이터 부재 시 적응 스코어링
- PPI 데이터가 없을 때 가중치가 재배분되어, Vina 단독 포켓이 과대 평가될 수 있다.
- Centroid spread > 15 A인 분산된 PPI 패치에 대해 별도 페널티가 없다.

### 4.4 포켓 도구 단일 의존
- 현재 fpocket만 사용한다. P2Rank 통합 시 druggability 축 판별력이 개선될 수 있다.
- 축 가중치(A1-A4)는 첫 결과 리뷰 후 전문가 조정이 필요하다 (`phase4/score_framework.py:421-432`).

### 4.5 사전 실행 친화도 한계
- Vina 실행 전 모든 리간드 지지 수준은 `pending_*` 상태이다. Vina 실행 후 재점수가 필수이다.

## 5. 교차 방법 검증 한계

- PyRosetta-LightDock 교차 검증은 잔기 수준 Jaccard 계수(중첩)만 사용한다.
- 에너지적 합의(energetic agreement)는 평가되지 않으며, 단위가 다르므로 불가능하다.
- Orientation filter의 ambiguous band 임계값(±0.15)은 경험적 검증이 없으며, 경계 구조가 잘못 분류될 수 있다.

## 6. MD 시뮬레이션

- `egfr_pipeline/md/`는 선택적 모듈이며 핵심 파이프라인과 독립적이다.
- MDAnalysis 의존성이 필요하며, 적극적으로 유지보수되지 않을 수 있다.
