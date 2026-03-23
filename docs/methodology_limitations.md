# Methodology Limitations

> 이 파이프라인의 계산적 예측에 내재한 방법론적 한계를 기술한다.
> 결과 해석 시 반드시 참조해야 한다.

## 1. Rigid-Body Docking

**한계:** AutoDock Vina와 PyRosetta global docking은 수용체-리간드 모두를 강체(rigid body)로 취급한다. 단백질의 induced-fit 효과(리간드 결합 시 구조 변화)가 반영되지 않는다.

**완화:** 3가지 EGFR 구조 상태(3GT8_raw, EGFR_160-185, EGFR_170-200)를 사용하여 MD 클러스터 대표 구조 간 구조적 다양성을 부분적으로 보완한다. 그러나 이는 리간드-유도 구조 변화를 직접 모델링하지 않으며, MD 스냅샷이 커버하지 못하는 구조 상태는 맹점으로 남는다.

**영향:** 깊은 포켓(induced-fit이 중요한 경우)의 affinity가 과소평가될 수 있다. 반대로, 유연한 루프 근처의 결합 사이트가 과대평가될 수 있다.

## 2. LightDock 독립성

**한계:** LightDock은 PyRosetta와 "다른 방법"이지만, 동일한 입력 구조(3개 PDB)를 사용한다. 따라서 cross_method_convergence는 "method-diverse" 증거이지 "method-independent" 증거가 아니다.

**완화:** cross_method_convergence.csv의 `method_agreement` 필드가 "strong_both"인 잔기도 공유 맹점에 의한 허위 일치일 수 있음을 인지해야 한다. `concordance_score`가 높더라도 입력 구조 편향을 배제할 수 없다.

**영향:** PyRosetta-LightDock 일치가 실제 생물학적 일치보다 과대하게 보일 수 있다. 특히 구조적 아티팩트(crystal contact 등)가 두 방법에 공통으로 영향을 미치는 경우.

## 3. 입력 구조 공유

**한계:** Vina blind docking, PyRosetta PPI docking, LightDock 검증 모두 동일한 3개 EGFR PDB 구조를 사용한다. 이 구조들은 하나의 X-ray 결정 구조(3GT8)에서 유래한 MD 클러스터 대표체이다.

**완화:** 3개 상태 간 cross-state robustness 분석으로 특정 구조에 의존적인 결과를 식별할 수 있다. 그러나 3개 구조가 공유하는 체계적 편향(예: 결정화 조건에 의한 특정 루프 위치)은 감지할 수 없다.

**영향:** 모든 워크플로우에서 공통으로 누락되는 결합 사이트가 존재할 수 있다 (공통 맹점). "robust" 분류가 구조적 다양성이 아닌 입력 공유에 의한 것일 수 있다.

## 4. Solvent 효과

**한계:** Vina는 implicit solvation만 사용하고, PyRosetta는 REF2015 에너지 함수(Lazaridis-Karplus implicit solvation)를 사용한다. 물 분자를 통한 수소 결합(water-mediated H-bond)이 명시적으로 모델링되지 않는다.

**완화:** 없음. 이 파이프라인에서는 explicit water를 사용하는 단계가 없다.

**영향:** PPI 인터페이스의 극성 잔기 근처 결합이 과소평가될 수 있다. 소수성 포켓의 결합은 상대적으로 과대평가된다. 이는 특히 PPI proximity scoring에 영향을 미칠 수 있다 — 실제로는 물 분자를 통해 PPI에 영향을 미치는 사이트가 "distant"로 분류될 수 있다.

## 5. Vina Scoring 편향

**한계:** AutoDock Vina의 scoring function은 소수성 상호작용을 과대평가하고, 극성/정전기 상호작용을 과소평가하는 경향이 있다. EGFR C-lobe surface의 MYO1D 결합면은 극성 잔기가 풍부하여, 이 편향이 특히 문제가 된다.

**완화:**
- pose_region_classifier.py로 포즈 분포를 영역별로 분류하여 소수성 편향 시각화
- bootstrap 안정성 분석으로 단일 고득점 포즈에 의존하는 결과 식별
- 3종 리간드 cross-chemical consensus로 특정 리간드 편향 감소

**영향:** C-lobe surface 포켓의 affinity가 실제보다 낮게 보고될 수 있다 (극성 인터페이스 과소평가). ATP binding site 근처 소수성 포켓은 과대평가될 수 있으나, 이는 `is_atp_site` 플래그로 필터링된다.

---

## 종합 권장 사항

1. STRONG 판정 포켓도 반드시 PyMOL 등으로 시각적 검토 필요
2. Consensus(Workflow A+B 일치) 포켓을 최우선 실험 후보로 선정
3. 단일 메서드에만 의존하는 결과(pyrosetta_only, lightdock_only)는 낮은 신뢰도로 취급
4. Allosteric 후보는 MD 시뮬레이션으로 구조 변화 기반 메커니즘 검증 필요
5. 최종 결과는 반드시 실험적 검증(SPR, ITC, cellular assay 등)을 거쳐야 함
