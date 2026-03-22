# Skill: Verdict 메트릭 교정

## 적용 범위
- Group: 4
- Features: F-4

## 도메인 지식

### Verdict 3축 체계 (현재)
| 축 | 만점 | 구성 |
|----|------|------|
| Vina | 50 | affinity_pts + consensus_pts + stability_pts |
| PPI | 20 | spatial_pts (centroid 거리 기반) |
| Cross-Receptor | 30 | n_receptors 기반 (현재: 2+ → 30, 동점) |
| **총점** | **100** | STRONG ≥ 55, MODERATE ≥ 35, WEAK < 35 |

- 단일 축만으로 STRONG 도달 불가능 → 최소 2축 기여 필요 (설계 의도)
- PPI 데이터 없을 시 적응적 재배분: 60/0/40

### 가중치 시뮬레이션 조합
| 조합 | Vina | PPI | Cross | 의도 |
|------|------|-----|-------|------|
| 현재 | 50 | 20 | 30 | baseline |
| PPI 강화 | 40 | 30 | 30 | PPI가 프로젝트 질문의 절반 |
| PPI 최대 | 35 | 35 | 30 | Vina-PPI 동등 |
| Cross 약화 | 40 | 30 | 30 | PPI 강화와 동일값이나 Cross 축 해석 관점 차이 |
| Vina 약화 | 30 | 40 | 30 | PPI 우선 |
| 균등 | 33 | 33 | 34 | 편향 없는 참고 |

- 핵심 관찰 대상: PPI 가중치 증가로 새로 STRONG에 진입하는 포켓 (생물학적 의미 검토 필요)

### Centroid 의미론적 오프셋
- Vina centroid: 리간드 원자 좌표 평균 → 포켓 "내부"
- PPI centroid: 수용체 Cα 좌표 평균 → 표면
- 체계적 오프셋: 3-5Å (깊은 포켓일수록 큼)
- pocket_depth = centroid에서 가장 가까운 비접촉 표면 Cα까지 거리
- 보정: corrected_distance = raw_distance - alpha × pocket_depth (alpha 0.5~1.0)

### 현재 임계값 체계
- **vina_affinity_pts**: -8.0 → 20점 / -6.5 → 15점 / -5.0 → 10점
  - C-lobe surface(-5~-8 범위)에서 대부분 15점 구간 → 차별력 부족 가능성
- **ppi_spatial_pts**: 8Å → 20점 / 15Å → 15점 / 25Å → 10점
  - 오프셋 3-5Å 감안 시 실질 범위 3-22Å
- **cross_receptor_pts**: 현재 2/3와 3/3이 동일(30점) → 차등 도입 대상

### Phase 4 A3 축 (Perturbation Relevance)
- A1(PPI interface) + A2(Druggability) + A3(Perturbation) + A4(State robustness)
- A1+A3 합산 가중치 = 60% → MYO1D 교란이 핵심 목적
- A3의 구체적 계산 로직이 코드에서 추출 필요 (D-5.1 바이브코딩)
- 가능한 입력: centroid 거리, 잔기 overlap 비율, 면적 비율, 메커니즘 분류(ortho/rim/allosteric)

## 주의 사항
- 가중치 변경 시 프로덕션 재실행 필요: `python run_production.py --from 5`
- 시뮬레이션은 기존 결과를 "읽기만" 하여 what-if 분석 — 코드 수정은 확정 후에만
- pocket_depth 계산 시 접촉 잔기 자체는 제외해야 올바른 "깊이"가 나옴
- cross_receptor_pts 차등 도입 시 기존 STRONG 경계(55점)에 걸친 포켓 확인 필수
