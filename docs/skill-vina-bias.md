# Skill: Vina Blind Docking 편향 정량화

## 적용 범위
- Group: 2
- Features: F-2

## 도메인 지식

### EGFR 5영역 분류 (B-1.1 확정)
| Region | 잔기 수 | 용도 |
|--------|---------|------|
| n_lobe | 107 | Vina 포즈 분류 (관심 대상 아님) |
| atp_site | 37 | False positive 배제 (F-1과 연동) |
| c_lobe_surface | 136 | PPI-relevant 핵심 관심 영역 |
| c_lobe_core | 29 | 매장 잔기, 도킹 결과에서 무시 |
- 총 309잔기, SASA > 10.0 Å² 기준, 3-state union 전략
- `region_definitions.py`의 `get_region(residue_number)` 사용
- 포즈의 "주 영역" = 접촉 잔기의 과반이 속하는 영역

### 편향의 핵심 메커니즘
- ATP binding site는 깊고 잘 정의된 포켓 → Vina 포즈가 과도하게 집중
- C-lobe surface 포켓은 얕고 작음 → min_pocket_size=3에 미달하여 탈락 가능
- WARNING 기준: C-lobe surface 포즈 < 전체의 10%
- 이 편향은 Workflow B의 focused docking에서 보완됨

### Bootstrap-Verdict 연동 규칙
- `pocket_exists_frac`: bootstrap 200회 중 해당 포켓이 재현되는 비율
- frac ≥ 0.8 → stability_pts = 만점 (bootstrap_confidence = "high")
- 0.5 ≤ frac < 0.8 → stability_pts = 만점 × 0.5 ("medium")
- frac < 0.5 → stability_pts = 0 ("low")
- bootstrap 미실행 → 기존 로직 유지 ("not_assessed"), 후방 호환 필수

### C-lobe Surface Affinity 해석
- C-lobe surface 포켓의 현실적 affinity 범위: -5 ~ -8 kcal/mol
- ATP site 포켓(-10 ~ -12)보다 낮은 것이 정상
- -5 ~ -7도 표면 포켓으로서 의미 있는 결합을 나타낼 수 있음
- Vina scoring은 소수성 접촉을 과대평가 → 극성 풍부한 PPI 근처 포켓 과소평가

## 주의 사항
- validate.py 통합 시 기존 exit code 체계 준수: 0=통과, 1=경고, 2=실패
- 잔기 번호 비교는 "같은 번호에서 다른 아미노산" 탐지 — 전체 오프셋이 아닌 부분 불일치에 주목
- affinity 임계값 조정 결정은 계획안 D(Group 4) 전에 완료해야 함 (의존 관계)
- report.py 수정 시 기존 섹션 순서를 깨지 않도록 주의
