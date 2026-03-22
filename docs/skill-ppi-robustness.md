# Skill: PPI Branch 강건성 검증

## 적용 범위
- Group: 3
- Features: F-3

## 도메인 지식

### Orientation Filter 메커니즘
- MYO1D beta-meander(두께 5-7Å)가 EGFR 표면을 향해야 유효한 접촉
- PCA normal + Cα→Cβ multi-probe consensus + dot product 판정
- 현재 threshold: |dot product| ≥ 0.15 → pass, < 0.15 → ambiguous (consensus에서 제외)
- sweep 대상: [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
- sheet 8/9 잔기(961-964, 968-972 — 965-967은 sheet 간 gap으로 포함하지 않음)가 모든 threshold에서 hotspot에 포함되면 "핵심 잔기 안정적"
- ambiguous 비율 < 10% → 현재 유지 / > 20% → 재검토

### Fragment 범위 Sensitivity
- 현재: 955-1006 (pilot에서 N-terminal charge artifact 해결 위해 조정)
- 범위 A: 945-1006 (N-terminal +10)
- 범위 B: 955-1006 (현재)
- 범위 C: 955-1015 (C-terminal +9)
- Pilot: 1 state(3GT8_raw) × 1 seed × 2,000 모델
- 판정: Jaccard > 0.7 → 범위 robust / < 0.5 → 재선정 / 0.5-0.7 → 부분 영향

### Sheet 12 Working Assumption
- 현재 가정: sheet 12(998-1004)는 structural support, direct contact 아님
- active face에서 제외하고 모니터링만
- sensitivity: 설정 A(8+9) vs 설정 B(8+9+12), 기존 모델 좌표로 orientation 재계산
- 겹침 ≥ 80% → 유지 / 50-80% → 검토 / < 50% → 병행 보고

### 핸드오프 검증 패턴
- `_validate_adv_handoff()`: 파일 존재 → 스키마 → (신규) 데이터 품질
- Phase 1: hotspot 0개 → FAIL / 전부 state-specific → WARNING
- Phase 2: 포켓 0개 → FAIL / 전부 irrelevant → WARNING / docking skip 전부 → FAIL
- Phase 3: 유효 포즈 < 5 → FAIL / 리간드 1종만 → WARNING
- 삽입 위치: 기존 스키마 체크 직후

### PPI 패치 과대추정
- orthosteric + rim > 전체 포켓의 80% → WARNING
- hotspot 잔기 > C-lobe 전체의 30% → WARNING
- 원인: occupancy threshold 너무 낮음, orientation-fail 혼입, 실제 넓은 interface

### Cross-method Agreement
- PyRosetta(ref2015) vs LightDock(DFIRE2): 같은 입력 구조 → 부분적 독립
- concordance_score = min(pyro_occ, light_occ) / max(pyro_occ, light_occ)
- "both" 내 세분화: > 0.5 → "strong_both" / < 0.5 → "weak_both"

## 주의 사항
- Group 3은 🔴(Red zone): 서버 pilot, 깊은 기존 코드 수정, 되돌리기 어려운 파라미터 결정 포함
- C-3.3(pilot 실행)은 서버 qsub 작업 — 에이전트가 직접 실행하지 않고 스크립트만 준비
- _validate_adv_handoff 수정 시 FAIL이 파이프라인을 완전히 중단시키므로 정상 케이스 PASS 테스트 필수
- conformational_selection_candidate 태그는 Phase 2 로직을 변경하지 않음 — 태그만 추가
