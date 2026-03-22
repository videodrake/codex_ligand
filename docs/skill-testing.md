# Skill: 테스트 전략

## 적용 범위
- Group: 전체 (0-7)
- Features: 모든 Feature

## 테스트 환경 제약

### HPC/PBS 의존성
- 실제 Vina/PyRosetta 도킹은 테스트에서 실행하지 않음 (시간: 수 시간~일)
- fpocket, P2Rank, LightDock도 설치 불확실 → 외부 도구 호출은 mock
- PBS qsub 관련 코드는 스크립트 생성만 테스트, 실제 제출 안함

### 데이터 크기
- 프로덕션: 4,500 Vina 포즈, 300K PPI 모델 → 테스트에서는 축소 fixture 사용
- valid_sites.csv: 프로덕션 수십 행 → fixture 5-10행

## Mock 전략

### Mock 대상
| 대상 | Mock 방법 | 이유 |
|------|-----------|------|
| Vina 실행 | 사전 생성 pose_table.csv fixture | 도킹 불필요 |
| PyRosetta scoring | scored_all_models.csv fixture | PPI 도킹 불필요 |
| PDB 파일 I/O | 최소 PDB fixture (10-20 잔기) | 구조 파싱만 테스트 |
| RDKit 분자 로드 | 테스트용 SDF fixture (간단한 분자) | 실제 리간드 파일 의존 제거 |
| PBS qsub | subprocess mock | 서버 의존 제거 |

### Mock하지 않을 것
- pandas CSV 읽기/쓰기 — 실제 I/O 테스트 필요
- region_definitions.py 상수 — 확정된 값으로 실제 사용
- 수학 계산 (Tanimoto, Jaccard, 점수 계산) — 정확성 검증 대상

## Fixture 전략

### 핵심 Fixture 파일
1. **fixture_pocket_table.csv**: 5개 포켓 (ATP site 1개, C-lobe surface 3개, borderline 1개)
2. **fixture_valid_sites.csv**: 기존 컬럼 + 기대 결과 포함
3. **fixture_hotspot_residues.csv**: sheet 8/9 잔기 포함/미포함 변형
4. **fixture_scored_models.csv**: 10-20개 모델, orientation pass/fail/ambiguous 포함
5. **fixture_phase4_table.csv**: 워크플로우 비교용, 5개 포켓 (매칭/미매칭 포함)

### Fixture 원칙
- 각 fixture에 "기대 결과"를 주석으로 기록
- edge case용 변형 fixture: 빈 CSV, 0개 hotspot, 전부 irrelevant 등

## 테스트 패턴

### 1. 후방 호환 (모든 그룹 필수)
- 기존 파서로 신규 컬럼 포함 CSV 읽기 가능, 기존 컬럼 값 불변, 구 버전 CSV도 처리 가능

### 2. 검증 로직 (Group 1, 3)
- PASS/WARNING/FAIL 3경로 각각 테스트. FAIL 시 파이프라인 중단 확인, WARNING 시 로그+계속 확인

### 3. 점수 계산 (Group 2, 4)
- 알려진 입력 → 기대 점수 수동 계산 비교. 경계값(frac=0.5, affinity=-6.5 등). 총점 합계 = 100

### 4. 분류 로직 (Group 5)
- 4가지 분류 정확성. allosteric 기준 (Vina ≥ 35 AND PPI ≤ 5) 경계값

### 5. Setup / Infra / E2E (Group 0, 6, 7)
- Group 0: region_definitions import + 상수 검증 + 기존 테스트 통과
- Group 6: precheck 플래그, import 호환. Group 7: 프로덕션 회귀, before/after, AC 커버리지

## 테스트 실행
- 프레임워크: pytest, 실행: `pytest tests/ -v`
- CI 없음 — 수동 실행. 회귀 테스트: 각 그룹 완료 시 전체 스위트 실행
