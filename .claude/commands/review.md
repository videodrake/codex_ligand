# /review — 시니어 코드 리뷰어

그룹 구현 완료 후 실행합니다. "PRD에 맞는가?"를 판정합니다.

## 프로세스
1. `docs/prd.md`에서 이 그룹의 Feature AC를 추출, 각 AC에 PASS/PARTIAL/FAIL 판정
2. Edge Case 커버리지 확인
3. 린트, 타입 체크, 기존 테스트를 돌린다. 오류가 있으면 수정한다
4. 셀프 체크: 에러 처리 누락, 엣지 케이스 미처리, 컨벤션 위반
5. Definition of Done 확인: `docs/tasks.md`에 명시된 그룹 DoD를 충족하는가?
6. 보고:

---
Group N 리뷰
PRD AC: {PASS}/{전체}, Edge Cases: {처리}/{전체}
DoD: {PASS | FAIL}
판정: {PASS | PARTIAL | FAIL}
수정 필요: {있으면 구체적 항목과 관련 AC 번호}
---

## 판정 기준
- **PASS:** 모든 Must-Have AC가 PASS이고, Edge Cases 80% 이상 처리, DoD 충족
- **PARTIAL:** Must-Have AC 중 일부가 PARTIAL. 수정 범위가 작아서 /review 내에서 직접 수정 가능하면 수정 후 재판정.
- **FAIL:** Must-Have AC 중 하나라도 FAIL이거나, DoD 미충족

PASS → /test.
PARTIAL (수정 후 PASS) → /test.
PARTIAL (수정 범위 큼) → /execute로 수정.
FAIL → /execute로 수정.

$ARGUMENTS