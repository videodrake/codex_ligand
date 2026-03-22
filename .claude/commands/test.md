# /test — 테스터

/review PASS 후 실행합니다. "어떻게 깨뜨릴 수 있는가?"를 따집니다.

## 프로세스
1. 이 그룹 코드의 공격 표면을 파악한다: PRD Edge Cases, 실패 경로, 잘못된 입력
2. `docs/skill-testing.md`가 있으면 mock 전략과 테스트 환경 제약을 참조한다
3. 테스트를 작성한다
4. 전체 테스트를 돌려 회귀 버그를 확인한다
5. 결과를 확인한다:
   - 모두 통과 → 6단계로
   - 실패 → 코드 버그면 /execute, 테스트 로직 버그면 여기서 수정
6. `docs/CONTEXT.md`에 이번 그룹의 결정 사항을 기록한다
7. 커밋하고 보고한다:

---
Group N 완료
테스트: {passed}/{total}
신규 테스트: {새로 추가된 테스트 수}개
커밋: {type}(group-N): {한줄 요약}
---

## 커밋 컨벤션
- feat(group-N): 기능 구현 완료
- fix(group-N): 리뷰 지적 수정
- test(group-N): 테스트 추가

$ARGUMENTS