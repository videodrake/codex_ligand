name: bug-history
description: PyRosetta/Vina 코드 수정 전 반드시 확인. 트리거 — PyRosetta 코드 수정, 도킹 엔진 변경, 스코어 계산 변경 시. 비트리거 — 문서만 수정할 때, PBS 스크립트만 바꿀 때, 테스트만 추가할 때.

## 역사적 버그 — 반복하면 안 되는 실수들

### BUG-001: DockingSlideIntoContact 누락 (V1.0, 심각도: Critical)
- 증상: 모든 dG 값이 0.0
- 원인: DockMCMProtocol 호출 전 DockingSlideIntoContact 누락
- 교훈: PyRosetta 도킹 프로토콜 수정 시 반드시 dG 값이 비-제로인지 확인
- 확장: SlideIntoContact와 DockMCMProtocol의 순서 변경, 제거, 조건부 스킵은
  모두 BUG-001과 동일한 실패(dG=0.0)를 유발할 수 있다.
  이 두 함수의 호출 순서를 변경하는 요청은 반드시 BUG-001을 인용하며 경고한다.

### BUG-002: FoldTree 미재설정 (심각도: High)
- 증상: 역직렬화된 pose에서 도킹 결과가 비정상
- 원인: setup_foldtree(pose, "A_B", movable_jumps) 호출 누락
- 교훈: pose 로드/역직렬화 후 반드시 FoldTree 재설정

### BUG-003: stdout/stderr 리다이렉트 깨짐 (심각도: Medium)
- 증상: PyRosetta 배너가 출력에 섞임
- 원인: pyrosetta_init.py의 리다이렉트 구조 변경
- 교훈: pyrosetta_init.py 수정 시 배너 억제 동작 확인

### BUG-004: Beta-meander 방향 미검증 (심각도: High)
- 증상: back face가 receptor를 향한 모델이 consensus에 포함
- 원인: orientation_filter 미적용 또는 통과 기준 오류
- 교훈: PPI 결과 처리 시 항상 orientation pass/fail 상태 확인

## 코드 수정 전 체크리스트
□ 이 파일의 관련 버그를 확인했는가?
□ 수정 후 dG 값이 합리적인가? (0.0이 아닌가?)
□ FoldTree가 올바르게 설정되어 있는가?
□ orientation_filter가 정상 동작하는가?

## 상세 참조
- 설계의도.md — PyRosetta 절대 주의사항 + 설계 판단 근거
