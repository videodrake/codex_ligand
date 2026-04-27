# 투두리스트

프로젝트의 미완료 작업, 개선 과제, 향후 계획을 추적한다.
완료된 항목은 체크 표시 후 날짜를 기록한다.

## 인프라/코드 정리

- [x] input/PPI/prepared/ 레거시 경로 참조 정리 (2026-04-06)
- [x] 하네스 구축 완료: 스킬 7 + 에이전트 3 + 훅 2 (2026-04-06)
- [x] CLAUDE.md / CONTEXT.md 재작성 (2026-04-06)
- [x] 레거시 파일 삭제 및 참조 정리 (2026-04-06)
- [x] smoke test 환경 의존성 해결 — `residue_utils.py`, `prepare_dimer_pdb.py` 복원 + pytest venv 의존성 설치 (2026-04-06)

## 문서

- [x] docs/methodology_limitations.md 생성 (2026-04-06)
- [x] docs/workflow_comparison_guide.md 생성 (2026-04-06)
- [x] docs/phase4_A3_axis_specification.md 생성 (2026-04-06)
- [x] 투두리스트.md 생성 (2026-04-06)

## 파이프라인 개선

- [x] P2Rank 통합 — 파서+하류 모듈 모두 동작 확인. volume 추정 heuristic 추가, druggability 임계값 호환성 문서화 (2026-04-06)
- [ ] Orientation filter ambiguous band 임계값(+-0.15) 경험적 검증 — 검증 계획 문서화 완료 (`docs/phase1_notes.md` 4.3a). seed 5-9 완료 후 dot product 분포 분석으로 실행
- [x] Sheet 12 역할 규명 — 구조적 지지로 판정 (MD 무접촉 + 도킹 관찰 일치). seed 5-9 결과로 재확인 예정 (2026-04-06)
- [ ] Vina 실행 후 `pending_*` 리간드 지지 수준 재점수 — Phase 4 A3 축의 잠정 평가를 확정 평가로 전환

## HPC 실행

- [ ] Workflow A 전체 production run (3 수용체 x 3 리간드)
- [ ] Workflow B Phase 1 PPI 도킹 실행 (3 상태 x 5 seeds x 20K models)
- [ ] Phase 3 focused docking 실행 (Phase 2 포켓 후보 확정 후)

## 검증/품질

- [ ] Ko et al. sheet 8/9 활성면 잔기 검증 — Workflow B Phase 1 결과에서 3개 이상 확인
- [ ] Cross-method 검증 (PyRosetta-LightDock Jaccard) 결과 리뷰
- [ ] Phase 4 축 가중치(A1-A4) 첫 결과 리뷰 후 전문가 조정 (사람 승인 필수)
- [ ] Centroid 거리 편향(~3-5 A 과대 추정) 임계값 적절성 리뷰

## 장기 과제

- [ ] 추가 수용체 상태 확보 — 현재 실험 결정 구조 1개(3GT8_raw) + MD 클러스터 2개
- [ ] 리간드 라이브러리 확장 — 현재 3종 고정 (Tanimoto < 0.4), SAR 구축에 부족
- [ ] MD 모듈(`egfr_pipeline/md/`) 활성화 여부 결정 — MDAnalysis 의존성, 핵심 파이프라인과 독립
