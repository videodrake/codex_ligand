# Test Suite Triage

현재 저장소의 테스트는 **marker 기반 경량 베이스라인**으로 운영합니다.

## 현재 테스트 레인

### 1) Pre-qsub blocking lane (기본 게이트)

- 대상: `tests/unit/`
- 마커: `unit and not slow`
- 목적: 경로 해석/호환 플래그/manifest·marker helper처럼
  최근 변경 위험이 높은 로직을 빠르게 검증

실행:

```bash
bash scripts/run_pre_qsub_checks.sh
```

또는 직접:

```bash
python -m pytest -m "unit and not slow" tests/unit -q
```

### 2) Reporting lane (선택적/비차단)

- 마커: `reporting`
- 목적: 사용자-facing 리포트/표현물 포맷 검증

실행:

```bash
python -m pytest -m reporting tests -q
```

## 마커 정책

- `unit`: 빠르고 결정적인 단위 테스트
- `integration`: 모듈 간 연동 검증 (파일시스템 wiring 포함)
- `e2e`: 시나리오 기반 end-to-end 검증
- `slow`: PR 기본 게이트에서 제외되는 장시간 테스트
- `reporting`: 보고서/표현물 검증

`pytest.ini`에 위 마커가 정의되어 있으며, 기본 pre-qsub는
`unit and not slow`만 실행합니다.

## 운영 원칙

1. PR 기본 게이트는 빠른 `unit` 중심으로 유지합니다.
2. 통합/시나리오 테스트(`integration`, `e2e`)는 별도 레인(예: nightly)에서 확장합니다.
3. 새 테스트 추가 시 파일명보다 **마커 계약**을 우선합니다.
