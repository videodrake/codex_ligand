# Precision Smoke Test (Production-like, Fast)

이 테스트는 `run_production.pbs`와 유사한 lane 경로를 사용하되,
`config/smoke-fast.yaml`로 입력/파라미터를 줄여 빠르게 검증합니다.

## What it does

`config/run_precision_smoke.pbs`(내부적으로 `scripts/run_precision_smoke.sh` 호출)는
아래 lane을 순서대로 실행합니다.

1. `vina-cpu`
2. `ppi` (`--state 3GT8_raw --seed 0`)
3. `ppi-post`
4. `vina-post`
5. `finalize`

## Important behavior

- 실행 중 `config/example-project.yaml`을 임시로 `config/smoke-fast.yaml`으로 교체합니다.
- 스크립트 종료 시(성공/실패 포함) 원래 `config/example-project.yaml`을 자동 복원합니다.

## Run (PBS 권장)

```bash
qsub config/run_precision_smoke.pbs
```

## Expected output location

- `output_smoke/workflow_a/...`
