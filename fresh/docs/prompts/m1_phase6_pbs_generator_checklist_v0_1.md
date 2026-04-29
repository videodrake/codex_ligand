# M1 Phase 6 Acceptance Checklist v0.1 — PBS Generator

Use this after the implementer applies M1 Phase 6.

## 1. Pre-Phase state preserved

```text
Old workflow files unchanged.
Phases 1-5 outputs/modules unchanged.
fresh/scripts/cleanup_run.py from Phase 1 unchanged.
```

## 2. New module importable

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -c "from egfr_myo1d.hpc.pbs import generate_pbs, render_pbs_content; print('OK')"
```

## 3. CLI registered

```bash
python -m egfr_myo1d.cli --help | grep prepare-pbs
python -m egfr_myo1d.cli prepare-pbs --help
```

Help text must include `--run-id`, `--job-name`, `--mode`, `--node`, `--ppn`, `--walltime`, `--output-path`, `--profile`.

## 4. Smoke-env PBS generation

```bash
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase6_local
python -m egfr_myo1d.cli prepare-pbs --run-id m1_phase6_local --job-name smoke_env_test --mode smoke_env --node node04
```

Expected output file: `fresh/runs/m1_phase6_local/scripts/smoke_env_test.pbs`

Inspect:

```bash
cat fresh/runs/m1_phase6_local/scripts/smoke_env_test.pbs
```

Expected lines:

```text
#!/bin/bash
#PBS -q workq
#PBS -N smoke_env_test
#PBS -l nodes=node04:ppn=4
#PBS -l walltime=02:00:00
#PBS -V
#PBS -o <ABS>/fresh/runs/m1_phase6_local/logs/jobs/smoke_env_test.stdout
#PBS -e <ABS>/fresh/runs/m1_phase6_local/logs/jobs/smoke_env_test.stderr

set -euo pipefail
cd <ABS>

source /usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh
conda activate pyrosetta

export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONIOENCODING=utf-8
export REPO_ROOT="<ABS>"
export PYTHONPATH="$REPO_ROOT/fresh/src:${PYTHONPATH:-}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

python -m egfr_myo1d.cli preflight --run-id m1_phase6_local --mode smoke_env --profile hpc_strict
python -m egfr_myo1d.cli cleanup --run-id m1_phase6_local --mode test --dry-run false
```

`<ABS>` should be an absolute path to repo root.

## 5. Production PBS uses ppn=32

```bash
python -m egfr_myo1d.cli prepare-pbs --run-id m1_phase6_local --job-name prod_test --mode production --node node05
grep "ppn=32" fresh/runs/m1_phase6_local/scripts/prod_test.pbs
grep "walltime=999:00:00" fresh/runs/m1_phase6_local/scripts/prod_test.pbs
```

Expected: matches.

## 6. Mini and scaling modes

```bash
python -m egfr_myo1d.cli prepare-pbs --run-id m1_phase6_local --job-name mini_test --mode mini --node node04
grep "ppn=16" fresh/runs/m1_phase6_local/scripts/mini_test.pbs
grep "walltime=12:00:00" fresh/runs/m1_phase6_local/scripts/mini_test.pbs

python -m egfr_myo1d.cli prepare-pbs --run-id m1_phase6_local --job-name scaling_test --mode scaling --node node06
grep "ppn=32" fresh/runs/m1_phase6_local/scripts/scaling_test.pbs
grep "walltime=02:00:00" fresh/runs/m1_phase6_local/scripts/scaling_test.pbs
```

## 7. Unknown node FAIL

```bash
python -m egfr_myo1d.cli prepare-pbs --run-id m1_phase6_local --job-name bad_node --mode smoke_env --node node99
```

Expected: nonzero exit, no .pbs file written for bad_node.

## 8. Path traversal

```bash
python -m egfr_myo1d.cli prepare-pbs --run-id ../bad_run --job-name x --mode smoke_env --node node04
```

Expected: nonzero exit, no outside writes.

## 9. submit_smoke_env.sh script

```bash
bash fresh/scripts/submit_smoke_env.sh m1_phase6_via_script
```

Expected:

```text
- prints PBS file path
- prints "qsub <path>" instruction
- DOES NOT actually call qsub
- file exists at fresh/runs/m1_phase6_via_script/scripts/<jobname>.pbs
```

## 10. submit_smoke_input.sh script

```bash
bash fresh/scripts/submit_smoke_input.sh m1_phase6_input_via_script
```

Same pattern as #9 but mode=smoke_input. PBS body should include `prepare-inputs` invocation (Phase 8 will fully validate this).

## 11. qsub never called automatically

```bash
grep -rn "^[[:space:]]*qsub" fresh/scripts/ fresh/src/egfr_myo1d/hpc/
```

Expected: empty (no auto-qsub). Comments and echo statements are fine.

```bash
grep -rn "qsub" fresh/scripts/ | grep -vE "echo|#|//"
```

Expected: empty (no executable qsub call).

## 12. Golden-file diffs

```bash
python -c "
import sys; sys.path.insert(0, 'fresh/src')
from egfr_myo1d.hpc.pbs import render_pbs_content
content = render_pbs_content(
    job_name='smoke_env_test',
    node='node04',
    ppn=4,
    walltime='02:00:00',
    queue='workq',
    command_lines=['python -m egfr_myo1d.cli preflight --run-id RUNID --mode smoke_env --profile hpc_strict', 'python -m egfr_myo1d.cli cleanup --run-id RUNID --mode test --dry-run false'],
    repo_root='/test/repo',
    conda_sh='/usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh',
    conda_env='pyrosetta',
    thread_limits={'OMP_NUM_THREADS':1,'OPENBLAS_NUM_THREADS':1,'MKL_NUM_THREADS':1,'NUMEXPR_NUM_THREADS':1,'VECLIB_MAXIMUM_THREADS':1},
    stdout_path='/test/repo/fresh/runs/RUNID/logs/jobs/smoke_env_test.stdout',
    stderr_path='/test/repo/fresh/runs/RUNID/logs/jobs/smoke_env_test.stderr',
    run_id='RUNID',
)
print(content)
" > /tmp/rendered_pbs.txt 2>&1 || true
diff /tmp/rendered_pbs.txt fresh/tests/fixtures/m1_phase6_pbs/golden_smoke_env_node04_ppn4.pbs
```

Expected: empty diff (or only trailing newline differences).

## 13. Tests

```bash
pytest -q fresh/tests/test_m1_phase6_pbs_generator.py
pytest -q fresh/tests
```

## 14. Old workflow protection

```bash
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

Empty.

## 15. What must not be in this phase

```text
- ligand work (Phase 7)
- prepare-inputs orchestrator (Phase 8 — Phase 6's smoke-input PBS may include the placeholder command line)
- Tasks 4-9 schema realignment (Phase 9)
- actual qsub submission
- modifying old workflow files
```

## 16. Phase 6 accepted if

```text
- hpc/__init__.py and hpc/pbs.py created.
- prepare-pbs CLI subcommand registered.
- generate_pbs.py / submit_smoke_env.sh / submit_smoke_input.sh real implementations replace placeholders.
- PBS files have absolute stdout/stderr, PYTHONPATH=fresh/src, 5 thread limits=1, conda activate pyrosetta.
- Mode → ppn / walltime mapping correct per hpc.yaml.
- Unknown node fails cleanly.
- qsub never auto-called.
- ≥10 phase tests + golden-file diffs pass; existing tests pass.
- M1 §23 #6 closed (file generation; qsub run is HPC-pending user step).
- Old workflow files unmodified.
```

## 17. Implementer final response must include

```text
M1 Phase 6 status: PASS / PASS WITH WARNINGS / FAIL
Files created
Files modified (cli.py + 3 scripts + module)
Commands run and results
Test summary
PBS structural validation evidence (per mode/node/ppn)
qsub auto-call audit: empty
HPC-pending note: user runs `bash fresh/scripts/submit_smoke_env.sh` on HPC for actual qsub
Acceptance closure: M1 §23 #6 closed
Old workflow protection
Known limitations
```
