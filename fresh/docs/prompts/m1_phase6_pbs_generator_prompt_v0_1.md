# Claude M1 Phase 6 Prompt — PBS Generator v0.1

Branch `claude/task10`. Phases 1-5 complete. This is **M1 Phase 6** — implements PBS job-file generation per `milestone1_foundation_codex_handoff_v0_5.md` §10, closing M1 §23 #6.

## 1. Project context

The HPC environment uses PBS/qsub on nodes node04, node05, node06 with conda env `pyrosetta` and Python 3.9.25. Currently `fresh/scripts/{generate_pbs.py, submit_smoke_env.sh, submit_smoke_input.sh}` are placeholders that print preview text without writing real PBS files.

Phase 6 implements concrete PBS file generation. Codex env has no `qsub`, so this phase **generates and validates files**, but does not run qsub. HPC validation is a separate user-side step.

## 2. Absolute rules

Do not modify the old workflow. Maintain Py2/3 syntax compatibility.

Source-of-truth values from `fresh/configs/hpc.yaml`:

```yaml
scheduler: PBS
queue: workq
nodes: [node04, node05, node06]
conda:
  env_name: pyrosetta
  conda_sh: /usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh
python:
  expected_path: /home/eunae/.conda/envs/pyrosetta/bin/python
ppn:
  smoke: 4
  mini: 16
  scaling: 32
  production: 32
walltime:
  smoke: "02:00:00"
  mini: "12:00:00"
  scaling: "02:00:00"
  production: "999:00:00"
thread_limits:
  OMP_NUM_THREADS: 1
  OPENBLAS_NUM_THREADS: 1
  MKL_NUM_THREADS: 1
  NUMEXPR_NUM_THREADS: 1
  VECLIB_MAXIMUM_THREADS: 1
```

No hardcoding inside the module.

## 3. Scope

In scope:
- Create `fresh/src/egfr_myo1d/hpc/__init__.py`
- Create `fresh/src/egfr_myo1d/hpc/pbs.py`
- Add `prepare-pbs` CLI subcommand
- Replace placeholder `fresh/scripts/generate_pbs.py` with real entry point
- Replace placeholder `fresh/scripts/submit_smoke_env.sh` and `submit_smoke_input.sh` with real PBS-emit + qsub-ready scripts (do NOT auto-call qsub)
- Tests under `fresh/tests/test_m1_phase6_pbs_generator.py` (≥10 tests)
- Fixtures under `fresh/tests/fixtures/m1_phase6_pbs/` — golden PBS files for diff comparison
- Docs `fresh/docs/m1_phase6_pbs_generator.md` and `m1_phase6_changes.md`

Out of scope:
- Calling `qsub` (HPC user-side step)
- Ligand work (Phase 7)
- Integration orchestrator (Phase 8)
- M2/M3 work
- Modifying old workflow files

## 4. Required CLI behavior

```bash
python -m egfr_myo1d.cli prepare-pbs \
  --run-id RUN \
  --job-name JOB_NAME \
  --mode smoke_env|smoke_input|mini|scaling|production \
  [--node node04|node05|node06] \
  [--ppn N]                         # override hpc.yaml default for the mode
  [--walltime HH:MM:SS]              # override hpc.yaml default
  [--output-path PATH]               # default: runs/<run_id>/scripts/<job_name>.pbs
  [--profile codex_dev|hpc_strict]
```

The script `fresh/scripts/submit_smoke_env.sh` calls this CLI to emit the smoke-env PBS, then prints the `qsub` command for the user to run manually:

```bash
#!/usr/bin/env bash
set -euo pipefail
RUN_ID="${1:-smoke_env_$(date +%Y%m%d_%H%M%S)}"
export PYTHONPATH="$(pwd)/fresh/src:${PYTHONPATH:-}"
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id "$RUN_ID"
PBS_FILE=$(python -m egfr_myo1d.cli prepare-pbs --run-id "$RUN_ID" --job-name "smoke_env_${RUN_ID}" --mode smoke_env --node node04 | tail -1)
echo "PBS file generated: $PBS_FILE"
echo "To submit on HPC, run:"
echo "  qsub $PBS_FILE"
```

`submit_smoke_input.sh` mirrors but uses `--mode smoke_input` and includes a `prepare-inputs` invocation in the PBS body (Phase 8 will provide that command).

## 5. Files to create / modify

Create:

```text
fresh/src/egfr_myo1d/hpc/__init__.py
fresh/src/egfr_myo1d/hpc/pbs.py
fresh/tests/test_m1_phase6_pbs_generator.py
fresh/tests/fixtures/m1_phase6_pbs/golden_smoke_env_node04_ppn4.pbs
fresh/tests/fixtures/m1_phase6_pbs/golden_smoke_input_node04_ppn4.pbs
fresh/tests/fixtures/m1_phase6_pbs/golden_production_node05_ppn32.pbs
fresh/docs/m1_phase6_pbs_generator.md
fresh/docs/m1_phase6_changes.md
```

Modify:

```text
fresh/src/egfr_myo1d/cli.py             # add prepare-pbs subparser + handler
fresh/scripts/generate_pbs.py           # replace placeholder with thin wrapper around CLI
fresh/scripts/submit_smoke_env.sh       # real qsub-ready emitter (does NOT call qsub)
fresh/scripts/submit_smoke_input.sh     # real qsub-ready emitter (does NOT call qsub)
```

## 6. Public API

`hpc/pbs.py`:

```python
def generate_pbs(ctx, job_name, node, ppn, walltime, command_lines, conda_env, mode, output_path):
    # type: (RunContext, str, str, int, str, list[str], str, str, Path) -> PBSFile
    """
    Render PBS file content, write to output_path under ctx.run_dir (or specified location).
    Refuses to write outside ctx.run_dir unless explicit absolute path supplied (in which case
    the path must be checked or warned).
    Returns PBSFile with path, sha256, content_preview.
    """

def render_pbs_content(job_name, node, ppn, walltime, queue, command_lines, repo_root, conda_sh, conda_env, thread_limits, stdout_path, stderr_path, run_id):
    # type: (...) -> str
    """Pure function returning PBS file text. Used by tests for golden-file diff."""

def derive_smoke_env_command_lines(run_id):
    # returns list of shell command lines for smoke-env PBS body

def derive_smoke_input_command_lines(run_id, input_root):
    # returns list of shell command lines for smoke-input PBS body
```

## 7. PBS file requirements (handoff §10.1)

Required content:

```bash
#!/bin/bash
#PBS -q workq
#PBS -N <job_name>
#PBS -l nodes=<node>:ppn=<ppn>
#PBS -l walltime=<walltime>
#PBS -V
#PBS -o <ABSOLUTE_REPO_ROOT>/fresh/runs/<run_id>/logs/jobs/<job_name>.stdout
#PBS -e <ABSOLUTE_REPO_ROOT>/fresh/runs/<run_id>/logs/jobs/<job_name>.stderr

set -euo pipefail
cd <ABSOLUTE_REPO_ROOT>

source /usr/local/anaconda/3/2023.09/etc/profile.d/conda.sh
conda activate pyrosetta

export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONIOENCODING=utf-8
export REPO_ROOT="<ABSOLUTE_REPO_ROOT>"
export PYTHONPATH="$REPO_ROOT/fresh/src:${PYTHONPATH:-}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

<command_lines from derive_*_command_lines>
```

Strict requirements:

```text
- stdout/stderr paths must be absolute (no PBS variable expansion like $PBS_JOBNAME)
- PYTHONPATH must include $REPO_ROOT/fresh/src
- Five thread limits must all be set to 1
- conda activate of pyrosetta env must be present
- queue must be workq (or hpc.yaml override)
- ppn matches mode mapping from hpc.yaml (smoke=4, mini=16, scaling=32, production=32)
- walltime matches mode mapping from hpc.yaml
```

Smoke-env command body (handoff §10.2 Smoke A):

```bash
python -m egfr_myo1d.cli preflight --run-id <run_id> --mode smoke_env --profile hpc_strict
python -m egfr_myo1d.cli cleanup --run-id <run_id> --mode test --dry-run false
```

Smoke-input command body (handoff §10.2 Smoke B):

```bash
python -m egfr_myo1d.cli prepare-inputs --run-id <run_id> --mode smoke_input --profile hpc_strict --input-root fresh/data/raw
python -m egfr_myo1d.cli cleanup --run-id <run_id> --mode test --dry-run false
```

(`prepare-inputs` is implemented in Phase 8. Phase 6's smoke-input PBS may include the placeholder command; Phase 8 verifies it works end-to-end.)

## 8. Behavior policy

```text
- Read all parameters from fresh/configs/hpc.yaml; allow CLI override.
- Validate node ∈ {node04, node05, node06}; reject others with FAIL.
- Validate ppn matches mode default unless explicitly overridden.
- Resolve ABSOLUTE_REPO_ROOT from RunContext (ctx.repo_root).
- Default output_path: fresh/runs/<run_id>/scripts/<job_name>.pbs (under run_dir).
- Refuse output_path outside run_dir unless --output-path is given as absolute path AND --profile=hpc_strict (allow strict-only override).
- Append phase status PASS / WARN / FAIL.
- Do NOT call qsub. Do NOT auto-submit.
- Append job_status.jsonl entry with status=GENERATED (not SUBMITTED) noting the PBS file path.
```

## 9. Severity rules

```text
PASS:  PBS file generated, golden-file equivalent for the mode/node/ppn combo, written under run_dir
WARN:  user-supplied --ppn or --walltime differs from hpc.yaml default (recorded but accepted)
FAIL:  unknown node, malformed walltime, output_path outside run_dir without strict override, write failure
```

## 10. Tests required (≥10)

```text
test_pbs_uses_concrete_stdout_stderr_absolute_paths
test_pbs_exports_pythonpath_to_fresh_src
test_pbs_sets_all_five_thread_limits_to_1
test_pbs_activates_pyrosetta_conda_env
test_pbs_smoke_env_uses_node04_ppn4_walltime_02h
test_pbs_smoke_input_uses_node04_ppn4
test_pbs_production_uses_ppn32
test_pbs_command_line_includes_run_id
test_pbs_does_not_call_qsub_automatically
test_pbs_walltime_matches_hpc_yaml_per_mode
test_pbs_unknown_node_fails
test_pbs_golden_file_diff_smoke_env
test_pbs_golden_file_diff_smoke_input
test_pbs_golden_file_diff_production
test_cli_help_includes_prepare_pbs
```

(15 tests; ≥10 required.)

Golden-file tests: render PBS content with fixed inputs (run_id="m1_phase6_test", ABSOLUTE_REPO_ROOT="/test/repo", etc.) and diff against committed golden file.

## 11. Acceptance commands

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"

pytest -q fresh/tests/test_m1_phase6_pbs_generator.py
pytest -q fresh/tests

python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase6_local

# Smoke-env PBS
python -m egfr_myo1d.cli prepare-pbs --run-id m1_phase6_local --job-name smoke_env_test --mode smoke_env --node node04
test -f fresh/runs/m1_phase6_local/scripts/smoke_env_test.pbs && echo OK

# Production PBS
python -m egfr_myo1d.cli prepare-pbs --run-id m1_phase6_local --job-name prod_test --mode production --node node05
grep "ppn=32" fresh/runs/m1_phase6_local/scripts/prod_test.pbs

# Unknown node FAIL
python -m egfr_myo1d.cli prepare-pbs --run-id m1_phase6_local --job-name bad_node --mode smoke_env --node node99 || echo "Expected FAIL"

# scripts/submit_smoke_env.sh emits but does NOT qsub
bash fresh/scripts/submit_smoke_env.sh m1_phase6_local_via_script

# Path traversal
python -m egfr_myo1d.cli prepare-pbs --run-id ../bad_run --job-name x --mode smoke_env --node node04

# Old workflow protection
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/

# qsub never called automatically — verify by absence of any qsub invocation in scripts
grep -rn "qsub " fresh/scripts/ | grep -v "echo\|#"
```

## 12. Final response format

```text
M1 Phase 6 status: PASS / PASS WITH WARNINGS / FAIL
Files created (incl. golden fixtures)
Files modified (cli.py + 3 scripts/ + module)
Commands run and results
Test summary
PBS golden-file diff results (smoke_env, smoke_input, production)
qsub auto-call audit: confirmed NOT auto-submitted
PBS structural validation:
- absolute stdout/stderr paths
- PYTHONPATH=fresh/src
- 5 thread limits = 1
- conda activate pyrosetta
- correct ppn/walltime per mode from hpc.yaml
HPC-pending: actual qsub run (user-side step)
Acceptance closure: M1 §23 #6 closed (file generation)
Old workflow protection: empty diff
Known limitations:
- No actual qsub run in Codex env; user runs `bash fresh/scripts/submit_smoke_env.sh` on HPC
```
