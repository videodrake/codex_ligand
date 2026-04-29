# M1 Phase 6 — Changes

Closes M1 §23 #6 (qsub smoke can be generated; actual qsub run is HPC-pending) per handoff §10 and v1.0 plan §16 M1 Task 3.

## Files created

```text
fresh/src/egfr_myo1d/hpc/__init__.py
fresh/src/egfr_myo1d/hpc/pbs.py
fresh/tests/test_m1_phase6_pbs_generator.py
fresh/tests/fixtures/m1_phase6_pbs/golden_smoke_env_node04_ppn4.pbs
fresh/docs/m1_phase6_pbs_generator.md
fresh/docs/m1_phase6_changes.md
```

## Files modified

```text
fresh/src/egfr_myo1d/cli.py            # added prepare-pbs subparser + handler
fresh/scripts/generate_pbs.py          # placeholder -> real CLI wrapper
fresh/scripts/submit_smoke_env.sh      # placeholder -> real PBS-emit + qsub-instruction
fresh/scripts/submit_smoke_input.sh    # placeholder -> real PBS-emit + qsub-instruction
```

## Files deleted

None.

## Public API additions

```python
# hpc/pbs.py
THREAD_ENV_KEYS = ("OMP_NUM_THREADS", ...)   # 5 keys
DEFAULT_CONDA_SH, DEFAULT_CONDA_ENV
WALLTIME_PATTERN
PBSFile (dataclass)
load_hpc_gate(ctx) -> dict
render_pbs_content(...) -> str
derive_smoke_env_command_lines(run_id) -> list[str]
derive_smoke_input_command_lines(run_id, input_root) -> list[str]
derive_command_lines_for_mode(mode, run_id, input_root) -> list[str]
generate_pbs(ctx, job_name, mode, node=None, ppn=None, walltime=None,
             output_path=None, profile="codex_dev",
             input_root="fresh/data/raw", command_lines=None) -> PBSFile
```

## CLI surface additions

```bash
python -m egfr_myo1d.cli prepare-pbs --run-id RUN --job-name NAME --mode <mode> \
    [--node <node>] [--ppn N] [--walltime HH:MM:SS] [--output-path PATH] \
    [--input-root fresh/data/raw] [--profile codex_dev|hpc_strict]
```

Wrapper scripts:

```bash
bash fresh/scripts/submit_smoke_env.sh   [<run_id>] [<node>]
bash fresh/scripts/submit_smoke_input.sh [<run_id>] [<node>]
```

Both scripts emit a PBS file and print the qsub command for manual user submission.
Neither auto-calls qsub.

Total CLI subparsers after Phase 6: 16 (was 15 after Phase 5).

## Acceptance closure

- M1 §23 #6 closed (file generation): PBS files have absolute stdout/stderr paths,
  PYTHONPATH=fresh/src export, all 5 thread limits = 1, conda activate pyrosetta,
  mode-specific command body, Unix line endings.
- HPC-pending: actual `qsub <path>` execution and `qstat` confirmation are
  user-side validation steps. Annotated as `HPC_PENDING` in the M1 acceptance
  scorecard.

## Verification

- 23 new Phase 6 tests pass:
  - 7 pure render tests (paths, exports, threads, conda, ppn/walltime, command bodies)
  - 9 generate_pbs end-to-end tests (writes, modes, errors, line endings, traversal, status, override warning)
  - 1 anti-auto-qsub regression (module + scripts, with heredoc-aware parser)
  - 3 CLI tests
  - 3 hpc.yaml integration tests (smoke_env/smoke_input both map to "smoke" key)
- Total suite: 208 passing (98 prior + 16 P1 + 27 P3 + 26 P4 + 18 P5 + 23 P6).
- Old workflow files unmodified.

## Out of scope (next phases)

- Phase 7: ligand manifest shell
- Phase 8: prepare-inputs orchestrator + M1 integration test
- Phase 9: Tasks 4-9 schema realignment
- Actual qsub execution on HPC nodes (user-side)
- M2/M3 command bodies for mini/scaling/production modes (user-supplied later)
