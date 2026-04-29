# M1 Phase 3 — Changes

Closes M1 §23 #13 (MYO1D 955-1006 construct QC) per handoff §16 and v1.0 plan §16 M1 Task 6 (MYO1D portion).

## Files created

```text
fresh/src/egfr_myo1d/myo1d/qc.py
fresh/tests/test_m1_phase3_myo1d_construct_qc.py
fresh/tests/fixtures/m1_phase3_myo1d/myo1d_955_1006_valid.pdb
fresh/tests/fixtures/m1_phase3_myo1d/myo1d_955_1001_short.pdb
fresh/tests/fixtures/m1_phase3_myo1d/myo1d_962_1006_terminal_bad.pdb
fresh/tests/fixtures/m1_phase3_myo1d/myo1d_with_ace_nme_caps.pdb
fresh/docs/m1_phase3_myo1d_construct_qc.md
fresh/docs/m1_phase3_changes.md
```

## Files modified

```text
fresh/src/egfr_myo1d/myo1d/construct.py    # added slice_myo1d_construct + emit_myo1d_construct_pdb
fresh/src/egfr_myo1d/cli.py                 # added prepare-myo1d subparser + _cmd_prepare_myo1d
```

## Files deleted

None.

## Public API additions

```python
# myo1d/construct.py (Phase 3)
slice_myo1d_construct(structure, start, end, include_caps=True) -> PDBStructure
emit_myo1d_construct_pdb(ctx, structure, output_path) -> dict

# myo1d/qc.py (new)
run_myo1d_qc(ctx, source_pdb, construct_range=None, profile="codex_dev") -> Myo1dQcReport
parse_construct_range(text_or_tuple) -> (int, int)
expand_residue_set(text) -> list[int]
load_myo1d_gate(ctx) -> dict
Myo1dQcReport (dataclass)
MYO1D_CONSTRUCT_QC_COLUMNS = [...] (15 columns)
```

## CLI surface additions

```bash
python -m egfr_myo1d.cli prepare-myo1d --run-id RUN --source PATH \
    [--construct 955-1006] [--profile codex_dev|hpc_strict] [--mode ...]
```

Total CLI subparsers after Phase 3: 13 (was 12 after Phase 1).

## Acceptance closure

- M1 §23 #13 closed: MYO1D 955-1006 construct QC works on synthetic fixture; output emitted under `runs/<run_id>/normalized/myo1d/`; QC CSV with the 15 columns from §16.3; manifest records sha256; key_residue_bonus_weight read from `gates.yaml` and asserted `0.0`; ACE/NME caps preserved; 962-start terminal artifact warns in codex_dev / fails in hpc_strict.

## Verification

- 27 new Phase 3 tests pass.
- Total suite: 141 passing (98 prior + 16 Phase 1 + 27 Phase 3).
- Old workflow files unmodified.
- All outputs land under `fresh/runs/<run_id>/`.

## Out of scope (next phases)

- Phase 4: receptor normalization (`prepare-receptor`)
- Phase 5: membrane frame generation (`compute-membrane-frame`)
- Phase 6: PBS generator
- Phase 7: ligand manifest shell
- Phase 8: `prepare-inputs` orchestrator + M1 integration test
- Phase 9: Tasks 4-9 schema realignment
