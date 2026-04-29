# M1 Phase 2 — MYO1D Module Relocation

Pure structural move aligning the package layout with `milestone1_foundation_codex_handoff_v0_5.md` §4. No logic changes; behavior is byte-identical to before the relocation. Phase 2 unblocks Phase 3 (MYO1D construct + QC extension).

## What moved

```text
fresh/src/egfr_myo1d/preparation/constructs.py  →  fresh/src/egfr_myo1d/myo1d/construct.py
fresh/src/egfr_myo1d/preparation/pdb_writer.py  →  fresh/src/egfr_myo1d/myo1d/pdb_writer.py
```

Both moves were performed via `git mv` so history follows. The two modules' content is unchanged except for one self-import inside `construct.py`:

```diff
-from egfr_myo1d.preparation.pdb_writer import select_atoms_by_residue_range, write_pdb_atoms
+from egfr_myo1d.myo1d.pdb_writer import select_atoms_by_residue_range, write_pdb_atoms
```

## What stays in `preparation/`

`preparation/{__init__.py, masks.py, restraints.py}` — these are EGFR-side restraint and mask helpers, not MYO1D-side, so they remain in `preparation/` per the spec module mapping.

## Import sites updated

```text
fresh/src/egfr_myo1d/validation/prepared_inputs.py    (3 imports)
fresh/src/egfr_myo1d/validation/real_inputs.py        (1 import)
fresh/tests/test_task4_ppi_input_preparation.py       (1 import — test fixture)
```

After the relocation, a recursive grep for `preparation.constructs` and `preparation.pdb_writer` in `fresh/src/` and `fresh/tests/` returns zero `import` statements (only the historical mention in `myo1d/__init__.py` docstring remains).

## Verification

- All 98 prior tests + 16 Phase 1 tests = 114 tests pass after the move.
- `python -c "from egfr_myo1d.myo1d import construct, pdb_writer"` succeeds.
- The two CLI commands that depend on these modules (`prepare-ppi-inputs`, `validate-real-inputs`) continue to produce identical outputs against the existing fixture inputs.

## Module tree progression

After Phase 2:

```text
fresh/src/egfr_myo1d/
├── analysis/                  (Tasks 7-9)
├── core/                      (Phase 1 added cleanup.py)
├── io/
├── myo1d/                     ← NEW (Phase 2)
│   ├── __init__.py
│   ├── construct.py           ← moved from preparation/
│   └── pdb_writer.py          ← moved from preparation/
├── planning/                  (Tasks 6)
├── preparation/               (EGFR-side helpers only now)
│   ├── __init__.py
│   ├── masks.py
│   └── restraints.py
├── structure/                 (Tasks 3)
└── validation/                (Tasks 3-9 + preflight)
```

Phase 3 will extend `myo1d/construct.py` with M1-spec functions (`slice_myo1d_construct`, `emit_myo1d_construct_pdb`) and add a new `myo1d/qc.py`.
