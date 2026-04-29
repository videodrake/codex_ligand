# M1 Phase 4 — Changes

Closes M1 §23 #10 (explicit A/B normalization) and #11 (+1000 runtime offset) per handoff §14 and v1.0 plan §16 M1 Task 5 (receptor portion).

## Files created

```text
fresh/src/egfr_myo1d/model/__init__.py
fresh/src/egfr_myo1d/model/receptor_normalize.py
fresh/src/egfr_myo1d/model/receptor_qc.py
fresh/src/egfr_myo1d/io/residue_mapping.py
fresh/tests/test_m1_phase4_receptor_normalization.py
fresh/tests/fixtures/m1_phase4_receptor/explicit_AB_dimer.pdb         (copy of mini_explicit_AB.pdb)
fresh/tests/fixtures/m1_phase4_receptor/duplicate_chain_X_dimer.pdb   (copy of mini_duplicate_chain_X.pdb)
fresh/tests/fixtures/m1_phase4_receptor/single_chain_monomer.pdb      (new)
fresh/tests/fixtures/m1_phase4_receptor/v924r_warn.pdb                (copy of task3_inputs/egfr_v924r_warn.pdb)
fresh/tests/fixtures/m1_phase4_receptor/dimer_with_TM_excluded_range.pdb (new; 634-1100 + CHL/HOH HETATM)
fresh/docs/m1_phase4_receptor_normalization.md
fresh/docs/m1_phase4_changes.md
```

## Files modified

```text
fresh/src/egfr_myo1d/cli.py    # added prepare-receptor subparser + _cmd_prepare_receptor handler
```

## Files deleted

None.

## Public API additions

```python
# model/receptor_normalize.py
normalize_receptor(ctx, source_pdb, state_id, profile, strict, ...) -> NormalizedReceptor
NormalizedReceptor (dataclass)
RECEPTOR_AUDIT_CSV_COLUMNS = [...]   # 15 columns
PRIMARY_STATES = ("EGFR_160-185", "EGFR_170-200")
REFERENCE_CONTROL_STATES = ("3GT8_raw",)
resolve_role(state_id) -> str
load_receptor_gate(ctx) -> dict

# model/receptor_qc.py
detect_normalization_case(structure) -> str
split_duplicate_chain(atoms) -> (a_atoms, b_atoms)
detect_warn_mutations(atoms, warn_mutations) -> list
compute_residue_audit_rows(...) -> list
is_receptor_atom(atom) -> bool
NON_RECEPTOR_HETATM = LIPID_RESNAMES | WATER_RESNAMES | ION_RESNAMES
CAP_RESNAMES = {"ACE", "NME"}
DEFAULT_WARN_MUTATIONS = (WarnMutation(924, "VAL", "ARG", "3GT8_V924R"),)

# io/residue_mapping.py
write_residue_mapping(path, rows, ctx=None) -> None
read_residue_mapping(path) -> list[dict]
MAPPING_CSV_COLUMNS = [...]   # 11 spec columns
```

## CLI surface additions

```bash
python -m egfr_myo1d.cli prepare-receptor --run-id RUN \
    --state EGFR_160-185|EGFR_170-200|3GT8_raw \
    --source PATH \
    [--profile codex_dev|hpc_strict] [--mode ...] [--strict]
```

Total CLI subparsers after Phase 4: 14 (was 13 after Phase 3).

## Acceptance closure

- M1 §23 #10 closed: explicit A/B chain normalization works on both Case A (passthrough) and Case B (duplicate-chain split). Mapping CSV round-trips identity.
- M1 §23 #11 closed: +1000 runtime offset applied only to protomer B in `runtime_offset_receptor_only.pdb`. Protomer A residue numbers unchanged. PDB column 22-25 rewritten to reflect the new residue number.

## Verification

- 26 new Phase 4 tests pass (7 helpers + 6 case A/B/C + 2 dockable crop + 2 runtime offset + 2 mapping + 1 3GT8_raw + 1 V924R + 4 outputs/manifest/log + 3 CLI).
- Total suite: 167 passing (98 prior + 16 Phase 1 + 27 Phase 3 + 26 Phase 4).
- Old workflow files unmodified.
- All outputs land under `fresh/runs/<run_id>/`.

## Out of scope (next phases)

- Phase 5: membrane frame generation
- Phase 6: PBS generator
- Phase 7: ligand manifest shell
- Phase 8: `prepare-inputs` orchestrator + M1 integration test
- Phase 9: Tasks 4-9 schema realignment
