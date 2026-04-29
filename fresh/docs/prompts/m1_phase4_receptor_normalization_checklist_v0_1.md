# M1 Phase 4 Acceptance Checklist v0.1 — Receptor Normalization

Use this after the implementer applies M1 Phase 4.

## 1. Pre-Phase state preserved

```text
Old workflow files unchanged.
Phases 1-3 outputs unchanged.
Tasks 1-9 logic unchanged (Phase 4 does not modify validation/prepared_inputs.py output paths; Phase 9 will).
preparation/{masks.py, restraints.py} unchanged.
```

## 2. New modules importable

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -c "from egfr_myo1d.model.receptor_normalize import normalize_receptor; print('OK')"
python -c "from egfr_myo1d.model.receptor_qc import audit_receptor_residues; print('OK')"
python -c "from egfr_myo1d.io.residue_mapping import write_residue_mapping, read_residue_mapping, MAPPING_CSV_COLUMNS; print('OK', len(MAPPING_CSV_COLUMNS))"
```

Expected: all OK.

## 3. CLI registered

```bash
python -m egfr_myo1d.cli --help | grep prepare-receptor
python -m egfr_myo1d.cli prepare-receptor --help
```

Help text must include `--state`, `--source`, `--profile`, `--mode`, `--strict`.

## 4. Case A — explicit AB receptor

```bash
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase4_local
python -m egfr_myo1d.cli prepare-receptor \
  --run-id m1_phase4_local \
  --state EGFR_160-185 \
  --source fresh/tests/fixtures/m1_phase4_receptor/explicit_AB_dimer.pdb \
  --profile codex_dev
```

Expected outputs (handoff §14.2):

```text
fresh/runs/m1_phase4_local/normalized/receptors/EGFR_160-185_full_frame_explicit_AB.pdb
fresh/runs/m1_phase4_local/normalized/receptors/EGFR_160-185_dockable_669_1014_explicit_AB.pdb
fresh/runs/m1_phase4_local/normalized/receptors/EGFR_160-185_runtime_offset_receptor_only.pdb
fresh/runs/m1_phase4_local/qc/EGFR_160-185_receptor_mapping.csv
fresh/runs/m1_phase4_local/qc/EGFR_160-185_receptor_normalization_audit.csv
fresh/runs/m1_phase4_local/manifest/EGFR_160-185_receptor_manifest.json
```

## 5. Case A — chain identity preservation

Inspect `EGFR_160-185_full_frame_explicit_AB.pdb`:

```text
- chain IDs A and B preserved
- protomer A residue numbers unchanged from source
- protomer B residue numbers unchanged from source (in full_frame)
- ATOM and HETATM records preserved (any caps, modified residues)
```

Inspect `EGFR_160-185_dockable_669_1014_explicit_AB.pdb`:

```text
- residues 634-668 absent
- residues 669-1014 present (or recorded missing in audit)
- chain IDs A and B preserved
- protomer A residue numbers unchanged
- protomer B residue numbers unchanged
- lipid HETATMs (POPC, POPS) absent
- water/ion HETATMs absent
```

Inspect `EGFR_160-185_runtime_offset_receptor_only.pdb`:

```text
- protomer A residue numbers unchanged
- protomer B residue numbers = source + 1000
- runtime_chain may be A and B or unified per implementation choice
```

## 6. Mapping CSV round-trip

```bash
head -3 fresh/runs/m1_phase4_local/qc/EGFR_160-185_receptor_mapping.csv
```

Expected header:

```text
state,source_file,protomer_id,source_chain,source_resseq,source_icode,source_resname,runtime_chain,runtime_resseq,atom_count,role
```

For protomer B rows: runtime_resseq = source_resseq + 1000. For protomer A rows: runtime_resseq = source_resseq.

```bash
python -c "
import sys; sys.path.insert(0, 'fresh/src')
from egfr_myo1d.io.residue_mapping import read_residue_mapping
rows = read_residue_mapping('fresh/runs/m1_phase4_local/qc/EGFR_160-185_receptor_mapping.csv')
print('rows:', len(rows))
print('protomer A samples:', sum(1 for r in rows if r['protomer_id']=='A'))
print('protomer B samples:', sum(1 for r in rows if r['protomer_id']=='B'))
print('offset valid:', all(int(r['runtime_resseq'])==int(r['source_resseq'])+1000 for r in rows if r['protomer_id']=='B'))
"
```

Expected: rows>0, both protomer counts >0, offset valid.

## 7. Case B — duplicate chain X

```bash
python -m egfr_myo1d.cli prepare-receptor \
  --run-id m1_phase4_local \
  --state EGFR_170-200 \
  --source fresh/tests/fixtures/m1_phase4_receptor/duplicate_chain_X_dimer.pdb \
  --profile codex_dev
```

Expected:

```text
- exit 0 (WARN)
- audit CSV records "case_B_duplicate_chain_split_into_A_B"
- normalized PDBs have chain IDs A and B (not X) after split
- manifest case == "B_duplicate_chain"
```

## 8. Case C — monomer in hpc_strict

```bash
python -m egfr_myo1d.cli prepare-receptor \
  --run-id m1_phase4_local \
  --state EGFR_160-185 \
  --source fresh/tests/fixtures/m1_phase4_receptor/single_chain_monomer.pdb \
  --profile hpc_strict
```

Expected:

```text
- exit 1 (FAIL)
- manifest status: FAIL
- audit records "single_chain_monomer_not_valid_for_dimer_analysis"
- no full_frame_explicit_AB output (or output marked blocked)
```

## 9. Case C — monomer in codex_dev

```bash
python -m egfr_myo1d.cli prepare-receptor \
  --run-id m1_phase4_local \
  --state EGFR_160-185 \
  --source fresh/tests/fixtures/m1_phase4_receptor/single_chain_monomer.pdb \
  --profile codex_dev
```

Expected: exit 0 (WARN), status WARN, manifest records `not_valid_for_dimer_only_main_analysis`.

## 10. V924R receptor

```bash
python -m egfr_myo1d.cli prepare-receptor \
  --run-id m1_phase4_local \
  --state EGFR_160-185 \
  --source fresh/tests/fixtures/m1_phase4_receptor/v924r_warn.pdb \
  --profile codex_dev
```

Expected:

```text
- exit 0 (WARN)
- v924r_warn: true in manifest
- v924r_handled: "warn_only_not_mutated"
- audit row for residue 924 records record_type, source_resname=ARG (or matching), warning_classification populated
- output PDB still contains ARG924 (NOT mutated to VAL)
```

## 11. 3GT8_raw reference control

```bash
python -m egfr_myo1d.cli prepare-receptor \
  --run-id m1_phase4_local \
  --state 3GT8_raw \
  --source fresh/tests/fixtures/m1_phase4_receptor/explicit_AB_dimer.pdb \
  --profile codex_dev
```

Expected:

```text
- manifest role: "crystallographic_reference_control_not_primary_membrane_state"
- normalized/receptors/3GT8_raw_explicit_AB.pdb exists
- normalized/receptors/3GT8_raw_runtime_offset_receptor_only.pdb exists
- qc/3GT8_raw_receptor_mapping.csv exists
- dockable_669_1014 may be omitted or marked optional per state config
```

## 12. Path traversal

```bash
python -m egfr_myo1d.cli prepare-receptor --run-id ../bad_run --state EGFR_160-185 --source fresh/tests/fixtures/m1_phase4_receptor/explicit_AB_dimer.pdb
```

Expected nonzero, no outside writes.

## 13. Tests

```bash
pytest -q fresh/tests/test_m1_phase4_receptor_normalization.py
pytest -q fresh/tests
```

Required tests (≥12) pass; existing tests still pass.

## 14. Old workflow protection

```bash
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

Empty.

## 15. What must not be in this phase

```text
- membrane frame computation (Phase 5)
- ligand work (Phase 7)
- integration orchestrator (Phase 8)
- modifying Task 4 prepared_inputs.py output paths (Phase 9)
- mutating V924R back to VAL
- promoting 3GT8_raw to primary membrane-validated state
- modifying old workflow files
```

## 16. Phase 4 accepted if

```text
- model/receptor_normalize.py, model/receptor_qc.py, io/residue_mapping.py created.
- prepare-receptor CLI subcommand registered.
- All 5 output artifacts per primary state emitted.
- Mapping CSV columns match spec; round-trip works.
- +1000 offset applied to protomer B only; protomer A unchanged.
- Case A passes; Case B splits into A/B with WARN; Case C is WARN(codex_dev)/FAIL(hpc_strict).
- 3GT8_raw marked reference_control_not_primary_membrane_state.
- V924R warned, never mutated.
- ATOM and HETATM (including caps) preserved; lipid/water/ion dropped from dockable.
- ≥12 phase tests pass; existing tests pass.
- M1 §23 #10 and #11 closed.
- Old workflow files untouched.
```

## 17. Implementer final response must include

```text
M1 Phase 4 status: PASS / PASS WITH WARNINGS / FAIL
Files created (3 src + 1 io + ~5 fixtures + 1 test + 2 docs)
Files modified (cli.py only)
Commands run and results
Test summary
Output verification per state per case
Mapping CSV round-trip evidence
V924R: warn-only, not mutated
3GT8_raw role: reference_control_not_primary
Old workflow protection
Acceptance closure: M1 §23 #10 and #11 closed
Known limitations: no membrane frame (Phase 5)
```
