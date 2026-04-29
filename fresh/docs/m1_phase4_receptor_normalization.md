# M1 Phase 4 — Receptor Normalization

Closes M1 §23 #10 (explicit A/B normalization) and #11 (+1000 runtime offset) per `milestone1_foundation_codex_handoff_v0_5.md` §14.

## What it does

For an EGFR receptor PDB and a known state (`EGFR_160-185`, `EGFR_170-200`, or `3GT8_raw`):

1. **Detect input case** (`model.receptor_qc.detect_normalization_case`):
   - **Case A** — explicit A/B chains
   - **Case B** — same-chain duplicate dimer (e.g., both copies stored as chain X)
   - **Case C** — single-chain monomer (NOT VALID for dimer-only main analysis)
   - **ambiguous** — multi-chain but neither A nor B explicitly present
2. **Apply receptor-only filter**: drop lipid (POPC/POPS/POPE/POPG/CHL/...), water (HOH/WAT/TIP3), and ion (NA/CL/K/MG/...) HETATM records globally; keep ACE/NME caps and standard amino acids written as HETATM.
3. **Split protomers** for Case B via duplicate-atom-identity tracking; rewrite chain IDs to `A` and `B` in both the dataclass field and the underlying raw PDB line.
4. **Detect V924R-like mutation** (handoff §13). Always **WARN only**, never mutate. Report records the hit and `v924r_handled = "warn_only_not_mutated"`.
5. **Emit three normalized PDBs** under `runs/<run_id>/normalized/receptors/`:
   - `<state>_full_frame_explicit_AB.pdb` (or `<state>_explicit_AB.pdb` for `3GT8_raw`)
   - `<state>_dockable_<crop_start>_<crop_end>_explicit_AB.pdb` (omitted for `3GT8_raw`; per handoff §14.2 it is a kinase-only crystal)
   - `<state>_runtime_offset_receptor_only.pdb` — protomer A unchanged, protomer B residue numbers + offset (default 1000 from gates.yaml). The raw PDB line is rewritten so cols [22:26] reflect the new residue number.
6. **Write the residue mapping CSV** at `runs/<run_id>/qc/<state>_receptor_mapping.csv` with the schema from handoff §14.3 (11 columns). Per residue: original chain/resseq/icode/resname + runtime chain/resseq + atom_count + role.
7. **Write the receptor normalization audit CSV** at `runs/<run_id>/qc/<state>_receptor_normalization_audit.csv` with one row per residue, recording crop membership, biological vs cap vs non-receptor HETATM classification, and any V924R hit.
8. **Write the receptor manifest JSON** at `runs/<run_id>/manifest/<state>_receptor_manifest.json` with source + output sha256, case, observed chains, V924R warn, dockable crop, runtime offset, status, warnings.
9. **Append a `prepare-receptor` phase status** entry to `phase_status.jsonl` and `master.log`.

## CLI

```bash
python -m egfr_myo1d.cli prepare-receptor \
    --run-id RUN \
    --state EGFR_160-185|EGFR_170-200|3GT8_raw \
    --source PATH/to/receptor.pdb \
    [--profile codex_dev|hpc_strict] \
    [--mode smoke_env|smoke_input] \
    [--strict]
```

## Module additions

```text
fresh/src/egfr_myo1d/model/__init__.py             new
fresh/src/egfr_myo1d/model/receptor_normalize.py   new   (main module)
fresh/src/egfr_myo1d/model/receptor_qc.py          new   (case detection + V924R + audit row helper)
fresh/src/egfr_myo1d/io/residue_mapping.py         new   (mapping CSV writer/reader)
```

Public API:

```python
# model/receptor_normalize.py
normalize_receptor(ctx, source_pdb, state_id, profile="codex_dev", strict=False, ...) -> NormalizedReceptor
NormalizedReceptor (dataclass)
RECEPTOR_AUDIT_CSV_COLUMNS  # list of 15 column names
PRIMARY_STATES = ("EGFR_160-185", "EGFR_170-200")
REFERENCE_CONTROL_STATES = ("3GT8_raw",)
resolve_role(state_id) -> str
load_receptor_gate(ctx) -> dict   # reads gates.yaml receptor section

# model/receptor_qc.py
detect_normalization_case(structure) -> "A_explicit_AB"|"B_duplicate_chain"|"C_monomer"|"ambiguous"
split_duplicate_chain(atoms) -> (a_atoms, b_atoms)
detect_warn_mutations(atoms, warn_mutations=DEFAULT_WARN_MUTATIONS) -> list[hit]
compute_residue_audit_rows(...) -> list[dict]
is_receptor_atom(atom) -> bool
NON_RECEPTOR_HETATM, CAP_RESNAMES, DEFAULT_WARN_MUTATIONS, WarnMutation

# io/residue_mapping.py
write_residue_mapping(path, rows, ctx=None)
read_residue_mapping(path) -> list[dict]
MAPPING_CSV_COLUMNS  # the 11 spec columns
```

## Severity

| Status | Conditions |
| --- | --- |
| `PASS` | Case A; required residues in dockable crop; no V924R; no missing chains |
| `WARN` | Case A with V924R; Case B successfully split; Case C in `codex_dev`; ambiguous chains in `codex_dev`; capped HETATMs preserved |
| `FAIL` | Case C in `hpc_strict` or `--strict`; ambiguous in strict; parse error; unknown profile; missing source in `hpc_strict`; write attempt outside `run_dir` |

## Outputs (handoff §14.2)

```text
For EGFR_160-185 / EGFR_170-200:
  fresh/runs/<run_id>/normalized/receptors/<state>_full_frame_explicit_AB.pdb
  fresh/runs/<run_id>/normalized/receptors/<state>_dockable_669_1014_explicit_AB.pdb
  fresh/runs/<run_id>/normalized/receptors/<state>_runtime_offset_receptor_only.pdb
  fresh/runs/<run_id>/qc/<state>_receptor_mapping.csv
  fresh/runs/<run_id>/qc/<state>_receptor_normalization_audit.csv
  fresh/runs/<run_id>/manifest/<state>_receptor_manifest.json

For 3GT8_raw:
  fresh/runs/<run_id>/normalized/receptors/3GT8_raw_explicit_AB.pdb
  fresh/runs/<run_id>/normalized/receptors/3GT8_raw_runtime_offset_receptor_only.pdb
  fresh/runs/<run_id>/qc/3GT8_raw_receptor_mapping.csv
  fresh/runs/<run_id>/qc/3GT8_raw_receptor_normalization_audit.csv
  fresh/runs/<run_id>/manifest/3GT8_raw_receptor_manifest.json
  (no dockable_*.pdb — 3GT8_raw is a kinase-only crystal reference/control)
```

Mapping CSV columns:

```csv
state,source_file,protomer_id,source_chain,source_resseq,source_icode,source_resname,runtime_chain,runtime_resseq,atom_count,role
```

Audit CSV columns:

```csv
state,fixture_role,chain_id,protomer_id,residue_number,insertion_code,residue_name,record_type,atom_count,in_dockable_crop,biological,is_cap,is_warn_mutation,classification,notes
```

## Behavior policy (handoff §14.3, §14.4 + §3 of phase prompt)

- Original residue numbers preserved on protomer A and on the full_frame/dockable outputs.
- Protomer B receives runtime residue number = source + offset (default 1000) only in the `runtime_offset_receptor_only.pdb`. The raw PDB line is rewritten so the column-22 chain ID and column-22-25 residue number both match the new values.
- Insertion codes preserved.
- V924R (or other warned mutations) is reported in the audit + manifest, never mutated to WT.
- Case B duplicate-chain split: first occurrence of each `(chain, resseq, icode, resname, atom_name)` becomes protomer A; second occurrence becomes protomer B. Chain IDs are rewritten in the raw line as well.
- Case C (monomer) is rejected as production-primary in `hpc_strict`/`--strict`.
- HETATM caps (ACE, NME, etc.) preserved; standard amino acids written as HETATM kept as biological residues; lipid (POPC/POPS/etc), water (HOH/WAT/TIP3), and ion (NA/CL/K/MG/...) HETATM records dropped from the dockable PDB.
- 3GT8_raw role is `crystallographic_reference_control_not_primary_membrane_state`. It does not produce the `dockable_<crop>_explicit_AB.pdb` artifact (M1→M2 transition gate item 4 must come from a primary state, not 3GT8_raw).

## Reusable fixtures

```text
fresh/tests/fixtures/m1_phase4_receptor/explicit_AB_dimer.pdb       (Case A; 2 chains × 2 residues)
fresh/tests/fixtures/m1_phase4_receptor/duplicate_chain_X_dimer.pdb (Case B; chain X with reset; 5+5 atoms)
fresh/tests/fixtures/m1_phase4_receptor/single_chain_monomer.pdb    (Case C; 8 ATOM records on chain A)
fresh/tests/fixtures/m1_phase4_receptor/v924r_warn.pdb              (V924R synthetic regression)
fresh/tests/fixtures/m1_phase4_receptor/dimer_with_TM_excluded_range.pdb (covers crop 634-1100 + lipid/water HETATM drop)
```

## What is intentionally not in this phase

- Membrane frame computation (Phase 5)
- PBS generation (Phase 6)
- Ligand manifest (Phase 7)
- prepare-inputs orchestrator (Phase 8)
- Tasks 4-9 schema realignment (Phase 9)
- Mutation repair (V924R is reported, never altered)
- 3GT8_raw promotion to primary membrane-validated state (forbidden)
- Modification of old workflow files
