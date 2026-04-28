# Task 4 PPI Input Preparation

Task 4 prepares deterministic, auditable PPI input packs for future EGFR-MYO1D sampling. It copies or range-selects validated fixture inputs, writes provenance manifests, and creates QC-only restraint/mask contracts.

Task 4 does not run docking, relaxation, pocket discovery, ligand docking, qsub/PBS, cleanup deletion, scoring, or candidate nomination.

## CLI

Bash:

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id test_task4_local
python -m egfr_myo1d.cli prepare-ppi-inputs --run-id test_task4_local --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task4_inputs
```

PowerShell:

```powershell
$env:PYTHONPATH = "$PWD\fresh\src;$env:PYTHONPATH"
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id test_task4_local
python -m egfr_myo1d.cli prepare-ppi-inputs --run-id test_task4_local --mode smoke_env --profile codex_dev --input-root fresh/tests/fixtures/task4_inputs
```

## Prepared Outputs

Task 4 writes only under `fresh/runs/<run_id>/`:

```text
prepared/egfr/egfr_receptor_normalized.pdb
prepared/myo1d/MYO1D_sheet8_9_12_core_955_1001.pdb
prepared/myo1d/MYO1D_ext_beta_meander_955_1006_tail_masked.pdb
prepared/restraints/myo1d_active_face_contract.json
prepared/restraints/egfr_membrane_exclusion_mask.json
prepared/restraints/egfr_non_target_region_mask.json
qc/egfr_receptor_normalization_audit.csv
qc/myo1d_construct_audit.csv
qc/terminal_artifact_audit.csv
qc/face_orientation_audit.csv
qc/membrane_exclusion_audit.csv
manifest/prepared_input_manifest.json
manifest/preparation_qc_report.json
```

## Terminal Artifacts

The primary MYO1D construct is residues `955-1001`. The `955-960` buffer is retained because starting at `962` exposes an artificial terminus and can create false contact behavior. The `962-1006` fixture is retained only as a negative regression fixture and is not emitted as a prepared production partner.

The `955-1006` construct is emitted only as a comparator/noise-monitor construct. Residues `998-1006` are marked as contact-monitoring or tail/noise residues, not primary binding evidence.

## HETATM Preservation

Task 4 preserves `HETATM` records when they represent caps or biological residues. Cap residues such as `ACE` and `NME` are retained and classified as caps. Standard amino acids encoded as `HETATM`, including the regression fixture `ILE1000`, are treated as biological residues and are not dropped by text filtering.

## Active-Face Annotation

MYO1D active-face residues are written to `myo1d_active_face_contract.json`:

```text
961-964
968-972
```

Sheet-12 support residues are:

```text
993-997
```

These are QC and future sampling annotations only. They do not create scoring bonuses, and `score_bonus_allowed` is always `false`.

## EGFR Masks

Task 4 writes membrane and non-target masks for future PPI sampling. The current synthetic fixtures validate membrane-frame metadata and ATP-site seed masks, but they do not infer production membrane-proximal residues. The mask policy is:

- membrane-proximal modes: exclude or warn in future PPI sampling
- ATP-site overlap: flag as non-target for PPI-modulator objectives
- dimer context: preserve and report accessibility

## Task 3 Warning Classification

Task 3 warnings are allowed only for explicit regression fixtures in `smoke_env/codex_dev`. The same classes become production blockers or quarantine signals for strict, real, or non-whitelisted inputs:

- ambiguous single-chain dimer
- duplicate atom identities
- missing expected EGFR chains `A/B`
- V924R-like ARG924 mutation
- MYO1D terminal/tail artifact fixtures

Task 4 serializes this policy in both `prepared_input_manifest.json` and `preparation_qc_report.json`.

## Not Implemented

Task 4 intentionally does not implement:

- PPI docking
- PyRosetta docking or relaxation
- LightDock
- Vina docking
- fpocket/P2Rank pocket discovery
- compound docking
- EGFR mutation repair
- full receptor modeling
- MYO1D structural prediction
- AI co-folding
- scoring or candidate nomination
- PBS/qsub generation
- cleanup deletion

