# M2 PyRosetta Input Sanitation Note

Date: 2026-05-02  
Code reference: commit `bd83048` (`Mirror legacy Rosetta input sanitation for M2 jobs`)  
Primary module: `fresh/src/egfr_myo1d/ppi/run_ppi_job.py`

## Purpose

This note records why M2.2 PyRosetta real execution sanitizes its input PDB immediately before `pose_from_pdb`, and how this differs from the earlier failed fresh input assembly.

The practical failure that triggered this note was:

```text
ERROR: too many tries in fill_missing_atoms!
RuntimeError: UtilityExitException
```

This happened while PyRosetta was importing the generated AB_C input PDB, before any docking move was executed. That means the failure was an input compatibility problem, not evidence that the PPI docking protocol had scientifically failed.

## Key Conclusion

The legacy pipeline did not feed raw or merely concatenated PDBs into PyRosetta. It first produced a Rosetta-friendly, protein-only PDB:

- only `ATOM` records were retained;
- water, ions, caps, and other `HETATM` records were removed;
- hydrogen atoms were removed;
- zero-occupancy atoms were removed;
- unsupported alternate locations were removed;
- CHARMM/PDB residue aliases were normalized;
- known atom-name aliases were normalized;
- residues missing backbone `N`, `CA`, `C`, or `O` were stripped;
- `TER` records were inserted between chains.

The fresh M2 real runner now mirrors that behavior at the final AB_C assembly step.

## Official Rosetta/PyRosetta Context

Relevant official documentation:

- [PyRosetta documentation](https://www.pyrosetta.org/documentation), especially the PDB preparation guidance.
- [RosettaCommons structure preparation](https://docs.rosettacommons.org/docs/latest/rosetta_basics/preparation/preparing-structures), which describes preparing/cleaning PDB files before Rosetta runs.
- [Rosetta docking protocol](https://docs.rosettacommons.org/docs/latest/application_documentation/docking/docking-protocol), which describes chain-based docking partners such as `A_B` or multi-chain partners.

The important operational reading for this project is:

- Rosetta is strict about residue/atom naming and PDB completeness.
- PPI docking inputs must use meaningful chain separation.
- Cleaning the PDB is preferable to hiding input problems behind permissive import flags.

## Legacy Behavior Used As Reference

The successful legacy path used input files like:

```text
output/workflow_a/phase2_ppi_docking/runtime_inputs/docking_EGFR_160-185_ext_beta_meander.pdb
output/workflow_a/phase2_ppi_docking/runtime_inputs/docking_EGFR_170-200_ext_beta_meander.pdb
```

These paths were referenced from:

```text
config/phase1/phase1_prod_EGFR_160-185_seed*.ini
config/phase1/phase1_prod_EGFR_170-200_seed*.ini
config/phase1/phase1_test_*.ini
```

The relevant legacy implementation is:

```text
egfr_pipeline/phase1/prepare_inputs.py
```

Notable legacy functions:

- `normalize_rosetta_input_lines`
- `check_backbone_completeness`
- `strip_incomplete_residues`
- `write_pdb`
- `prepare_docking_pair`

The legacy code also used clear chain separation and `TER` records when writing docking PDBs.

## Fresh M2.2 Current Behavior

Fresh M2.2 uses the biologically intended dimer convention:

```text
EGFR receptor protomer A + EGFR receptor protomer B docked against MYO1D chain C
Rosetta partners: AB_C
```

The M2.1 receptor packs remain the source of truth for chain and residue mapping:

```text
fresh/runs/<run_id>/prepared/m2_1_ppi_inputs/<state>/receptor/<state>_runtime_offset_receptor_only.pdb
fresh/runs/<run_id>/qc/<state>_receptor_mapping.csv
```

The M2.2 real runner then writes a run-local sanitized AB_C input:

```text
fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/outputs/<job_name>/inputs/<job_name>_AB_C_input.pdb
```

This AB_C input is the file passed to:

```python
pyrosetta.pose_from_pdb(...)
```

## Sanitation Rules

The current fresh runner applies these rules immediately before PyRosetta import:

| Input issue | Current action | Why |
| --- | --- | --- |
| `HETATM` records | Drop | Prevent water/ion/cap/unknown residue import failures |
| Hydrogen atoms | Drop | Match legacy Rosetta-friendly protein-only input |
| Zero occupancy atoms | Drop | Avoid invalid or placeholder coordinates |
| Altloc not blank or `A` | Drop | Avoid duplicate/ambiguous atom identities |
| Altloc `A` | Keep and clear altloc field | Collapse to a single concrete atom position |
| `HSD`, `HSE`, `HSP`, `HID`, `HIE`, `HIP` | Rename to `HIS` | Rosetta histidine naming compatibility |
| `CYX` | Rename to `CYS` | Rosetta cysteine naming compatibility |
| `ILE CD` | Rename to `ILE CD1` | Rosetta atom naming compatibility |
| Non-standard residue names | Drop | Avoid unrecognized residue import failures |
| Duplicate atom identity | Drop duplicate | Avoid atom identity collisions |
| Residue missing `N`, `CA`, `C`, or `O` | Drop entire residue | Avoid `fill_missing_atoms` import failure |
| Receptor chain transition A to B | Insert `TER` | Preserve chain separation |
| Partner chain | Rewrite as `C` | Maintain `AB_C` docking convention |

## What This Fix Does Not Yet Prove

Passing `pose_from_pdb` proves only that PyRosetta can import the cleaned AB_C input.

It does not prove:

- the docking protocol has produced biologically meaningful poses;
- the receptor and MYO1D relative starting geometry is optimal;
- side-chain packing has been fully prepared;
- the input is prepacked according to the full RosettaDock recommendation;
- the resulting contacts are scientifically accepted.

RosettaDock documentation commonly expects a cleaned/prepared PDB and, for production-quality work, a prepacking step. If real docking proceeds but scores or clashes look suspicious, the next reviewed improvement should be an explicit prepack or constrained-relax preparation step, not arbitrary score filtering.

## HPC Reproduction Command

After pulling a commit at or after `bd83048`, run one real model as an import/docking smoke:

```bash
cd /work4/hwang/onepack/new/codex_ligand
git pull --ff-only origin main
git rev-parse --short HEAD

conda activate pyrosetta
export PYTHONPATH="$PWD/fresh/src:$PYTHONPATH"

RUN_ID=input_clean_20260501_211242
JOB=$(python - <<'PY'
import csv
from pathlib import Path
run = Path("fresh/runs") / "input_clean_20260501_211242"
with (run / "phase1_ppi/pyrosetta_adapter/pyrosetta_job_manifest.csv").open() as f:
    print(next(csv.DictReader(f))["job_name"])
PY
)

python -m egfr_myo1d.ppi.run_ppi_job \
  --run-id "$RUN_ID" \
  --job-name "$JOB" \
  --dry-run false \
  --allow-real true \
  --max-models 1
```

## 32-Core Production Planning

The old workflow used PBS allocation as the local worker limit:

```text
config/run_ppi_state_seed.pbs
run_production.py --allocated-cpus "${PBS_NP:-16}"
egfr_pipeline/pyrosetta_docking/pipeline_manager.py -> multiprocessing.Pool(n_cpus)
```

That pattern is useful, but it is not copied blindly. In fresh M2.2, the safer production plan is:

- PBS requests the intended node and core count, e.g. `nodes=node04:ppn=32`;
- the PBS body reads `PBS_NP` and caps local concurrent workers at that value;
- each worker is an isolated `python -m egfr_myo1d.ppi.run_ppi_job` subprocess with its own `--chunk-id`;
- every chunk writes its own `chunks/<chunk_id>/real_run_status.json` and `chunks/<chunk_id>/pose_scores.csv`;
- all successful pose PDBs still land in the job-level `poses/` directory with unique model indices.

This keeps the useful legacy resource contract while avoiding shared PyRosetta objects across a Python multiprocessing pool. The practical advantage is cleaner failure isolation: if one chunk fails, the completed chunk statuses and pose scores remain inspectable and the failed model range can be rerun directly.

For final-scale sampling, regenerate the dry harness with an explicit model count before making the PBS plan. The legacy production shape used 20,000 models per state/seed; use a smaller value only for smoke, mini, or calibration runs:

```bash
RUN_ID=input_clean_20260501_211242

python -m egfr_myo1d.cli prepare-m2-pyrosetta-harness \
  --run-id "$RUN_ID" \
  --mode production \
  --profile hpc_strict \
  --states EGFR_160-185,EGFR_170-200 \
  --models-per-seed 20000
```

Then generate the 3-node, 32-core plan:

```bash
RUN_ID=input_clean_20260501_211242

python fresh/scripts/submit_m2_pyrosetta_real_jobs.py \
  --run-id "$RUN_ID" \
  --profile production \
  --nodes node04,node05,node06 \
  --ppn 32
```

The script does not call `qsub` by default. It writes:

```text
fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/pyrosetta_real_chunk_manifest.csv
fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/pyrosetta_real_pbs_manifest.csv
fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/qsub/m2_pyrosetta_node04.pbs
fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/qsub/m2_pyrosetta_node05.pbs
fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/qsub/m2_pyrosetta_node06.pbs
fresh/runs/<run_id>/reports/m2_2_pyrosetta_real_pbs_plan.md
```

By default, chunk sizing is automatic: each state/seed job is split into roughly one 32-worker wave. With `--models-per-seed 20000 --ppn 32`, that means about 625 models per subprocess chunk. You can override this with `--models-per-chunk N` if a calibration run shows that chunks are too short or too long.

Submit manually only after inspecting the PBS manifest:

```bash
cat fresh/runs/$RUN_ID/phase1_ppi/pyrosetta_adapter/pyrosetta_real_pbs_manifest.csv

qsub fresh/runs/$RUN_ID/phase1_ppi/pyrosetta_adapter/qsub/m2_pyrosetta_node04.pbs
qsub fresh/runs/$RUN_ID/phase1_ppi/pyrosetta_adapter/qsub/m2_pyrosetta_node05.pbs
qsub fresh/runs/$RUN_ID/phase1_ppi/pyrosetta_adapter/qsub/m2_pyrosetta_node06.pbs
```

If import still fails, inspect:

```bash
cat fresh/runs/$RUN_ID/phase1_ppi/pyrosetta_adapter/outputs/$JOB/real_run_status.json
ls fresh/runs/$RUN_ID/phase1_ppi/pyrosetta_adapter/outputs/$JOB/inputs/
```

The status JSON records:

- `status`;
- `docking_executed`;
- `input_summary.combined_input_pdb`;
- `input_summary.receptor_rosetta_sanitize`;
- `input_summary.partner_rosetta_sanitize`;
- failure text if `pose_from_pdb` failed.

## Expected Debug Interpretation

### `REAL_RUN_INPUT_IMPORT_FAIL`

PyRosetta still could not import the sanitized AB_C input. Check:

- whether the sanitized input has any non-standard residue names left;
- whether chain A, B, and C are present;
- whether any chain has zero atoms after sanitation;
- whether residues near the failure have missing backbone atoms;
- whether Rosetta requires a residue patch or parameter file not currently supplied.

### `REAL_RUN_FAIL`

PyRosetta imported the input, but one or more docking models failed after import. This is no longer the same problem as the original `fill_missing_atoms` failure. Inspect pose-level failures in:

```text
fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/outputs/<job_name>/real_run_status.json
fresh/runs/<run_id>/phase1_ppi/pyrosetta_adapter/outputs/<job_name>/pose_scores.csv
```

### `REAL_RUN_PASS`

At least the requested model count completed. This is still not the final scientific gate. Continue with contact extraction, residue mapping restoration, clustering/consensus, and manual sanity checks.

## Guardrails

- Do not bypass this sanitation by passing raw receptor or MYO1D PDBs directly into PyRosetta.
- Do not use import flags to hide unrecognized residues unless the project explicitly adds matching params files and tests.
- Do not commit private ligand structures or private compound maps unless explicitly approved.
- Do not treat dry-run success as docking success.
- Do not treat `pose_from_pdb` success as biological acceptance.
- Before larger qsub submission, always run a single `--max-models 1` real job and inspect the generated AB_C input.

## Related Files

```text
fresh/src/egfr_myo1d/ppi/run_ppi_job.py
fresh/src/egfr_myo1d/ppi/pyrosetta_real_jobs.py
fresh/scripts/submit_m2_pyrosetta_real_jobs.py
fresh/tests/test_m2_phase2_pyrosetta_adapter.py
fresh/tests/test_m2_phase2c_pyrosetta_real_pbs.py
fresh/tests/test_m2_phase2b_pose_contacts.py
egfr_pipeline/phase1/prepare_inputs.py
config/phase1/phase1_prod_EGFR_160-185_seed*.ini
config/phase1/phase1_prod_EGFR_170-200_seed*.ini
```
