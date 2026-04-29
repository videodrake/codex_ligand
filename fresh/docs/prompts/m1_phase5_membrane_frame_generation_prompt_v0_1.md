# Claude M1 Phase 5 Prompt — Membrane Frame Generation v0.1

Branch `claude/task10`. Phases 1-4 complete. This is **M1 Phase 5** — implements state-aware membrane frame generation per `milestone1_foundation_codex_handoff_v0_5.md` §15, closing M1 §23 #12.

## 1. Project context

The project's primary receptor model is the +10° rotation membrane-validated EGFR TM-JM-kinase dimer (provenance: `EGFR_MYO1D_additional_MD_clustering_provenance.md`). The downstream geometry frame is inherited from this +10° model.

Coordinate convention (handoff §15.1):

```text
Z+ = C2 symmetry axis = membrane normal (extracellular -> intracellular/cytosolic)
X+ = chain A -> chain B
lower    = TM/JM/membrane-proximal kinase-domain side
lateral  = outward-facing surface away from central dimer interface
```

Currently `structure/contracts.py` (Task 3) **validates** an existing `membrane_frame.json` schema, but the workflow has no module that **computes** the frame from coordinates. Phase 5 fills this gap.

## 2. Absolute rules

Do not modify the old workflow. Maintain Py2/3 syntax compatibility.

Coordinate vectors must be **computed from input coordinates**, not hardcoded. Hardcoded `[0,0,1]` or `[1,0,0]` may appear only in documentation/example schema text — never as a fallback computation result.

Source-of-truth values from `fresh/configs/receptor_states.yaml`:

```yaml
states:
  EGFR_160-185:  primary_membrane_validated_state
  EGFR_170-200:  primary_membrane_validated_state
  3GT8_raw:      crystallographic_reference_control
plus10_full_frame:
  expected_input: fresh/data/raw/receptors/plus10_full_frame.pdb
  role: model_derived_membrane_frame_source
```

## 3. Scope

In scope:
- Create `fresh/src/egfr_myo1d/model/membrane_frame.py`
- Add `compute-membrane-frame` CLI subcommand
- Tests under `fresh/tests/test_m1_phase5_membrane_frame_generation.py` (≥9 tests)
- Fixtures under `fresh/tests/fixtures/m1_phase5_membrane_frame/` (synthetic full-frame PDBs with TM/JM coordinates that produce computable frames)
- Docs `fresh/docs/m1_phase5_membrane_frame_generation.md` and `m1_phase5_changes.md`

Out of scope:
- PBS generation (Phase 6)
- Ligand work (Phase 7)
- Validating an externally-supplied membrane_frame.json — that is `structure/contracts.py` (Task 3) territory and stays unchanged
- Modifying any Task 1-9 modules

## 4. Required CLI behavior

```bash
python -m egfr_myo1d.cli compute-membrane-frame \
  --run-id RUN \
  [--state EGFR_160-185|EGFR_170-200|3GT8_raw|all] \
  [--full-frame-source fresh/data/raw/receptors/plus10_full_frame.pdb] \
  [--profile codex_dev|hpc_strict] \
  [--mode smoke_env|smoke_input]
```

Behavior:
- `--state all` (default): compute frames for all three states
- `--state X`: compute only that state's frame
- `--full-frame-source`: fallback PDB containing TM/JM coordinates (default: plus10_full_frame.pdb)
- Process exit: 0 PASS/WARN, 1 FAIL

## 5. Files to create / modify

Create:

```text
fresh/src/egfr_myo1d/model/membrane_frame.py
fresh/tests/test_m1_phase5_membrane_frame_generation.py
fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_dimer_with_TM_JM.pdb
fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_dimer_kinase_only.pdb
fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_3gt8_raw_kinase.pdb
fresh/docs/m1_phase5_membrane_frame_generation.md
fresh/docs/m1_phase5_changes.md
```

Modify:

```text
fresh/src/egfr_myo1d/cli.py   # add compute-membrane-frame subparser + handler
```

## 6. Public API

`model/membrane_frame.py`:

```python
def compute_membrane_frame(ctx, full_frame_pdb, plus10_full_frame_pdb, state_id, profile):
    # type: (RunContext, Path | None, Path | None, str, str) -> StateMembraneFrame
    """
    Compute (or inherit) membrane frame for a state.
    Strategy:
      1. If full_frame_pdb has TM/JM atoms (residues 634-674), compute frame from those coords.
      2. Else if plus10_full_frame_pdb is supplied and has TM/JM, inherit frame from it.
      3. Else mark status=missing_frame_source; do NOT invent vectors.
    Compute:
      protomer_a_centroid, protomer_b_centroid (from C-alpha within respective chains)
      x_dimer_axis = normalize(protomer_b_centroid - protomer_a_centroid)
      n_membrane = principal axis of TM/JM atoms (e.g., PCA of CA coords in 634-674)
      p_jm_anchor = midpoint of TM/JM CA cloud
    Returns StateMembraneFrame with status PASS / WARN / FAIL_MISSING_SOURCE.
    """

def write_state_aware_membrane_frame_json(ctx, frames):
    # frames: list of StateMembraneFrame for each state
    # writes runs/<run_id>/manifest/membrane_frame.json (state-aware schema, handoff §15.2)

def write_membrane_frame_qc_csv(ctx, frames):
    # writes runs/<run_id>/qc/membrane_frame_qc.csv
```

`StateMembraneFrame` dataclass:

```text
state_id: str
role: "primary_membrane_validated_state" | "crystallographic_reference_control"
frame_source: "state_full_frame" | "plus10_inherited" | "missing"
n_membrane: tuple[float, float, float] | None
x_dimer_axis: tuple[float, float, float] | None
protomer_a_centroid: tuple[float, float, float] | None
protomer_b_centroid: tuple[float, float, float] | None
p_jm_anchor: tuple[float, float, float] | None
status: "PASS" | "WARN" | "FAIL"
warnings: list[str]
notes: str
```

## 7. Required output files (handoff §15.2)

After `compute-membrane-frame --state all` runs:

```text
fresh/runs/<run_id>/manifest/membrane_frame.json     # state-aware schema below
fresh/runs/<run_id>/qc/membrane_frame_qc.csv
fresh/runs/<run_id>/logs/phase_status.jsonl          # appended
fresh/runs/<run_id>/logs/master.log                  # appended
```

`membrane_frame.json` schema (handoff §15.2):

```json
{
  "coordinate_convention": "Z+ is membrane normal from extracellular to intracellular/cytosolic; X+ is chain A to chain B",
  "frame_source_policy": {
    "primary":         "state full-frame coordinates when available",
    "fallback":        "plus10_full_frame coordinate frame or TM/JM 634-674 axis",
    "reference_only":  "3GT8_raw is not a primary membrane-frame source"
  },
  "states": {
    "EGFR_160-185": {
      "role": "primary_membrane_validated_state",
      "frame_source": "state_full_frame",
      "n_membrane": [<computed_x>, <computed_y>, <computed_z>],
      "x_dimer_axis": [<computed_x>, <computed_y>, <computed_z>],
      "protomer_a_centroid": [...],
      "protomer_b_centroid": [...],
      "p_jm_anchor": [...],
      "status": "PASS"
    },
    "EGFR_170-200": { ... },
    "3GT8_raw": {
      "role": "crystallographic_reference_control",
      "frame_source": "not_primary; optional alignment-derived only",
      "status": "reference_control"
    }
  }
}
```

`membrane_frame_qc.csv` columns:

```csv
state,role,frame_source,status,n_membrane_norm,x_dimer_axis_norm,centroid_distance,warnings,notes
```

## 8. Behavior policy

```text
- Compute n_membrane from TM/JM (residues 634-674) C-alpha coords:
    select CA atoms in residues 634-674 across both chains
    if >=4 atoms, compute principal axis (largest variance) via 3x3 covariance / PCA
    orient toward extracellular side (per convention; if unsure, document and warn)
- Compute x_dimer_axis from chain A and B centroids:
    centroid = mean(CA positions per chain)
    x_dimer_axis = normalize(centroid_B - centroid_A)
- Compute p_jm_anchor as midpoint of TM/JM CA cloud (mean of selected CAs).
- If state full_frame lacks 634-674: try plus10_full_frame fallback.
- If plus10_full_frame missing or also lacks 634-674: status=FAIL_MISSING_SOURCE; vectors null; do NOT invent [0,0,1].
- 3GT8_raw: mark role=crystallographic_reference_control, frame_source not_primary; vectors omitted unless explicit alignment-derived value provided (out of scope this phase).
- Vectors must be computed values, not literal [0,0,1] or [1,0,0]. The schema example may show literal vectors; the actual output for valid input must show numerical values derived from coordinates.
- Warn (status=WARN) when:
    - n_membrane and x_dimer_axis are nearly parallel (|dot| > 0.95) — suggests bad geometry
    - protomer centroids are identical (suggests Case C monomer — should have been blocked by Phase 4 already)
    - TM/JM residue count below threshold (e.g., <4 CAs)
- FAIL when:
    - missing source AND no fallback
    - malformed PDB
    - vectors have zero norm
```

## 9. Severity rules

```text
PASS: vectors computed from state full_frame TM/JM coords with non-zero norms, sensible geometry
WARN: vectors computed from plus10 fallback; nearly parallel n/x; partial TM/JM residues
FAIL: missing both state full_frame and plus10 fallback; malformed PDB; zero norm
```

## 10. Tests required (≥9)

```text
test_membrane_normal_computed_from_TM_JM_coords_when_available
test_x_dimer_axis_computed_from_chain_A_to_B_centroids
test_plus10_fallback_used_when_state_lacks_TM_JM
test_3gt8_raw_marked_reference_control_not_primary_source
test_missing_frame_source_reported_cleanly_no_invented_vectors
test_membrane_frame_json_state_aware_schema
test_membrane_frame_qc_csv_columns_match_spec
test_protomer_centroids_recorded
test_nearly_parallel_n_x_warns
test_no_hardcoded_001_or_100_vectors_in_module_source
test_cli_help_includes_compute_membrane_frame
```

(11 tests; ≥9 required.)

The `test_no_hardcoded_001_or_100_vectors_in_module_source` test reads `model/membrane_frame.py` source and greps for literal `[0, 0, 1]` / `[0,0,1]` / `[1, 0, 0]` etc. as actual computation results (allow them only inside docstrings/example schema).

## 11. Acceptance commands

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"

pytest -q fresh/tests/test_m1_phase5_membrane_frame_generation.py
pytest -q fresh/tests

python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase5_local

# All-state computation
python -m egfr_myo1d.cli compute-membrane-frame --run-id m1_phase5_local --state all --full-frame-source fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_dimer_with_TM_JM.pdb --profile codex_dev

# 3GT8_raw alone (kinase-only fallback)
python -m egfr_myo1d.cli compute-membrane-frame --run-id m1_phase5_local --state 3GT8_raw --full-frame-source fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_3gt8_raw_kinase.pdb --profile codex_dev

# Missing source case
python -m egfr_myo1d.cli compute-membrane-frame --run-id m1_phase5_local --state EGFR_160-185 --full-frame-source fresh/data/raw/receptors/does_not_exist.pdb --profile hpc_strict || echo "Expected FAIL"

# Path traversal
python -m egfr_myo1d.cli compute-membrane-frame --run-id ../bad_run --state all

# Old workflow protection
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/

# Inspect output
python -c "import json; print(json.dumps(json.load(open('fresh/runs/m1_phase5_local/manifest/membrane_frame.json')), indent=2))" | head -60
```

## 12. Final response format

```text
M1 Phase 5 status: PASS / PASS WITH WARNINGS / FAIL
Files created
Files modified
Commands run and results
Test summary
Vector computation evidence:
- n_membrane: numerical, non-trivial, derived from TM/JM
- x_dimer_axis: numerical, derived from chain A->B centroids
- p_jm_anchor: numerical
- no hardcoded [0,0,1] / [1,0,0] in actual computation paths
Per-state status:
- EGFR_160-185: PASS or WARN
- EGFR_170-200: PASS or WARN
- 3GT8_raw: reference_control (not primary)
Missing-source behavior: FAIL in hpc_strict, no invented vectors
Acceptance closure: M1 §23 #12 closed
Old workflow protection: empty diff
Known limitations:
- No alignment-derived 3GT8 frame in this phase
- No PBS generation (Phase 6)
```
