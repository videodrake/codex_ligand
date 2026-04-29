# M1 Phase 5 — Membrane Frame Generation

Closes M1 §23 #12 (state-aware `membrane_frame.json` schema is implemented and populated from coordinates) per `milestone1_foundation_codex_handoff_v0_5.md` §15.

## What it does

For each known state (`EGFR_160-185`, `EGFR_170-200`, `3GT8_raw`):

1. Resolve the source PDB path:
   - First try the state's own `input_files[<state>]` from `fresh/configs/receptor_states.yaml`
   - If absent, fall back to `--full-frame-source PATH` (or the receptor_states.yaml `plus10_full_frame` entry)
   - If both missing → `frame_source = "missing"`, status = `WARN` (codex_dev) or `FAIL` (hpc_strict). No vectors fabricated.
2. Parse the source. If chains A and B aren't both present → status = `WARN`/`FAIL`, no vectors.
3. Compute protomer centroids from C-α atoms; `centroid_distance = ||centroid_B - centroid_A||`.
4. Compute `x_dimer_axis = unit(centroid_B - centroid_A)`.
5. Compute `n_membrane` from C-α atoms in residues 634-674 (TM/JM segment) via SVD principal-axis. Sign-orient so the Z component is positive (extracellular → cytosolic per project convention).
6. Compute `p_jm_anchor = mean(TM/JM C-α coords)`.
7. For `3GT8_raw` (reference/control), return a marker entry without computing any vectors. The status field is `"reference_control"` (not PASS/WARN/FAIL).
8. Write the consolidated state-aware `manifest/membrane_frame.json` and `qc/membrane_frame_qc.csv`. Append a `compute-membrane-frame` phase status entry.

Hardcoded fallback axis vectors (`[0, 0, 1]`, `[1, 0, 0]`, etc.) are forbidden in module source outside docstrings — enforced programmatically by `test_no_hardcoded_001_or_100_vectors_in_module_source`.

## CLI

```bash
python -m egfr_myo1d.cli compute-membrane-frame \
    --run-id RUN \
    [--state EGFR_160-185|EGFR_170-200|3GT8_raw|all] \
    [--full-frame-source PATH] \
    [--profile codex_dev|hpc_strict] \
    [--mode smoke_env|smoke_input]
```

## Module additions

```text
fresh/src/egfr_myo1d/model/membrane_frame.py    new (Phase 5)
fresh/src/egfr_myo1d/cli.py                     +compute-membrane-frame subparser + handler
```

Public API:

```python
TM_JM_RESIDUE_START = 634
TM_JM_RESIDUE_END = 674
PRIMARY_STATES = ("EGFR_160-185", "EGFR_170-200")
REFERENCE_CONTROL_STATES = ("3GT8_raw",)
ALL_STATES = PRIMARY_STATES + REFERENCE_CONTROL_STATES

MEMBRANE_FRAME_QC_COLUMNS = [...]   # 10 columns
COORDINATE_CONVENTION = "..."
FRAME_SOURCE_POLICY = {...}

StateMembraneFrame (dataclass)
compute_membrane_frame(state_full_frame_pdb, plus10_full_frame_pdb, state_id, profile) -> StateMembraneFrame
write_state_aware_membrane_frame_json(ctx, frames) -> Path
write_membrane_frame_qc_csv(ctx, frames) -> Path
run_membrane_frame_computation(ctx, state_ids, full_frame_source, profile) -> (frames, overall_status)
```

## Severity

| Status | Conditions |
| --- | --- |
| `PASS` | vectors computed from state full_frame TM/JM with non-zero norms; sensible geometry (n and x not nearly parallel) |
| `WARN` | plus10 fallback used; n and x nearly parallel (`|dot| > 0.95`); chains A/B not both present in codex_dev; missing source in codex_dev |
| `FAIL` | malformed PDB; missing source in hpc_strict; chain A or B C-α atoms missing entirely; zero-norm vector |
| `reference_control` | special-case status for 3GT8_raw — not promoted as a primary frame source; vectors are null |

## Outputs (handoff §15.2)

```text
fresh/runs/<run_id>/manifest/membrane_frame.json
fresh/runs/<run_id>/qc/membrane_frame_qc.csv
fresh/runs/<run_id>/logs/phase_status.jsonl     (appended)
fresh/runs/<run_id>/logs/master.log             (appended)
```

`membrane_frame.json` schema:

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
      "n_membrane": [0.0, 0.0, 1.0],          // computed unit vector
      "x_dimer_axis": [1.0, 0.0, 0.0],         // computed unit vector
      "protomer_a_centroid": [...],
      "protomer_b_centroid": [...],
      "p_jm_anchor": [...],
      "n_membrane_norm": 1.0,                  // post-normalization sanity (1 for non-zero)
      "x_dimer_axis_norm": 1.0,
      "centroid_distance": 16.43,              // raw inter-protomer distance
      "n_tm_jm_residues": 10,
      "status": "PASS",
      "warnings": [],
      "notes": ""
    },
    "EGFR_170-200": { ... },
    "3GT8_raw": {
      "role": "crystallographic_reference_control",
      "frame_source": "not_primary; optional alignment-derived only",
      "n_membrane": null,
      "x_dimer_axis": null,
      "...": null,
      "status": "reference_control",
      "notes": "3GT8_raw is not a primary membrane-frame source"
    }
  },
  "timestamp": "..."
}
```

QC CSV columns:

```csv
state,role,frame_source,status,n_membrane_norm,x_dimer_axis_norm,centroid_distance,n_tm_jm_residues,warnings,notes
```

## Behavior policy (handoff §15)

- Z+ = C2 symmetry axis = membrane normal (extracellular → cytosolic).
- X+ = chain A → chain B (centroid difference).
- TM/JM principal axis from residues 634-674 C-α atoms via SVD.
- Sign convention: principal direction with positive Z component.
- 3GT8_raw is `crystallographic_reference_control`. Vectors are not computed regardless of input.
- Missing source in `codex_dev` → status WARN, vectors null. In `hpc_strict` → FAIL.
- No literal axis vectors in computation paths. The schema may show literal vectors in docstrings/examples only.

## Reusable fixtures

```text
fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_dimer_with_TM_JM.pdb
    Chain A residues 640, 650, 660, 670, 674 along z = -50..0
    Chain A residues 700, 800 (KD) above
    Chain B mirrored at x = 15
    Designed so n_membrane ~ +z and x_dimer_axis ~ +x

fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_dimer_kinase_only.pdb
    Chains A/B with residues 700, 750, 800, 900 only
    No 634-674; tests the missing-TM/JM warning path

fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_3gt8_raw_kinase.pdb
    Synthetic kinase-only crystal-like dimer; for the reference_control case
```

## What is intentionally not in this phase

- Alignment-derived 3GT8_raw frame (out of scope; possible future extension)
- PBS generation (Phase 6)
- Ligand manifest (Phase 7)
- prepare-inputs orchestrator (Phase 8)
- Tasks 4-9 schema realignment (Phase 9)
- Real EGFR/MYO1D PDB placement (user-side, not required for Phase 5 work)
