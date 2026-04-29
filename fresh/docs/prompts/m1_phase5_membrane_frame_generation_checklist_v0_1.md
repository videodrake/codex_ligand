# M1 Phase 5 Acceptance Checklist v0.1 — Membrane Frame Generation

Use this after the implementer applies M1 Phase 5.

## 1. Pre-Phase state preserved

```text
Old workflow files unchanged.
Phases 1-4 outputs/modules unchanged.
structure/contracts.py (Task 3 membrane_frame validation) unchanged.
```

## 2. New module importable

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"
python -c "from egfr_myo1d.model.membrane_frame import compute_membrane_frame, write_state_aware_membrane_frame_json, write_membrane_frame_qc_csv; print('OK')"
```

## 3. CLI registered

```bash
python -m egfr_myo1d.cli --help | grep compute-membrane-frame
python -m egfr_myo1d.cli compute-membrane-frame --help
```

## 4. All-state computation produces correct outputs

```bash
python -m egfr_myo1d.cli init-run --mode smoke_env --run-id m1_phase5_local
python -m egfr_myo1d.cli compute-membrane-frame \
  --run-id m1_phase5_local \
  --state all \
  --full-frame-source fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_dimer_with_TM_JM.pdb \
  --profile codex_dev
```

Expected outputs:

```text
fresh/runs/m1_phase5_local/manifest/membrane_frame.json
fresh/runs/m1_phase5_local/qc/membrane_frame_qc.csv
```

## 5. Schema validation

`membrane_frame.json` must contain (handoff §15.2):

```text
- "coordinate_convention" string
- "frame_source_policy" object with primary / fallback / reference_only entries
- "states" object with three keys: EGFR_160-185, EGFR_170-200, 3GT8_raw
- For each primary state:
  - role: "primary_membrane_validated_state"
  - frame_source: "state_full_frame" or "plus10_inherited"
  - n_membrane: array length 3, numeric, non-zero norm
  - x_dimer_axis: array length 3, numeric, non-zero norm
  - protomer_a_centroid: array length 3, numeric
  - protomer_b_centroid: array length 3, numeric
  - p_jm_anchor: array length 3, numeric
  - status: "PASS" or "WARN"
- For 3GT8_raw:
  - role: "crystallographic_reference_control"
  - frame_source: "not_primary" or similar, NOT "state_full_frame"
  - vectors may be null
  - status: "reference_control"
```

```bash
python -c "
import json
data = json.load(open('fresh/runs/m1_phase5_local/manifest/membrane_frame.json'))
states = data['states']
print('keys:', sorted(states.keys()))
import math
for sid in ['EGFR_160-185', 'EGFR_170-200']:
    s = states[sid]
    n = s.get('n_membrane'); x = s.get('x_dimer_axis')
    nn = math.sqrt(sum(c*c for c in n)) if n else 0
    xn = math.sqrt(sum(c*c for c in x)) if x else 0
    print(sid, 'n_norm=', round(nn,4), 'x_norm=', round(xn,4), 'status=', s['status'])
print('3GT8_raw:', states['3GT8_raw'])
"
```

Expected:
- both n_norm and x_norm close to 1.0 (normalized) for primary states
- 3GT8_raw marked reference_control

## 6. Vectors are computed, not hardcoded

Inspect module source for literal `[0, 0, 1]` / `[1, 0, 0]` in computation paths:

```bash
grep -nE "\\[\\s*0\\s*,\\s*0\\s*,\\s*1\\s*\\]|\\[\\s*1\\s*,\\s*0\\s*,\\s*0\\s*\\]" fresh/src/egfr_myo1d/model/membrane_frame.py
```

If matches exist, they must be inside docstrings or example schema strings, not as fallback computation results. The corresponding test `test_no_hardcoded_001_or_100_vectors_in_module_source` enforces this programmatically.

## 7. Plus10 fallback path

```bash
python -m egfr_myo1d.cli compute-membrane-frame \
  --run-id m1_phase5_local \
  --state EGFR_160-185 \
  --full-frame-source fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_dimer_kinase_only.pdb \
  --profile codex_dev
```

If state's full_frame lacks 634-674, the implementation should fall back to `plus10_full_frame` (default expected at `fresh/data/raw/receptors/plus10_full_frame.pdb`). If that's also unavailable in the test env, frame_source should reflect the actual chain of fallback or `missing`.

## 8. Missing-source behavior

```bash
python -m egfr_myo1d.cli compute-membrane-frame \
  --run-id m1_phase5_local \
  --state EGFR_160-185 \
  --full-frame-source fresh/data/raw/receptors/does_not_exist.pdb \
  --profile hpc_strict
```

Expected:

```text
- exit 1 (FAIL)
- For EGFR_160-185 in membrane_frame.json (if file written): status=FAIL_MISSING_SOURCE, vectors null
- No invented vectors
```

## 9. 3GT8_raw alone

```bash
python -m egfr_myo1d.cli compute-membrane-frame \
  --run-id m1_phase5_local \
  --state 3GT8_raw \
  --full-frame-source fresh/tests/fixtures/m1_phase5_membrane_frame/synthetic_3gt8_raw_kinase.pdb \
  --profile codex_dev
```

Expected:

```text
- 3GT8_raw entry: role=crystallographic_reference_control, frame_source not_primary
- No primary membrane-validated frame promotion from 3GT8_raw
```

## 10. QC CSV

```bash
head -5 fresh/runs/m1_phase5_local/qc/membrane_frame_qc.csv
```

Expected header:

```text
state,role,frame_source,status,n_membrane_norm,x_dimer_axis_norm,centroid_distance,warnings,notes
```

## 11. Path traversal

```bash
python -m egfr_myo1d.cli compute-membrane-frame --run-id ../bad_run --state all
```

Nonzero exit, no outside writes.

## 12. Tests

```bash
pytest -q fresh/tests/test_m1_phase5_membrane_frame_generation.py
pytest -q fresh/tests
```

## 13. Old workflow protection

```bash
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/
```

Empty.

## 14. What must not be in this phase

```text
- PBS generation (Phase 6)
- ligand work (Phase 7)
- prepare-inputs orchestrator (Phase 8)
- Tasks 4-9 schema realignment (Phase 9)
- alignment-based 3GT8 frame derivation
- modifying structure/contracts.py validation logic
- modifying old workflow files
```

## 15. Phase 5 accepted if

```text
- model/membrane_frame.py created with three public functions.
- compute-membrane-frame CLI subcommand registered.
- manifest/membrane_frame.json state-aware schema written.
- qc/membrane_frame_qc.csv columns match spec.
- Vectors computed from coordinates, not hardcoded.
- 3GT8_raw marked reference_control, not primary frame source.
- Missing-source path produces FAIL or WARN cleanly with vectors=null (no invented).
- ≥9 phase tests pass; existing tests pass.
- M1 §23 #12 closed.
- Old workflow files unmodified.
```

## 16. Implementer final response must include

```text
M1 Phase 5 status: PASS / PASS WITH WARNINGS / FAIL
Files created
Files modified
Commands run and results
Vector computation evidence (numerical values, non-trivial)
Per-state status
Missing-source FAIL behavior verified
3GT8_raw not promoted
Old workflow protection
Acceptance closure: M1 §23 #12 closed
Known limitations
```
