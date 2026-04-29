# Claude M1 Phase 2 Prompt — MYO1D Module Relocation v0.1

You are working in the `fresh/` workflow on branch `claude/task10` (or successor). M1 Phase 1 (cleanup manager) is complete. This is **M1 Phase 2** — a pure relocation/restructure step that aligns the package layout with `milestone1_foundation_codex_handoff_v0_5.md` §4.

## 1. Project context

The current package places MYO1D-specific construct logic under `preparation/`:

```text
fresh/src/egfr_myo1d/preparation/constructs.py    # MYO1D residue annotation, construct prep, terminal-artifact detection
fresh/src/egfr_myo1d/preparation/pdb_writer.py    # MYO1D-aware PDB writer
```

The handoff §4 spec puts MYO1D-specific code under `myo1d/`:

```text
fresh/src/egfr_myo1d/myo1d/__init__.py
fresh/src/egfr_myo1d/myo1d/construct.py
fresh/src/egfr_myo1d/myo1d/qc.py            # Phase 3 will create
fresh/src/egfr_myo1d/myo1d/pdb_writer.py    # implicit in spec; receptor-side writer is separate
```

This phase is a structural move only. **No logic changes.** It unblocks Phase 3 (MYO1D construct + QC).

## 2. Absolute rules

Do not modify the old workflow:

```text
run_production.py
main.py
egfr_pipeline/**
config/**
docs/runbook.md
output/**
results_export/**
```

Maintain Python 2.7.11 / 3.9 syntax compatibility (commit 40336dd pattern).

The two files being moved must arrive at their new location with **byte-identical content** modulo only the module's own self-import statements. Logic, function signatures, constants, docstrings: unchanged.

## 3. Scope

In scope:
- Create `fresh/src/egfr_myo1d/myo1d/__init__.py`
- Move `fresh/src/egfr_myo1d/preparation/constructs.py` → `fresh/src/egfr_myo1d/myo1d/construct.py`
- Move `fresh/src/egfr_myo1d/preparation/pdb_writer.py` → `fresh/src/egfr_myo1d/myo1d/pdb_writer.py`
- Update import statements in `validation/prepared_inputs.py` and `validation/real_inputs.py` from `preparation.constructs` → `myo1d.construct`, `preparation.pdb_writer` → `myo1d.pdb_writer`
- If `preparation/__init__.py` re-exports any of the moved symbols, remove those re-exports
- Add docs `fresh/docs/m1_phase2_myo1d_relocation.md` and `fresh/docs/m1_phase2_changes.md`

Out of scope:
- Any logic change in the moved files
- Adding new functions to `myo1d/construct.py` (Phase 3 does that)
- Adding `myo1d/qc.py` (Phase 3)
- Modifying `preparation/{masks.py, restraints.py}` — those are EGFR-side, stay in preparation/
- New tests (existing 98 tests are the regression net)

## 4. Required CLI behavior

No CLI changes. Existing CLI commands must continue to work identically.

## 5. Files to create / modify / delete

Create:

```text
fresh/src/egfr_myo1d/myo1d/__init__.py       # may be empty or re-export the public symbols moved
fresh/src/egfr_myo1d/myo1d/construct.py      # moved from preparation/constructs.py
fresh/src/egfr_myo1d/myo1d/pdb_writer.py     # moved from preparation/pdb_writer.py
fresh/docs/m1_phase2_myo1d_relocation.md
fresh/docs/m1_phase2_changes.md
```

Modify:

```text
fresh/src/egfr_myo1d/validation/prepared_inputs.py     # update imports only (lines ~14-23)
fresh/src/egfr_myo1d/validation/real_inputs.py         # update imports only
fresh/src/egfr_myo1d/preparation/__init__.py           # remove re-exports of moved symbols if any
```

Delete (after move + import updates verified):

```text
fresh/src/egfr_myo1d/preparation/constructs.py
fresh/src/egfr_myo1d/preparation/pdb_writer.py
```

## 6. Verification protocol

Order of operations:
1. Create `myo1d/__init__.py`, `myo1d/construct.py`, `myo1d/pdb_writer.py` (copies of preparation/ originals).
2. Run `pytest -q fresh/tests` — should still pass (old import path still resolves).
3. Update `validation/prepared_inputs.py` imports to point to `myo1d.*`.
4. Update `validation/real_inputs.py` imports to point to `myo1d.*`.
5. Run `pytest -q fresh/tests` — should still pass.
6. Search for any remaining `preparation.constructs` or `preparation.pdb_writer` references:
   ```bash
   grep -rn "preparation.constructs\|preparation.pdb_writer" fresh/src fresh/tests
   ```
   Should return zero matches.
7. Delete `preparation/constructs.py` and `preparation/pdb_writer.py`.
8. Update `preparation/__init__.py` if it re-exports moved symbols.
9. Run `pytest -q fresh/tests` — should still pass.

If any step fails, stop and revert that step before continuing.

## 7. Tests required

No new test file in this phase. Pass criterion is the **existing 98 tests** remain green throughout the relocation.

If a test fails after the import update, the failure indicates an import path was missed; fix the import, do not modify the test.

## 8. Acceptance commands

```bash
export PYTHONPATH="$PWD/fresh/src:${PYTHONPATH:-}"

# Full suite must pass
pytest -q fresh/tests

# Verify no orphan imports
grep -rn "preparation.constructs" fresh/src fresh/tests   # expect: empty
grep -rn "preparation.pdb_writer" fresh/src fresh/tests   # expect: empty

# Verify new module is importable
python -c "from egfr_myo1d.myo1d import construct, pdb_writer; print('OK', construct.__name__, pdb_writer.__name__)"

# Verify old paths are gone
test ! -f fresh/src/egfr_myo1d/preparation/constructs.py && echo OK
test ! -f fresh/src/egfr_myo1d/preparation/pdb_writer.py && echo OK

# Old workflow protection
git diff --name-only -- run_production.py main.py egfr_pipeline/ config/ docs/runbook.md output/ results_export/   # expect: empty
```

## 9. Final response format

```text
M1 Phase 2 status: PASS / FAIL

Files created:
- fresh/src/egfr_myo1d/myo1d/__init__.py
- fresh/src/egfr_myo1d/myo1d/construct.py
- fresh/src/egfr_myo1d/myo1d/pdb_writer.py
- fresh/docs/m1_phase2_myo1d_relocation.md
- fresh/docs/m1_phase2_changes.md

Files modified:
- fresh/src/egfr_myo1d/validation/prepared_inputs.py (imports only)
- fresh/src/egfr_myo1d/validation/real_inputs.py (imports only)
- fresh/src/egfr_myo1d/preparation/__init__.py (if needed)

Files deleted:
- fresh/src/egfr_myo1d/preparation/constructs.py
- fresh/src/egfr_myo1d/preparation/pdb_writer.py

Test result:
- prior 98 tests still pass after move
- no new tests in this phase

Import audit:
- grep for preparation.constructs / preparation.pdb_writer returns zero matches
- new myo1d.construct / myo1d.pdb_writer importable

Old workflow protection:
- git diff prints nothing for protected paths

Acceptance closure:
- Phase 3 unblocked (myo1d/construct.py now exists for extension)
- Module tree advances toward handoff §4 spec

Known limitations / not implemented by design:
- No new MYO1D functions yet (Phase 3)
- No myo1d/qc.py yet (Phase 3)
- No logic changes
```
