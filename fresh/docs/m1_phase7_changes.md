# M1 Phase 7 — Changes

Closes M1 §23 #14 (ligand manifest shell exists without exposing private IDs by default) per handoff §17.

## Files created

```text
fresh/src/egfr_myo1d/ligand/__init__.py
fresh/src/egfr_myo1d/ligand/manifest.py
fresh/tests/test_m1_phase7_ligand_manifest.py
fresh/tests/fixtures/m1_phase7_ligand/Cpd-A.sdf
fresh/tests/fixtures/m1_phase7_ligand/Cpd-B.sdf
fresh/tests/fixtures/m1_phase7_ligand/Cpd-C.sdf
fresh/tests/fixtures/m1_phase7_ligand/compound_id_map.csv  (synthetic INTERNAL_TEST_PLACEHOLDER_*)
fresh/docs/m1_phase7_ligand_manifest.md
fresh/docs/m1_phase7_changes.md
```

## Files modified

```text
fresh/src/egfr_myo1d/cli.py    # added manifest-ligands subparser + handler
```

## Files deleted

None.

## Public API additions

```python
# ligand/manifest.py
SUPPORTED_FORMATS = ("sdf", "mol", "mol2", "pdb")
PRIVATE_MAPPING_REQUIRED_COLUMNS = ("public_id", "internal_id", "notes")
LIGAND_MANIFEST_QC_COLUMNS = [...]   # 8 columns
LigandFileRecord (dataclass)
LigandManifest (dataclass)
load_public_ids(ctx) -> list[str]
load_default_paths(ctx) -> (raw_ligands, private_mapping)
build_ligand_manifest(ctx, public_ids=None, ligands_dir=None, private_mapping_path=None,
                      profile="codex_dev", compound_stage_enabled=False) -> LigandManifest
```

## CLI surface additions

```bash
python -m egfr_myo1d.cli manifest-ligands --run-id RUN \
    [--ligands-dir PATH] [--private-mapping PATH] \
    [--profile codex_dev|hpc_strict] [--mode smoke_env|smoke_input] \
    [--compound-stage-enabled true|false]
```

Total CLI subparsers after Phase 7: 17 (was 16 after Phase 6).

## Acceptance closure

- M1 §23 #14 closed: ligand manifest shell exists, public IDs only in run outputs,
  internal-ID leak detection escalates to FAIL on any substring match,
  4-cell severity matrix (profile × compound_stage_enabled) implemented correctly,
  .gitignore protection of `fresh/data/private/compound_id_map.csv` validated by
  regression test.

## Verification

- 21 new Phase 7 tests pass:
  - 2 config loaders
  - 4 primary cases (public-IDs-only, no leak, sha256, status PASS)
  - 4 severity-matrix cells
  - 3 private-mapping tests (schema valid, absent, invalid)
  - 1 leak-detection FAIL path
  - 2 output containment / phase status / invariant
  - 1 .gitignore protection
  - 3 CLI tests (help, subcommand help, path traversal)
- Total suite: 229 passing (98 prior + 16 P1 + 27 P3 + 26 P4 + 18 P5 + 23 P6 + 21 P7).
- Old workflow files unmodified.

## Out of scope (next phases)

- Phase 8: prepare-inputs orchestrator + M1 integration test (final §15 closure)
- Phase 9: Tasks 4-9 schema realignment
- M3 ligand prep (RDKit/OpenBabel/PDBQT) and docking
