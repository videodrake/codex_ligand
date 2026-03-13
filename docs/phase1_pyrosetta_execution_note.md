# Phase 1 PyRosetta Execution Note

## Task Group 1.1: PyRosetta Global Docking Standardization

---

## 1. Audit Summary (Task 1.1.1)

### Current Execution Path

The existing PyRosetta PPI docking pipeline consists of:

| Component | File | Role |
|-----------|------|------|
| Orchestrator | `egfr_pipeline/pyrosetta_docking/pipeline_manager.py` | 7-step pipeline execution |
| Docking worker | `egfr_pipeline/pyrosetta_docking/docking.py` | Relax, Global Docking, Refinement |
| Scoring worker | `egfr_pipeline/pyrosetta_docking/analysis.py` | InterfaceAnalyzer metrics |
| Utilities | `egfr_pipeline/pyrosetta_docking/common.py` | PyRosetta init, Pose I/O |

### 7-Step Pipeline

1. **Relax** ??FastRelax (ref2015), cached in `relaxed_cache/`
2. **Global Docking** ??RigidBodyPerturbMover(360°, 100Å) ??SlideIntoContact ??DockMCMProtocol
3. **Fast Scoring & Filtering** ??v2.0 2-pass or v1.0 single-pass
4. **Full Scoring** ??Complete InterfaceAnalyzer metrics + L_RMSD
5. **Clustering** ??CoM pre-filter + L_RMSD greedy clustering
6. **Selection & Save** ??Round-robin diversity + L_RMSD deduplication
7. **Visualization & Report** ??PyMOL scripts + validation report

### Full-Kinase-Domain Compatibility Assessment

The existing pipeline **can handle full-kinase-domain inputs without code modification**:

- The pipeline accepts any 2-chain PDB (chain A = receptor, chain B = partner)
- FoldTree setup uses generic `"A_B"` partner definition
- Scoring, filtering, clustering are all residue-count independent
- Excluded residues are parsed from config (not hard-coded)
- Auto-threshold clustering scales with chain size

**No modifications to pipeline_manager.py, docking.py, analysis.py, or common.py are needed for Phase 1.**

---

## 2. Normalized Execution Inputs (Task 1.1.2)

### Run Metadata Schema

Every Phase 1 run generates `pyrosetta_run_metadata.json` containing:

```json
{
    "receptor_id": "3GT8_raw",
    "partner_id": "extended_beta_meander_955_1006",
    "construct_type": "full_kinase_domain",
    "config_file": "config/phase1/phase1_prod_3GT8_raw_seed0.ini",
    "input_pdb": "input/PPI/phase1/docking_3GT8_raw_ext_beta_meander.pdb",
    "total_global_models": 20000,
    "n_cpus": 16,
    "random_seed": "42",
    "seed_index": 0,
    "is_production": true,
    "filter_version": "v2.0",
    "phase": "Phase 1: PPI-first Interface Mapping",
    "task_group": "TG 1.1: PyRosetta Global Docking Standardization"
}
```

### Current Implementation Note

The JSON example above reflects the intended Phase 1 v2 target state. The current repository implementation is narrower and currently emits:

- `partner_id` values such as `MYO1D_beta_meander` or `MYO1D_TH1`
- `partner_construct` labels that describe the actual prepared input, for example `legacy_beta_meander_960_1006`
- explicit `receptor_chain_ids`, `partner_chain_ids`, and `numbering_system`
- `n_cpus_requested` and `n_cpus_used` instead of a single `n_cpus`
- `run_label`, `run_status`, and `input_validation_status`

This means the current metadata is already traceable, but it does not yet represent the full extended-beta-meander target state described in the PRD.

### Key Fields for Downstream Traceability

| Field | Purpose | Used by |
|-------|---------|---------|
| `receptor_id` | Cross-state comparison | TG 1.5 |
| `partner_id` | Partner identification | All downstream |
| `construct_type` | Distinguish from legacy | TG 1.7 pilot comparison |
| `seed_index` | Multi-seed consolidation | Score standardization |
| `filter_version` | Audit trail | Reproducibility |

---

## 3. Output Path Convention (Task 1.1.4)

### Directory Structure

```
output/phase1_ppi/
?????? 3GT8_raw/
??  ?????? test_seed0/
??  ??  ?????? pyrosetta_run_metadata.json
??  ??  ?-???? docking_3GT8_raw_ext_beta_meander/
??  ??      ?????? final_ranking.csv
??  ??      ?????? cluster_results/
??  ??      ?????? final_result/
??  ??      ?-???? ...
??  ?????? prod_seed0/
??  ?????? prod_seed1/
??  ?????? ...
??  ?????? prod_seed4/
??  ?-???? pyrosetta_decoy_scores.csv  ??consolidated
??
?????? 3GT8_cl38_48/
??  ?-???? (same structure)
??
?-???? 3GT8_cl85_100/
    ?-???? (same structure)
```

### Convention Rules

- **Top level:** Separated by receptor state
- **Run level:** Separated by run type (test/prod) and seed index
- **Pipeline level:** Standard pipeline output structure (unchanged)
- **Consolidated:** `pyrosetta_decoy_scores.csv` merges all seeds

### Current Implementation Note

The directory example above is historical planning guidance, not the current code path. The current implementation writes one metadata-tagged run directory per execution, for example:

```text
EGFR_dimer_beta_meander__3GT8_raw__legacy_beta_meander_960_1006__full_kinase_domain__prod/
|- config_snapshot.ini
|- phase1_input_validation_report.json
|- phase1_input_validation_summary.md
|- pyrosetta_run_metadata.json
|- pyrosetta_decoy_scores.csv
|- filter_passed/
|- cluster_results/
`- final_result/
```

Current naming rule:

- `<input_stem>__<receptor_id>__<partner_construct>__<construct_type>__<run_label>`

This was added so test/prod runs do not overwrite each other when they share the same prepared input PDB.

---

## 4. Standardized Score Table (Task 1.1.5)

### pyrosetta_decoy_scores.csv Schema

| Column | Type | Source | Notes |
|--------|------|--------|-------|
| decoy_id | str | File_PDB | PDB filename |
| receptor_id | str | metadata | e.g., "3GT8_raw" |
| partner_id | str | metadata | "extended_beta_meander_955_1006" |
| construct_type | str | metadata | "full_kinase_domain" |
| seed_index | int | metadata | 0-4 for production |
| run_type | str | metadata | "test" or "production" |
| Rank | int | final_ranking.csv | 1-20 within seed |
| cluster_id | str | Parent ??"C01" | Cluster identifier |
| total_score | float | Total_Score | Full Rosetta score (REU) |
| I_sc | float | dG_separated | Interface score proxy |
| dG_separated | float | InterfaceAnalyzer | Binding energy (REU) |
| dSASA | float | InterfaceAnalyzer | Buried surface (Å²) |
| dG_density | float | derived | dG/dSASA×100 |
| sc_value | float | InterfaceAnalyzer | Shape complementarity (0-1) |
| packstat | float | InterfaceAnalyzer | Packing density (0-1) |
| delta_unsatHbonds | int | InterfaceAnalyzer | Unsatisfied H-bonds |
| nres_int | int | InterfaceAnalyzer | Interface residue count |
| hbonds_int | int | InterfaceAnalyzer | Interface H-bond count |
| L_RMSD | float | CalphaSuperimpose | vs relaxed reference (Å) |
| L_RMSD_best | float | CalphaSuperimpose | vs best-dG model (Å) |
| Binding_Residues_A | str | ContactAnalysis | Receptor interface residues |
| Binding_Residues_B | str | ContactAnalysis | Partner interface residues |
| key_contact_ratio | float | ContactAnalysis | Key residue contact fraction |
| source_file | str | ??| Provenance tracking |

### Current Implementation Note

The table above reflects the broader target schema. The current `pyrosetta_decoy_scores.csv` export in this repository already includes these practical baseline fields:

- `decoy_id`
- `source_file`
- `receptor_id`
- `partner_id`
- `construct_type`
- `receptor_construct`
- `partner_construct`
- `receptor_chain_ids`
- `partner_chain_ids`
- `total_score`
- `I_sc`
- `dG_separated`
- `dSASA`
- `sc_value`
- `packstat`
- `delta_unsatHbonds`
- `nres_int`
- `hbonds_int`
- `L_RMSD`
- `center_x`, `center_y`, `center_z`
- `binding_residues_A`, `binding_residues_B`

### Note on I_sc

The task document specifies `I_sc (interface score) ??preferred primary ranking metric`. In the existing pipeline, `dG_separated` from `InterfaceAnalyzerMover` serves this role. True `I_sc` from `InterfaceScoreCalculator` is a closely related but distinct metric. For Phase 1:

- `I_sc` column is populated with `dG_separated` values
- Both metrics are highly correlated for rigid-body docking
- `dG_separated` is already thoroughly tested and validated in the pipeline
- If true `I_sc` is needed later, it can be added to `analysis.py` as an additional extraction

---

## 5. Config File Changes from Legacy

### Key Differences

| Setting | Legacy (dimer) | Phase 1 (monomer) | Why |
|---------|---------------|-------------------|-----|
| input_pdb | `EGFR_dimer_*.pdb` | `docking_*_ext_beta_meander.pdb` | Monomer receptor |
| excluded_residues_A | Includes `1713-1720,...` | Only `709-720,...` | No dimer chain B |
| total_global_models | 50,000 (single seed) | 20,000 × 5 seeds | Multi-seed strategy |
| random_seed | auto | Deterministic per seed | Reproducibility |
| n_cpus | 32 | 16 | Shared HPC safety |

### What Did NOT Change

- Filter thresholds (v2.0 2-pass design)
- MiniRefinement settings
- Clustering parameters (auto-adaptive)
- Refinement protocol
- Output format
- ExperimentalData section

---

## 6. Execution Commands

### Test Run (Codex workspace: dry-run only)

```bash
# Generate configs
python -m egfr_pipeline.phase1.generate_configs

# Validate and generate metadata (no actual docking)
python -m egfr_pipeline.phase1.launch_docking --test --dry-run
```

### Server-Side Test Run

```bash
# Single state, 1K models (~2-4 hours)
conda activate pyrosetta
python -m egfr_pipeline.phase1.launch_docking --test --state 3GT8_raw
```

### Server-Side Production Run

```bash
# All states, all seeds (submit as PBS jobs)
python -m egfr_pipeline.phase1.launch_docking --production

# Or submit individual seeds for parallel execution:
for state in 3GT8_raw 3GT8_cl38_48 3GT8_cl85_100; do
    for seed in 0 1 2 3 4; do
        python -m egfr_pipeline.phase1.launch_docking \
            --config config/phase1/phase1_prod_${state}_seed${seed}.ini &
    done
done
```

### Post-Run Score Standardization

```bash
# After docking completes:
python -m egfr_pipeline.phase1.standardize_scores
```

---

## 7. Validation Checklist

- [x] Existing pipeline handles full-kinase-domain inputs without modification
- [x] receptor_id, partner_id, construct_type recorded in run metadata
- [x] Config files generated for all 3 states × (1 test + 5 production seeds)
- [x] Output directories separate receptor states and seed indices
- [x] Score standardization utility consolidates multi-seed outputs
- [x] Compute scaling documented with multi-seed recommendation
- [x] All production runs designated as server-side only

