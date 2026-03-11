# PRD — Phase 1→2 Gate: PPI Complex MD Validation

## Version: 1.0
## Date: 2026-03
## Status: Design complete, pending full-kinase-domain docking results

---

## 1. Why This Gate Exists

Rosetta docking produces static rigid-body models. A high-scoring, orientation-validated docking pose does not guarantee that the predicted PPI interface is dynamically stable. Specific failure modes that static docking cannot detect include:

1. **Interface collapse:** The beta-meander detaches from the receptor within nanoseconds because the static pose occupied a local energy minimum that disappears under thermal fluctuation.

2. **Hotspot contact loss:** Key contacts identified by alanine scanning (e.g., VAL962, ILE925) break during MD, revealing that the static contact geometry was strained and unsustainable.

3. **Face reorientation:** The beta-meander rotates on the receptor surface during simulation, adopting a face-flipped or edge-on orientation that the static orientation filter could not predict.

4. **Artifactual stability:** Rosetta's ref2015 scoring function overestimates certain interactions (e.g., due to rigid backbone) that relax away under explicit-solvent MD with backbone flexibility.

5. **Frustrated interface instability:** C04's "frustrated interface" (hotspot sum > WT dG) suggests competing attractive and repulsive forces. MD reveals whether this frustration causes oscillation, partial detachment, or stable frustrated binding.

Without MD validation, committing to Phase 2 pocket proposal based solely on static docking evidence carries unacceptable risk of building the entire downstream pipeline on an unstable interface definition.

---

## 2. Scope and Boundaries

### In scope
- All-atom explicit-solvent MD of the top 2 PPI complex representatives from the new full-kinase-domain docking
- Quantitative assessment of interface dynamic stability
- Go/no-go decision for Phase 2 entry

### Out of scope
- Membrane-embedded MD (that is a separate project track; see EGFR_Symmetric_Inactive_Dimer_MD_Guideline.md)
- Free energy perturbation or enhanced sampling methods
- MD of all cluster representatives (only top 2 are required at this gate)
- Ligand-bound or small-molecule co-simulation

### Relationship to existing MD work
- The existing GROMACS MD pipeline (md_173940_2/, charmm_test_02/) was designed for Vina ligand–EGFR complexes, not PPI complexes. The system setup, topology generation, and analysis scripts require adaptation for a protein–protein complex.
- The symmetric dimer membrane MD (130 ns, ongoing) is a separate track that provides membrane context but does not validate specific PPI poses.

---

## 3. System Definition

### 3.1 Complex structures to simulate

Simulate the top 2 cluster representative poses from the new full-kinase-domain × extended-beta-meander docking. These should be selected by:

1. **Rank 1:** Highest-confidence cluster — best combination of population, I_sc, dG_separated, and orientation score among the orientation-validated poses.
2. **Rank 2:** Second-highest cluster that occupies a **biologically distinct receptor surface region** from Rank 1. This ensures that two independent binding hypotheses are tested, not two variants of the same site.

Selection must be based entirely on the new docking data. Legacy pilot site names (C02, C04, C07) must not influence cluster selection or bias expectations about which receptor surface regions should appear.

### 3.2 System components

| Component | Description |
|-----------|-------------|
| Receptor | EGFR full kinase domain (~280 residues, chain A) |
| Partner | Extended beta-meander (~955–1006, chain B) |
| Solvent | Explicit water (TIP3P) |
| Ions | 0.15 M NaCl for physiological ionic strength |
| Box | Cubic or dodecahedral, minimum 12 Å padding from any protein atom to box edge |
| Force field | AMBER ff19SB (protein) or CHARMM36m — both are suitable for PPI MD; choose based on existing infrastructure |

### 3.3 Why not the C-lobe fragment system

The C-lobe fragment (45 residues) is too small for meaningful MD. Without the N-lobe providing structural context, the fragment would undergo unrealistic unfolding or backbone fluctuations at its truncated termini. The full kinase domain provides the necessary structural framework for stable simulation.

---

## 4. Simulation Protocol

### 4.1 System preparation

1. **Topology generation:** Use `pdb2gmx` (GROMACS) or `tleap` (AmberTools) to generate topology for the two-chain complex. Ensure both chains are treated as separate molecules with correct chain termini.

2. **Solvation and ionization:**
   - Solvate in a cubic box with ≥ 12 Å padding
   - Add 0.15 M NaCl
   - Verify net charge neutralization

3. **Energy minimization:**
   - Steepest descent, ≤ 50,000 steps or until max force < 1000 kJ/mol/nm
   - No position restraints during minimization

### 4.2 Equilibration

| Phase | Ensemble | Duration | Restraints | Thermostat | Barostat |
|-------|----------|----------|------------|------------|----------|
| EQ-1 | NVT | 1 ns | Heavy atoms 1000 kJ/mol/nm² | V-rescale (τ=0.1 ps) | — |
| EQ-2 | NPT | 1 ns | Heavy atoms 500 kJ/mol/nm² | V-rescale (τ=0.1 ps) | Parrinello-Rahman (τ=2 ps) |
| EQ-3 | NPT | 1 ns | Backbone only 200 kJ/mol/nm² | V-rescale | Parrinello-Rahman |
| EQ-4 | NPT | 1 ns | Interface Cα only 100 kJ/mol/nm² | V-rescale | Parrinello-Rahman |

**Rationale for 4-stage equilibration:** Gradual release of restraints prevents interface disruption artifacts that can occur when all restraints are removed simultaneously. The final stage (interface Cα restraints only) allows side chains and non-interface regions to relax while preserving the docked interface geometry during the transition to unrestrained production.

### 4.3 Production MD

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Duration | **200 ns per complex** (minimum) | PPI interface equilibration requires ≥100 ns; 200 ns provides sufficient post-equilibration sampling window |
| dt | 2 fs | Standard with LINCS H-bond constraints |
| Temperature | 300 K | Physiological |
| Pressure | 1 bar | Standard NPT |
| Thermostat | V-rescale (τ=0.1 ps) | Correct ensemble, well-tested |
| Barostat | Parrinello-Rahman (τ=2 ps) | Correct ensemble for production |
| PME | rcoulomb = 1.2 nm | Standard for AMBER/CHARMM |
| Coordinate output | Every 10 ps (100 frames/ns) | Sufficient for contact analysis |
| Energy output | Every 10 ps | For convergence monitoring |

**Why 200 ns, not 100 ns:** The existing pilot Handoff document suggested 100–200 ns. However, the symmetric dimer membrane MD showed that EGFR kinase domain conformational transitions can require >100 ns (RMSD still rising at 130 ns). For a PPI complex where interface stability is the primary question, 200 ns provides a better chance of observing whether the interface reaches a stable plateau or exhibits progressive deterioration.

### 4.4 Replicate strategy

Run each complex with **3 independent replicas** using different random velocity seeds. This is essential because:
- A single trajectory can be trapped in a kinetic artifact (stable by accident or unstable by accident)
- 3 replicas provide minimum statistical confidence for convergence assessment
- If 2/3 replicas show interface disruption, the pose is classified as unstable regardless of the third

Total simulations: 2 complexes × 3 replicas × 200 ns = **1.2 μs aggregate simulation time**

---

## 5. Analysis Protocol

### 5.1 Global stability metrics

| Metric | Tool | Pass criterion |
|--------|------|----------------|
| Backbone RMSD (each chain separately) | `gmx rms` | Plateau within ±0.15 nm over last 50 ns in ≥2/3 replicas |
| Complex RMSD (both chains aligned to receptor) | `gmx rms` | Plateau within ±0.2 nm over last 50 ns |
| Radius of gyration | `gmx gyrate` | No systematic drift |
| Potential energy | `gmx energy` | No drift; stable fluctuation |

### 5.2 Interface-specific stability metrics (critical)

| Metric | Definition | Pass criterion | Fail criterion |
|--------|------------|----------------|----------------|
| **Interface contact survival** | Fraction of initial docking interface contacts (Cα < 10 Å) maintained at each frame | ≥ 60% contacts maintained on average over last 100 ns in ≥2/3 replicas | < 40% contacts in ≥2/3 replicas |

**Important caveat on contact survival interpretation:** Pilot alanine scanning showed that several high-frequency contact residues (ASN963, SER971, CYS970) are actually energetically unfavorable (stabilizing = ALA substitution improves binding). If MD naturally breaks these unfavorable contacts while maintaining the energetically important ones (VAL962, ILE925 etc.), contact survival may drop below 60% even though the interface is becoming more stable. Therefore, contact survival should be interpreted jointly with hotspot contact persistence. If overall contact survival drops but hotspot persistence remains high (≥70%), the complex should be classified as STABLE (interface optimization), not METASTABLE.
| **Hotspot contact persistence** | Per-residue contact occupancy for alanine-scanning-identified hotspot residues | Top 3 hotspot residues maintain contact ≥ 70% of last 100 ns | Top 3 hotspot residues contact < 40% |
| **Interface dSASA stability** | Buried surface area at interface over time | dSASA remains ≥ 60% of initial value over last 100 ns | dSASA drops to < 40% of initial |
| **Inter-chain minimum distance** | Closest heavy-atom distance between chains | < 5 Å continuously | > 8 Å for > 10 ns (partial detachment) |
| **Interface H-bond count** | Number of inter-chain hydrogen bonds over time | Stable or fluctuating around initial count | Progressive decline to < 30% of initial |

### 5.3 Orientation stability

| Metric | Definition | Pass criterion |
|--------|------------|----------------|
| **Sheet 8/9 orientation score over time** | Apply the same PCA-based orientation metric from the static filter at each frame | Score remains positive (active face toward receptor) in ≥ 90% of frames in last 100 ns |
| **Beta-meander structural integrity** | RMSD of beta-meander Cα relative to docking starting structure, aligned to itself | < 0.3 nm (sheet structure preserved) |

### 5.4 Energetic assessment (informative, not gating)

| Metric | Tool | Notes |
|--------|------|-------|
| **MM-PBSA binding energy** | `gmx_MMPBSA` or `g_mmpbsa` | Computed over last 50 ns with 100-frame sampling. Not used as a pass/fail criterion because MM-PBSA has known systematic biases for PPI, but provides a quantitative energy comparison between sites. |
| **Per-residue energy decomposition** | MM-PBSA per-residue | Compare against alanine scanning hotspot rankings. Convergence between static Rosetta ΔΔG and dynamic MM-PBSA per-residue contributions strengthens confidence. |

**Why MM-PBSA is informative but not gating:** MM-PBSA for protein–protein complexes has larger systematic errors (~5–10 kcal/mol) than for protein–ligand systems, and is sensitive to snapshot selection, dielectric constants, and entropy treatment. It is useful for relative ranking between sites (C02-equivalent vs C04-equivalent) but should not be used as an absolute stability criterion.

---

## 6. Go/No-Go Decision Framework

### 6.1 Classification per complex

Based on the analysis metrics above, each complex is classified as:

| Classification | Criteria | Consequence |
|----------------|----------|-------------|
| **STABLE** | All critical interface metrics PASS in ≥2/3 replicas | Proceed to Phase 2 using this patch as receptor-side reference |
| **METASTABLE** | Mixed results: some metrics pass, some fail, or 1/3 replicas unstable | Proceed with caution; flag reduced confidence in Phase 2 patch reference |
| **UNSTABLE** | Critical interface metrics FAIL in ≥2/3 replicas | Do NOT use this patch for Phase 2. Investigate why (wrong site? need refinement?) |

### 6.2 Overall Phase 2 gate decision

| Scenario | Decision |
|----------|----------|
| At least 1 complex is STABLE | Phase 2 can proceed with that patch |
| Both complexes are METASTABLE | Phase 2 can proceed with reduced confidence; consider additional docking/refinement |
| Both complexes are UNSTABLE | Phase 2 blocked. Return to Phase 1 for re-evaluation (different cluster, refined docking, or reassessment of receptor state) |

### 6.3 What "proceed with reduced confidence" means

If the Phase 2 patch reference is based on a METASTABLE complex, downstream Phase 2–4 outputs must:
- Carry a confidence flag indicating MD support is partial
- Not be used for strong mechanistic claims without additional validation
- Be flagged for priority MD re-evaluation if Phase 2 pocket candidates are found

---

## 7. GROMACS Implementation Notes

### 7.1 Topology for two-chain PPI complex

The PPI complex differs from the existing ligand–EGFR setup because both partners are proteins. Key differences:

- No GAFF/CGenFF ligand parameterization needed
- `pdb2gmx` handles both chains with standard protein force field
- Use `gmx make_ndx` to create custom groups: `Protein_A`, `Protein_B`, `Protein_AB`, `Water_and_ions`
- `tc-grps` for temperature coupling: `Protein_AB Water_and_ions` (couple both protein chains together)

### 7.2 Interface contact analysis script

A custom analysis script is needed to compute:
- Per-frame contact map between chains (heavy atom distance < 5 Å)
- Contact survival relative to frame 0
- Per-residue contact occupancy
- Interface dSASA per frame

This can be implemented using MDAnalysis or GROMACS built-in tools (`gmx mindist`, `gmx hbond`, `gmx sasa`).

### 7.3 Orientation score time series

The orientation filter from `orientation_filter.py` (Task 1.2A) must be adapted to work on trajectory frames. Extract PDB snapshots at each analysis frame and compute the orientation score, or implement the equivalent vector calculation in MDAnalysis directly.

### 7.4 Server resource estimate

| Item | Estimate |
|------|----------|
| System size | ~80,000 atoms (protein + water + ions) |
| Time per ns (GPU) | ~2–4 hours on single GPU (GROMACS 2024, RTX 3090 class) |
| Time per 200 ns replica | ~17–33 days wall clock |
| Total for 6 replicas | ~100–200 days single-GPU, or ~35–65 days with 3 GPUs |

This is substantial but feasible on the project server within 1–2 months. Parallel execution across replicas and complexes is straightforward.

---

## 8. Deliverables

| Deliverable | Description |
|-------------|-------------|
| `md_ppi_complex_1_rep{1,2,3}/` | Trajectory directories for complex 1 |
| `md_ppi_complex_2_rep{1,2,3}/` | Trajectory directories for complex 2 |
| `md_stability_metrics.csv` | All quantitative metrics per complex per replica |
| `md_contact_survival.csv` | Time-resolved contact survival data |
| `md_hotspot_persistence.csv` | Per-hotspot-residue contact occupancy |
| `md_orientation_timeseries.csv` | Orientation score at each analysis frame |
| `md_mmpbsa_summary.csv` | MM-PBSA results (informative) |
| `md_gate_decision.md` | Formal go/no-go report with classification and rationale |
| `phase1_md_validation_report.md` | Full analysis report |

---

## 9. Known Limitations and Caveats

1. **200 ns may be insufficient for slow conformational transitions.** If the beta-meander needs to undergo a large-scale reorientation to reach the true binding mode, 200 ns may not capture this. Mitigation: 3 replicas increase the chance of sampling; if all 3 show the same stable interface, confidence is high.

2. **Force field dependence.** AMBER ff19SB and CHARMM36m give somewhat different results for PPI dynamics. The project should use one force field consistently and note this as a systematic uncertainty.

3. **No membrane context.** The PPI MD is in bulk solvent, not at the membrane surface where the interaction naturally occurs. Membrane proximity could stabilize or destabilize certain poses (e.g., by restricting approach angles). This is a known simplification; the membrane MD track addresses this separately.

4. **Backbone rigidity in Rosetta vs flexibility in MD.** Rosetta docking uses limited backbone flexibility, so the starting pose may have strained backbone angles that relax during MD equilibration. The 4-stage equilibration protocol with gradual restraint release is designed to minimize this issue, but some interface restructuring during early equilibration is expected and should not be interpreted as instability.

5. **MM-PBSA entropy estimation.** The standard MM-PBSA approach does not properly account for conformational entropy loss upon binding. -TΔS corrections from normal mode analysis are computationally expensive and noisy. Results should be interpreted as relative comparisons, not absolute binding free energies.

---

## 10. Korean Summary (간단 요약)

이 문서는 Phase 1에서 Phase 2로 넘어가기 전에 **top PPI 포즈의 동적 안정성을 MD로 검증**하는 게이트를 정의한다.

핵심:
- 전장 kinase domain × 확장 beta-meander 복합체의 top 2 클러스터를 각 3 replica × 200 ns 시뮬레이션
- **Interface contact survival**, **hotspot contact persistence**, **dSASA 안정성**, **orientation score 시계열** 등 정량적 기준으로 STABLE / METASTABLE / UNSTABLE 판정
- MM-PBSA는 참고용 (gating 기준 아님)
- STABLE이면 Phase 2 진행, 둘 다 UNSTABLE이면 Phase 1 재평가
- 기존 GROMACS 인프라를 PPI 복합체용으로 적응하여 사용
- 예상 총 시뮬레이션: 1.2 μs (2 복합체 × 3 replica × 200 ns)

