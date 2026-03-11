# Compute Scaling Estimate: Full Kinase Domain vs C-lobe Fragment

## Phase 1 Task Group 1.1.3

---

## System Comparison

| Parameter | Legacy (C-lobe fragment) | Phase 1 (Full kinase domain) |
|-----------|-------------------------|------------------------------|
| Receptor | ~45 residues (C-lobe only) | ~309 residues (N+C-lobe) |
| Partner | 47 residues (960-1006) | 52 residues (955-1006) |
| Total residues | ~90 | ~361 |
| Receptor chain | Dimer (A+B merged, ~600 res) | Monomer (A only, 309 res) |
| System size | ~645 residues | ~361 residues |

**Note:** The legacy system used a dimer (both EGFR chains merged ~600 res) even though only ~45 res were the "active" C-lobe fragment. The full kinase domain monomer is actually **smaller** than the legacy dimer system in total residue count, though the receptor chain itself is larger.

## Scaling Analysis

### RosettaDock Per-Decoy Cost

RosettaDock scoring scales with the number of residue pairs at the interface. The dominant cost comes from:

1. **Rigid-body perturbation:** O(1) — independent of system size
2. **SlideIntoContact:** O(N) — linear scan along docking axis
3. **DockMCMProtocol:** O(N_interface × M_rotamers) — the expensive step
   - Interface residue count scales with receptor surface area
   - Rotamer sampling scales with interface residues
4. **Score evaluation:** O(N × M) for neighbor-based terms

### Estimated Per-Decoy Time

| Component | Legacy dimer | Phase 1 monomer | Ratio |
|-----------|-------------|-----------------|-------|
| SlideIntoContact | ~0.5s | ~0.3s | 0.6× |
| DockMCMProtocol | ~3-5s | ~2-4s | ~0.8× |
| Fast scoring | ~0.5s | ~0.5s | 1.0× |
| **Total per decoy** | **~4-6s** | **~3-5s** | **~0.7-0.8×** |

**Key insight:** Despite the receptor being larger (309 vs 45 res), the total system is **smaller** than the legacy dimer system (361 vs 645 res). The per-decoy cost may actually be **comparable or slightly lower** for the Phase 1 monomer system.

However, the full kinase domain has a much larger receptor surface to explore, which means:
- More of the docking landscape is "empty" (partner lands far from viable binding sites)
- Early rejection rate may be higher (more membrane-proximal surface to avoid)
- Need more total decoys to adequately sample the larger receptor surface

### Surface Area Scaling

| System | Receptor surface area (approx) | Viable binding surface | Sampling efficiency |
|--------|-------------------------------|----------------------|-------------------|
| C-lobe fragment (45 res) | ~3,000 Å² | ~2,000 Å² (67%) | High |
| Full kinase domain (309 res) | ~18,000 Å² | ~4,000 Å² (22%) | Lower |

The full kinase domain has ~6× more surface area but only ~2× more viable binding surface. This means ~3× more decoys are "wasted" on non-productive surface sampling.

## Target Decoy Count

### Recommended: Multi-Seed Strategy

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Seeds per state | 5 | Independent sampling diversity |
| Models per seed | 20,000 | Sufficient per-seed sampling |
| Total per state | 100,000 | Adequate for global blind docking |
| Total (3 states) | 300,000 | Full Phase 1 campaign |

### Why Multi-Seed Over Single Mega-Run

1. **Sampling diversity:** Different random seeds explore different regions of conformation space. A single 100K run may oversample some regions and miss others.

2. **Fault tolerance:** If one run fails or produces poor results, only 20K models are lost. The remaining 4 seeds still provide 80K models.

3. **Parallelization:** Seeds can run simultaneously on different HPC nodes (5 × 16 cores = 80 cores), reducing wall-clock time by ~5×.

4. **Incremental analysis:** Results from each seed can be analyzed as they complete, providing early signal about binding site preferences.

5. **Literature support:** Multi-seed approaches are standard practice in global docking (Gray et al. RosettaDock recommendations: "multiple independent trajectories provide better sampling than a single long trajectory").

## Estimated Compute Time

### Per-Seed Run (20K models, 16 cores)

| Step | Est. time | Notes |
|------|-----------|-------|
| Relax | 10-30 min | Cached after first run |
| Global Docking | 8-16 hours | ~3-5s per decoy × 20K / 16 cores |
| Fast Scoring | 30-60 min | ~1s per survivor |
| Mini Refinement | 30-60 min | Stage 1 survivors only |
| Expensive Scoring | 15-30 min | Stage 2 survivors only |
| Clustering | 5-15 min | Depends on survivor count |
| Refinement | 30-60 min | Top clusters only |
| Final Scoring | 10-20 min | Top 20 models |
| **Total per seed** | **~12-20 hours** | |

### Full Phase 1 Campaign

| Scenario | Wall clock | Total CPU-hours |
|----------|-----------|-----------------|
| Sequential (1 node, 16 cores) | ~15 days | ~3,600 |
| 3-way parallel (3 nodes) | ~5 days | ~3,600 |
| 5-way parallel (5 nodes per state) | ~20 hours | ~3,600 |

**Recommendation:** Run 5 seeds per state in parallel (if HPC resources allow), then analyze results incrementally. Total campaign can complete in ~1-2 days with adequate node availability.

## Comparison with Literature

| Study | System size | Decoy count | Method |
|-------|------------|-------------|--------|
| Gray & Moughon 2006 | ~300 res × 80 res | 100,000 | RosettaDock |
| Chaudhury & Gray 2008 | ~200 res × 150 res | 50,000 | RosettaDock 3.0 |
| Marze et al. 2018 | Various | 10,000-100,000 | RosettaDock 4.0 |
| **This work** | **309 res × 52 res** | **100,000 (5×20K)** | **RosettaDock + v2.0 filter** |

Our target of 100K total decoys per state is well within the standard range for global rigid-body protein-protein docking.

## Risk Factors

1. **Early rejection rate:** Full kinase domain has more membrane-proximal surface, so early rejection may discard 30-50% of models. The 100K target accounts for this.

2. **Cluster diversity:** Larger receptor surface may produce more distinct clusters (binding sites). The auto-threshold clustering should handle this automatically (cluster_top_n = auto scales with chain size).

3. **I/O bottleneck:** 20K PDB strings in memory (~10-20 GB) may approach node memory limits. The pipeline already implements `gc.collect()` memory management.

4. **Activation loop flexibility:** If the activation loop is poorly modeled, it may attract spurious docking poses. This is a structural input issue, not a compute scaling issue.

---

## Summary

The Phase 1 full-kinase-domain monomer system is comparable in computational cost to the legacy dimer system on a per-decoy basis. The main scaling concern is the need for more decoys to adequately sample the larger receptor surface. A multi-seed strategy (5 × 20K = 100K per state) provides adequate sampling with built-in fault tolerance and parallelization benefits.

**All production runs are server-side only.** The Codex workspace cannot validate compute performance.
