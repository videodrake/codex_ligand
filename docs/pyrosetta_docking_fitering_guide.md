# Comprehensive filtering strategies for PyRosetta global docking

**Moving beyond total_score + RMSD is essential for accurate PPI prediction from global docking.** Interface-specific metrics—particularly I_sc (interface score), dSASA, shape complementarity, and unsatisfied buried polar atoms—discriminate true binding modes far more effectively than total energy alone. The validated approach combines multi-stage energy filtering, InterfaceAnalyzerMover-based interface quality assessment, RMSD-based clustering of top decoys, and funnel analysis, applied across **10,000–100,000 decoys** from RosettaDock's global search. This report provides the complete filtering pipeline with specific metric cutoffs, clustering strategies, code patterns, and literature-backed best practices.

---

## The essential scoring metrics beyond total_score

Total_score in Rosetta sums all energy terms—van der Waals, solvation, hydrogen bonding, electrostatics, rotamer probabilities—across the entire complex. For docking, this is problematic because **internal protein energies dominate the signal**, burying interface-specific information in noise. The Rosetta documentation explicitly states that total_score "is often a bad indicator of near-native conformations" for docking. It remains useful only as a coarse filter to discard structurally poor models.

**I_sc (interface score)** is the primary discriminating metric. Defined as the total score of the complex minus the total scores of each partner in isolation, I_sc isolates the energetic contribution of the interface. Good docking decoys typically show **I_sc values of −5 to −10 REU**, though values are system-dependent—large obligate complexes can reach −50 to −100 REU. Positive I_sc values indicate repulsive interfaces and should be immediately discarded.

The full suite of interface metrics available through Rosetta's InterfaceAnalyzerMover provides far richer discrimination:

| Metric | Rosetta term | What it measures | Good values |
|--------|-------------|-----------------|-------------|
| Binding energy | `dG_separated` | ΔG by physical chain separation ± repacking | < −5 REU (system-dependent) |
| Binding energy density | `dG_separated/dSASAx100` | Energy per unit interface area | **< −1.5** |
| Buried surface area | `dSASA_int` | SASA buried upon complexation | **800–2000 Å²** (biological) |
| Shape complementarity | `sc_value` | Lawrence & Colman Sc statistic | **> 0.65** |
| Packing quality | `packstat` | RosettaHoles interface packing | **> 0.65** |
| Unsatisfied H-bonds | `delta_unsatHbonds` | Buried polars lacking H-bond partners | **< 5** |
| Interface H-bonds | `hbonds_int` | Cross-interface hydrogen bonds | ≥ 1; ~1 per 170 Å² interface |
| Interface residues | `nres_int` | Residue count at interface | > 15 for biological PPI |

Among these, the **three most discriminative metrics** for separating biological interfaces from artifacts are dSASA_int (crystal contacts typically < 400–500 Å² vs. > 800 Å² for biological), dG_separated/dSASAx100 (normalizes for interface size), and sc_value (biological interfaces 0.65–0.76 vs. < 0.55 for crystal contacts). Delta_unsatHbonds provides a critical negative filter—each buried unsatisfied polar costs ~1.5 REU and signals a physically unrealistic interface.

---

## The four-stage filtering pipeline

The recommended workflow proceeds through four stages, each progressively enriching for true binding modes. Start by generating **at minimum 10,000 decoys, ideally 100,000** for global docking, as Rosetta's stochastic search requires extensive sampling to adequately cover the rotational and translational space around a protein surface.

**Stage 1: Coarse energy filtering.** Retain the **top 5–10% by total_score** to remove structurally broken models (those with severe clashes, chain breaks, or poor internal geometry). This step is purely defensive—it eliminates garbage without attempting to identify good models. For 100,000 decoys, this yields 5,000–10,000 survivors. Additionally, discard any decoy with I_sc > 0 REU.

**Stage 2: Interface quality filtering.** Run InterfaceAnalyzerMover on surviving decoys and apply multi-metric thresholds. The critical filters are: `dSASA_int` ≥ 800 Å² (eliminates sub-biological interfaces), `dG_separated/dSASAx100` < −1.0 to −1.5 (ensures dense favorable contacts, not just large contact area), `sc_value` > 0.55–0.65 (adequate geometric complementarity), `packstat` > 0.65 (well-packed interface), and `delta_unsatHbonds` < 5 (few energetic penalties from buried unsatisfied polars). The Rosetta community emphasizes that **all thresholds are system-dependent** and should be calibrated against known benchmark complexes of similar size and type before applying to novel targets.

**Stage 3: Clustering.** Cluster the **top 200–500 surviving decoys** by pairwise RMSD. Clustering identifies convergent binding modes—regions of conformational space where independent docking trajectories repeatedly find similar solutions. Large, low-energy clusters represent broad energy basins and carry the highest prediction confidence. The ClusPro team demonstrated that the **30 largest clusters contain at least one near-native structure for 93% of benchmark complexes**.

**Stage 4: Local refinement and re-scoring.** Take 5–10 cluster center structures and perform local refinement docking (`-docking_local_refine`), generating 500–1,000 decoys per cluster center. Re-run InterfaceAnalyzerMover on refined models. Generate score-vs-RMSD funnel plots for each cluster center's refinement ensemble. The cluster that produces the deepest, most well-defined energy funnel is the most confident prediction.

---

## How clustering separates true binding modes from noise

Clustering effectiveness depends critically on the RMSD metric and algorithm chosen. **Interface RMSD (I_rmsd)** is preferred over full Cα RMSD for docking because it focuses on the biologically relevant region—the contact surface—rather than distant loop regions that contribute noise. Ligand RMSD (L_rmsd), computed as the Cα RMSD of the smaller chain after superimposing the larger chain, provides a complementary view of overall binding mode accuracy. Full Cα RMSD can mask important interface differences and is generally inferior for docking analysis.

For clustering algorithms, Rosetta's built-in `energy_based_clustering` application works well for moderate datasets. It uses a greedy leader-follower approach, sorting structures by energy first and assigning the lowest-energy structure as each cluster center. For larger datasets or more flexibility, hierarchical agglomerative clustering with a precomputed pairwise RMSD matrix (using scikit-learn's `AgglomerativeClustering` with `metric='precomputed'`) provides better control. Typical **cluster radii of 5–10 Å** for Cα RMSD or 2–4 Å for I_rmsd work well for protein-protein docking, though these should be adjusted based on system size.

The key interpretive principle is that **high population combined with low energy yields the highest confidence**. A large cluster indicates a broad, accessible energy basin—multiple independent trajectories converge on the same solution, suggesting thermodynamic favorability. Conversely, a single or few low-energy structures that are structurally isolated (small cluster) indicate a narrow energy minimum that is likely an artifact. The Rosetta developers note: "A single or just a few low energy structures that are structurally close indicates a narrow minimum...unlikely to be native-like."

---

## Funnel analysis reveals true binding landscapes

An energy funnel plot (score vs. RMSD) is the single most informative diagnostic for docking success. A **good funnel** shows low-energy decoys tightly clustered at low RMSD, with energy rising progressively as RMSD increases—resembling a funnel or "V" shape. This pattern indicates that the scoring function successfully discriminates the native-like binding mode from alternatives and that sampling has converged. RosettaDock defines success using the **N5 metric**: at least 3 of the 5 lowest-I_sc decoys must have I_rmsd ≤ 4.0 Å to the native structure.

Plot **I_sc (not total_score) on the y-axis** and **I_rmsd on the x-axis** for the most informative funnel. A flat landscape (no funnel) indicates either insufficient sampling or that the scoring function cannot discriminate this particular system—common when significant backbone conformational changes occur between unbound and bound forms. Multiple competing funnels suggest alternative binding modes that may require experimental data to distinguish.

For **blind docking without a known native structure**, use the lowest-energy structure as a pseudo-native reference. Look for structural convergence: large populations of low-energy decoys at similar RMSD values. The **P_near metric** (Bhardwaj et al., Nature, 2016) quantifies funnel quality as a continuous value from 0 (no funnel) to 1 (perfect funnel), computed as a Boltzmann-weighted fraction of near-native structures. For docking with the ref2015 score function, use **k_BT = 1.0** and **λ = 2.0–4.0 Å**. P_near > 0.5 generally indicates good funnel quality.

---

## InterfaceAnalyzerMover extracts the metrics that matter most

InterfaceAnalyzerMover is the workhorse for post-docking interface assessment. It physically separates chains, optionally repacks exposed residues, and computes binding energy as the difference between bound and separated states. Critical setup details: always set `pack_separated=True` for accurate dG_separated calculation, set `pack_input=True` for non-Rosetta input structures, and use `compute_packstat=True` and `compute_interface_sc=True` to obtain packing and shape complementarity metrics (these are off by default due to computational cost).

The practical PyRosetta implementation for batch analysis of docking decoys follows this pattern:

```python
import pyrosetta, glob, pandas as pd
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

pyrosetta.init("-ex1 -ex2aro")
scorefxn = pyrosetta.create_score_function("ref2015")

results = []
for pdb_file in sorted(glob.glob("docking_output_*.pdb")):
    pose = pyrosetta.pose_from_pdb(pdb_file)
    
    ia = InterfaceAnalyzerMover("A_B")
    ia.set_pack_separated(True)
    ia.set_pack_input(True)
    ia.apply(pose)
    
    results.append({
        'filename':           pdb_file,
        'total_score':        scorefxn(pose),
        'interface_dG':       ia.get_interface_dG(),
        'dSASA':              ia.get_interface_delta_sasa(),
        'unsat_hbonds':       ia.get_interface_delta_hbond_unsat(),
        'packstat':           ia.get_interface_packstat(),
        'nres_int':           ia.get_num_interface_residues(),
    })

df = pd.DataFrame(results)
```

For clustering implementation, build a pairwise RMSD matrix and apply hierarchical clustering:

```python
import numpy as np
from pyrosetta.rosetta.core.scoring import CA_rmsd
from sklearn.cluster import AgglomerativeClustering

rmsd_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(i+1, n):
        r = CA_rmsd(poses[i], poses[j])
        rmsd_matrix[i, j] = rmsd_matrix[j, i] = r

clustering = AgglomerativeClustering(
    n_clusters=None, distance_threshold=5.0,
    metric='precomputed', linkage='average'
)
labels = clustering.fit_predict(rmsd_matrix)
```

Key practical tips for large-scale runs: **never store all Pose objects in memory**—write each decoy to disk after scoring and discard it. PyRosetta is inherently serial at the C++ level due to the GIL, so use **process-based parallelism** (GNU Parallel, multiple Python processes, or `pyrosetta.distributed` with Dask) rather than threading. Use `pyrosetta.init("-constant_seed -jran SEED")` with different seeds per process for reproducibility.

---

## Avoiding the most common false positive traps

**Crystal packing artifacts** are the most insidious false positives in docking-based PPI prediction. Crystal contacts typically show dSASA < 400–500 Å², Sc < 0.55, loose packing, and more planar/polar interfaces compared to biological contacts. Tools like PISA, EPPIC, and the HADDOCK interface-classifier can help discriminate these. Evolutionary conservation mapping provides an orthogonal signal—biological interfaces tend to show higher residue conservation than crystal contacts, and coevolution signals between interface residues (detectable via EVcomplex2 or InterEvDock) are strong discriminators.

**Over-reliance on total_score** remains the most common methodological error. The internal energy of each protein partner typically dwarfs the interface energy, meaning that a model with slightly better backbone energetics but a poor interface can outscore a model with an excellent interface. Always use I_sc or dG_separated as the primary ranking metric.

**Insufficient sampling** causes failures when the native-like binding mode is never adequately explored. RosettaDock benchmarks show that global docking with < 10,000 decoys frequently misses the correct solution. For complexes > 450 total residues, the conformational space becomes prohibitively large for adequate global sampling. Backbone conformational changes > 2 Å between unbound and bound forms cause the majority of docking failures—consider RosettaDock 4.0's Adaptive Conformer Selection or ReplicaDock 2.0 for flexible targets.

**Non-specific hydrophobic contacts** can produce artificially favorable I_sc values. The dG_separated/dSASAx100 metric helps detect these—a large hydrophobic patch may have favorable total binding energy but poor energy density. Cross-validate with delta_unsatHbonds and hbonds_int: real biological interfaces almost always feature hydrogen bonds and salt bridges alongside hydrophobic contacts.

---

## What the benchmarks and literature recommend

The **CAPRI quality criteria** (Méndez et al., Proteins, 2003) remain the gold standard for evaluating docking predictions, classifying models as high quality (Fnat ≥ 0.5, I_rmsd ≤ 1.0 Å), medium (Fnat ≥ 0.3, I_rmsd ≤ 2.0 Å), acceptable (Fnat ≥ 0.1, I_rmsd ≤ 4.0 Å), or incorrect. The **DockQ score** (Basu & Wallner, PLoS ONE, 2016) combines these into a single continuous [0, 1] metric that reproduces CAPRI classification with ~94% precision.

RosettaDock performance on the Protein Docking Benchmark has improved substantially across versions. Chaudhury et al. (PLoS ONE, 2011) reported **48% success** on Benchmark 3.0 with RosettaDock v3.2, while Marze et al. (Bioinformatics, 2018) achieved **8-fold higher enrichment** of near-native structures with RosettaDock 4.0's Motif Dock Score (MDS), reaching 77% success on rigid-body targets. The most recent AlphaRED pipeline (Harmalkar et al., eLife, 2024) combines AlphaFold-Multimer with ReplicaDock 2.0, achieving **63% overall success** on Benchmark 5.5 and 43% on the notoriously difficult antibody-antigen category.

Machine learning-based rescoring of docking decoys represents a frontier approach. GNN-DOVE (Wang et al., 2021) selected acceptable models for 49/58 benchmark targets. DeepRank-GNN (Réau et al., Bioinformatics, 2023) achieved AUC 0.71 on the CAPRI Score_set using only geometric and physicochemical features. The most recent comprehensive survey (BMC Bioinformatics, 2024) found that **PIsToN and dMaSIF outperformed** all other scoring functions, including both classical and earlier deep learning approaches.

For users performing general PPI studies, the most impactful single improvement is incorporating **AlphaFold-Multimer as a complementary approach**. When interface-pLDDT > 85, AFm predictions are often sufficient; when interface-pLDDT < 85, this triggers physics-based global docking with Rosetta, as implemented in the AlphaRED protocol.

---

## Conclusion

The transition from simple total_score + RMSD filtering to a multi-metric pipeline represents a qualitative leap in docking discrimination power. The most critical upgrades are: adopting **I_sc as the primary ranking metric** instead of total_score, applying **interface quality filters** (dSASA > 800 Å², Sc > 0.65, packstat > 0.65, delta_unsatHbonds < 5), and using **energy-weighted clustering** to identify convergent binding modes. The binding energy density metric (dG_separated/dSASAx100 < −1.5) is particularly valuable for eliminating large but non-specific interfaces. For blind predictions without a native structure, cluster population combined with funnel analysis using the lowest-energy structure as reference provides the strongest signal. Finally, integrating AlphaFold-Multimer predictions—using interface-pLDDT as a confidence gauge—with physics-based Rosetta refinement now represents the state of the art for PPI structure prediction.