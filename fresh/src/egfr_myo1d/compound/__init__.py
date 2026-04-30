"""Milestone 3 compound-stage helpers.

Implemented tasks stay separated by scope:
M3-T0 readiness, M3-T1 ligand QC, M3-T2 ligand preparation, and M3-T3
receptor PDBQT plus docking-box validation. M3-T4 is limited to one technical
Vina smoke command. M3-T5 plans focused docking jobs and PBS wrappers without
submitting qsub or executing Vina. M3-T6 collects expected Vina outputs into
raw pose/completion tables. M3-T7 performs pose-level scalar geometry
attribution. M3-T8 clusters hard-gate-pass poses for convergence readiness.
M3-T9 aggregates cluster support across compounds/states/pockets without
scoring, ranking, tiering, broad scans, or candidate claims. M3-T10 integrates
M2/M3 evidence into cautious candidate-hypothesis tiers while keeping hard gates
separate from soft scores and affinity descriptive-only. No module in this
package should silently run production docking jobs or expose private compound
data.
"""
