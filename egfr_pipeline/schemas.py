"""Single source of truth for CSV column schemas.

Both the producing modules (Phase outputs) and validate.py import from here.
Until now, the canonical column lists lived only inside validate.py:218-336,
which detected drift after the fact rather than preventing it. By centralizing
the definitions, future producers can also import these lists when writing
headers, so absolute rule #5 ("preserve CSV schemas") is enforced by
construction rather than by post-hoc validation.

This module deliberately depends only on the standard library — keep it that
way to avoid import cycles with csv_utils, paths, etc.
"""
from __future__ import annotations

from typing import Dict, List


VINA_POSE_TABLE: List[str] = [
    "receptor_id", "ligand_id", "pose_rank", "affinity",
    "rmsd_lb", "rmsd_ub", "centroid_x", "centroid_y", "centroid_z",
    "raw_pose_file", "docking_mode", "exhaustiveness", "requested_n_poses",
    "energy_range", "cpu_per_job", "vina_seed", "pose_source_status",
    "pocket_id", "contact_residues", "n_contact_residues", "contact_distances",
]

VINA_POCKET_TABLE: List[str] = [
    "receptor_id", "pocket_id", "centroid_x", "centroid_y", "centroid_z",
    "n_pose", "n_ligand", "best_affinity", "mean_affinity",
    "union_contact_residues", "top_residues",
    "centroid_spread_A", "affinity_std", "affinity_iqr",
    "dominant_ligand_fraction", "ligand_pose_entropy",
]

VINA_DRUG_POCKET_MAP: List[str] = [
    "receptor_id", "ligand_id", "dominant_pocket_id",
    "dominant_pocket_pose_count", "dominant_pocket_fraction",
    "best_affinity", "best_pose_rank", "top_pose_residues",
    "alternative_pockets", "is_multimodal_binding",
]

VINA_POCKET_COMPARISON: List[str] = [
    "receptor_a", "pocket_a", "receptor_b", "pocket_b",
    "centroid_dist", "residue_jaccard", "residue_overlap_coeff",
    "shared_residues", "n_shared_residues",
    "residues_only_a", "residues_only_b",
    "n_residues_a", "n_residues_b",
    "shared_ligands", "n_shared_ligands", "n_ligands_a", "n_ligands_b",
    "affinity_a", "affinity_b", "n_pose_a", "n_pose_b",
    "same_patch_candidate",
    "centroid_dist_bootstrap_ci",
]

PPI_PYROSETTA_RESIDUES: List[str] = [
    "receptor_id", "partner_id", "source", "chain", "residue_id", "residue_num",
    "residue_name", "lobe_label", "construct_type", "orientation_validation_status",
    "n_runs_total", "n_runs_supporting", "frac_runs_supporting", "supporting_seed_indices",
    "frequency_final_ranking", "frequency_cluster_summary",
    "n_models_final_ranking", "occupancy",
    "mean_interface_delta_e", "best_interface_delta_e",
]

PPI_PYROSETTA_SUMMARY: List[str] = [
    "receptor_id", "partner_id", "source", "construct_type", "orientation_validation_status",
    "n_runs_total", "n_runs_completed", "seed_indices",
    "n_final_models", "n_clusters", "n_interface_residues",
    "n_nlobe_interface_residues", "n_clobe_interface_residues",
    "top_residues", "best_dg", "mean_dg", "best_dsasa",
]

PPI_PYROSETTA_RESIDUE_LONG: List[str] = [
    "model_id", "receptor_id", "partner_id", "source", "construct_type",
    "orientation_validation_status", "seed_index", "run_label", "run_dir",
    "rank", "cluster_id", "chain", "residue_id", "residue_num",
    "residue_name", "lobe_label", "delta_e_total", "delta_e_fa_atr",
    "delta_e_fa_rep", "delta_e_fa_sol", "delta_e_fa_elec", "source_file",
]

PPI_PYROSETTA_MODEL_TABLE: List[str] = [
    "model_id", "receptor_id", "partner_id", "source", "construct_type",
    "orientation_validation_status", "seed_index", "run_label", "run_dir",
    "rank", "cluster_id", "dG_separated", "dSASA", "sc_value", "packstat",
    "nres_int", "n_receptor_interface_residues", "n_partner_interface_residues",
    "n_nlobe_interface_residues", "n_clobe_interface_residues",
    "receptor_interface_residues", "partner_interface_residues", "source_file",
]

PPI_AFM_RESIDUES: List[str] = [
    "receptor_id", "source", "residue_id", "residue_num",
    "min_ca_distance",
]

COMBINED_RESIDUE_EVIDENCE: List[str] = [
    "receptor_id", "residue_id", "vina_pockets", "n_vina_pockets",
    "vina_best_affinity", "vina_best_pocket_stability",
    "vina_best_confidence_score", "vina_cross_receptor_support",
    "vina_evidence_profile",
    "ppi_occupancy", "ppi_frequency", "ppi_delta_e",
    "ppi_frac_runs_supporting", "ppi_partners", "evidence_sources",
]

CROSS_METHOD_AGREEMENT: List[str] = [
    "receptor_id", "pocket_id", "n_vina_residues", "n_ppi_residues",
    "n_shared_residues", "jaccard", "overlap_coeff", "shared_residue_list",
    "ppi_mean_occupancy_of_shared", "spatial_dist_A", "spatial_proximity",
    "closest_ppi_partner", "n_ppi_partners_near",
    "ppi_frac_runs_supporting", "ppi_best_frac_runs_supporting",
    "ppi_best_interface_delta_e", "ppi_shared_partner_count",
    "vina_best_affinity_kcal", "ppi_best_dg_REU", "agreement_level",
]

VALID_SITES: List[str] = [
    "receptor_id", "pocket_id", "verdict", "confidence_score",
    "vina_quality_score", "ppi_proximity_score", "cross_receptor_score",
    "vina_affinity_pts", "vina_convergence_pts", "vina_stability_pts",
    "vina_diversity_pts", "vina_consensus_pts",
    "ppi_spatial_pts", "ppi_overlap_pts", "ppi_reproducibility_pts",
    "cross_receptor_pts", "cross_receptor_support_pts",
    "exp_sensitivity_pts", "exp_enrichment_pts", "exp_specificity_pts",
    "exp_score",
    "score_denominator",
    "ppi_data_available", "best_affinity", "n_pose", "n_ligand",
    "dominant_ligand_fraction", "ligand_pose_entropy",
    "spatial_dist_to_ppi", "closest_ppi_partner", "n_ppi_partners_near",
    "n_shared_with_ppi", "ppi_frac_runs_supporting", "ppi_best_interface_delta_e",
    "cross_receptor_matches", "cross_receptor_support", "consensus_site_id",
    "exp_sensitivity", "exp_specificity", "exp_enrichment", "exp_rank_impact",
    "pocket_stability", "evidence_profile", "reason_tags", "decision_trace",
    "reasons",
    "is_atp_site", "exclusion_reason",
]

VINA_CONSENSUS_SITES: List[str] = [
    "consensus_site_id", "n_receptors", "receptor_list", "pocket_list",
    "centroid_x", "centroid_y", "centroid_z",
    "best_affinity", "total_n_ligand", "total_n_pose", "consensus_residues",
]

VINA_POCKET_BOOTSTRAP: List[str] = [
    "receptor_id", "pocket_id", "pocket_exists_frac", "centroid_std_A",
    "affinity_mean", "affinity_std", "affinity_iqr",
    "n_pose_mean", "n_pose_std", "n_replicates", "sample_fraction",
    "stability_scope",
]

VINA_POSTPROCESS_COVERAGE: List[str] = [
    "receptor_id", "ligand_id", "raw_pose_file", "status",
    "parsed_n_poses", "requested_n_poses", "coverage_fraction",
    "docking_mode", "exhaustiveness", "energy_range",
    "cpu_per_job", "base_seed", "vina_seed",
]


ALL_SCHEMAS: Dict[str, List[str]] = {
    "vina_pose_table.csv": VINA_POSE_TABLE,
    "vina_pocket_table.csv": VINA_POCKET_TABLE,
    "vina_drug_pocket_map.csv": VINA_DRUG_POCKET_MAP,
    "vina_pocket_comparison.csv": VINA_POCKET_COMPARISON,
    "ppi_pyrosetta_residues.csv": PPI_PYROSETTA_RESIDUES,
    "ppi_pyrosetta_summary.csv": PPI_PYROSETTA_SUMMARY,
    "ppi_pyrosetta_residue_long.csv": PPI_PYROSETTA_RESIDUE_LONG,
    "ppi_pyrosetta_model_table.csv": PPI_PYROSETTA_MODEL_TABLE,
    "ppi_afm_residues.csv": PPI_AFM_RESIDUES,
    "combined_residue_evidence.csv": COMBINED_RESIDUE_EVIDENCE,
    "cross_method_agreement.csv": CROSS_METHOD_AGREEMENT,
    "valid_sites.csv": VALID_SITES,
    "vina_consensus_sites.csv": VINA_CONSENSUS_SITES,
    "vina_pocket_bootstrap.csv": VINA_POCKET_BOOTSTRAP,
    "vina_postprocess_coverage.csv": VINA_POSTPROCESS_COVERAGE,
}
