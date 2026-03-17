import csv
import json
from pathlib import Path

from egfr_pipeline.report import generate_report


def _write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def test_generate_report_includes_reproducibility_aware_verdict_details(tmp_path: Path) -> None:
    output_root = tmp_path / "output"
    vina_post = output_root / "workflow_a" / "phase4_vina_postprocess"
    ppi_post = output_root / "workflow_a" / "phase3_ppi_postprocess"
    verdict_dir = output_root / "workflow_a" / "phase5_verdict"
    config = {
        "project_name": "report_project",
        "output_root": str(output_root),
        "receptors": [{"id": "3GT8_raw", "pdb": str(tmp_path / "input" / "receptors" / "3GT8_raw.pdb")}],
        "ligands": [{"id": "lig_001"}],
    }
    config_path = tmp_path / "report_config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")

    _write_csv(
        vina_post / "vina_pocket_table.csv",
        [
            "receptor_id", "pocket_id", "centroid_x", "centroid_y", "centroid_z",
            "n_pose", "n_ligand", "best_affinity", "mean_affinity",
            "union_contact_residues", "top_residues",
            "centroid_spread_A", "affinity_std", "affinity_iqr",
            "dominant_ligand_fraction", "ligand_pose_entropy",
        ],
        [[
            "3GT8_raw", "P001", "10.0", "11.0", "12.0",
            "9", "3", "-7.8", "-7.2",
            "A:ALA699;A:GLY700", "A:ALA699;A:GLY700",
            "1.1", "0.2", "0.1", "0.60", "1.10",
        ]],
    )
    _write_csv(
        vina_post / "vina_drug_pocket_map.csv",
        [
            "receptor_id", "ligand_id", "dominant_pocket_id", "dominant_pocket_pose_count",
            "dominant_pocket_fraction", "best_affinity", "best_pose_rank",
            "top_pose_residues", "alternative_pockets", "is_multimodal_binding",
        ],
        [[
            "3GT8_raw", "lig_001", "P001", "9", "1.0", "-7.8", "1",
            "A:ALA699;A:GLY700", "", "False",
        ]],
    )
    _write_csv(
        vina_post / "vina_pocket_comparison.csv",
        [
            "receptor_a", "pocket_a", "receptor_b", "pocket_b", "centroid_dist",
            "residue_jaccard", "residue_overlap_coeff", "shared_residues",
            "n_shared_residues", "residues_only_a", "residues_only_b",
            "n_residues_a", "n_residues_b", "shared_ligands", "n_shared_ligands",
            "n_ligands_a", "n_ligands_b", "affinity_a", "affinity_b", "n_pose_a",
            "n_pose_b", "same_patch_candidate", "centroid_dist_bootstrap_ci",
        ],
        [[
            "3GT8_raw", "P001", "EGFR_160-185", "P010", "6.5",
            "0.55", "0.60", "ALA699;GLY700", "2", "", "",
            "2", "2", "lig_001", "1", "1", "1", "-7.8", "-7.5", "9",
            "8", "true", "6.0-7.0",
        ]],
    )
    _write_csv(
        ppi_post / "ppi_pyrosetta_summary.csv",
        [
            "receptor_id", "partner_id", "source", "construct_type", "orientation_validation_status",
            "n_runs_total", "n_runs_completed", "seed_indices", "n_final_models", "n_clusters",
            "n_interface_residues", "n_nlobe_interface_residues", "n_clobe_interface_residues",
            "top_residues", "best_dg", "mean_dg", "best_dsasa",
        ],
        [[
            "3GT8_raw", "TH1", "pyrosetta_ppi:TH1", "full_kinase_domain", "passed",
            "5", "5", "0;1;2;3;4", "100", "4", "8", "2", "6",
            "ALA699;GLY700", "-12.0", "-10.0", "450.0",
        ]],
    )
    _write_csv(
        ppi_post / "ppi_pyrosetta_residues.csv",
        [
            "receptor_id", "partner_id", "source", "chain", "residue_id", "residue_num",
            "residue_name", "lobe_label", "construct_type", "orientation_validation_status",
            "n_runs_total", "n_runs_supporting", "frac_runs_supporting", "supporting_seed_indices",
            "frequency_final_ranking", "frequency_cluster_summary", "n_models_final_ranking",
            "occupancy", "mean_interface_delta_e", "best_interface_delta_e",
        ],
        [[
            "3GT8_raw", "TH1", "pyrosetta_ppi:TH1", "A", "ALA699", "699",
            "ALA", "c_lobe", "full_kinase_domain", "passed",
            "5", "4", "0.80", "0;1;2;4", "0.50", "0.40", "80",
            "0.65", "-2.5", "-3.1",
        ]],
    )
    _write_csv(
        verdict_dir / "cross_method_agreement.csv",
        [
            "receptor_id", "pocket_id", "n_vina_residues", "n_ppi_residues", "n_shared_residues",
            "jaccard", "overlap_coeff", "shared_residue_list", "ppi_mean_occupancy_of_shared",
            "spatial_dist_A", "spatial_proximity", "closest_ppi_partner", "n_ppi_partners_near",
            "ppi_frac_runs_supporting", "ppi_best_frac_runs_supporting", "ppi_best_interface_delta_e",
            "ppi_shared_partner_count", "vina_best_affinity_kcal", "ppi_best_dg_REU", "agreement_level",
        ],
        [[
            "3GT8_raw", "P001", "2", "2", "2", "0.50", "0.60", "ALA699;GLY700", "0.65",
            "5.5", "adjacent", "TH1", "1", "0.80", "0.80", "-3.1", "1", "-7.8", "-12.0", "strong",
        ]],
    )
    _write_csv(
        verdict_dir / "valid_sites.csv",
        [
            "receptor_id", "pocket_id", "verdict", "confidence_score",
            "vina_quality_score", "ppi_proximity_score", "cross_receptor_score",
            "vina_affinity_pts", "vina_convergence_pts", "vina_stability_pts",
            "vina_diversity_pts", "vina_consensus_pts", "ppi_spatial_pts",
            "ppi_overlap_pts", "ppi_reproducibility_pts", "cross_receptor_pts",
            "cross_receptor_support_pts", "score_denominator", "ppi_data_available",
            "best_affinity", "n_pose", "n_ligand", "dominant_ligand_fraction",
            "ligand_pose_entropy", "spatial_dist_to_ppi", "closest_ppi_partner",
            "n_ppi_partners_near", "n_shared_with_ppi", "ppi_frac_runs_supporting",
            "ppi_best_interface_delta_e", "cross_receptor_support", "cross_receptor_matches",
            "consensus_site_id", "exp_sensitivity", "exp_specificity", "exp_enrichment",
            "exp_rank_impact", "pocket_stability", "evidence_profile", "reason_tags",
            "decision_trace", "reasons",
        ],
        [[
            "3GT8_raw", "P001", "STRONG", "81.5",
            "44.0", "18.0", "19.5",
            "10.0", "12.0", "10.0",
            "12.0", "12.0", "10.0",
            "4.0", "4.0", "19.5",
            "9.5", "100.0", "yes",
            "-7.8", "9", "3", "0.60",
            "1.10", "5.5", "TH1",
            "1", "2", "0.80",
            "-3.1", "0.82", "EGFR_160-185",
            "CS001", "", "", "",
            "", "0.78", "cross_state+multi_ligand+ppi_supported+recurrent+stable", "multi_ligand;ppi_multi_seed;cross_state_recurrent",
            "VINA[aff=10.0,conv=12.0,stab=10.0,div=12.0] PPI[spatial=10.0,overlap=4.0,repro=4.0] CROSS[count=10.0,support=9.5]",
            "affinity=-7.8; n_ligand=3; ppi_repro=0.80; cross_support=0.82",
        ]],
    )

    report_path, combined_path = generate_report(str(config_path))

    report_text = report_path.read_text(encoding="utf-8")
    assert "Vina = affinity + within-run convergence + pose-resampling stability + ligand diversity" in report_text
    assert "Recurrent pockets across receptor states:" in report_text
    assert "Recurrent verdict-supported pockets: 1" in report_text

    with open(combined_path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows
    ala699 = next(row for row in rows if row["residue_id"] == "ALA699")
    assert ala699["vina_best_pocket_stability"] == "0.78"
    assert ala699["vina_cross_receptor_support"] == "0.82"
    assert ala699["vina_evidence_profile"] == "cross_state+multi_ligand+ppi_supported+recurrent+stable"
    assert ala699["ppi_partners"] == "TH1"
