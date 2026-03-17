import csv
import json
from pathlib import Path

from egfr_pipeline.report import format_combined_residue_table
from egfr_pipeline.verdict import DEFAULT_THRESHOLDS, generate_verdict, score_pocket


def test_score_pocket_dampens_convergence_when_stability_is_low():
    pocket = {
        "best_affinity": "-8.2",
        "n_pose": "12",
        "n_ligand": "2",
        "pocket_exists_frac": "0.45",
        "dominant_ligand_fraction": "0.95",
        "ligand_pose_entropy": "0.20",
    }

    total, verdict, reasons, vina_score, _ppi, cross, raw = score_pocket(
        pocket,
        ppi_agreement=None,
        cross_receptor_matches=[],
        cross_support=None,
        thresholds=dict(DEFAULT_THRESHOLDS),
        has_ppi_data=False,
    )

    assert verdict in {"MODERATE", "STRONG", "WEAK"}
    assert total > 0
    assert vina_score > 0
    assert cross == 0
    assert raw["vina_convergence_pts"] == 7.5
    assert raw["vina_stability_pts"] == 3.0
    assert raw["vina_diversity_pts"] == 4.0
    assert "ligand_dominance=0.95" in reasons
    assert raw["evidence_profile"] == "exploratory+multi_ligand"


def test_score_pocket_adds_ppi_reproducibility_and_cross_support():
    pocket = {
        "best_affinity": "-7.6",
        "n_pose": "9",
        "n_ligand": "3",
        "pocket_exists_frac": "0.85",
        "dominant_ligand_fraction": "0.60",
        "ligand_pose_entropy": "1.10",
    }
    ppi_agreement = {
        "spatial_proximity": "adjacent",
        "n_shared_residues": "4",
        "spatial_dist_A": "5.5",
        "ppi_mean_occupancy_of_shared": "0.7",
        "n_ppi_partners_near": "2",
        "ppi_frac_runs_supporting": "0.80",
        "ppi_best_interface_delta_e": "-3.5",
    }

    total, verdict, reasons, vina_score, ppi_score, cross_score, raw = score_pocket(
        pocket,
        ppi_agreement=ppi_agreement,
        cross_receptor_matches=["EGFR_160-185", "EGFR_170-200"],
        cross_support={"support_frac": 0.82},
        thresholds=dict(DEFAULT_THRESHOLDS),
        has_ppi_data=True,
    )

    assert verdict == "STRONG"
    assert total > 70
    assert vina_score > 0
    assert ppi_score == 20.0
    assert cross_score == 30.0
    assert raw["ppi_reproducibility_pts"] == 4.0
    assert raw["cross_receptor_support_pts"] == 10.0
    assert raw["evidence_profile"] == "cross_state+exploratory+multi_ligand+ppi_supported+recurrent+stable"
    assert "ppi_repro=0.80" in reasons
    assert "cross_support=0.82" in reasons


def test_format_combined_residue_table_carries_verdict_annotations():
    pocket_rows = [
        {
            "receptor_id": "3GT8_raw",
            "pocket_id": "P001",
            "best_affinity": "-7.8",
            "union_contact_residues": "A:ALA699;A:GLY700",
            "pocket_exists_frac": "0.75",
        }
    ]
    pyrosetta_residues = [
        {
            "receptor_id": "3GT8_raw",
            "residue_id": "ALA699",
            "occupancy": "0.60",
            "frequency_final_ranking": "0.40",
            "mean_interface_delta_e": "-2.1",
            "frac_runs_supporting": "0.66",
            "partner_id": "TH1",
            "source": "pyrosetta_ppi:TH1",
        }
    ]
    verdict_rows = [
        {
            "receptor_id": "3GT8_raw",
            "pocket_id": "P001",
            "confidence_score": "74.2",
            "cross_receptor_support": "0.80",
            "pocket_stability": "0.75",
            "evidence_profile": "exploratory+stable+cross_state",
        }
    ]

    combined = format_combined_residue_table(
        pocket_rows,
        pyrosetta_residues,
        verdict_rows=verdict_rows,
    )

    assert len(combined) == 2
    ala699 = next(row for row in combined if row["residue_id"] == "ALA699")
    assert ala699["vina_best_affinity"] == -7.8
    assert ala699["vina_best_pocket_stability"] == 0.75
    assert ala699["vina_best_confidence_score"] == 74.2
    assert ala699["vina_cross_receptor_support"] == "0.80"
    assert ala699["vina_evidence_profile"] == "exploratory+stable+cross_state"
    assert ala699["ppi_frac_runs_supporting"] == "0.66"
    assert ala699["ppi_partners"] == "TH1"


def _write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def _write_simple_pdb(path: Path, chain: str = "A") -> None:
    lines = []
    coords = [
        (699, 10.0, 10.0, 10.0),
        (700, 12.0, 10.0, 10.0),
        (701, 30.0, 10.0, 10.0),
        (702, 32.0, 10.0, 10.0),
    ]
    for idx, (resnum, x, y, z) in enumerate(coords, 1):
        lines.append(
            f"ATOM  {idx:5d}  CA  ALA {chain}{resnum:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
        )
    lines.append("END")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_generate_verdict_prefers_stable_recurrent_site_over_pose_heavy_site(tmp_path: Path):
    output_root = tmp_path / "output"
    project_root = output_root / "verdict_project"
    receptors_dir = tmp_path / "input" / "receptors"

    for name in ["3GT8_raw", "EGFR_160-185", "EGFR_170-200"]:
        _write_simple_pdb(receptors_dir / f"{name}.pdb")

    config = {
        "project_name": "verdict_project",
        "output_root": str(output_root),
        "receptors": [
            {"id": "3GT8_raw", "pdb": str(receptors_dir / "3GT8_raw.pdb")},
            {"id": "EGFR_160-185", "pdb": str(receptors_dir / "EGFR_160-185.pdb")},
            {"id": "EGFR_170-200", "pdb": str(receptors_dir / "EGFR_170-200.pdb")},
        ],
        "ligands": [{"id": "lig_001"}, {"id": "lig_002"}, {"id": "lig_003"}],
    }
    config_path = tmp_path / "verdict_config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")

    _write_csv(
        project_root / "vina_pocket_table.csv",
        [
            "receptor_id", "pocket_id", "centroid_x", "centroid_y", "centroid_z",
            "n_pose", "n_ligand", "best_affinity", "mean_affinity",
            "union_contact_residues", "top_residues",
            "dominant_ligand_fraction", "ligand_pose_entropy",
        ],
        [
            ["3GT8_raw", "P001", "11.0", "10.0", "10.0", "8", "3", "-7.4", "-7.0", "A:ALA699;A:ALA700", "A:ALA699;A:ALA700", "0.55", "1.10"],
            ["3GT8_raw", "P002", "31.0", "10.0", "10.0", "20", "1", "-8.3", "-8.0", "A:ALA701;A:ALA702", "A:ALA701;A:ALA702", "1.00", "0.00"],
            ["EGFR_160-185", "P101", "11.5", "10.0", "10.0", "7", "2", "-7.2", "-7.0", "A:ALA699;A:ALA700", "A:ALA699;A:ALA700", "0.60", "0.90"],
            ["EGFR_170-200", "P201", "10.5", "10.0", "10.0", "6", "2", "-7.1", "-6.9", "A:ALA699;A:ALA700", "A:ALA699;A:ALA700", "0.65", "0.80"],
        ],
    )
    _write_csv(
        project_root / "vina_drug_pocket_map.csv",
        [
            "receptor_id", "ligand_id", "dominant_pocket_id", "dominant_pocket_pose_count",
            "dominant_pocket_fraction", "best_affinity", "best_pose_rank",
            "top_pose_residues", "alternative_pockets", "is_multimodal_binding",
        ],
        [
            ["3GT8_raw", "lig_001", "P001", "3", "1.0", "-7.4", "1", "A:ALA699", "", "False"],
            ["3GT8_raw", "lig_002", "P001", "3", "1.0", "-7.2", "1", "A:ALA699", "", "False"],
            ["3GT8_raw", "lig_003", "P001", "2", "1.0", "-7.1", "1", "A:ALA700", "", "False"],
            ["EGFR_160-185", "lig_001", "P101", "4", "1.0", "-7.2", "1", "A:ALA699", "", "False"],
            ["EGFR_160-185", "lig_002", "P101", "3", "1.0", "-7.0", "1", "A:ALA700", "", "False"],
            ["EGFR_170-200", "lig_001", "P201", "3", "1.0", "-7.1", "1", "A:ALA699", "", "False"],
            ["EGFR_170-200", "lig_003", "P201", "3", "1.0", "-6.9", "1", "A:ALA700", "", "False"],
        ],
    )
    _write_csv(
        project_root / "vina_pocket_comparison.csv",
        [
            "receptor_a", "pocket_a", "receptor_b", "pocket_b", "centroid_dist",
            "residue_jaccard", "residue_overlap_coeff", "shared_residues",
            "n_shared_residues", "residues_only_a", "residues_only_b",
            "n_residues_a", "n_residues_b", "shared_ligands", "n_shared_ligands",
            "n_ligands_a", "n_ligands_b", "affinity_a", "affinity_b", "n_pose_a",
            "n_pose_b", "same_patch_candidate", "centroid_dist_bootstrap_ci",
        ],
        [
            ["3GT8_raw", "P001", "EGFR_160-185", "P101", "5.0", "0.60", "0.70", "ALA699;ALA700", "2", "", "", "2", "2", "lig_001;lig_002", "2", "3", "2", "-7.4", "-7.2", "8", "7", "true", "4.5-5.5"],
            ["3GT8_raw", "P001", "EGFR_170-200", "P201", "5.5", "0.55", "0.65", "ALA699;ALA700", "2", "", "", "2", "2", "lig_001;lig_003", "2", "3", "2", "-7.4", "-7.1", "8", "6", "true", "5.0-6.0"],
        ],
    )
    _write_csv(
        project_root / "ppi_pyrosetta_residues.csv",
        [
            "receptor_id", "partner_id", "source", "chain", "residue_id", "residue_num",
            "residue_name", "lobe_label", "construct_type", "orientation_validation_status",
            "n_runs_total", "n_runs_supporting", "frac_runs_supporting", "supporting_seed_indices",
            "frequency_final_ranking", "frequency_cluster_summary", "n_models_final_ranking",
            "occupancy", "mean_interface_delta_e", "best_interface_delta_e",
        ],
        [
            ["3GT8_raw", "TH1", "pyrosetta_ppi:TH1", "A", "ALA699", "699", "ALA", "c_lobe", "full_kinase_domain", "passed", "5", "4", "0.80", "0;1;2;4", "0.50", "0.40", "80", "0.70", "-2.5", "-3.0"],
            ["3GT8_raw", "TH1", "pyrosetta_ppi:TH1", "A", "ALA700", "700", "ALA", "c_lobe", "full_kinase_domain", "passed", "5", "4", "0.75", "0;1;3;4", "0.45", "0.35", "76", "0.60", "-2.3", "-2.8"],
        ],
    )
    _write_csv(
        project_root / "ppi_pyrosetta_summary.csv",
        [
            "receptor_id", "partner_id", "source", "construct_type", "orientation_validation_status",
            "n_runs_total", "n_runs_completed", "seed_indices",
            "n_final_models", "n_clusters", "n_interface_residues",
            "n_nlobe_interface_residues", "n_clobe_interface_residues",
            "top_residues", "best_dg", "mean_dg", "best_dsasa",
        ],
        [[
            "3GT8_raw", "TH1", "pyrosetta_ppi:TH1", "full_kinase_domain", "passed",
            "5", "5", "0;1;2;3;4", "100", "4", "8", "2", "6",
            "ALA699;ALA700", "-12.0", "-10.0", "450.0",
        ]],
    )
    _write_csv(
        project_root / "vina_pocket_bootstrap.csv",
        [
            "receptor_id", "pocket_id", "pocket_exists_frac", "centroid_std_A",
            "affinity_mean", "affinity_std", "affinity_iqr",
            "n_pose_mean", "n_pose_std", "n_replicates", "sample_fraction",
            "stability_scope",
        ],
        [
            ["3GT8_raw", "P001", "0.85", "1.2", "-7.2", "0.3", "0.2", "7.5", "1.0", "200", "0.80", "pose_resampling"],
            ["3GT8_raw", "P002", "0.25", "3.5", "-8.0", "0.5", "0.4", "18.0", "3.0", "200", "0.80", "pose_resampling"],
            ["EGFR_160-185", "P101", "0.80", "1.0", "-7.1", "0.2", "0.2", "7.0", "1.0", "200", "0.80", "pose_resampling"],
            ["EGFR_170-200", "P201", "0.75", "1.1", "-7.0", "0.2", "0.2", "6.0", "1.0", "200", "0.80", "pose_resampling"],
        ],
    )

    _, verdict_path = generate_verdict(str(config_path))

    with open(verdict_path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    rows_3gt8 = [row for row in rows if row["receptor_id"] == "3GT8_raw"]
    assert rows_3gt8[0]["pocket_id"] == "P001"
    assert rows_3gt8[0]["verdict"] == "STRONG"
    assert float(rows_3gt8[0]["confidence_score"]) > float(rows_3gt8[1]["confidence_score"])
    assert rows_3gt8[0]["evidence_profile"] == "cross_state+exploratory+multi_ligand+ppi_supported+recurrent+stable"
    assert rows_3gt8[1]["evidence_profile"] == "exploratory"
    assert float(rows_3gt8[0]["cross_receptor_support"]) > float(rows_3gt8[1]["cross_receptor_support"] or 0.0)
    assert float(rows_3gt8[0]["ppi_frac_runs_supporting"]) > 0.7
