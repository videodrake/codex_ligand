from __future__ import annotations

import csv
import json
from pathlib import Path

from egfr_pipeline.vina import bootstrap as bootstrap_module
from egfr_pipeline.vina.compare import compare_from_config
from egfr_pipeline.vina.parse_poses import build_pose_table_from_config


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _make_pose_block(affinity: float, x: float, y: float, z: float) -> str:
    return "\n".join(
        [
            "MODEL 1",
            f"REMARK VINA RESULT: {affinity:.3f} 0.000 0.000",
            f"ATOM      1  C1  LIG A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C",
            "ENDMDL",
        ]
    )


def test_build_pose_table_writes_coverage_and_pose_provenance(tmp_path: Path) -> None:
    project_root = tmp_path / "output" / "postprocess_project"
    raw_pose = project_root / "3GT8_raw" / "173940_blind.pdbqt"
    _write_text(
        raw_pose,
        "\n".join(
            [
                _make_pose_block(-7.5, 10.0, 11.0, 12.0),
                _make_pose_block(-7.0, 11.0, 12.0, 13.0),
                "",
            ]
        ),
    )

    config = {
        "project_name": "postprocess_project",
        "output_root": str(tmp_path / "output"),
        "mode": "blind",
        "vina": {
            "mode": "blind",
            "exhaustiveness": 384,
            "n_poses": 100,
            "energy_range": 5.0,
            "cpu": 8,
            "seed": 20260316,
        },
        "receptors": [{"id": "3GT8_raw", "pdbqt": "input/receptors/3GT8_raw_receptor.pdbqt"}],
        "ligands": [
            {"id": "173940", "pdbqt": "input/ligands/173940_ligand.pdbqt"},
            {"id": "97806", "pdbqt": "input/ligands/97806_ligand.pdbqt"},
        ],
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")

    pose_csv = build_pose_table_from_config(str(config_path))
    coverage_csv = project_root / "vina_postprocess_coverage.csv"

    pose_rows = _read_csv(pose_csv)
    coverage_rows = _read_csv(coverage_csv)

    assert len(pose_rows) == 2
    assert pose_rows[0]["docking_mode"] == "blind"
    assert pose_rows[0]["exhaustiveness"] == "384"
    assert pose_rows[0]["requested_n_poses"] == "100"
    assert pose_rows[0]["energy_range"] == "5.0"
    assert pose_rows[0]["cpu_per_job"] == "8"
    assert pose_rows[0]["pose_source_status"] == "parsed"
    assert pose_rows[0]["vina_seed"] == pose_rows[1]["vina_seed"]

    assert len(coverage_rows) == 2
    by_ligand = {row["ligand_id"]: row for row in coverage_rows}
    assert by_ligand["173940"]["status"] == "parsed"
    assert by_ligand["173940"]["parsed_n_poses"] == "2"
    assert by_ligand["173940"]["coverage_fraction"] == "0.02"
    assert by_ligand["97806"]["status"] == "missing"
    assert by_ligand["97806"]["parsed_n_poses"] == "0"


def test_compare_from_config_uses_configured_cutoffs(tmp_path: Path) -> None:
    project_root = tmp_path / "output" / "compare_project"
    _write_csv(
        project_root / "vina_pocket_table.csv",
        [
            "receptor_id",
            "pocket_id",
            "centroid_x",
            "centroid_y",
            "centroid_z",
            "n_pose",
            "n_ligand",
            "best_affinity",
            "mean_affinity",
            "union_contact_residues",
            "top_residues",
            "centroid_spread_A",
            "affinity_std",
            "affinity_iqr",
            "dominant_ligand_fraction",
            "ligand_pose_entropy",
        ],
        [
            {
                "receptor_id": "3GT8_raw",
                "pocket_id": "P001",
                "centroid_x": "0.0",
                "centroid_y": "0.0",
                "centroid_z": "0.0",
                "n_pose": "5",
                "n_ligand": "1",
                "best_affinity": "-7.5",
                "mean_affinity": "-7.2",
                "union_contact_residues": "A:ALA699;A:ALA700",
                "top_residues": "A:ALA699",
                "centroid_spread_A": "1.0",
                "affinity_std": "0.2",
                "affinity_iqr": "0.1",
                "dominant_ligand_fraction": "1.0",
                "ligand_pose_entropy": "0.0",
            },
            {
                "receptor_id": "EGFR_160-185",
                "pocket_id": "P001",
                "centroid_x": "10.0",
                "centroid_y": "0.0",
                "centroid_z": "0.0",
                "n_pose": "4",
                "n_ligand": "1",
                "best_affinity": "-7.1",
                "mean_affinity": "-6.9",
                "union_contact_residues": "A:ALA699;A:ALA700",
                "top_residues": "A:ALA699",
                "centroid_spread_A": "1.2",
                "affinity_std": "0.3",
                "affinity_iqr": "0.2",
                "dominant_ligand_fraction": "1.0",
                "ligand_pose_entropy": "0.0",
            },
        ],
    )
    _write_csv(
        project_root / "vina_drug_pocket_map.csv",
        [
            "receptor_id",
            "ligand_id",
            "dominant_pocket_id",
            "dominant_pocket_pose_count",
            "dominant_pocket_fraction",
            "best_affinity",
            "best_pose_rank",
            "top_pose_residues",
            "alternative_pockets",
            "is_multimodal_binding",
        ],
        [
            {
                "receptor_id": "3GT8_raw",
                "ligand_id": "173940",
                "dominant_pocket_id": "P001",
                "dominant_pocket_pose_count": "5",
                "dominant_pocket_fraction": "1.0",
                "best_affinity": "-7.5",
                "best_pose_rank": "1",
                "top_pose_residues": "A:ALA699",
                "alternative_pockets": "",
                "is_multimodal_binding": "False",
            },
            {
                "receptor_id": "EGFR_160-185",
                "ligand_id": "173940",
                "dominant_pocket_id": "P001",
                "dominant_pocket_pose_count": "4",
                "dominant_pocket_fraction": "1.0",
                "best_affinity": "-7.1",
                "best_pose_rank": "1",
                "top_pose_residues": "A:ALA699",
                "alternative_pockets": "",
                "is_multimodal_binding": "False",
            },
        ],
    )

    config = {
        "project_name": "compare_project",
        "output_root": str(tmp_path / "output"),
        "postprocess": {
            "comparison_centroid_cutoff": 5.0,
            "same_patch_centroid_cutoff": 8.0,
            "same_patch_jaccard": 0.3,
            "same_patch_overlap": 0.5,
            "keep_chain": False,
        },
    }
    config_path = tmp_path / "compare_config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")

    out_csv = compare_from_config(str(config_path))
    assert _read_csv(out_csv) == []

    config["postprocess"]["comparison_centroid_cutoff"] = 15.0
    config["postprocess"]["same_patch_centroid_cutoff"] = 12.0
    config_path.write_text(json.dumps(config), encoding="utf-8")

    rows = _read_csv(compare_from_config(str(config_path)))
    assert len(rows) == 1
    assert rows[0]["same_patch_candidate"] == "True"


def test_bootstrap_from_config_respects_cluster_params_and_labels_scope(
    tmp_path: Path,
    monkeypatch,
) -> None:
    project_root = tmp_path / "output" / "bootstrap_project"
    _write_csv(
        project_root / "vina_pose_table.csv",
        [
            "receptor_id",
            "ligand_id",
            "pose_rank",
            "affinity",
            "rmsd_lb",
            "rmsd_ub",
            "centroid_x",
            "centroid_y",
            "centroid_z",
            "raw_pose_file",
            "docking_mode",
            "exhaustiveness",
            "requested_n_poses",
            "energy_range",
            "cpu_per_job",
            "vina_seed",
            "pose_source_status",
            "pocket_id",
            "contact_residues",
            "n_contact_residues",
            "contact_distances",
        ],
        [
            {
                "receptor_id": "3GT8_raw",
                "ligand_id": "173940",
                "pose_rank": "1",
                "affinity": "-7.5",
                "rmsd_lb": "0.0",
                "rmsd_ub": "0.0",
                "centroid_x": "0.0",
                "centroid_y": "0.0",
                "centroid_z": "0.0",
                "raw_pose_file": "raw1.pdbqt",
                "docking_mode": "blind",
                "exhaustiveness": "384",
                "requested_n_poses": "100",
                "energy_range": "5.0",
                "cpu_per_job": "8",
                "vina_seed": "1",
                "pose_source_status": "parsed",
                "pocket_id": "",
                "contact_residues": "A:ALA699",
                "n_contact_residues": "1",
                "contact_distances": "A:ALA699:3.5",
            },
            {
                "receptor_id": "3GT8_raw",
                "ligand_id": "173940",
                "pose_rank": "2",
                "affinity": "-7.0",
                "rmsd_lb": "1.0",
                "rmsd_ub": "2.0",
                "centroid_x": "1.0",
                "centroid_y": "1.0",
                "centroid_z": "1.0",
                "raw_pose_file": "raw1.pdbqt",
                "docking_mode": "blind",
                "exhaustiveness": "384",
                "requested_n_poses": "100",
                "energy_range": "5.0",
                "cpu_per_job": "8",
                "vina_seed": "1",
                "pose_source_status": "parsed",
                "pocket_id": "",
                "contact_residues": "A:ALA700",
                "n_contact_residues": "1",
                "contact_distances": "A:ALA700:3.7",
            },
        ],
    )

    calls: list[tuple[float, int, int]] = []

    def fake_assign_pockets(rows, cutoff=8.0, max_iterations=10, min_pocket_size=2):
        calls.append((float(cutoff), int(max_iterations), int(min_pocket_size)))
        for row in rows:
            row["pocket_id"] = "P001"
        return rows

    monkeypatch.setattr(bootstrap_module, "assign_pockets", fake_assign_pockets)

    config = {
        "project_name": "bootstrap_project",
        "output_root": str(tmp_path / "output"),
        "bootstrap": {"n_replicates": 3, "sample_fraction": 0.5, "seed": 7},
        "postprocess": {
            "pocket_cutoff": 6.0,
            "cluster_max_iterations": 15,
            "min_pocket_size": 3,
            "merge_by_residue": False,
        },
    }
    config_path = tmp_path / "bootstrap_config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")

    out_csv = bootstrap_module.bootstrap_from_config(str(config_path))
    rows = _read_csv(out_csv)

    assert calls
    assert all(call == (6.0, 15, 3) for call in calls)
    assert rows[0]["sample_fraction"] == "0.5"
    assert rows[0]["stability_scope"] == "pose_resampling"
