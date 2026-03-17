"""Fixture-based validation smoke test for pre-qsub checks."""

import json
from pathlib import Path

from egfr_pipeline.validate import run_validation


POSE_HEADER = [
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
]

POCKET_HEADER = [
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
]

MAP_HEADER = [
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
]


def _write(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def _make_pdb(path: Path, chain: str, start: int = 699, count: int = 4) -> None:
    lines = []
    for idx in range(count):
        resnum = start + idx
        x = 10.0 + idx
        lines.append(
            f"ATOM  {idx + 1:5d}  CA  ALA {chain}{resnum:4d}    "
            f"{x:8.3f}{10.0:8.3f}{10.0:8.3f}  1.00 20.00           C"
        )
    lines.append("END")
    _write(path, lines)


def test_run_validation_smoke(tmp_path: Path) -> None:
    output_root = tmp_path / "output"
    project_name = "smoke_project"
    vina_post = output_root / "workflow_a" / "phase4_vina_postprocess"
    receptors_dir = tmp_path / "input" / "receptors"
    raw_pose = vina_post / "raw_pose_1.pdbqt"

    _make_pdb(receptors_dir / "3GT8_raw.pdb", "A")
    _make_pdb(receptors_dir / "EGFR_160-185.pdb", "A")
    _make_pdb(receptors_dir / "EGFR_170-200.pdb", "A")
    _write(raw_pose, ["MODEL 1", "ENDMDL"])

    _write(
        vina_post / "vina_pose_table.csv",
        [
            ",".join(POSE_HEADER),
            ",".join(
                [
                    "3GT8_raw",
                    "lig_001",
                    "1",
                    "-7.2",
                    "0.0",
                    "0.0",
                    "10.0",
                    "11.0",
                    "12.0",
                    str(raw_pose),
                    "blind",
                    "384",
                    "100",
                    "5.0",
                    "8",
                    "20260316",
                    "parsed",
                    "3GT8_raw_PKT01",
                    "A:ALA699;A:ALA700",
                    "2",
                    "3.5;4.1",
                ]
            ),
        ],
    )
    _write(
        vina_post / "vina_pocket_table.csv",
        [
            ",".join(POCKET_HEADER),
            ",".join(
                [
                    "3GT8_raw",
                    "3GT8_raw_PKT01",
                    "10.0",
                    "11.0",
                    "12.0",
                    "1",
                    "1",
                    "-7.2",
                    "-7.2",
                    "A:ALA699;A:ALA700",
                    "A:ALA699;A:ALA700",
                    "0.0",
                    "0.0",
                    "0.0",
                    "1.0",
                    "0.0",
                ]
            ),
        ],
    )
    _write(
        vina_post / "vina_drug_pocket_map.csv",
        [
            ",".join(MAP_HEADER),
            ",".join(
                [
                    "3GT8_raw",
                    "lig_001",
                    "3GT8_raw_PKT01",
                    "1",
                    "1.0",
                    "-7.2",
                    "1",
                    "A:ALA699;A:ALA700",
                    "",
                    "False",
                ]
            ),
        ],
    )

    config = {
        "project_name": project_name,
        "output_root": str(output_root),
        "receptors": [
            {"id": "3GT8_raw", "pdb": str(receptors_dir / "3GT8_raw.pdb")},
            {"id": "EGFR_160-185", "pdb": str(receptors_dir / "EGFR_160-185.pdb")},
            {"id": "EGFR_170-200", "pdb": str(receptors_dir / "EGFR_170-200.pdb")},
        ],
        "ligands": [{"id": "lig_001"}],
    }
    config_path = tmp_path / "smoke_config.json"
    with open(config_path, "w", encoding="utf-8") as handle:
        json.dump(config, handle)

    result = run_validation(str(config_path), repo_root=str(Path(__file__).resolve().parents[1]))

    assert result.exit_code in (0, 1), result.summary()
    assert not result.failures, result.summary()
