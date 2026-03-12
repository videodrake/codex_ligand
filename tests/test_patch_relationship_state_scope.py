import csv
from pathlib import Path

from egfr_pipeline.phase2.patch_relationship import run_patch_relationship


def _write_csv(path: Path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_receptor_pdb(path: Path, residue_num: int, x: float, y: float, z: float):
    path.write_text(
        (
            f"ATOM      1  CA  ALA A{residue_num:4d}"
            f"    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
            "END\n"
        ),
        encoding="utf-8",
    )


def test_state_specific_hotspots_do_not_leak_between_receptors(tmp_path):
    output_dir = tmp_path / "phase2"
    receptor_dir = tmp_path / "receptors"
    output_dir.mkdir()
    receptor_dir.mkdir()

    _write_csv(
        output_dir / "candidate_pockets.csv",
        ["pocket_id", "receptor_id", "residue_ids", "centroid_x", "centroid_y", "centroid_z"],
        [
            {
                "pocket_id": "PKT_01",
                "receptor_id": "3GT8_cl38_48",
                "residue_ids": "LEU838",
                "centroid_x": "1.000",
                "centroid_y": "1.000",
                "centroid_z": "1.000",
            }
        ],
    )
    _write_csv(
        output_dir / "phase2_patch_reference_normalized.csv",
        [
            "patch_residue_id",
            "residue_num",
            "residue_name",
            "chain",
            "lobe",
            "construct_type",
            "orientation_validated",
            "robustness",
            "n_states",
            "states",
            "max_occupancy",
            "is_hotspot",
            "pyrosetta_support",
            "lightdock_support",
            "method_agreement",
            "confidence",
            "phase1_evidence_source",
            "notes",
        ],
        [
            {
                "patch_residue_id": "LEU838",
                "residue_num": "838",
                "residue_name": "LEU",
                "chain": "A",
                "lobe": "C-lobe",
                "construct_type": "full_kinase_domain",
                "orientation_validated": "True",
                "robustness": "state_specific",
                "n_states": "1",
                "states": "3GT8_raw",
                "max_occupancy": "0.90",
                "is_hotspot": "True",
                "pyrosetta_support": "True",
                "lightdock_support": "False",
                "method_agreement": "pyrosetta_only",
                "confidence": "medium",
                "phase1_evidence_source": "primary",
                "notes": "",
            }
        ],
    )

    _write_receptor_pdb(receptor_dir / "receptor_3GT8_raw.pdb", 838, 1.0, 1.0, 1.0)
    _write_receptor_pdb(receptor_dir / "receptor_3GT8_cl38_48.pdb", 838, 1.0, 1.0, 1.0)

    rel_path, metrics_path, _ = run_patch_relationship(output_dir, receptor_dir)

    with open(rel_path, encoding="utf-8") as handle:
        relationships = list(csv.DictReader(handle))
    with open(metrics_path, encoding="utf-8") as handle:
        metrics = list(csv.DictReader(handle))

    assert relationships[0]["relationship_class"] == "low_relevance_candidate"
    assert relationships[0]["hotspot_overlap_count"] == "0"
    assert metrics[0]["patch_hotspot_ids"] == ""
    assert metrics[0]["patch_centroid_x"] == ""
