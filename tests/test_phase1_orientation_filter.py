import csv

from egfr_pipeline.phase1.orientation_filter import merge_orientation_into_models


def _write_csv(path, rows, fieldnames=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def test_merge_orientation_into_models_uses_seed_index_and_basename_fallback(tmp_path):
    output_dir = tmp_path / "phase1_ppi"
    state_dir = output_dir / "3GT8_raw"

    _write_csv(
        state_dir / "orientation_filter_log.csv",
        [
            {
                "model_id": "Rank01_C01_M01_S-15.20.pdb",
                "receptor_id": "3GT8_raw",
                "seed_index": "0",
                "orientation_score": "0.7832",
                "orientation_class": "pass",
            },
            {
                "model_id": "Rank01_C01_M01_S-15.20.pdb",
                "receptor_id": "3GT8_raw",
                "seed_index": "1",
                "orientation_score": "-0.4521",
                "orientation_class": "fail",
            },
            {
                "model_id": "Rank02_C02_M01_S-12.80.pdb",
                "receptor_id": "3GT8_raw",
                "seed_index": "0",
                "orientation_score": "0.0812",
                "orientation_class": "ambiguous",
            },
        ],
    )
    _write_csv(
        state_dir / "pyrosetta_interface_models.csv",
        [
            {
                "model_id": "Rank01_C01_M01_S-15.20.pdb",
                "seed_index": "0",
                "cluster_id": "C01",
            },
            {
                "model_id": "Rank01_C01_M01_S-15.20.pdb",
                "seed_index": "1",
                "cluster_id": "C02",
            },
            {
                "model_id": "output/phase1_ppi/3GT8_raw/test_seed0/final_result/Rank02_C02_M01_S-12.80.pdb",
                "seed_index": "0",
                "cluster_id": "C03",
            },
        ],
    )

    assert merge_orientation_into_models("3GT8_raw", output_dir) is True

    with open(state_dir / "pyrosetta_interface_models.csv", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert rows[0]["orientation_class"] == "pass"
    assert rows[0]["orientation_score"] == "0.7832"
    assert rows[1]["orientation_class"] == "fail"
    assert rows[1]["orientation_score"] == "-0.4521"
    assert rows[2]["orientation_class"] == "ambiguous"
    assert rows[2]["orientation_score"] == "0.0812"
