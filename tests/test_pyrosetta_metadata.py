"""Unit tests for PyRosetta run metadata helpers."""

import configparser

from egfr_pipeline.pyrosetta_docking.metadata import (
    build_input_validation_report,
    build_input_validation_markdown,
    build_output_root_name,
    build_decoy_score_rows,
    build_run_metadata,
    infer_construct_type,
    infer_partner_id,
    infer_receptor_id,
)


def _make_config() -> configparser.ConfigParser:
    config = configparser.ConfigParser()
    config.read_dict(
        {
            "Docking": {"total_global_models": "1000"},
            "Output": {"save_top_n": "20", "save_filter_max": "200"},
            "Cluster": {"threshold": "auto", "cluster_top_n": "auto"},
        }
    )
    return config


def test_infer_ids_from_input_name() -> None:
    config = _make_config()
    input_pdb = "input/PPI/prepared/3GT8_raw_beta_meander_model.pdb"
    assert infer_receptor_id(config, input_pdb) == "3GT8_raw"
    assert infer_partner_id(config, input_pdb) == "MYO1D_beta_meander"


def test_construct_type_respects_explicit_metadata() -> None:
    config = _make_config()
    config.read_dict({"Metadata": {"construct_type": "full_kinase_domain"}})
    assert (
        infer_construct_type(config, "input/PPI/prepared/EGFR_dimer_beta_meander_wt.pdb")
        == "full_kinase_domain"
    )


def test_build_run_metadata_contains_traceability_fields() -> None:
    config = _make_config()
    config.read_dict(
        {
            "Metadata": {
                "receptor_id": "EGFR_160-185",
                "partner_id": "MYO1D_beta_meander",
                "construct_type": "full_kinase_domain",
                "receptor_construct": "full_kinase_domain",
                "partner_construct": "extended_beta_meander",
                "numbering_system": "pdb",
            }
        }
    )
    metadata = build_run_metadata(
        config=config,
        config_file="config/ppi_test_beta_meander.ini",
        input_pdb="input/PPI/prepared/EGFR_dimer_beta_meander_wt.pdb",
        root_dir="EGFR_dimer_beta_meander_wt",
        filename="EGFR_dimer_beta_meander_wt.pdb",
        requested_cpus=16,
        used_cpus=8,
        master_seed=1234,
        dir_filter="EGFR_dimer_beta_meander_wt/filter_passed",
        dir_cluster="EGFR_dimer_beta_meander_wt/cluster_results",
        dir_final="EGFR_dimer_beta_meander_wt/final_result",
    )
    assert metadata["receptor_id"] == "EGFR_160-185"
    assert metadata["partner_id"] == "MYO1D_beta_meander"
    assert metadata["construct_type"] == "full_kinase_domain"
    assert metadata["numbering_system"] == "pdb"
    assert metadata["receptor_chain_ids"] == "A"
    assert metadata["partner_chain_ids"] == "B"
    assert metadata["root_dir_mode"] == "legacy_stem"
    assert metadata["run_label"] == "run"
    assert metadata["n_cpus_requested"] == 16
    assert metadata["n_cpus_used"] == 8


def test_build_output_root_name_metadata_tagged() -> None:
    config = _make_config()
    config.read_dict(
        {
            "Metadata": {
                "receptor_id": "3GT8_raw",
                "construct_type": "full_kinase_domain",
                "partner_construct": "legacy_beta_meander_960_1006",
            },
            "Output": {
                "root_dir_mode": "metadata_tagged",
                "run_label": "prod",
            },
        }
    )
    root_name = build_output_root_name(
        config,
        "input/PPI/prepared/EGFR_dimer_beta_meander.pdb",
        "EGFR_dimer_beta_meander",
    )
    assert (
        root_name
        == "EGFR_dimer_beta_meander__3GT8_raw__legacy_beta_meander_960_1006__full_kinase_domain__prod"
    )


def test_build_decoy_score_rows_adds_run_context() -> None:
    run_metadata = {
        "receptor_id": "3GT8_raw",
        "partner_id": "MYO1D_TH1",
        "construct_type": "unknown_construct_type",
        "receptor_construct": "unknown_construct_type",
        "partner_construct": "TH1_domain",
        "receptor_chain_ids": "A",
        "partner_chain_ids": "B",
    }
    fully_scored = [
        {
            "id": "decoy_001",
            "File": "global_0001.pdb",
            "Total_Score": -123.4,
            "dG_separated": -12.3,
            "dSASA": 800.0,
            "sc_value": 0.66,
            "packstat": 0.52,
            "delta_unsatHbonds": 3,
            "nres_int": 18,
            "hbonds_int": 2,
            "L_RMSD": 4.2,
            "center_x": 10.0,
            "center_y": 11.0,
            "center_z": 12.0,
            "Binding_Residues_A": "A:LEU838;A:ASP855",
            "Binding_Residues_B": "B:VAL962;B:LEU964",
        }
    ]
    rows = build_decoy_score_rows(fully_scored, run_metadata)
    assert len(rows) == 1
    assert rows[0]["decoy_id"] == "decoy_001"
    assert rows[0]["receptor_id"] == "3GT8_raw"
    assert rows[0]["partner_id"] == "MYO1D_TH1"
    assert rows[0]["I_sc"] == -12.3
    assert rows[0]["receptor_chain_ids"] == "A"
    assert rows[0]["partner_chain_ids"] == "B"


def test_build_input_validation_report_surfaces_legacy_partner_warning(tmp_path) -> None:
    input_pdb = tmp_path / "mock_beta_meander.pdb"
    input_pdb.write_text(
        "\n".join(
            [
                "ATOM      1  N   GLY A 699       0.000   0.000   0.000  1.00  0.00           N",
                "ATOM      2  CA  GLY A 700       0.000   0.000   0.000  1.00  0.00           C",
                "ATOM      3  N   VAL B 960       0.000   0.000   0.000  1.00  0.00           N",
                "ATOM      4  CA  VAL B 961       0.000   0.000   0.000  1.00  0.00           C",
                "END",
            ]
        ),
        encoding="utf-8",
    )
    (tmp_path / "mock_beta_meander_info.json").write_text(
        '{"partner_residues": "960-1006"}', encoding="utf-8"
    )
    (tmp_path / "mock_beta_meander_mapping.csv").write_text(
        "\n".join(
            [
                "new_chain,new_resnum,original_chain,original_resnum",
                "A,699,A,699",
                "A,700,A,700",
                "B,960,B,960",
                "B,961,B,961",
            ]
        ),
        encoding="utf-8",
    )

    config = _make_config()
    config.read_dict(
        {
            "Metadata": {
                "receptor_id": "3GT8_raw",
                "partner_id": "MYO1D_beta_meander",
                "construct_type": "full_kinase_domain",
                "receptor_construct": "full_kinase_domain",
                "partner_construct": "legacy_beta_meander_960_1006",
                "receptor_chain_ids": "A",
                "partner_chain_ids": "B",
                "numbering_system": "prepared_model_numbering_699_1007_chainB_offset_1000",
            }
        }
    )
    run_metadata = build_run_metadata(
        config=config,
        config_file="config/ppi_test_beta_meander.ini",
        input_pdb=str(input_pdb),
        root_dir=str(tmp_path / "output"),
        filename=input_pdb.name,
        requested_cpus=16,
        used_cpus=8,
        master_seed=1234,
        dir_filter=str(tmp_path / "output" / "filter_passed"),
        dir_cluster=str(tmp_path / "output" / "cluster_results"),
        dir_final=str(tmp_path / "output" / "final_result"),
    )

    report = build_input_validation_report(
        config=config,
        input_pdb=str(input_pdb),
        run_metadata=run_metadata,
    )
    assert report["status"] == "warn"
    assert report["checks"]["expected_receptor_chains_present"] is True
    assert report["checks"]["expected_partner_chains_present"] is True
    assert report["checks"]["mapping_csv_exists"] is True
    assert report["observed_chains"]["B"]["min_resnum"] == 960
    assert any("legacy beta-meander" in warning for warning in report["warnings"])


def test_build_input_validation_markdown_contains_key_sections() -> None:
    report = {
        "status": "warn",
        "input_pdb": "/tmp/mock_beta_meander.pdb",
        "receptor_id": "3GT8_raw",
        "partner_id": "MYO1D_beta_meander",
        "receptor_construct": "full_kinase_domain",
        "partner_construct": "legacy_beta_meander_960_1006",
        "numbering_system": "prepared_model_numbering_699_1007_chainB_offset_1000",
        "checks": {
            "input_pdb_exists": True,
            "mapping_csv_exists": True,
        },
        "observed_chains": {
            "A": {"min_resnum": 699, "max_resnum": 1007, "n_unique_residues": 309},
            "B": {"min_resnum": 960, "max_resnum": 1006, "n_unique_residues": 47},
        },
        "sidecar": {
            "info_json_exists": True,
            "mapping_csv_exists": True,
            "mapping_summary": {
                "B": {
                    "new_min_resnum": 960,
                    "new_max_resnum": 1006,
                    "original_min_resnum": 960,
                    "original_max_resnum": 1006,
                    "original_chain_ids": ["B"],
                }
            },
        },
        "warnings": ["Partner input is still a legacy beta-meander construct starting at residue 960."],
    }

    markdown = build_input_validation_markdown(report)
    assert "# Phase 1 Input Validation Summary" in markdown
    assert "## Checks" in markdown
    assert "## Observed Chains" in markdown
    assert "## Sidecar Files" in markdown
    assert "## Warnings" in markdown
    assert "legacy beta-meander" in markdown
    assert "Chain B: 960-1006 (47 residues)" in markdown
