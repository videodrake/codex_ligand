"""Checks for PyRosetta config traceability and input integrity."""

from __future__ import annotations

import configparser
from pathlib import Path

import pytest

from run_production import PPI_TARGETS, _ppi_docking_dir


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    ("config_name", "partner_id", "partner_construct"),
    [
        ("ppi_test_beta_meander.ini", "MYO1D_beta_meander", "legacy_beta_meander_960_1006"),
        ("ppi_prod_beta_meander.ini", "MYO1D_beta_meander", "legacy_beta_meander_960_1006"),
        ("ppi_test_TH1.ini", "MYO1D_TH1", "TH1_domain_801_1006"),
        ("ppi_prod_TH1.ini", "MYO1D_TH1", "TH1_domain_801_1006"),
    ],
)
def test_pyrosetta_configs_have_explicit_metadata(
    config_name: str, partner_id: str, partner_construct: str
) -> None:
    config = configparser.ConfigParser()
    config_path = REPO_ROOT / "config" / config_name
    config.read(config_path, encoding="utf-8")

    assert config.has_section("Metadata")
    assert config.get("Metadata", "receptor_id") == "3GT8_raw"
    assert config.get("Metadata", "partner_id") == partner_id
    assert config.get("Metadata", "construct_type") == "full_kinase_domain"
    assert config.get("Metadata", "receptor_construct") == "full_kinase_domain"
    assert config.get("Metadata", "partner_construct") == partner_construct
    assert config.get("Metadata", "receptor_chain_ids") == "A"
    assert config.get("Metadata", "partner_chain_ids") == "B"
    assert (
        config.get("Metadata", "numbering_system")
        == "prepared_model_numbering_699_1007_chainB_offset_1000"
    )
    assert config.get("Output", "root_dir_mode") == "metadata_tagged"
    assert config.get("Output", "run_label") in {"test", "prod"}


@pytest.mark.parametrize(
    "config_name",
    [
        "ppi_test_beta_meander.ini",
        "ppi_prod_beta_meander.ini",
        "ppi_test_TH1.ini",
        "ppi_prod_TH1.ini",
    ],
)
def test_pyrosetta_config_inputs_exist(config_name: str) -> None:
    config = configparser.ConfigParser()
    config_path = REPO_ROOT / "config" / config_name
    config.read(config_path, encoding="utf-8")

    input_path = REPO_ROOT / config.get("Path", "input_pdb_name")
    assert input_path.exists(), f"Missing configured PyRosetta input: {input_path}"


def test_run_production_resolves_metadata_tagged_ppi_output_dir() -> None:
    beta_target = next(target for target in PPI_TARGETS if target["name"] == "beta-meander")
    docking_dir = _ppi_docking_dir(beta_target)
    assert docking_dir is not None
    assert docking_dir.name.endswith(
        "__3GT8_raw__legacy_beta_meander_960_1006__full_kinase_domain__prod"
    )
