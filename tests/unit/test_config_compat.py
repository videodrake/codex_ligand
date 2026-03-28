import pytest
from egfr_pipeline.config import compat_allow_legacy_ppi_paths

pytestmark = pytest.mark.unit


def test_compat_allow_legacy_ppi_paths_defaults_false() -> None:
    assert compat_allow_legacy_ppi_paths({}) is False
    assert compat_allow_legacy_ppi_paths({"compat": {}}) is False


def test_compat_allow_legacy_ppi_paths_truthy_strings() -> None:
    for value in ("true", "TRUE", "yes", "on", "1"):
        assert compat_allow_legacy_ppi_paths({"compat": {"allow_legacy_ppi_paths": value}}) is True


def test_compat_allow_legacy_ppi_paths_boolean() -> None:
    assert compat_allow_legacy_ppi_paths({"compat": {"allow_legacy_ppi_paths": True}}) is True
    assert compat_allow_legacy_ppi_paths({"compat": {"allow_legacy_ppi_paths": False}}) is False
